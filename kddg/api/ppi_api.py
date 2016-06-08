#!/usr/bin/python2.4
# encoding: utf-8
"""
ppi_api.py
High-level functions for interacting with the protein-protein interaction sections of the ddG database.

Classes:
BindingAffinityDDGInterface - an class used to interface with the database. Call get_interface to get a user API based on this class.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

import pprint
from io import BytesIO
import os
import sys
import copy
import json
import zipfile
import re
import random
import traceback
import StringIO
import gzip
import shutil
import sqlite3
import cPickle as pickle
import datetime
import time
import getpass

import numpy
from sqlalchemy import and_, or_, func

from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation, residue_type_1to3_map
from klab.fs.fsio import read_file, write_temp_file
from klab.benchmarking.analysis.ddg_binding_affinity_analysis import DBBenchmarkRun as BindingAffinityBenchmarkRun
from klab.bio.alignment import ScaffoldModelChainMapper, DecoyChainMapper
from klab.db.sqlalchemy_interface import row_to_dict, get_or_create_in_transaction, get_single_record_from_query
from klab.stats.misc import get_xy_dataset_statistics_pandas

import kddg.api.schema as dbmodel
from kddg.api.layers import *
from kddg.api.db import ddG, PartialDataException, SanityCheckException
from import_api import json_dumps

import settings # from ddg.ddglib import settings
sys_settings = settings.load()

DeclarativeBase = dbmodel.DeclarativeBase


def get_interface(passwd, username = sys_settings.database.username, hostname = sys_settings.database.hostname, rosetta_scripts_path = None, rosetta_database_path = None, port = sys_settings.database.port):
    '''This is the function that should be used to get a BindingAffinityDDGInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(BindingAffinityDDGInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port)


def get_interface_with_config_file(host_config_name = sys_settings.database.host_config_name, rosetta_scripts_path = None, rosetta_database_path = None, get_interface_factory = get_interface, passed_port = None):
    # Uses ~/.my.cnf to get authentication information
    ### Example .my.cnf (host_config_name will equal myserver):
    ### [clientmyserver]
    ### user=username
    ### password=notmyrealpass
    ### host=server.domain.com
    my_cnf_path = os.path.expanduser(os.path.join('~', '.my.cnf'))
    if not os.path.isfile( os.path.expanduser(my_cnf_path) ):
        raise Exception("A .my.cnf file must exist at: " + my_cnf_path)

    # These three variables must be set in a section of .my.cnf named host_config_name
    user = None
    password = None
    host = None
    port = None
    with open(my_cnf_path, 'r') as f:
        parsing_config_section = False
        for line in f:
            if line.strip() == '[client%s]' % host_config_name:
                parsing_config_section = True
            elif line.strip() == '':
                parsing_config_section = False
            elif parsing_config_section:
                if '=' in line:
                    tokens = line.strip().split('=')
                    key, val = tokens[0], '='.join(tokens[1:]) # values may contain '=' signs
                    key, val = key.strip(), val.strip()
                    if key == 'user':
                        user = val
                    elif key == 'password':
                        password = val
                    elif key == 'host':
                        host = val
                    elif key == 'port':
                        port = int(val)
                else:
                    parsing_config_section = False
    port = passed_port or port or 3306

    if not user or not password or not host:
        raise Exception("Couldn't find host(%s), username(%s), or password in section %s in %s" % (host, user, host_config_name, my_cnf_path) )

    return get_interface_factory(password, username = user, hostname = host, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port)


class BindingAffinityDDGInterface(ddG):
    '''This is the internal API class that should be NOT used to interface with the database.'''


    def __init__(self, passwd = None, username = sys_settings.database.username, hostname = sys_settings.database.hostname, rosetta_scripts_path = None, rosetta_database_path = None, port = sys_settings.database.port, file_content_buffer_size = None):
        super(BindingAffinityDDGInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port, file_content_buffer_size = file_content_buffer_size)
        self.prediction_data_path = self.DDG_db.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionPPIDataPath"')[0]['Value']
        self.unfinished_prediction_ids_cache = {}

    def get_prediction_ids_with_scores(self, prediction_set_id, score_method_id = None):
        '''Returns a set of all prediction_ids that already have an associated score in prediction_set_id
        '''
        score_table = self._get_sqa_prediction_structure_scores_table()
        prediction_table = self.PredictionTable

        if score_method_id != None:
            return set([r['ID'] for r in self.DDG_db.execute_select('''
                SELECT DISTINCT PredictionPPI.ID FROM PredictionPPIStructureScore
                INNER JOIN PredictionPPI
                ON PredictionPPI.ID=PredictionPPIStructureScore.PredictionPPIID
                WHERE PredictionPPI.PredictionSet=%s AND PredictionPPIStructureScore.ScoreMethodID=%s''', parameters=(prediction_set_id, score_method_id))])
        else:
            return set([r['ID'] for r in self.DDG_db.execute_select('''
                SELECT DISTINCT PredictionPPI.ID FROM PredictionPPIStructureScore
                INNER JOIN PredictionPPI
                ON PredictionPPI.ID=PredictionPPIStructureScore.PredictionPPIID
                WHERE PredictionPPI.PredictionSet=%s''', parameters=(prediction_set_id,))])


    def get_unfinished_prediction_ids(self, prediction_set_id):
        '''Returns a set of all prediction_ids that have Status != "done"
        '''
        if prediction_set_id in self.unfinished_prediction_ids_cache:
            return self.unfinished_prediction_ids_cache[prediction_set_id]
        else:
            unfinished_ids = [r.ID for r in self.get_session().query(self.PredictionTable).filter(and_(self.PredictionTable.PredictionSet == prediction_set_id, self.PredictionTable.Status != 'done'))]
            self.unfinished_prediction_ids_cache[prediction_set_id] = unfinished_ids
            return unfinished_ids


    def get_prediction_ids_without_scores(self, prediction_set_id, score_method_id = None):
        all_prediction_ids = [x for x in self.get_prediction_ids(prediction_set_id)]
        all_prediction_ids_set = set()
        for prediction_id in all_prediction_ids:
            all_prediction_ids_set.add( prediction_id )
        scored_prediction_ids_set = self.get_prediction_ids_with_scores(prediction_set_id, score_method_id = score_method_id)
        return [x for x in all_prediction_ids_set.difference(scored_prediction_ids_set)]


    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================


    @informational_pdb
    def get_pdb_chains_for_prediction(self, prediction_id):
        # look up the complex associated with the dataset record for the list of chains
        raise Exception('This needs to be implemented.')


    @informational_pdb
    def get_chain_sets_for_mutatagenesis(self, mutagenesis_id, complex_id = None):
        '''Gets a list of possibilities for the associated complex and calls get_chains_for_mutatagenesis on each.
             e.g. returns {('1KI1', 0) : {'L' : ['A','B'], 'R' : ['C']}, ('12AB', 2) : {'L' : ['L','H'], 'R' : ['A']}, ...}
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.

           This is a useful method for listing the possible complexes to use in a prediction or to determine whether one
           may be missing. and we need to update the database.'''

        pp_mutagenesis = self.DDG_db.execute_select("SELECT * FROM PPMutagenesis WHERE ID=%s", parameters = (mutagenesis_id,))
        # Sanity checks
        assert(len(pp_mutagenesis) == 1)
        if complex_id:
            assert(pp_mutagenesis[0]['PPComplexID'] == complex_id)
        else:
            complex_id = pp_mutagenesis[0]['PPComplexID']

        d = {}
        for pdb_set in self.DDG_db.execute_select("SELECT * FROM PPIPDBSet WHERE PPComplexID=%s AND IsComplex=1", parameters = (complex_id,)):
            pdb_set_number = pdb_set['SetNumber']
            pdb_file_ids = self.DDG_db.execute_select("SELECT DISTINCT PDBFileID FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND SetNumber=%s", parameters = (complex_id, pdb_set_number))
            assert(len(pdb_file_ids) == 1)
            pdb_file_id = pdb_file_ids[0]['PDBFileID']
            d[(pdb_file_id, pdb_set_number)] = self.get_chains_for_mutatagenesis(mutagenesis_id, pdb_file_id, pdb_set_number)
        return d


    @informational_pdb
    def get_chains_for_mutatagenesis(self, mutagenesis_id, pdb_file_id, pdb_set_number, complex_id = None, tsession = None):
        '''Returns a dictionary mapping 'L' to the list of left chains and 'R' to the list of right chains.
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.
        '''

        tsession = tsession or self.get_session() # do not create a new session

        pp_mutagenesis = None
        for r in tsession.execute('''SELECT * FROM PPMutagenesis WHERE ID=:mutagenesis_id''', dict(mutagenesis_id = mutagenesis_id)):
            assert(pp_mutagenesis == None)
            pp_mutagenesis = r

        # Sanity checks
        if complex_id:
            assert(pp_mutagenesis['PPComplexID'] == complex_id)
            pdb_set = None
            for r in tsession.execute('''SELECT * FROM PPIPDBSet WHERE PPComplexID=:complex_id AND SetNumber=:pdb_set_number''', dict(complex_id = complex_id, pdb_set_number = pdb_set_number)):
                assert(pdb_set == None)
                pdb_set = r
            assert(pdb_set['IsComplex'] == 1) # complex structure check
        else:
            complex_id = pp_mutagenesis['PPComplexID']

        pdb_file_id, complex_chains = self.get_bound_pdb_set_details(complex_id, pdb_set_number, pdb_file_id = pdb_file_id, tsession = tsession)
        return complex_chains


    def get_bound_pdb_set_details(self, complex_id, pdb_set_number, pdb_file_id = None, tsession = None):
        '''Returns the pdb_id and complex partner definitions (left PDB chains, right PDB chains) for complexes where all chains share the same PDB structure.'''

        tsession = tsession or self.get_session() # do not create a new session

        assert(complex_id != None and pdb_set_number != None)
        complex_chains = dict(L = [], R = [])
        for c in tsession.execute('''SELECT * FROM PPIPDBPartnerChain WHERE PPComplexID=:complex_id AND SetNumber=:pdb_set_number ORDER BY ChainIndex''', dict(complex_id = complex_id, pdb_set_number = pdb_set_number)):
            if pdb_file_id:
                assert(c['PDBFileID'] == pdb_file_id) # complex structure check
            else:
                pdb_file_id = c['PDBFileID']
            complex_chains[c['Side']].append(c['Chain'])
        assert(complex_chains['L'] and complex_chains['R'])
        assert(len(set(complex_chains['L']).intersection(set(complex_chains['R']))) == 0) # in one unbound case, the same chain appears twice on one side (2CLR_DE|1CD8_AA, may be an error since this was published as 1CD8_AB but 1CD8 has no chain B) but it seems reasonable to assume that a chain should only appear on one side
        return pdb_file_id, complex_chains


    @informational_pdb
    def get_pdb_mutations_for_mutagenesis(self, mutagenesis_id, pdb_file_id, set_number, complex_id = None):
        '''Returns the PDB mutations for a mutagenesis experiment as well as the PDB residue information.'''
        pdb_mutations = []
        for pdb_mutation in self.DDG_db.execute_select('''
            SELECT PPMutagenesisPDBMutation.*, PDBResidue.ResidueType,
            PDBResidue.BFactorMean, PDBResidue.BFactorDeviation,
            PDBResidue.ComplexExposure, PDBResidue.ComplexDSSP, PDBResidue.MonomericExposure, PDBResidue.MonomericDSSP
            FROM
            PPMutagenesisPDBMutation
            INNER JOIN
            PDBResidue ON PPMutagenesisPDBMutation.PDBFileID = PDBResidue.PDBFileID AND PPMutagenesisPDBMutation.Chain = PDBResidue.Chain AND PPMutagenesisPDBMutation.ResidueID = PDBResidue.ResidueID AND PPMutagenesisPDBMutation.WildTypeAA = PDBResidue.ResidueAA
            WHERE PPMutagenesisID=%s AND PDBResidue.PDBFileID=%s AND SetNumber=%s ORDER BY Chain, ResidueID''', parameters=(mutagenesis_id, pdb_file_id, set_number)):
                if complex_id:
                    assert(pdb_mutation['PPComplexID'] == complex_id)
                pdb_mutations.append(pdb_mutation)
        return pdb_mutations


    @sanity_check
    def find_pdb_files_involved_in_multiple_complexes(self):

        known_exceptions = {
            # These need to be checked - 1OYV only has 1 chain besides Subtilisin Carlsberg
            '1OYV' : 2, # Subtilisin Carlsberg bound to: i) domain 1 of its inhibitor; and ii) domain 2 of its inhibitor.
            '1QFW' : 2, # Human chorionic gonadotropin (chains A, B) bound to: i) Fv anti-alpha (chains L, H); and ii) Fv anti-beta (chain M, I).
        }

        d = {}
        for r in self.DDG_db.execute_select('SELECT ID FROM PDBFile ORDER BY ID'):
            pdb_id = r['ID']
            complex_ids = self.search_complexes_by_pdb_id(pdb_id)
            if pdb_id.upper() in known_exceptions:
                assert(len(complex_ids) == known_exceptions[pdb_id])
            else:
                if len(complex_ids) > 1:
                    d[pdb_id] = {'complex_ids' : complex_ids, 'complexes' : {}}
                    for complex_id in complex_ids:
                        d[pdb_id]['complexes'][complex_id] = self.get_complex_details(complex_id)
        if d:
            raise SanityCheckException('Some PDB files are associated with multiple complexes:\n{0}'.format(pprint.pformat(d)))


    @informational_complex
    def search_complexes_by_pdb_id(self, pdb_id):
        '''Returns the list of PPComplexIDs which are related to the PDB ID. Typically this list will be empty or have one
           ID. In rarer cases, the same structure may be used as a structural basis for multiple complexes.'''
        results = self.DDG_db_utf.execute_select('''
            SELECT DISTINCT PPIPDBSet.PPComplexID FROM PPIPDBPartnerChain
            INNER JOIN PPIPDBSet ON PPIPDBPartnerChain.PPComplexID=PPIPDBSet.PPComplexID AND PPIPDBPartnerChain.SetNumber=PPIPDBSet.SetNumber
            WHERE PDBFileID=%s AND IsComplex=1
            ''', parameters=(pdb_id,))
        return [r['PPComplexID'] for r in results]


    @informational_job
    def get_complex_details(self, complex_id):
        results = self.DDG_db_utf.execute_select('SELECT * FROM PPComplex WHERE ID=%s', parameters=(complex_id, ))
        if len(results) == 1:
            return results[0]
        return None


    def _get_dataset_record_with_checks(self, dataset_experiment_id, dataset_id = None):
        if dataset_id:
            de = self.DDG_db_utf.execute_select('SELECT * FROM PPIDataSetDDG WHERE ID=%s AND DataSetID=%s', parameters=(dataset_experiment_id, dataset_id))
            if len(de) != 1:
                raise colortext.Exception('Dataset record #%d does not exist for/correspond to the dataset %s.' % (dataset_experiment_id, dataset_id))
        else:
            de = self.DDG_db_utf.execute_select('SELECT * FROM PPIDataSetDDG WHERE ID=%s', parameters=(dataset_experiment_id,))
            if len(de) != 1:
                raise colortext.Exception('Dataset record #%d does not exist.' % (dataset_experiment_id, ))
        return de[0]


    @informational_job
    def get_job_details(self, prediction_id, include_files = True, truncate_content = None):
        try:
            prediction_record = self.get_session().query(self.PredictionTable).filter(self.PredictionTable.ID == prediction_id).one()
        except Exception, e:
            raise colortext.Exception('No details could be found for prediction #{0} in the database.\n{1}\n{2}'.format(prediction_id, str(e), traceback.format_exc()))

        # mutfile_content = self.create_mutfile(prediction_id)

        # Read the UserPPDataSetExperiment details
        user_dataset_experiment_id = prediction_record.UserPPDataSetExperimentID
        ude_details = self.get_user_dataset_experiment_details(user_dataset_experiment_id)
        assert(ude_details['Mutagenesis']['PPMutagenesisID'] == prediction_record.PPMutagenesisID)

        # Convert the record to dict
        prediction_record = row_to_dict(prediction_record)
        prediction_record['Files'] = {}
        if include_files:
            prediction_record['Files'] = self.get_job_files(prediction_id, truncate_content = truncate_content)

        for k, v in ude_details.iteritems():
            assert(k not in prediction_record)
            prediction_record[k] = v
        return prediction_record


    @informational_job
    def get_dataset_experiment_details(self, dataset_experiment_id, dataset_id = None):
        de = self._get_dataset_record_with_checks(dataset_experiment_id, dataset_id = dataset_id)
        PDBFileID = de['PDBFileID']
        PPMutagenesisID = de['PPMutagenesisID']
        ComplexID = self.DDG_db.execute_select('SELECT PPComplexID FROM PPMutagenesis WHERE ID=%s', parameters=(PPMutagenesisID,))[0]['PPComplexID']
        SetNumber = None

        # todo: this is a nasty hack due to the fact that we do not currently store the SetNumber and PPComplexID in the PPIDataSetDDG table. See ticket:1457.
        pdb_sets = self.DDG_db.execute_select('SELECT * FROM PPIPDBSet WHERE PPComplexID=%s AND IsComplex=1', parameters=(ComplexID,))
        if len(pdb_sets) > 1:
            probable_sets = self.DDG_db.execute_select('SELECT DatabaseKey FROM PPIDatabaseComplex WHERE DatabaseName LIKE "%%SKEMPI%%" AND DatabaseKey LIKE "%%%s%%" AND PPComplexID=%s' % (PDBFileID, ComplexID))
            assert(len(probable_sets) == 1)
            match_pdb_chains = sorted(list(''.join(probable_sets[0]['DatabaseKey'].split('_')[1:])))

            pdb_sets = {}
            for set_record in self.DDG_db.execute_select('SELECT * FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND PDBFileID=%s', parameters=(ComplexID, PDBFileID)):
                pdb_sets[set_record['SetNumber']] = pdb_sets.get(set_record['SetNumber'], [])
                pdb_sets[set_record['SetNumber']].append(set_record['Chain'])
                pdb_sets[set_record['SetNumber']] = sorted(pdb_sets[set_record['SetNumber']])

            hits = []
            for k, v in pdb_sets.iteritems():
                if v == match_pdb_chains:
                    hits.append(k)
            if not len(hits) == 1:
                raise Exception('Error: multiple possible PDB sets for dataset record #%d and PPMutagenesisID=%s.' % (dataset_experiment_id, PPMutagenesisID))
            SetNumber = hits[0]
        elif len(pdb_sets) == 0:
            raise Exception('Error: no possible PDB sets for dataset record #%d and PPMutagenesisID=%s.' % (dataset_experiment_id, PPMutagenesisID))
        else:
            SetNumber = pdb_sets[0]['SetNumber']

        pdb_mutations = self.get_pdb_mutations_for_mutagenesis(PPMutagenesisID, PDBFileID, SetNumber, complex_id = ComplexID)

        d = dict(
            _DataSetID = de['ID'],
            RecordID = de['RecordNumber'],
            PublishedDDG = de['PublishedDDG'],
            PDBFileID = PDBFileID,
            DerivedMutation = de['RecordIsDerivative'] == 1,
            PossiblyBadRecord = de['PossibleError'] == 1,
            Notes = [de['Remark'], de['CorrectionRemark']],
            Mutagenesis = dict(
                PPMutagenesisID = PPMutagenesisID,
            ),
            Complex = self.get_complex_details(ComplexID),
            Structure = dict(
                PDBFileID = PDBFileID,
                SetNumber = SetNumber,
                Partners = self.get_chains_for_mutatagenesis(PPMutagenesisID, PDBFileID, SetNumber, complex_id = ComplexID),
            ),
            PDBMutations = pdb_mutations,
        )
        if de['PublishedPDBFileID'] != PDBFileID:
            d['Notes'].append("The PDB ID was changed by Shane O'Connor from %s to %s." % (de['PublishedPDBFileID'], PDBFileID))
        d['Notes'] = '. '.join([x for x in d['Notes'] if x])
        d['ExperimentalDDGs'] = self.get_ddg_values_for_dataset_record(dataset_experiment_id, dataset_id = dataset_id)
        d['DDG'] = sum([((e.get('Positive') or {}).get('DDG', 0) - (e.get('Negative') or {}).get('DDG', 0)) for e in d['ExperimentalDDGs']])
        # todo: add SCOPe class, Pfam domain
        return d


    def _get_ddg_values_for_dataset_record(self, dataset_experiment_id, dataset_id = None):
        de = self._get_dataset_record_with_checks(dataset_experiment_id, dataset_id = dataset_id)
        ddg_pairs = self.DDG_db.execute_select('SELECT PositiveDependentPPIDDGID, NegativeDependentPPIDDGID FROM PPIDataSetDDGSource WHERE PPIDataSetDDGID=%s', parameters=(dataset_experiment_id,))
        assert(ddg_pairs)
        ddgs = []
        for ddg_pair in ddg_pairs:
            paired_record = {'Positive' : None, 'Negative' : None}
            if ddg_pair['PositiveDependentPPIDDGID']:
                positive_record = self.DDG_db.execute_select('SELECT * FROM PPIDDG WHERE ID=%s', parameters=(ddg_pair['PositiveDependentPPIDDGID'],))[0]
                paired_record['Positive'] = dict(
                    DDG = positive_record['DDG'],
                    LocationOfValueInPublication = positive_record['LocationOfValueInPublication'],
                    Publication = positive_record['Publication'],
                    Temperature = positive_record['Temperature'],
                    pH = positive_record['pH'],
                )
            if ddg_pair['NegativeDependentPPIDDGID']:
                negative_record = self.DDG_db.execute_select('SELECT * FROM PPIDDG WHERE ID=%s', parameters=(ddg_pair['NegativeDependentPPIDDGID'],))[0]
                paired_record['Negative'] = dict(
                    DDG = negative_record['DDG'],
                    LocationOfValueInPublication = negative_record['LocationOfValueInPublication'],
                    Publication = negative_record['Publication'],
                    Temperature = negative_record['Temperature'],
                    pH = negative_record['pH'],
                )
            ddgs.append(paired_record)
        return ddgs


    @informational_job
    def get_user_dataset_experiment_details(self, user_dataset_experiment_id, user_dataset_id = None):
        if user_dataset_id:
            colortext.ppurple('PRE-SELECT')
            ude = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(user_dataset_experiment_id, user_dataset_id))
            colortext.ppurple('POST-SELECT')
            if len(ude) != 1:
                raise colortext.Exception('User dataset experiment %d does not exist for/correspond to the user dataset %s.' % (user_dataset_experiment_id, user_dataset_id))
        else:
            ude = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s', parameters=(user_dataset_experiment_id,))
            if len(ude) != 1:
                raise colortext.Exception('User dataset experiment %d does not exist.' % (user_dataset_experiment_id, ))
        ude = ude[0]
        user_dataset_id = ude['UserDataSetID']
        assert(ude['IsComplex'] == 1)

        pdb_mutations = self.get_pdb_mutations_for_mutagenesis(ude['PPMutagenesisID'], ude['PDBFileID'], ude['SetNumber'], complex_id = ude['PPComplexID'])
        return dict(
            Mutagenesis = dict(
                PPMutagenesisID = ude['PPMutagenesisID'],
            ),
            Complex = self.get_complex_details(ude['PPComplexID']),
            Structure = dict(
                PDBFileID = ude['PDBFileID'],
                SetNumber = ude['SetNumber'],
                Partners = self.get_chains_for_mutatagenesis(ude['PPMutagenesisID'], ude['PDBFileID'], ude['SetNumber'], complex_id = ude['PPComplexID']),
            ),
            PDBMutations = pdb_mutations,
        )


    def _export_dataset(self, dataset_id):
        '''Returns a dict containing the dataset information.'''
        dataset_record = self.DDG_db.execute_select('SELECT * FROM DataSet WHERE ID=%s', parameters=(dataset_id,))
        if not dataset_record:
            raise Exception('Dataset %s does not exist in the database.' % dataset_id)
        dataset_record = dataset_record[0]
        if dataset_record['DatasetType'] != 'Binding affinity' and dataset_record['DatasetType'] != 'Protein stability and binding affinity':
            raise Exception('The dataset %s does not contain any binding affinity data..' % dataset_id)

        # Read the UserPPDataSetExperiment details
        data = []
        ref_ids = set()
        for dataset_ddg in self.DDG_db.execute_select('SELECT * FROM PPIDataSetDDG WHERE DataSetID=%s ORDER BY Section, RecordNumber', parameters=(dataset_id,)):
            de_details = self.get_dataset_experiment_details(dataset_ddg['ID'], dataset_id)
            for ddg_pair in de_details['ExperimentalDDGs']:
                if ddg_pair['Positive']: ref_ids.add(ddg_pair['Positive']['Publication'])
                if ddg_pair['Negative']: ref_ids.add(ddg_pair['Negative']['Publication'])
            data.append(de_details)

        references = {}
        for ref_id in sorted(ref_ids):
            references[ref_id] = self.get_publication(ref_id)

        return dict(
            Data = data,
            References = references
            )


    @informational_job
    def export_dataset_to_csv(self, dataset_id):
        '''Returns the dataset information in CSV format.'''
        dataset_set = self._export_dataset(dataset_id)['Data']
        lines = ['\t'.join(['Record #', 'Mutagenesis #', 'Partner 1', 'Partner 2', 'PDB ID', 'Partner 1 chains', 'Partner 2 chains', 'Mutations', 'DDG', 'PublishedDDG', 'IsDerivedMutation'])]
        for record in dataset_set:
            line = '\t'.join([
                str(record['RecordID']),
                str(record['Mutagenesis']['PPMutagenesisID']),
                record['Complex']['LShortName'],
                record['Complex']['RShortName'],
                record['PDBFileID'],
                ','.join(sorted(record['Structure']['Partners']['L'])),
                ','.join(sorted(record['Structure']['Partners']['R'])),
                ','.join(['%s:%s%s%s' % (m['Chain'], m['WildTypeAA'], m['ResidueID'], m['MutantAA']) for m in record['PDBMutations']]),
                str(record['DDG']),
                str(record['PublishedDDG']),
                str(int(record['DerivedMutation'])),
            ])
            lines.append(line)
        return ('\n'.join(lines)).encode('utf8', 'replace')


    @informational_job
    def get_predictions_experimental_details(self, prediction_id, userdatset_experiment_ids_to_subset_ddgs = None, include_files = False, reference_ids = set(), include_experimental_data = True):

        details = self.get_job_details(prediction_id, include_files = include_files)

        # Sanity checks and redundancy removal
        PPMutagenesisID = details['PPMutagenesisID']
        ComplexID = details['Complex']['ID']
        chains = set([item for sublist in [v for k, v in details['Structure']['Partners'].iteritems()] for item in sublist])
        PDBFileID = details['Structure']['PDBFileID']
        SetNumber = details['Structure']['SetNumber']
        for m in details['PDBMutations']:
            assert(m['PPMutagenesisID'] == PPMutagenesisID)
            del m['PPMutagenesisID']
            assert(ComplexID == m['PPComplexID'])
            del m['PPComplexID']
            assert(PDBFileID == m['PDBFileID'])
            del m['PDBFileID']
            assert(SetNumber == m['SetNumber'])
            del m['SetNumber']
            assert(m['Chain'] in chains)
        assert(details['Mutagenesis']['PPMutagenesisID'] == PPMutagenesisID)
        del details['Mutagenesis']

        # Add the DDG values for the related analysis sets
        user_dataset_experiment_id = details['UserPPDataSetExperimentID']
        if include_experimental_data:
            userdatset_experiment_ids_to_subset_ddgs = userdatset_experiment_ids_to_subset_ddgs or self.get_experimental_ddgs_by_analysis_set(user_dataset_experiment_id, reference_ids = reference_ids)
            assert('DDG' not in details)
            details['DDG'] = userdatset_experiment_ids_to_subset_ddgs[user_dataset_experiment_id]
        else:
            details['DDG'] = None

        return details


    @informational_job
    def get_experimental_ddgs_by_analysis_set(self, user_dataset_experiment_id = None, reference_ids = set()):

        # Determine the set of analysis sets
        userdatset_experiment_ids_to_subset_ddgs = {}
        analysis_sets = [r['Subset'] for r in self.DDG_db.execute_select('SELECT DISTINCT Subset FROM UserPPAnalysisSet')]

        # Query the database, restricting to one user_dataset_experiment_id if passed
        parameters = None
        qry = '''
            SELECT UserPPAnalysisSet.*,
            (IFNULL(PositiveDDG.DDG, 0) - IFNULL(NegativeDDG.DDG, 0)) AS ExperimentalDDG,
            IF(ISNULL(NegativeDDG.DDG), 0, 1) AS DerivedMutation,
            PositiveDDG.PPMutagenesisID, PositiveDDG.Publication AS PositiveDDGPublication, PositiveDDG.DDG as PositiveDDGValue,
            NegativeDDG.PPMutagenesisID, NegativeDDG.Publication AS NegativeDDGPublication, NegativeDDG.DDG as NegativeDDGValue
            FROM UserPPAnalysisSet
            LEFT JOIN PPIDDG AS PositiveDDG ON PositiveDependentPPIDDGID=PositiveDDG.ID
            LEFT JOIN PPIDDG AS NegativeDDG ON NegativeDependentPPIDDGID=NegativeDDG.ID'''
        if user_dataset_experiment_id != None:
            qry += ' WHERE UserPPAnalysisSet.UserPPDataSetExperimentID=%s'
            parameters = (user_dataset_experiment_id,)
        results = self.DDG_db.execute_select(qry, parameters)

        # Return the mapping
        for r in results:
            if not userdatset_experiment_ids_to_subset_ddgs.get(r['UserPPDataSetExperimentID']):
                d = dict.fromkeys(analysis_sets, None)
                for analysis_set in analysis_sets:
                    d[analysis_set] = {}
                userdatset_experiment_ids_to_subset_ddgs[r['UserPPDataSetExperimentID']] = d

            userdatset_experiment_ids_to_subset_ddgs[r['UserPPDataSetExperimentID']][r['Subset']] = userdatset_experiment_ids_to_subset_ddgs[r['UserPPDataSetExperimentID']][r['Subset']] or dict(
                Cases = set(),
                DDGs = [],
                IsDerivedValue = False,
                MeanDDG = None
            )

            # Store the references IDs
            reference = None
            if r['PositiveDDGPublication'] and r['NegativeDDGPublication']:
                reference = r['PositiveDDGPublication'] + ', ' + r['NegativeDDGPublication']
                reference_ids.add(r['PositiveDDGPublication'])
                reference_ids.add(r['NegativeDDGPublication'])
            elif r['PositiveDDGPublication']:
                reference = r['PositiveDDGPublication']
                reference_ids.add(r['PositiveDDGPublication'])
            elif r['NegativeDDGPublication']:
                reference = r['NegativeDDGPublication']
                reference_ids.add(r['NegativeDDGPublication'])

            record_d = userdatset_experiment_ids_to_subset_ddgs[r['UserPPDataSetExperimentID']][r['Subset']]
            record_d['Cases'].add((r['Subset'], r['Section'], r['RecordNumber']))
            record_d['DDGs'].append({'Value' : r['ExperimentalDDG'], 'IsDerivedValue' : r['DerivedMutation'], 'Reference' : reference})
            record_d['IsDerivedValue'] = record_d['IsDerivedValue'] or r['DerivedMutation']

        # Calculate the mean of the DDG values
        # Note: Based on experience, summing in Python over small lists can be faster than creating temporary numpy arrays due to the array creation overhead
        for k, v in userdatset_experiment_ids_to_subset_ddgs.iteritems():
            for subset, subset_ddgs in v.iteritems():
                if subset_ddgs:
                    num_points = len(subset_ddgs['DDGs'])
                    if num_points > 1:
                        subset_ddgs['MeanDDG'] = sum([float(ddg['Value'])for ddg in subset_ddgs['DDGs']]) / float(num_points)
                    else:
                        # Avoid unnecessary garbage creation and division
                        subset_ddgs['MeanDDG'] = subset_ddgs['DDGs'][0]['Value']

        return userdatset_experiment_ids_to_subset_ddgs


    @informational_job
    def export_prediction_cases_to_json(self, prediction_set_id, retrieve_references = True):
        print('This will probably break - I need to dump datetime.datetime objects to ISO strings.')
        return json_dumps(self.get_prediction_set_case_details(prediction_set_id, retrieve_references = retrieve_references))

    @informational_job
    def export_prediction_cases_to_pickle(self, prediction_set_id, retrieve_references = True):
        return pickle.dumps(self.get_prediction_set_case_details(prediction_set_id, retrieve_references = retrieve_references))

    ##### Public API: Rosetta-related functions


    @job_input
    def create_resfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    @job_input
    def create_mutfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    ###########################################################################################
    ## Prediction layer
    ##
    ## This part of the API is responsible for inserting prediction jobs in the database via
    ## the trickle-down proteomics paradigm.
    ###########################################################################################


    #== Job creation/management API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via the
    # trickle-down proteomics paradigm.


    #   PredictionSet interface

    @job_creator
    def add_prediction_set(self, prediction_set_id, halted = True, priority = 5, batch_size = 40, allow_existing_prediction_set = False,
                                 series_name = None, series_color = 'ff0000', series_alpha = 1.0, description = None):
        return super(BindingAffinityDDGInterface, self).add_prediction_set(prediction_set_id, halted = halted, priority = priority, batch_size = batch_size, allow_existing_prediction_set = allow_existing_prediction_set, contains_protein_stability_predictions = False, contains_binding_affinity_predictions = True, series_name = series_name, series_color = series_color, series_alpha = series_alpha, description = description)


    @job_creator
    def add_development_protocol_command_lines(self, prediction_set_id, protocol_name, application, template_command_line, rosetta_script_file = None):
        dev_protocol_id = self._get_dev_protocol_id(protocol_name)
        if not dev_protocol_id:
            dev_protocol_id = self._create_dev_protocol(protocol_name, application, template_command_line)

        rosetta_script = None
        if rosetta_script_file:
            with open(rosetta_script_file, 'r') as f:
                rosetta_script = f.read()

        prediction_ids = self.get_prediction_ids(prediction_set_id)

        # All functions within the next with block should use the same database cursor.
        # The commands then function as parts of a transaction which is rolled back if errors occur within the block
        # or else is committed.
        file_content_id = None

        tsession = self.get_session(new_session = True)
        try:
            for prediction_id in prediction_ids:
                prediction_record = tsession.query(dbmodel.PredictionPPI).filter(dbmodel.PredictionPPI.ID == prediction_id)
                prediction_record.DevelopmentProtocolID = dev_protocol_id
                tsession.flush()
                if rosetta_script:
                    # Add the Rosetta script to the database, not using cursor
                    file_content_id = self._add_prediction_file(tsession, prediction_id, rosetta_script, os.path.basename(rosetta_script_file), 'RosettaScript', 'RosettaScript', 'Input', rm_trailing_line_whitespace = True, forced_mime_type = 'text/xml', file_content_id = file_content_id)
            tsession.commit()
            tsession.close()
        except Exception, e:
            colortext.error('Failure: {0}.'.format(str(e)))
            colortext.error(traceback.format_exc())
            tsession.rollback()
            tsession.close()

    @job_creator
    def add_job(self, tsession, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, keep_all_lines = False, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = None, suppress_warnings = False):
        '''This function inserts a prediction into the database.
            The parameters define:
                - the prediction set id used to group this prediction with other predictions for analysis;
                - the protocol to be used to run the prediction;
                - the set of mutations and PDB complex associated with the mutagenesis experiment;
                - whether HETATM lines are to be kept or not.
                - additional Rosetta flags e.g. "-ignore_zero_occupancy false" used to determine the mapping from PDB to Rosetta numbering. These flags should correspond to those used in the protocol otherwise errors could occur.
            We strip the PDB based on the chains defined by the complex and keep_all_lines and keep_hetatm_lines and store the PDB in the database.
            Next, the mapping from Rosetta numbering to PDB numbering is determined and stored in the database.
            Then, the appropriate input files e.g. resfiles or mutfiles are generated and stored in the database.
            Finally, we add the prediction record and associate it with the generated files.'''
        return self._add_job(tsession, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_all_lines = keep_all_lines, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache, suppress_warnings = suppress_warnings)


    @job_creator
    def add_job_by_user_dataset_record(self, prediction_set_id, user_dataset_name, user_dataset_experiment_id, protocol_id, extra_rosetta_command_flags = None, keep_all_lines = False, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = None, suppress_warnings = False, tsession = None, allowed_user_datasets = None):
        '''Add a prediction job based on a user dataset record. This is typically called during add_prediction_run rather than directly by the user.
           user_dataset_name is implied by user_dataset_experiment_id but we include it for sanity checking errors in data-entry.

           The extra_rosetta_command_flags variable is used to add additional flags e.g. "-ignore_zero_occupancy false". These should be added if they are used in the protocol.'''

        new_session = False
        if not tsession:
            new_session = True
            tsession = self.get_session(new_session = True)

        if not allowed_user_datasets:
            allowed_user_datasets = self.get_defined_user_datasets(tsession)

        try:
            user_dataset_id = allowed_user_datasets[user_dataset_name]['ID']
        except:
            raise colortext.Exception('The user dataset "%s" does not exist for this API.' % user_dataset_name)

        udse_table = self._get_sqa_user_dataset_experiment_table()
        ude = None
        for r in tsession.execute('''SELECT * FROM UserPPDataSetExperiment WHERE ID=:udse AND UserDataSetID=:uds''', dict(udse = user_dataset_experiment_id, uds = user_dataset_id)):
            assert(not ude)
            ude = r
        if not ude:
            raise colortext.Exception('User dataset experiment {0} does not exist for/correspond to this user dataset.'.format(user_dataset_experiment_id))

        prediction_id = self._add_job(tsession, prediction_set_id, protocol_id, ude.PPMutagenesisID, ude.PPComplexID, ude.PDBFileID, ude.SetNumber, extra_rosetta_command_flags = extra_rosetta_command_flags, user_dataset_experiment_id = user_dataset_experiment_id, keep_all_lines = keep_all_lines, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache, suppress_warnings = suppress_warnings)
        if new_session:
            tsession.close()
        return prediction_id


    @job_creator
    def merge_prediction_run(self, from_prediction_set_id, to_prediction_set_id, create_if_does_not_exist = True, series_color = 'ff0000', description = None):

        # Start a new transaction
        tsession = self.get_session(new_session = True)
        try:
            # Look up the source prediction set details
            try:
                from_prediction_set = self.get_session().query(dbmodel.PredictionSet).filter(dbmodel.PredictionSet.ID == from_prediction_set_id).one()
            except Exception, e:
                print(str(e))
                print(traceback.format_exc())
                raise Exception('Could not retrieve details for source PredictionSet "{0}".'.format(from_prediction_set_id))

            # Look up or create the target prediction set details
            try:
                to_prediction_set_details = self.get_session().query(dbmodel.PredictionSet).filter(dbmodel.PredictionSet.ID == to_prediction_set_id).one()
            except:
                if create_if_does_not_exist:
                    prediction_set_dict = row_to_dict(from_prediction_set)
                    prediction_set_dict['ID'] = to_prediction_set_id
                    prediction_set_dict['EntryDate'] = datetime.datetime.now()
                    prediction_set_dict['Description'] = description or 'Clone of {0}'.format(from_prediction_set_id)
                    db_ligand_synonym = get_or_create_in_transaction(tsession, dbmodel.PredictionSet, prediction_set_dict)
                else:
                    raise Exception('Could not retrieve details for target PredictionSet "{0}". To create a new PredictionSet, set create_if_does_not_exist to True.'.format(to_prediction_set_id))

            # Create prediction records
            num_predictions = len(from_prediction_set.ppi_predictions)
            colortext.message('Merging/cloning prediction set.'.format())
            c = 1
            for prediction in from_prediction_set.ppi_predictions:
                colortext.wyellow('{0}/{1}: Prediction #{2}\r'.format(c, num_predictions, str(prediction.ID).ljust(15)))
                c += 1

                # Add a prediction record if it does not already exist
                new_prediction_id = None
                if self.get_session().query(self.PredictionTable).filter(and_(
                        self.PredictionTable.PredictionSet == to_prediction_set_id,
                        self.PredictionTable.UserPPDataSetExperimentID == prediction.UserPPDataSetExperimentID,
                        self.PredictionTable.ProtocolID == prediction.ProtocolID)).count() > 0:
                    continue
                else:
                    new_prediction = prediction.clone(to_prediction_set_id)
                    tsession.add(new_prediction)
                    tsession.flush()
                    new_prediction_id = new_prediction.ID

                # Add the prediction file records. The underlying FileContent tables will already exist.
                for prediction_file in prediction.files:
                    new_prediction_file = prediction_file.clone(new_prediction_id)
                    tsession.add(new_prediction_file)
                    tsession.flush()

            print('\nSuccess.\n')
            tsession.commit()
            tsession.close()
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise


    @job_creator
    def add_prediction_run(self, prediction_set_id, user_dataset_name, extra_rosetta_command_flags = None, protocol_id = None, tagged_subset = None, keep_all_lines = False, keep_hetatm_lines = False, input_files = {}, quiet = False, test_only = False, only_single_mutations = False, short_run = False, test_run_first = True, show_full_errors = False, suppress_warnings = False):
        '''Adds all jobs corresponding to a user dataset e.g. add_prediction_run("my first run", "AllBindingAffinityData", tagged_subset = "ZEMu").
           If keep_hetatm_lines is False then all HETATM records for the PDB prediction chains will be removed. Otherwise, they are kept.
           input_files is a global parameter for the run which is generally empty. Any files added here will be associated to all predictions in the run.

           The extra_rosetta_command_flags parameter e.g. "-ignore_zero_occupancy false" is used to determine the mapping
           from PDB to Rosetta numbering. These flags should correspond to those used in the protocol otherwise errors could occur.

           Returns False if no predictions were added to the run else return True if all predictions (and there were some) were added to the run.'''

        # For test runs, this number of predictions will be created
        short_run_limit = 100

        # Create a new session
        tsession = self.get_session(new_session = True)
        try:
            # Check preconditions
            assert(not(input_files)) # todo: do something with input_files when we use that here - call self._add_file_content, associate the filenames with the FileContent IDs, and pass that dict to add_job which will create PredictionPPIFile records
            assert(only_single_mutations == False) # todo: support this later? it may make more sense to just define new UserDataSets
            allowed_user_datasets = self._add_prediction_run_preconditions(tsession, prediction_set_id, user_dataset_name, tagged_subset)

            # Get the list of user dataset experiment records
            user_dataset_experiments = self.get_user_dataset_experiments(tsession, user_dataset_name, tagged_subset = tagged_subset)
            assert(set([u.IsComplex for u in user_dataset_experiments]) == set([1,]))
            num_user_dataset_experiments = user_dataset_experiments.count()
            if not user_dataset_experiments:
                tsession.close()
                return False

            # Count the number of individual PDB files
            pdb_file_ids = set([u.PDBFileID for u in user_dataset_experiments])
            tagged_subset_str = ''
            if not quiet:
                if tagged_subset:
                    tagged_subset_str = 'subset "%s" of ' % tagged_subset

            # Create a cache to speed up job insertion
            pdb_residues_to_rosetta_cache = {}

            t1 = time.time()

            # Run one query over the PredictionSet
            result_set = None
            if protocol_id:
                result_set = tsession.execute('''SELECT * FROM PredictionPPI WHERE PredictionSet=:prediction_set AND ProtocolID=:protocol_id''', dict(prediction_set = prediction_set_id, protocol_id = protocol_id))
            else:
                result_set = tsession.execute('''SELECT * FROM PredictionPPI WHERE PredictionSet=:prediction_set AND ProtocolID IS NULL''', dict(prediction_set = prediction_set_id))
            existing_results = set()
            for r in result_set:
                existing_results.add(r['UserPPDataSetExperimentID'])

            # Test all predictions before creating records
            if test_only or test_run_first:
                if not quiet:
                    colortext.message('Testing %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (num_user_dataset_experiments, len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))

                # Progress counter setup
                count, records_per_dot = 0, 50
                showprogress = not(quiet) and num_user_dataset_experiments > 300
                if showprogress: print("|" + ("*" * (int(num_user_dataset_experiments/records_per_dot)-2)) + "|")
                for ude in user_dataset_experiments:
                    # If the mutagenesis already exists in the prediction set, do not test it again
                    if not(ude.ID in existing_results):
                        # Test the prediction setup
                        prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude.ID, protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_all_lines = keep_all_lines, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = True, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache, suppress_warnings = suppress_warnings, tsession = tsession, allowed_user_datasets = allowed_user_datasets)

                    # Progress counter
                    count += 1
                    if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
                    if short_run and count >= short_run_limit: break
                if not quiet: print('')

            t2 = time.time()
            print('Time taken for dry run: {0}s.'.format(t2 - t1))

            if test_only:
                tsession.rollback()
                tsession.close()
                return

            # Progress counter setup
            failed_jobs = {}
            if not quiet:
                colortext.message('Adding %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (num_user_dataset_experiments, len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
            count, records_per_dot = 0, 50
            showprogress = not(quiet) and num_user_dataset_experiments > 300
            if showprogress: print("|" + ("*" * (int(num_user_dataset_experiments/records_per_dot)-2)) + "|")

            t1 = time.time()
            time_to_ignore = 0

            # Add the individual predictions
            for ude in user_dataset_experiments:

                # If the mutagenesis already exists in the prediction set, do not add it again
                if not(ude.ID in existing_results):
                    t3 = time.time()
                    try:
                        # Add the prediction
                        user_dataset_id = allowed_user_datasets[user_dataset_name]['ID']
                        prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude.ID, protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_all_lines = keep_all_lines, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = False, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache, suppress_warnings = suppress_warnings, tsession = tsession, allowed_user_datasets = allowed_user_datasets)
                    except Exception, e:
                        time_to_ignore += time.time() - t3
                        user_dataset_id = allowed_user_datasets[user_dataset_name]['ID']

                        ude_record = None
                        for r in tsession.execute('SELECT * FROM UserPPDataSetExperiment WHERE ID=:ude_id AND UserDataSetID=:uds_id', dict(ude_id = ude.ID, uds_id = user_dataset_id)):
                            assert(ude_record == None)
                            ude_record = r
                        assert(ude_record['ID'] == ude.ID)
                        colortext.error('Adding the prediction for UserPPDataSetExperimentID %(ID)d failed (%(PDBFileID)s).' % ude_record)
                        failed_jobs[ude_record['PDBFileID']] = failed_jobs.get(ude_record['PDBFileID'], 0)
                        failed_jobs[ude_record['PDBFileID']] += 1
                        if show_full_errors:
                            print(e)
                            print(traceback.format_exc())

                # Progress counter
                count += 1
                if showprogress and count % records_per_dot == 0: colortext.write(".", "green", flush = True)
                if short_run and count >= short_run_limit: break
            t2 = time.time()
            print('Time taken for actual run: {0}s.'.format(t2 - t1 - time_to_ignore))

            if failed_jobs:
                colortext.error('Some jobs failed to run:\n%s' % pprint.pformat(failed_jobs))
            if not quiet: print('')

            print('Success')
            tsession.commit()
            tsession.close()
            return True
        except Exception, e:
            print(str(e))
            print(traceback.format_exc())
            tsession.rollback()
            tsession.close()
            raise

    def _create_pdb_residues_to_rosetta_cache_mp(self, pdb_residues_to_rosetta_cache, pdb_file_id, pdb_chains_to_keep, extra_rosetta_command_flags, keep_hetatm_lines):
        # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
        raise Exception('Shane should finish this and add keep_all_lines')
        assert(type(pdb_residues_to_rosetta_cache) == None)# use the manager dictproxy)
        pdb_file = self.DDG_db.execute_select("SELECT * FROM PDBFile WHERE ID=%s", parameters = (pdb_file_id,))
        p = PDB(pdb_file[0]['Content'])
        p.strip_to_chains(list(pdb_chains_to_keep))
        if not keep_hetatm_lines:
            p.strip_HETATMs()
        stripped_p = PDB('\n'.join(p.lines))

        stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, self.rosetta_database_path, extra_command_flags = extra_rosetta_command_flags)
        atom_to_rosetta_residue_map = stripped_p.get_atom_sequence_to_rosetta_json_map()
        rosetta_to_atom_residue_map = stripped_p.get_rosetta_sequence_to_atom_json_map()
        cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, self.rosetta_database_path, extra_rosetta_command_flags)
        pdb_residues_to_rosetta_cache[cache_key] = dict(
            stripped_p = stripped_p,
            atom_to_rosetta_residue_map = atom_to_rosetta_residue_map,
            rosetta_to_atom_residue_map = rosetta_to_atom_residue_map)


    @job_creator
    def add_prediction_run_mp(self, prediction_set_id, user_dataset_name, extra_rosetta_command_flags = None, protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, input_files = {}, quiet = False, only_single_mutations = False, short_run = False, show_full_errors = False):
        '''This is a multiprocessing version of add_prediction_run and should be used in favor of that function as it runs faster.

           It takes advantage of parallelism at two points - creating the stripped PDB files and mutfiles for input and
           inserting the jobs (MD5 is run multiple times for each job).
           It was simple/quicker to write this as a 2-step method with a bottleneck in the middle i.e. it waits until all
           stripped PDB files are generated before adding the jobs.
           This could be made even more parallel by removing the bottleneck i.e. the process which strips the PDBs could
           then call _add_job immediately rather than waiting for the other calls to _create_pdb_residues_to_rosetta_cache_mp
           to complete.
           '''

        # Check preconditions
        assert(keep_all_lines)
        assert(suppress_warnings)
        assert(tsession)
        assert(not(input_files)) # todo: do something with input_files when we use that here - call self._add_file_content, associate the filenames with the FileContent IDs, and pass that dict to add_job which will create PredictionPPIFile records
        assert(only_single_mutations == False) # todo: support this later? it may make more sense to just define new UserDataSets
        self._add_prediction_run_preconditions(tsession, prediction_set_id, user_dataset_name, tagged_subset)

        # Get the list of user dataset experiment records
        user_dataset_experiments = self.get_user_dataset_experiments(tsession, user_dataset_name, tagged_subset = tagged_subset)
        assert(set([u['IsComplex'] for u in user_dataset_experiments]) == set([1,]))
        if not user_dataset_experiments:
            return False

        # Count the number of individual PDB files
        pdb_file_ids = set([u['PDBFileID'] for u in user_dataset_experiments])
        tagged_subset_str = ''
        if not quiet:
            if tagged_subset:
                tagged_subset_str = 'subset "%s" of ' % tagged_subset

        # Create a cache to speed up job insertion
        #todo: start back here  pdb_residues_to_rosetta_cache = manager dictproxy

        # Create the stripped PDBs and residue maps in parallel using the multiprocessing module
    #todo: write this function on Monday - get_user_dataset_pdb_partner_chains should return a set (<list of {'id' : pdb_file_id, 'L' : <list of chain ids>, , 'R' : <list of chain ids>} dicts>)
        pdb_partner_chains = self.get_user_dataset_pdb_partner_chains(user_dataset_name, tagged_subset = tagged_subset)
        #todo: start back here for ppc in pdb_partner_chains:

        #todo: start back here    apply_async self._create_pdb_residues_to_rosetta_cache_mp(pdb_residues_to_rosetta_cache, ppc['id'], set(ppc['L'] + ppc['R']), extra_rosetta_command_flags, keep_hetatm_lines)
        #todo: start back here .join()

        # Progress counter setup
        failed_jobs = {}
        if not quiet:
            colortext.message('Adding %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
        count, records_per_dot = 0, 50
        showprogress = not(quiet) and len(user_dataset_experiments) > 300
        if showprogress: print("|" + ("*" * (int(len(user_dataset_experiments)/records_per_dot)-2)) + "|")

        # Add the individual predictions
        for ude in user_dataset_experiments:

            # If the mutagenesis already exists in the prediction set, do not add it again
            if protocol_id:
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
            else:
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID IS NULL", parameters=(prediction_set_id, ude['ID']))
            if len(existing_results) == 0:
                # Add the prediction
                try:
                    user_dataset_id = self.get_defined_user_datasets(tsession)[user_dataset_name]['ID']
                    prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_all_lines = keep_all_lines, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = False, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache, suppress_warnings = suppress_warnings)
                except Exception, e:
                    user_dataset_id = self.get_defined_user_datasets(tsession)[user_dataset_name]['ID']
                    ude_record = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(ude['ID'], user_dataset_id))
                    ude_record = ude_record[0]
                    assert(ude_record['ID'] == ude['ID'])
                    colortext.error('Adding the prediction for UserPPDataSetExperimentID %(ID)d failed (%(PDBFileID)s).' % ude_record)
                    failed_jobs[ude_record['PDBFileID']] = failed_jobs.get(ude_record['PDBFileID'], 0)
                    failed_jobs[ude_record['PDBFileID']] += 1
                    if show_full_errors:
                        print(e)
                        print(traceback.format_exc())

            # Progress counter
            count += 1
            if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
            if short_run and count > 4: break

        if failed_jobs:
            colortext.error('Some jobs failed to run:\n%s' % pprint.pformat(failed_jobs))
        if not quiet: print('')
        return True


    @job_creator
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set):
        raise Exception('not implemented yet')
        #assert(existing_prediction_set exists and has records)
        #assert(new_prediction_set is empty)
        #for each prediction record, add the record and all associated predictionfile records,


    def _add_job(self, tsession, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, user_dataset_experiment_id = None, keep_all_lines = False, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = {}, suppress_warnings = False):
        '''This is the internal function which adds a prediction job to the database. We distinguish it from add_job as
           prediction jobs added using that function should have no associated user dataset experiment ID.

           pdb_residues_to_rosetta_cache can be used to speed up job insertion. When the same PDB/chains combination is used again, this cache uses the old mapping rather than run RosettaScripts again.

           The extra_rosetta_command_flags variable is used to add additional flags e.g. "-ignore_zero_occupancy false".
           These are used to generate a mapping from PDB to Rosetta numbering so they should be set according to how they
           are set in the protocol. In particular, include any flags which have an effect on what residues are present.
           '-ignore_zero_occupancy false' and '-ignore_unrecognized_res' are typically used.
           '''

        # todo: do something with input_files when we use that here - see add_prediction_run
        assert(not(input_files))

        # Preliminaries
        if not self.rosetta_scripts_path or not os.path.exists(self.rosetta_scripts_path):
            raise Exception('The path "%s" to the RosettaScripts executable does not exist.' % self.rosetta_scripts_path)

        cache_maps = False
        if isinstance(pdb_residues_to_rosetta_cache, dict):
            cache_maps = True

        # Information for debugging
        pp_complex = None
        for r in tsession.execute('''SELECT * FROM PPComplex WHERE ID=:pp_complex_id''', dict(pp_complex_id = pp_complex_id)):
            assert(pp_complex == None)
            pp_complex = r

        # Determine the list of PDB chains that will be kept
        pdb_chains = self.get_chains_for_mutatagenesis(pp_mutagenesis_id, pdb_file_id, pp_complex_pdb_set_number, complex_id = pp_complex_id, tsession = tsession)

        pdb_chains_to_keep = set(pdb_chains['L'] + pdb_chains['R'])
        if self.rosetta_database_path:
            cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, self.rosetta_database_path, extra_rosetta_command_flags)
        else:
            cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, extra_rosetta_command_flags)

        if cache_maps and pdb_residues_to_rosetta_cache.get(cache_key):
            stripped_p = pdb_residues_to_rosetta_cache[cache_key]['stripped_p']
        else:
            # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
            p = PDB(tsession.query(dbmodel.PDBFile).filter(dbmodel.PDBFile.ID == pdb_file_id).one().Content)
            stripped_p = p
            if not keep_all_lines:
                p.strip_to_chains(list(pdb_chains_to_keep))
                if not keep_hetatm_lines:
                    p.strip_HETATMs()
                stripped_p = PDB('\n'.join(p.lines))

        # Determine PDB chains to move
        pdb_chains_to_move_str = ','.join(sorted(set(pdb_chains['R'])))

        # Check for CSE and MSE
        try:
            if 'CSE' in stripped_p.residue_types:
                raise Exception('This case contains a CSE residue which may (or may not) cause an issue.')
            elif 'MSE' in stripped_p.residue_types:
                raise Exception('This case contains an MSE residue which may (or may not) cause an issue.')
                # It looks like MSE (and CSE?) may now be handled - https://www.rosettacommons.org/content/pdb-files-rosetta-format
        except Exception, e:
            if not suppress_warnings:
                colortext.error('%s: %s, chains %s' % (str(e), stripped_p.pdb_id or pdb_file_id, str(pdb_chains_to_keep)))

        # Assert that there are no empty sequences
        assert(sorted(stripped_p.atom_sequences.keys()) == sorted(pdb_chains_to_keep))
        for chain_id, sequence in stripped_p.atom_sequences.iteritems():
            assert(len(sequence) > 0)

        # Get the PDB mutations and check that they make sense in the context of the stripped PDB file
        # Note: the schema assumes that at most one set of mutations can be specified per PDB file per complex per mutagenesis. We may want to relax that in future by adding the SetNumber to the PPMutagenesisPDBMutation table
        complex_mutations = [m for m in tsession.execute('SELECT * FROM PPMutagenesisMutation WHERE PPMutagenesisID=:pp_mutagenesis_id', dict(pp_mutagenesis_id = pp_mutagenesis_id))]
        pdb_complex_mutations = [m for m in tsession.execute('SELECT * FROM PPMutagenesisPDBMutation WHERE PPMutagenesisID=:pp_mutagenesis_id AND PPComplexID=:pp_complex_id AND PDBFileID=:pdb_file_id', dict(pp_mutagenesis_id = pp_mutagenesis_id, pp_complex_id = pp_complex_id, pdb_file_id = pdb_file_id))]

        assert(len(complex_mutations) == len(pdb_complex_mutations))
        mutations = [ChainMutation(m['WildTypeAA'], m['ResidueID'], m['MutantAA'], Chain = m['Chain']) for m in pdb_complex_mutations]
        try:
            stripped_p.validate_mutations(mutations)
        except Exception, e:
            colortext.error('%s: %s' % (str(e), str(mutations)))
            #colortext.warning('PPMutagenesisID=%d, ComplexID=%d, PDBFileID=%s, SetNumber=%d, UserDatasetExperimentID=%d' % (pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, user_dataset_experiment_id))
            #colortext.warning('SKEMPI record: %s' % self.DDG_db.execute_select('SELECT * FROM PPMutagenesis WHERE ID=%s', parameters=(pp_mutagenesis_id,))[0]['SKEMPI_KEY'])
            #colortext.warning('PDB chains to keep: %s' % str(pdb_chains_to_keep))
            #colortext.warning('PPIPDBPartnerChain records: %s' % pprint.pformat(self.DDG_db.execute_select('SELECT PPIPDBPartnerChain.* FROM PPIPDBPartnerChain INNER JOIN PPIPDBSet ON PPIPDBSet.PPComplexID=PPIPDBPartnerChain.PPComplexID AND PPIPDBSet.SetNumber=PPIPDBPartnerChain.SetNumber WHERE PPIPDBPartnerChain.PPComplexID=%s AND IsComplex=1 ORDER BY PPIPDBPartnerChain.SetNumber, PPIPDBPartnerChain.ChainIndex', parameters=(pp_complex_id,))))

        # Determine the mapping from the stripped PDB to Rosetta numbering
        # Note: we assume that this stripped PDB will be the input to the Rosetta protocol and that

        # Make JSON mappings
        if cache_maps and pdb_residues_to_rosetta_cache.get(cache_key):
            atom_to_rosetta_residue_map = pdb_residues_to_rosetta_cache[cache_key]['atom_to_rosetta_residue_map']
            rosetta_to_atom_residue_map = pdb_residues_to_rosetta_cache[cache_key]['rosetta_to_atom_residue_map']
        else:
            if self.rosetta_database_path:
                stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, self.rosetta_database_path, extra_command_flags = extra_rosetta_command_flags)
            else:
                stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, extra_command_flags = extra_rosetta_command_flags)
            atom_to_rosetta_residue_map = stripped_p.get_atom_sequence_to_rosetta_json_map()
            rosetta_to_atom_residue_map = stripped_p.get_rosetta_sequence_to_atom_json_map()

        if cache_maps and (not pdb_residues_to_rosetta_cache.get(cache_key)):
            pdb_residues_to_rosetta_cache[cache_key] = dict(
                stripped_p = stripped_p,
                atom_to_rosetta_residue_map = atom_to_rosetta_residue_map,
                rosetta_to_atom_residue_map = rosetta_to_atom_residue_map)

        # Assert that there are no empty sequences in the Rosetta-processed PDB file
        total_num_residues = 0
        d = json.loads(rosetta_to_atom_residue_map)
        stripped_p_chains = stripped_p.atom_sequences.keys()
        for chain_id in stripped_p_chains:
            num_chain_residues = len([z for z in d.values() if z[0] == chain_id])
            total_num_residues += num_chain_residues
            assert(num_chain_residues > 0)

        pdb_filename = '%s_%s.pdb' % (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)))

        # Create parameter substitution dictionary
        mutfile_name = 'mutations.mutfile'
        resfile_name = 'mutations.resfile'
        parameter_sub_dict = {
            '%%input_pdb%%' : pdb_filename,
            '%%chainstomove%%' : pdb_chains_to_move_str,
            '%%pathtoresfile%%' : resfile_name,
            '%%pathtomutfile%%' : mutfile_name,
        }

        if test_only:
            return

        # All functions below use tsession which allows us to use transactions which can be rolled back if errors occur
        if protocol_id:
            existing_records = [r for r in tsession.execute('SELECT * FROM {0} WHERE PredictionSet=:prediction_set AND UserPPDataSetExperimentID=:user_dataset_experiment_id AND ProtocolID=:protocol_id'.format(self._get_prediction_table()), dict(prediction_set = prediction_set_id, user_dataset_experiment_id = user_dataset_experiment_id, protocol_id = protocol_id))]
        else:
            existing_records = [r for r in tsession.execute('SELECT * FROM {0} WHERE PredictionSet=:prediction_set AND UserPPDataSetExperimentID=:user_dataset_experiment_id AND ProtocolID IS NULL'.format(self._get_prediction_table()), dict(prediction_set = prediction_set_id, user_dataset_experiment_id = user_dataset_experiment_id))]
        assert(len(existing_records) == 0)

        prediction_record = dict(
            PredictionSet = prediction_set_id,
            PPMutagenesisID = pp_mutagenesis_id,
            UserPPDataSetExperimentID = user_dataset_experiment_id,
            ProtocolID = protocol_id,
            JSONParameters = json_dumps(parameter_sub_dict),
            DevelopmentProtocolID = None,
            ExtraParameters = extra_rosetta_command_flags,
            Status = 'queued',
            Cost = total_num_residues,
            KeptHETATMLines = keep_hetatm_lines,
        )
        prediction_ppi = get_or_create_in_transaction(tsession, self._get_sqa_prediction_table(), dict(
            PredictionSet = prediction_set_id,
            PPMutagenesisID = pp_mutagenesis_id,
            UserPPDataSetExperimentID = user_dataset_experiment_id,
            ProtocolID = protocol_id,
            JSONParameters = json_dumps(parameter_sub_dict),
            DevelopmentProtocolID = None,
            ExtraParameters = extra_rosetta_command_flags,
            Status = 'queued',
            Cost = total_num_residues,
            KeptHETATMLines = keep_hetatm_lines,
        ), missing_columns = ['ID', 'EntryDate', 'StartDate', 'EndDate', 'Errors', 'AdminCommand', 'maxvmem', 'DDGTime', 'NumberOfMeasurements'])
        #sql, params, record_exists = self.DDG_db.create_insert_dict_string(self._get_prediction_table(), prediction_record, ['PredictionSet', 'UserPPDataSetExperimentID', 'ProtocolID'])
        #cur.execute(sql, params)
        #prediction_id = cur.lastrowid
        prediction_id = prediction_ppi.ID

        # Add the stripped PDB file
        self._add_prediction_file(tsession, prediction_id, '\n'.join(stripped_p.lines), pdb_filename, 'PDB', 'StrippedPDB', 'Input', rm_trailing_line_whitespace = True, forced_mime_type = 'chemical/x-pdb')

        # Make and add the mutfile
        rosetta_mutations = stripped_p.map_pdb_residues_to_rosetta_residues(mutations)
        self._add_mutfile_to_prediction(tsession, prediction_id, rosetta_mutations, mutfile_name)

        # Make and add the resfile
        self._add_resfile_to_prediction(tsession, prediction_id, mutations, resfile_name)

        # Add the residue mappings
        self._add_residue_map_json_to_prediction(tsession, prediction_id, rosetta_to_atom_residue_map, 'Rosetta residue->PDB residue map')
        self._add_residue_map_json_to_prediction(tsession, prediction_id, atom_to_rosetta_residue_map, 'PDB residue->Rosetta residue map')

        # Add the params files
        self._add_ligand_params_files_to_prediction(tsession, prediction_id, pdb_file_id)

        if protocol_id:
            existing_records = [r for r in tsession.execute('SELECT * FROM {0} WHERE PredictionSet=:prediction_set_id AND UserPPDataSetExperimentID=:user_dataset_experiment_id AND ProtocolID=:protocol_id'.format(self._get_prediction_table()),
                             dict(prediction_set_id = prediction_set_id, user_dataset_experiment_id = user_dataset_experiment_id, protocol_id = protocol_id))]
        else:
            existing_records = [r for r in tsession.execute('SELECT * FROM {0} WHERE PredictionSet=:prediction_set_id AND UserPPDataSetExperimentID=:user_dataset_experiment_id AND ProtocolID IS NULL'.format(self._get_prediction_table()),
                             dict(prediction_set_id = prediction_set_id, user_dataset_experiment_id = user_dataset_experiment_id))]

        assert(len(existing_records) == 1)
        prediction_id = existing_records[0]['ID']
        return prediction_id


    #== Job execution/completion API ===========================================================
    #
    # This part of the API is responsible for starting jobs and setting them as failed or
    # completed


    @job_execution
    def set_job_temporary_protocol_field(self, prediction_id, prediction_set_id, temporary_protocol_field):
        raise Exception('not implemented yet')


    @job_execution
    def start_job(self, prediction_id, prediction_set_id):
        '''Sets the job status to "active". prediction_set must be passed and is used as a sanity check.'''
        prediction_record = self.DDG_db.execute_select('SELECT * FROM PredictionPPI WHERE ID=%s AND PredictionSet=%s', parameters=(prediction_id, prediction_set_id))
        if prediction_record['Protocol'] == None:
            print('empty Protocol')
            if prediction_record['DevelopmentProtocolID'] == None:
                raise Exception('Neither the Protocol nor the DevelopmentProtocolID is set for this job - it cannot be started without this information.')

        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_execution
    def get_max_number_of_cluster_jobs(self, prediction_set_id, priority):
        return self.DDG_db.execute_select('SELECT Value FROM _DBCONSTANTS WHERE VariableName="MaxStabilityClusterJobs"')['Value']


    @job_completion
    def complete_job(self, prediction_id, prediction_set, scores, maxvmem, ddgtime, files = []):
        '''Sets a job to 'completed' and stores scores. prediction_set must be passed and is used as a sanity check.'''

        raise Exception('This function needs to be implemented by subclasses of the API.')


    ###########################################################################################
    ## Analysis layer
    ##
    ## This part of the API is responsible for running analysis on completed predictions
    ###########################################################################################

    @analysis_api
    def determine_best_pairs(self, prediction_id, score_method_id = None, expectn = None, top_x = 3):
        '''This returns the top_x lowest-scoring wildtype/mutants for a prediction given a scoring method.
           The results are returned as a dict:
               "wildtype" -> list(tuple(score, structure_id))
               "mutant" -> list(tuple(score, structure_id))
           If no scoring method is supplied then the first (i.e. random) top_x structures are returned (with scores set
           to zero) as we have no method of scoring or discerning them.
        .'''

        scores = self.get_prediction_scores(prediction_id, expectn = expectn)
        if score_method_id != None:
            assert(isinstance(top_x, int) and top_x > 0)
            scores = scores.get(score_method_id)
            mutant_complexes = []
            wildtype_complexes = []
            for structure_id, scores in scores.iteritems():
                if scores.get('MutantComplex'):
                    mutant_complexes.append((scores['MutantComplex']['total'], structure_id))
                if scores.get('WildTypeComplex'):
                    wildtype_complexes.append((scores['WildTypeComplex']['total'], structure_id))
            wildtype_complexes = sorted(wildtype_complexes)[:top_x]
            mutant_complexes = sorted(mutant_complexes)[:top_x]
        else:
            wt_structure_ids = set()
            mut_structure_ids = set()
            for score_method_id, scores in scores.iteritems():
                for structure_id, scores in scores.iteritems():
                    if scores.get('WildTypeComplex'):
                        wt_structure_ids.add(structure_id)
                    if scores.get('MutantComplex'):
                        mut_structure_ids.add(structure_id)
            wildtype_complexes = [(0, i) for i in sorted(wt_structure_ids)]
            mutant_complexes = [(0, i) for i in sorted(mut_structure_ids)]
            if top_x != None:
                # If no score method is specified then we cannot choose the top X so we arbitrarily choose X structures
                assert(isinstance(top_x, int) and top_x > 0)
                wildtype_complexes = wildtype_complexes[:top_x]
                mutant_complexes = mutant_complexes[:top_x]

        # Truncate so that we have an equal number of both types
        max_len = min(len(wildtype_complexes), len(mutant_complexes))
        wildtype_complexes, mutant_complexes = wildtype_complexes[:max_len], mutant_complexes[:max_len]

        if wildtype_complexes and mutant_complexes:
            return {'wildtype' : wildtype_complexes, 'mutant' : mutant_complexes}
        return {}


    @app_pymol
    def create_pymol_session_in_memory(self, prediction_id, wt_task_number, mutant_task_number, pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol'):

        # Retrieve and unzip results
        archive = self.get_job_data(prediction_id)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            # Get the name of the files from the zip
            wildtype_filename = 'repacked_wt_round_%d.pdb.gz' % wt_task_number
            mutant_filename = None
            for filepath in sorted(zipped_content.namelist()):
                filename = os.path.split(filepath)[1]
                if filename.startswith('mut_') and filename.endswith('_round_%d.pdb.gz' % mutant_task_number):
                    mutant_filename = filename
                    break
            print(wildtype_filename, mutant_filename)
            PyMOL_session = None
            file_list = zipped_content.namelist()
            print(file_list)

            # If both files exist in the zip, extract their contents in memory and create a PyMOL session pair (PSE, script)
            if (mutant_filename in file_list) and (wildtype_filename in file_list):
                wildtype_pdb = zipped_content.open(wildtype_filename, 'r').read()
                mutant_pdb = zipped_content.open(mutant_filename, 'U').read()

                wildtype_pdb = read_file(write_temp_file('/tmp', wildtype_pdb, ftype = 'w', suffix = '.gz', prefix = ''))
                mutant_pdb = read_file(write_temp_file('/tmp', mutant_pdb, ftype = 'w', suffix = '.gz', prefix = ''))

                # todo: this should be structure_1_name = 'Wildtype', structure_2_name = 'Mutant' but the underlying PyMOL script needs to be parameterized
                chain_mapper = ScaffoldModelChainMapper.from_file_contents(wildtype_pdb, mutant_pdb, structure_1_name = 'Scaffold', structure_2_name = 'Model')
                PyMOL_session = chain_mapper.generate_pymol_session(pymol_executable = pymol_executable)

            zipped_content.close()
            return PyMOL_session

        except Exception, e:
            zipped_content.close()
            raise Exception(str(e))


    @app_pymol
    def create_full_pymol_session_in_memory(self, prediction_id, score_method_id = None, top_x = 3, mutation_string = None, settings = {}, pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol', wt_chain_seed = None, mutant_chain_seed = None):

        wt_chain_seed = wt_chain_seed or 'blue'
        mutant_chain_seed = mutant_chain_seed or 'yellow'

        best_pairs = self.determine_best_pairs(prediction_id, score_method_id = score_method_id, expectn = None, top_x = top_x)

        # Retrieve and unzip results
        archive = self.get_job_data(prediction_id)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            file_paths = {'wildtype' : {}, 'mutant' : {}}

            # Get the name of the files from the zip
            zip_filenames = set([os.path.split(filepath)[1] for filepath in zipped_content.namelist()])

            # Retrieve the input structure
            input_pdb_contents = None
            try:
                file_content_id = self.get_session().query(dbmodel.PredictionPPIFile).filter(and_(dbmodel.PredictionPPIFile.PredictionPPIID == prediction_id, dbmodel.PredictionPPIFile.FileRole == 'StrippedPDB')).one().FileContentID
                input_pdb_contents = self.importer.get_file_content_from_cache(file_content_id)
            except Exception, e:
                # Report the error but continue
                colortext.error(str(e))
                colortext.error(traceback.format_exc())

            # Find all wildtype structures
            for p in best_pairs['wildtype']:
                structure_id = p[1]
                expected_filename = 'repacked_wt_round_{0}.pdb.gz'.format(structure_id)
                if expected_filename in zip_filenames:
                    file_paths['wildtype'][structure_id] = expected_filename

            # Find all mutant structures
            mutant_ids = [p[1] for p in best_pairs['mutant']]
            for filename in zip_filenames:
                if filename.startswith('mut_'):
                    mtch = re.match('^mut_(.*?)_round_(\d+).pdb.*$', filename)
                    if mtch:
                        structure_id = int(mtch.group(2))
                        if structure_id in mutant_ids:
                            if not mutation_string:
                                mutation_string = mtch.group(1)
                            file_paths['mutant'][structure_id] = filename

            PyMOL_session = None
            file_list = zipped_content.namelist()

            # If both files exist in the zip, extract their contents in memory and create a PyMOL session pair (PSE, script)
            chain_mapper = DecoyChainMapper()

            for stypep in [('wildtype', 'wt', wt_chain_seed, 'white'), ('mutant', mutation_string or 'mutant', mutant_chain_seed, 'red')]:
                stype = stypep[0]
                prefix = stypep[1].replace(' ', '_')
                for structure_id, filename in file_paths[stype].iteritems():
                    if filename in file_list:
                        if filename.endswith('.gz'):
                            wildtype_pdb_stream = StringIO.StringIO(zipped_content.open(filename, 'r').read())
                            wildtype_pdb = gzip.GzipFile(fileobj=wildtype_pdb_stream).read()
                        else:
                            wildtype_pdb = zipped_content.open(filename, 'r').read()
                        pdb_object = PDB(wildtype_pdb)
                        chain_mapper.add(pdb_object, '{0}_n{1}'.format(prefix, structure_id), chain_seed_color = stypep[2], backbone_color = stypep[2], sidechain_color = stypep[3])
            if input_pdb_contents:
                chain_mapper.add(PDB(input_pdb_contents), 'input', backbone_color = 'grey50', sidechain_color = 'grey50')
            zipped_content.close()

            PyMOL_session = chain_mapper.generate_pymol_session(settings = settings, pymol_executable = pymol_executable)
            return PyMOL_session

        except Exception, e:
            zipped_content.close()
            raise Exception('{0}\n{1}'.format(str(e), traceback.format_exc()))


    def _get_prediction_data(self, prediction_id, score_method_id, main_ddg_analysis_type, expectn = None, extract_data_for_case_if_missing = False, root_directory = None, dataframe_type = "Binding affinity", prediction_data = {}):
        assert( main_ddg_analysis_type.startswith('DDG_') )
        analysis_type = main_ddg_analysis_type[4:]
        top_x = 3
        if analysis_type.startswith('Top'):
            analysis_function = self.get_top_x_ddg
            analysis_parameter = int( analysis_type[3:] )
            top_x = analysis_parameter
        elif analysis_type.startswith('Random'):
            analysis_function = self.get_random_pairing_ddg
            if len(analysis_type) > len('Random'):
                analysis_parameter = int( analysis_type[len('Random'):] )
            else:
                analysis_parameter = None
        elif analysis_type == 'AvgAllPairs':
            analysis_function = self.get_avg_all_pairings_ddg
            analysis_parameter = None
        elif analysis_type == 'MatchPairs':
            analysis_function = self.get_match_pairs_ddg
            analysis_parameter = None
        elif analysis_type.startswith('CplxBoltzWT'):
            assert( len(analysis_type) > len('CplxBoltzWT') )
            analysis_function = self.get_wt_complex_weighted_boltzmann_ddg
            analysis_parameter = float( analysis_type[len('CplxBoltzWT'):] )
        elif analysis_type.startswith('CplxBoltzMut'):
            assert( len(analysis_type) > len('CplxBoltzMut') )
            analysis_function = self.get_mut_complex_weighted_boltzmann_ddg
            analysis_parameter = float( analysis_type[len('CplxBoltzMut'):] )
        elif analysis_type.startswith('CplxBoltzBoth'):
            assert( len(analysis_type) > len('CplxBoltzBoth') )
            analysis_function = self.get_both_complex_weighted_boltzmann_ddg
            analysis_parameter = float( analysis_type[len('CplxBoltzBoth'):] )
        else:
            raise Exception("Didn't recognize analysis type: " + str(main_ddg_analysis_type))


        try:
            predicted_ddg = analysis_function(prediction_id, score_method_id, analysis_parameter, expectn = expectn)
        except Exception, e:
            colortext.pcyan(str(e))
            colortext.warning(traceback.format_exc())
            if extract_data_for_case_if_missing:
                self.extract_data_for_case(prediction_id, root_directory = root_directory, force = True, score_method_id = score_method_id)
            try:
                predicted_ddg = analysis_function(prediction_id, score_method_id, analysis_parameter, expectn = expectn)
            except PartialDataException, e:
                raise
            except Exception, e:
                raise
        top_x_ddg_stability = self.get_top_x_ddg_stability(prediction_id, score_method_id, top_x = top_x, expectn = expectn)

        prediction_data[main_ddg_analysis_type] = predicted_ddg
        prediction_data['DDGStability_Top%d' % top_x] = top_x_ddg_stability
        return prediction_data

    @analysis_api
    def get_wt_complex_weighted_boltzmann_ddg(self, prediction_id, score_method_id, temperature, expectn = None):
        return self.get_complex_weighted_boltzmann_ddg(prediction_id, score_method_id, temperature, expectn = expectn, scores_to_weight = 'wt_complex')

    @analysis_api
    def get_mut_complex_weighted_boltzmann_ddg(self, prediction_id, score_method_id, temperature, expectn = None):
        return self.get_complex_weighted_boltzmann_ddg(prediction_id, score_method_id, temperature, expectn = expectn, scores_to_weight = 'mut_complex')

    @analysis_api
    def get_both_complex_weighted_boltzmann_ddg(self, prediction_id, score_method_id, temperature, expectn = None):
        return self.get_complex_weighted_boltzmann_ddg(prediction_id, score_method_id, temperature, expectn = expectn, scores_to_weight = 'both_complexes')

    @analysis_api
    def get_complex_weighted_boltzmann_ddg(self, prediction_id, score_method_id, temperature, expectn = None, scores_to_weight = 'wt_complex'):
        '''
        Returns DDG for this prediction by averaging all values for paired output structures
        '''

        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
            return None

        if self.scores_contains_ddg_score(scores):
            raise Exception("This scoring analysis doesn't make sense to use without complex scores")

        def boltz_exponent(x, t):
            return numpy.exp( -1.0 * x / t )

        try:
            np_type = numpy.float64

            struct_nums = scores.keys()
            mut_complex = numpy.array( [np_type( scores[struct_num]['MutantComplex']['total'] ) for struct_num in struct_nums]  )
            mut_lpartner = numpy.array( [np_type( scores[struct_num]['MutantLPartner']['total'] ) for struct_num in struct_nums] )
            mut_rpartner = numpy.array( [np_type( scores[struct_num]['MutantRPartner']['total'] ) for struct_num in struct_nums] )
            wt_complex = numpy.array( [np_type( scores[struct_num]['WildTypeComplex']['total'] ) for struct_num in struct_nums] )
            wt_lpartner = numpy.array( [np_type( scores[struct_num]['WildTypeLPartner']['total'] ) for struct_num in struct_nums] )
            wt_rpartner = numpy.array( [np_type( scores[struct_num]['WildTypeRPartner']['total'] ) for struct_num in struct_nums] )

            matched_ddgs = (mut_complex - mut_lpartner - mut_rpartner) - (wt_complex - wt_lpartner - wt_rpartner)

            if scores_to_weight == 'wt_complex':
                scores_for_weighting = wt_complex
            elif scores_to_weight == 'mut_complex':
                scores_for_weighting = mut_complex
            elif scores_to_weight == 'both_complexes':
                scores_for_weighting = mut_complex + wt_complex
            else:
                raise Exception('Unrecognized scores_to_weight argument: ' + str(scores_to_weight) )

            max_scores_for_weighting = numpy.max(scores_for_weighting)
            normalized_scores_for_weighting = scores_for_weighting - max_scores_for_weighting

            exponented_scores = numpy.exp( np_type(-1.0) * normalized_scores_for_weighting / np_type(temperature) )

            weighted_ddg = numpy.divide(
                numpy.sum( numpy.multiply(matched_ddgs, exponented_scores) ),
                numpy.sum( exponented_scores )
            )

            return weighted_ddg
        except PartialDataException:
            sys.exit(0)
            raise PartialDataException('The case is missing some data.')


    @analysis_api
    def get_match_pairs_ddg(self, prediction_id, score_method_id, structs_to_use, expectn = None):
        '''
        Returns DDG for this prediction by averaging all values for paired output structures
        '''

        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
            return None

        if self.scores_contains_ddg_score(scores):
            raise Exception("This scoring analysis doesn't make sense to use without complex scores")

        try:
            structs_to_use_score = numpy.average([
                (scores[struct_num]['MutantComplex']['total'] - scores[struct_num]['MutantLPartner']['total'] - scores[struct_num]['MutantRPartner']['total']) -
                (scores[struct_num]['WildTypeComplex']['total'] - scores[struct_num]['WildTypeLPartner']['total'] - scores[struct_num]['WildTypeRPartner']['total'])
                for struct_num in scores
            ])
            return structs_to_use_score
        except PartialDataException:
            sys.exit(0)
            raise PartialDataException('The case is missing some data.')


    @analysis_api
    def get_avg_all_pairings_ddg(self, prediction_id, score_method_id, structs_to_use, expectn = None):
        '''
        Returns DDG for this prediction by averaging together all possible pairings
        '''

        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
            return None

        if self.scores_contains_ddg_score(scores):
            raise Exception("This scoring analysis doesn't make sense to use without complex scores")

        try:
            all_struct_num_pairs = []
            for wt_struct_num in scores:
                if 'WildTypeComplex' in scores[wt_struct_num]:
                    for mut_struct_num in scores:
                        if 'MutantComplex' in scores[mut_struct_num]:
                            all_struct_num_pairs.append( (wt_struct_num, mut_struct_num) )

            structs_to_use_score = numpy.average([
                (scores[mut_struct_num]['MutantComplex']['total'] - scores[mut_struct_num]['MutantLPartner']['total'] - scores[mut_struct_num]['MutantRPartner']['total']) -
                (scores[wt_struct_num]['WildTypeComplex']['total'] - scores[wt_struct_num]['WildTypeLPartner']['total'] - scores[wt_struct_num]['WildTypeRPartner']['total'])
                for wt_struct_num, mut_struct_num in all_struct_num_pairs
            ])
            return structs_to_use_score
        except PartialDataException:
            sys.exit(0)
            raise PartialDataException('The case is missing some data.')


    @analysis_api
    def get_random_pairing_ddg(self, prediction_id, score_method_id, structs_to_use, expectn = None):
        '''
        Returns DDG for this prediction by randomly pairing mutant structures with wildtype structures
        '''

        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
            return None

        if self.scores_contains_ddg_score(scores):
            try:
                total_scores = [scores[struct_num]['DDG']['total'] for struct_num in scores]
                if structs_to_use == None:
                    structs_to_use = len(total_scores)
                structs_to_use_score = numpy.average(
                    random.sample(total_scores, structs_to_use)
                )
                return structs_to_use_score
            except:
                raise PartialDataException('The case is missing some data.')

        try:
            if structs_to_use == None:
                structs_to_use = len(scores)
            else:
                structs_to_use = min(structs_to_use, len(scores))

            structs_to_use_wt_struct_nums = random.sample(scores.keys(), structs_to_use)
            structs_to_use_mut_struct_nums = random.sample(scores.keys(), structs_to_use)

            structs_to_use_score = numpy.average([
                (scores[mut_struct_num]['MutantComplex']['total'] - scores[mut_struct_num]['MutantLPartner']['total'] - scores[mut_struct_num]['MutantRPartner']['total']) -
                (scores[wt_struct_num]['WildTypeComplex']['total'] - scores[wt_struct_num]['WildTypeLPartner']['total'] - scores[wt_struct_num]['WildTypeRPartner']['total'])
                for wt_struct_num, mut_struct_num in zip(structs_to_use_wt_struct_nums, structs_to_use_mut_struct_nums)
            ])
            return structs_to_use_score
        except PartialDataException:
            raise PartialDataException('The case is missing some data.')


    @analysis_api
    def get_top_x_ddg(self, prediction_id, score_method_id, top_x , expectn = None):
        '''Returns the TopX value for the prediction. Typically, this is the mean value of the top X predictions for a
           case computed using the associated Score records in the database.'''

        # scores is a mapping from nstruct -> ScoreType -> score record where ScoreType is one of 'DDG', 'WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex'
        # if we do the calculation in Python, pull scores out to the top level first
        # otherwise, we can add a stored procedure to determine the TopX
        # if we go the Python route, we can implement different variations on TopX (including a stored procedure) and pass the function pointers as an argument to the main analysis function

        # Make sure that we have as many cases as we expect
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)

        if scores == None:
            return None

        if self.scores_contains_ddg_score(scores):
            try:
                total_scores = [(scores[struct_num]['DDG']['total'], struct_num) for struct_num in scores]
                total_scores.sort()
                top_x_struct_nums = [t[1] for t in total_scores[:top_x]]
                top_x_score = numpy.average([
                    scores[struct_num]['DDG']['total']
                    for struct_num in top_x_struct_nums
                ])
                return top_x_score
            except:
                print scores[struct_num]
                raise PartialDataException('The case is missing some data.')

        try:
            wt_total_scores = [(scores[struct_num]['WildTypeComplex']['total'], struct_num) for struct_num in scores]
            wt_total_scores.sort()
            top_x_wt_struct_nums = [t[1] for t in wt_total_scores[:top_x]]

            mut_total_scores = [(scores[struct_num]['MutantComplex']['total'], struct_num) for struct_num in scores]
            mut_total_scores.sort()
            top_x_mut_struct_nums = [t[1] for t in mut_total_scores[:top_x]]

            top_x_score = numpy.average([
                (scores[mut_struct_num]['MutantComplex']['total'] - scores[mut_struct_num]['MutantLPartner']['total'] - scores[mut_struct_num]['MutantRPartner']['total']) -
                (scores[wt_struct_num]['WildTypeComplex']['total'] - scores[wt_struct_num]['WildTypeLPartner']['total'] - scores[wt_struct_num]['WildTypeRPartner']['total'])
                for wt_struct_num, mut_struct_num in zip(top_x_wt_struct_nums, top_x_mut_struct_nums)
            ])
            return top_x_score
        except:
            raise PartialDataException('The case is missing some data.')


    def scores_contains_ddg_score(self, scores):
        for struct_num, score_dict in scores.iteritems():
            if 'DDG' not in score_dict:
                return False
        return True


    def scores_contains_complex_scores(self, scores):
        for struct_num, score_dict in scores.iteritems():
            if 'WildTypeComplex' not in score_dict or 'MutantComplex' not in score_dict:
                return False
        return True


    @analysis_api
    def get_top_x_ddg_stability(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''Returns the TopX value for the prediction only considering the complex scores. This computation may work as a
           measure of a stability DDG value.'''
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None or not self.scores_contains_complex_scores(scores):
            return None

        wt_total_scores = [(scores[struct_num]['WildTypeComplex']['total'], struct_num) for struct_num in scores]
        wt_total_scores.sort()
        top_x_wt_struct_nums = [t[1] for t in wt_total_scores[:top_x]]

        mut_total_scores = [(scores[struct_num]['MutantComplex']['total'], struct_num) for struct_num in scores]
        mut_total_scores.sort()
        top_x_mut_struct_nums = [t[1] for t in mut_total_scores[:top_x]]

        return numpy.average([scores[mut_struct_num]['MutantComplex']['total'] - scores[wt_struct_num]['WildTypeComplex']['total']
                           for wt_struct_num, mut_struct_num in zip(top_x_wt_struct_nums, top_x_mut_struct_nums)])


    @analysis_api
    def get_analysis_dataframe(self, prediction_set_id,
            experimental_data_exists = True,
            prediction_set_series_name = None, prediction_set_description = None, prediction_set_credit = None,
            prediction_set_color = None, prediction_set_alpha = None,
            use_existing_benchmark_data = True,
            include_derived_mutations = False,
            use_single_reported_value = False,
            ddg_analysis_type = 'DDG_Top3',
            take_lowest = None,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            score_method_id = None,
            expectn = None,
            allow_failures = False,
            extract_data_for_case_if_missing = False,
            debug = False,
            restrict_to = set(),
            remove_cases = set(),
        ):

        #todo: rename function since we return BenchmarkRun objects

        assert(score_method_id)
        dataframe_type = 'Binding affinity'

        parameters = copy.copy(locals())
        del parameters['self']
        return super(BindingAffinityDDGInterface, self)._get_analysis_dataframe(BindingAffinityBenchmarkRun, **parameters)


    @analysis_api
    def get_existing_analysis(self, prediction_set_id = None, analysis_dataframe_id = None, return_dataframe = True):
        '''Returns a list of the summary statistics for any existing dataframes in the database.
           Each item in the list is a dict corresponding to a dataframe. These dicts are structured as e.g.

            {
                'AnalysisDataFrameID': 185L,
                'analysis_sets': ['SKEMPI', 'BeAtMuSiC', 'ZEMu'],
                'analysis_type': 'DDG_Top3',
                'analysis_type_description': '...',
                'dataframe': <pandas dataframe>,
                'scalar_adjustments': {
                    'BeAtMuSiC': 2.383437079488905,
                    'SKEMPI': 2.206268329703589,
                    'ZEMu': 2.2046199780552374
                },
                'stats': {
                    'BeAtMuSiC': {
                        'MAE': nan,
                        'fraction_correct': 0.7308900523560209,
                        'fraction_correct_fuzzy_linear': 0.74128683025321573,
                        'gamma_CC': 0.4047074501135616,
                        'ks_2samp': (0.24269480519480513, 2.9466866316296972e-32),
                        'kstestx': (nan, nan),
                        'kstesty': (nan, nan),
                        'normaltestx': (nan, nan),
                        'normaltesty': (nan, nan),
                        'pearsonr': (nan, 1.0),
                        'spearmanr': (0.41841534629950339, 2.1365219255798831e-53)
                    },
                    'SKEMPI': {...},
                    'ZEMu': {...},
                }
            }
        '''
        ### KAB TODO: this function is not adjusted for new changes in top_x
        if analysis_dataframe_id == None:
            # Get a valid PredictionSet record if one exists
            assert(prediction_set_id != None)
            try:
                prediction_set = self.get_session().query(dbmodel.PredictionSet).filter(and_(dbmodel.PredictionSet.ID == prediction_set_id, dbmodel.PredictionSet.BindingAffinity == 1)).one()
            except:
                return None
            dataframes = self.get_session().query(dbmodel.AnalysisDataFrame).filter(and_(dbmodel.AnalysisDataFrame.PredictionSet == prediction_set_id, dbmodel.AnalysisDataFrame.DataFrameType == 'Binding affinity')).order_by(dbmodel.AnalysisDataFrame.ScoreMethodID, dbmodel.AnalysisDataFrame.TopX, dbmodel.AnalysisDataFrame.StabilityClassicationExperimentalCutoff, dbmodel.AnalysisDataFrame.StabilityClassicationPredictedCutoff)
        else:
            try:
                dataframe = self.get_session().query(dbmodel.AnalysisDataFrame).filter(dbmodel.AnalysisDataFrame.ID == analysis_dataframe_id).one()
                assert(dataframe.DataFrameType == 'Binding affinity')
                dataframes = [dataframe]
            except Exception, e:
                colortext.error(str(e))
                colortext.error(traceback.format_exc())
                return None

        analysis_results = []
        dataframes = [dfr for dfr in dataframes]
        for dfr in dataframes:
            # The dict to return
            dfi = dfr.get_dataframe_info()
            dfi['stats'] = {}

            # Compute the stats per analysis set
            df = dfi['dataframe']
            if dfi['analysis_sets']:
                # Case where there are analysis sets
                for analysis_set in dfi['analysis_sets']:
                    dfi['stats'][analysis_set] = get_xy_dataset_statistics_pandas(
                        df,
                        BindingAffinityBenchmarkRun.get_analysis_set_fieldname('Experimental', analysis_set),
                        BindingAffinityBenchmarkRun.get_analysis_set_fieldname('Predicted_adj', analysis_set),
                        fcorrect_x_cutoff = float(dfr.StabilityClassicationExperimentalCutoff),
                        fcorrect_y_cutoff = float(dfr.StabilityClassicationPredictedCutoff),
                        ignore_null_values = True)
            elif 'Experimental' in df.columns:
                # Case where there are no analysis sets
                dfi['stats']['Global'] = get_xy_dataset_statistics_pandas(
                    df,
                    'Experimental',
                    'Predicted_adj',
                    fcorrect_x_cutoff = float(dfr.StabilityClassicationExperimentalCutoff),
                    fcorrect_y_cutoff = float(dfr.StabilityClassicationPredictedCutoff),
                    ignore_null_values = True)
            else:
                # Case where there are no experimental data
                dfi['stats'] = None

            if not return_dataframe:
                # May be useful if we are keeping a lot of these in memory and the dataframe is not useful
                dfi['dataframe'] = None
            analysis_results.append(dfi)

        return analysis_results


    @analysis_api
    def analyze(self, prediction_set_ids, score_method_ids,
            experimental_data_exists = True,
            analysis_set_ids = [],
            prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_credits = {}, prediction_set_colors = {}, prediction_set_alphas = {},
            use_published_data = False,
            allow_failures = False,
            use_existing_benchmark_data = True, recreate_graphs = False,
            include_derived_mutations = False,
            expectn = 50,
            use_single_reported_value = False,
            take_lowests = [],
            ddg_analysis_types = [],
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            output_directory = None,
            output_directory_root = None,
            generate_plots = True,
            generate_matplotlib_plots = False,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            debug = False,
            restrict_to = set(),
            remove_cases = set(),
            call_analysis = True,
            ):
        '''Runs the analyses for the specified PredictionSets and cross-analyzes the sets against each other if appropriate.

           * Analysis setup arguments *

           PredictionSets is a list of PredictionSet IDs. Each PredictionSet will be analyzed separately and appropriate
           pairs will be cross-analyzed.
           PredictionSetSeriesNames, PredictionSetDescriptions, and PredictionSetCredits are mappings from PredictionSet IDs
           to series names (in plots), descriptions, and credits respectively. The details are stored in PredictionSet so
           they are not necessary. The mappings can be used to override the database values to customize the analysis
           reports. Likewise, PredictionSetColors and PredictionSetAlphas are mappings to series colors and transparency values
           for use in the plots.
           use_published_data. todo: implement later. This should include any published data e.g. the Kellogg et al. data for protein stability.
           use_existing_benchmark_data and recreate_graphs are data creation arguments i.e. "should we use existing data or create it from scratch?"
           include_derived_mutations is used to filter out dataset cases with derived mutations.
           expectn declares how many predictions we expect to see per dataset case. If the actual number is less than expectn
           then a warning will be included in the analysis.

           * Dataframe arguments *

           use_single_reported_value is specific to ddg_monomer. If this is True then the DDG value reported by the application is used and take_lowest is ignored. This is inadvisable - take_lowest = 3 is a better default.
           take_lowest AKA Top_X. Specifies how many of the best-scoring groups of structures to consider when calculating the predicted DDG value.
           analysis_types defines if other analysis methods other than TopX/take_lowest will be used. Not mutually exclusive.
           burial_cutoff defines what should be considered buried (DSSPExposure field). Values around 1.0 are fully exposed, values of 0.0 are fully buried. For technical reasons, the DSSP value can exceed 1.0 but usually not by much.
           stability_classication_experimental_cutoff AKA x_cutoff. This defines the neutral mutation range for experimental values in kcal/mol i.e. values between -1.0 and 1.0 kcal/mol are considered neutral by default.
           stability_classication_predicted_cutoff AKA y_cutoff. This defines the neutral mutation range for predicted values in energy units.

           * Reporting arguments *

           output_directory : The directory in which to save plots and reports.
           output_directory_root : A place to create an autogenerated output directory.
           generate_plots   : if plots are not needed, setting this to False can shorten the analysis time.
           report_analysis  : Whether or not to print analysis to stdout.
           silent = False   : Whether or not anything should be printed to stdout (True is useful for webserver interaction).
        '''
        for ddg_analysis_type in ddg_analysis_types:
            assert( ddg_analysis_type.startswith('DDG_') )
        for take_lowest in take_lowests:
            assert(take_lowest > 0 and (int(take_lowest) == take_lowest))
            ddg_analysis_types.append( 'DDG_Top%d' % take_lowest )

        # Remove duplicate analysis types
        ddg_analysis_types = set( ddg_analysis_types )
        ddg_analysis_types = sorted( list(ddg_analysis_types) )

        assert(0 <= burial_cutoff <= 2.0)
        assert(stability_classication_experimental_cutoff > 0)
        assert(stability_classication_predicted_cutoff > 0)
        assert(expectn > 0 and (int(expectn) == expectn))

        # Can't specify both output_directory and output_directory_root
        if output_directory_root != None:
            assert( output_directory == None )
            if not os.path.isdir( output_directory_root ):
                os.makedirs( output_directory_root )
        if output_directory != None:
            assert( output_directory_root == None )

        benchmark_runs = []

        for prediction_set_id in prediction_set_ids:
            if len(prediction_set_ids) > 1:
                print 'Generating benchmark run for prediction set: %s' % prediction_set_id
            for score_method_id in score_method_ids:
                if len(score_method_ids) > 1:
                    print 'Generating benchmark run for score method ID: %d' % score_method_id
                for ddg_analysis_type in ddg_analysis_types:
                    if len(ddg_analysis_types) > 1:
                        print 'Generating benchmark run for DDG analysis type: %s' % ddg_analysis_type

                    benchmark_run = self.get_analysis_dataframe(prediction_set_id,
                        experimental_data_exists = experimental_data_exists,
                        prediction_set_series_name = prediction_set_series_names.get(prediction_set_id),
                        prediction_set_description = prediction_set_descriptions.get(prediction_set_id),
                        prediction_set_color = prediction_set_colors.get(prediction_set_id),
                        prediction_set_alpha = prediction_set_alphas.get(prediction_set_id),
                        prediction_set_credit = prediction_set_credits[prediction_set_id],
                        use_existing_benchmark_data = use_existing_benchmark_data,
                        include_derived_mutations = include_derived_mutations,
                        use_single_reported_value = use_single_reported_value,
                        ddg_analysis_type = ddg_analysis_type,
                        burial_cutoff = burial_cutoff,
                        stability_classication_experimental_cutoff = 1.0,
                        stability_classication_predicted_cutoff = 1.0,
                        report_analysis = report_analysis,
                        silent = silent,
                        root_directory = root_directory, # where to find the
                        score_method_id = score_method_id,
                        expectn = expectn,
                        allow_failures = allow_failures,
                        debug = debug,
                        restrict_to = restrict_to,
                        remove_cases = remove_cases,
                    )

                    # The keys of scalar_adjustments are the stored analysis sets
                    analysis_sets_to_run = benchmark_run.scalar_adjustments.keys()
                    if analysis_set_ids:
                        analysis_sets_to_run = set(analysis_sets_to_run).intersection(set(analysis_set_ids))
                    benchmark_runs.append(benchmark_run)

        analysis_sets_to_run = sorted(analysis_sets_to_run)
        if experimental_data_exists:
            #todo: hack. this currently seems to expect all datapoints to be present. handle the case when we are missing data e.g. prediction set "ZEMu run 1"
            analysis_sets_to_run = ['ZEMu'] # ['BeAtMuSiC', 'SKEMPI', 'ZEMu']

        if call_analysis:
            if len(benchmark_runs) == 1 and len(analysis_sets_to_run) == 1:
                if output_directory_root:
                    # Create output directory inside output_directory_root
                    output_directory = os.path.join(output_directory_root, '%s-%s-%s_n-%d_topx-%d_score_method_%d-analysis_%s' % (time.strftime("%y%m%d"), getpass.getuser(), prediction_set_id, expectn, take_lowest, score_method_id, analysis_set_id))

                colortext.message(analysis_set_id)

                benchmark_run.full_analysis(analysis_set_id, output_directory)
            else:
                if output_directory or not output_directory_root:
                    raise Exception("Multiple benchmark run objects will be analyzed and output created; this requires setting output_directory_root instead of output_directory")

                BindingAffinityBenchmarkRun.analyze_multiple(
                    benchmark_runs,
                    analysis_sets = analysis_sets_to_run,
                    analysis_directory = output_directory_root,
                )
        else:
            return (benchmark_runs, analysis_sets_to_run)


    ################################################################################################
    ## Private API layer
    ## These are helper functions used internally by the class but which are not intended for export
    ################################################################################################


    ###########################################################################################
    ## Subclass layer
    ##
    ## These functions need to be implemented by subclasses
    ###########################################################################################

    # Concrete functions


    def _get_sqa_prediction_table(self): return dbmodel.PredictionPPI
    def _get_sqa_prediction_structure_scores_table(self): return dbmodel.PredictionPPIStructureScore

    def _get_sqa_user_dataset_experiment_table(self): return dbmodel.UserPPDataSetExperiment
    def _get_sqa_user_dataset_experiment_tag_table(self): return dbmodel.UserPPDataSetExperimentTag
    def _get_sqa_user_dataset_experiment_tag_table_udsid(self): return dbmodel.UserPPDataSetExperimentTag.UserPPDataSetExperimentID
    def _get_sqa_predictions_user_dataset_experiment_id(self, p): return p.UserPPDataSetExperimentID
    def _get_sqa_prediction_type(self): return dbmodel.PredictionSet.BindingAffinity

    prediction_table = 'PredictionPPI'
    def _get_prediction_table(self): return self.prediction_table
    prediction_structure_scores_table = 'PredictionPPIStructureScore'
    def _get_prediction_structure_scores_table(self): return self.prediction_structure_scores_table
    def _get_prediction_type(self): return 'BindingAffinity'
    def _get_prediction_dataset_type(self): return 'Binding affinity'
    def _get_prediction_type_description(self): return 'binding affinity'
    def _get_user_dataset_experiment_table(self): return 'UserPPDataSetExperiment'
    def _get_user_dataset_experiment_tag_table(self): return 'UserPPDataSetExperimentTag'
    def _get_allowed_score_types(self): return set(['DDG', 'WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex'])


    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================


    @informational_job
    def get_development_protocol(self, development_protocol_id):
        results = self.DDG_db.execute_select('SELECT * FROM DevelopmentProtocol WHERE ID = %s', parameters=(development_protocol_id,) )
        assert( len(results) == 1 )
        return results[0]


    @informational_pdb
    def get_complex_ids_matching_protein_name(self, partial_name, tsession = None):
        '''Returns a list of PPComplex IDs where at least one of the partner names matches partial_name.'''

        tsession = self.importer.get_session(utf = True)
        tsession_utf = self.importer.get_session()

        results = []
        partial_name_ascii = partial_name.encode('ascii', errors='ignore').decode('ascii') # ugh
        if len(partial_name.split()) == 1 and len(partial_name) <= 4:
            results += [c.ID for c in tsession_utf.query(dbmodel.PPComplex).filter(or_(
                            dbmodel.PPComplex.LName.like(u'^' + partial_name),
                            dbmodel.PPComplex.LShortName.like(u'^' + partial_name),
                            dbmodel.PPComplex.RName.like(u'^' + partial_name),
                            dbmodel.PPComplex.RShortName.like(u'^' + partial_name)))]
            results += [c.ID for c in tsession.query(dbmodel.PPComplex).filter(or_(
                            dbmodel.PPComplex.LHTMLName.like('^' + partial_name_ascii),
                            dbmodel.PPComplex.RHTMLName.like('^' + partial_name_ascii)))]
            results += [c.ID for c in tsession_utf.query(dbmodel.PPComplex).filter(or_(
                            dbmodel.PPComplex.LName.like(partial_name + u'$'),
                            dbmodel.PPComplex.LShortName.like(partial_name + u'$'),
                            dbmodel.PPComplex.RName.like(partial_name + u'$'),
                            dbmodel.PPComplex.RShortName.like(partial_name + u'$')))]
            results += [c.ID for c in tsession.query(dbmodel.PPComplex).filter(or_(
                            dbmodel.PPComplex.LHTMLName.like(partial_name_ascii + '$'),
                            dbmodel.PPComplex.RHTMLName.like(partial_name_ascii + '$')))]
        else:
            results += [c.ID for c in tsession_utf.query(dbmodel.PPComplex).filter(or_(
                            dbmodel.PPComplex.LName.like(u'%' + partial_name + u'%'),
                            dbmodel.PPComplex.LShortName.like(u'%' + partial_name + u'%'),
                            dbmodel.PPComplex.RName.like(u'%' + partial_name + u'%'),
                            dbmodel.PPComplex.RShortName.like(u'%' + partial_name + u'%')))]
            results += [c.ID for c in tsession.query(dbmodel.PPComplex).filter(or_(
                            dbmodel.PPComplex.LHTMLName.like('%' + partial_name_ascii + '%'),
                            dbmodel.PPComplex.RHTMLName.like('%' + partial_name_ascii + '%')))]
        return results

        qry = '''SELECT ID FROM PPComplex
                 WHERE
                 LName LIKE %s
                 OR LShortName LIKE %s
                 OR LHTMLName LIKE %s
                 OR RName LIKE %s
                 OR RShortName LIKE %s
                 OR RHTMLName LIKE %s ORDER BY ID'''


        if len(partial_name.split()) == 1 and len(partial_name) <= 4:
            # for short names, we require that any matches have the string as a prefix or suffix as otherwise we may get many matches e.g. 'RAN' matches 'transferase', 'membrane', etc.
            partial_name_ascii = partial_name.encode('ascii', errors='ignore').decode('ascii') # ugh
            results += self.DDG_db_utf.execute_select(qry, parameters=(u'%{0}'.format(partial_name), u'%{0}'.format(partial_name), '%{0}'.format(partial_name_ascii), u'%{0}'.format(partial_name), u'%{0}'.format(partial_name), '%{0}'.format(partial_name_ascii)))
            results += self.DDG_db_utf.execute_select(qry, parameters=(u'{0}%'.format(partial_name), u'{0}%'.format(partial_name), '{0}%'.format(partial_name_ascii), u'{0}%'.format(partial_name), u'{0}%'.format(partial_name), '{0}%'.format(partial_name_ascii)))
        else:
            partial_name_ascii = partial_name.encode('ascii', errors='ignore').decode('ascii') # ugh
            results += self.DDG_db_utf.execute_select(qry, parameters=(u'%{0}%'.format(partial_name), u'%{0}%'.format(partial_name), '%{0}%'.format(partial_name_ascii), u'%{0}%'.format(partial_name), u'%{0}%'.format(partial_name), '%{0}%'.format(partial_name_ascii)))
        return [r['ID'] for r in results]


    @informational_pdb
    def _get_pdb_chains_used_for_prediction_set(self, prediction_set):
        raise Exception('not implemented yet')
        return self.DDG_db.execute_select('''
            SELECT Prediction.ID, Experiment.PDBFileID, Chain
            FROM Prediction
            INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID
            INNER JOIN ExperimentChain ON ExperimentChain.ExperimentID=Prediction.ExperimentID
            WHERE PredictionSet=%s''', parameters=(prediction_set,))


    ###########################################################################################
    ## Prediction layer
    ##
    ## This part of the API is responsible for inserting prediction jobs in the database via
    ## the trickle-down proteomics paradigm.
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    def _charge_prediction_set_by_residue_count(self, PredictionSet):
        '''This function assigns a cost for a prediction equal to the number of residues in the chains.'''
        raise Exception('This function needs to be rewritten.')
        from klab.bio.rcsb import parseFASTAs

        DDG_db = self.DDG_db
        predictions = DDG_db.execute_select("SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet=%s", parameters=(PredictionSet,))

        PDB_chain_lengths ={}
        for prediction in predictions:
            chain_records = DDG_db.execute_select('SELECT PDBFileID, Chain FROM Experiment INNER JOIN ExperimentChain ON ExperimentID=Experiment.ID WHERE ExperimentID=%s', parameters=(prediction['ExperimentID']))
            num_residues = 0
            for chain_record in chain_records:
                key = (chain_record['PDBFileID'], chain_record['Chain'])

                if PDB_chain_lengths.get(key) == None:
                    fasta = DDG_db.execute_select("SELECT FASTA FROM PDBFile WHERE ID=%s", parameters = (chain_record['PDBFileID'],))
                    assert(len(fasta) == 1)
                    fasta = fasta[0]['FASTA']
                    f = parseFASTAs(fasta)
                    PDB_chain_lengths[key] = len(f[chain_record['PDBFileID']][chain_record['Chain']])
                chain_length = PDB_chain_lengths[key]
                num_residues += chain_length

            print("UPDATE Prediction SET Cost=%0.2f WHERE ID=%d" % (num_residues, prediction['ID']))

            predictions = DDG_db.execute("UPDATE Prediction SET Cost=%s WHERE ID=%s", parameters=(num_residues, prediction['ID'],))


    def _get_dev_protocol_id(self, name):
        dev_protocol_ids = self.DDG_db.execute_select("SELECT ID FROM DevelopmentProtocol WHERE Name=%s", parameters = (name,))
        if len(dev_protocol_ids) == 0:
            return None
        elif len(dev_protocol_ids) == 1:
            return int(dev_protocol_ids[0]['ID'])
        else:
            raise Exception("DevelopmentProtocol table was originally set up so that names are unique; this has obviously changed")


    def _create_dev_protocol(self, name, application, template_command_line):
        dev_prot_record = {
            'Name' : name,
            'Application' : application,
            'TemplateCommandLine' : template_command_line,
        }
        sql, params, record_exists = self.DDG_db.create_insert_dict_string('DevelopmentProtocol', dev_prot_record)
        self.DDG_db.execute(sql, params)


    ###########################################################################################
    ## Data entry layer
    ##
    ## This part of the API is responsible for data entry (e.g. complex definitions)
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    #######################################
    #                                     #
    #  Protein-protein complex data entry #
    #                          public API #
    #                                     #
    #                                     #
    #   PPComplex                         #
    #   PPIPDBPartnerChain                #
    #   PPIPDBSet                         #
    #                                     #
    #  Missing tables:                    #
    #       PPIConformationalChange       #
    #       PPIDatabaseComplex            #
    #       PPIDataSetCrossmap            #
    #                                     #
    #######################################


    @ppi_data_entry
    def find_complex(self, pdb_ids, keywords = [], tsession = None, quiet = True):
        possible_match_ids = []
        for pdb_id in pdb_ids:
            existing_records = self.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
            if existing_records and not quiet:
                colortext.warning('The PDB file {0} exists in the database.'.format(pdb_id))
            complex_ids = self.search_complexes_by_pdb_id(pdb_id)
            if complex_ids:
                if existing_records and not quiet:
                    colortext.warning('The PDB file {0} has associated complexes: {1}'.format(pdb_id, ', '.join(map(str, complex_ids))))
                assert(len(complex_ids) == 1)
                complex_id = complex_ids[0]
                #colortext.warning('Complex #{0}'.format(complex_id))
                #pprint.pprint(self.get_complex_details(complex_id))

            assert(type(keywords) == list)
            keywords = set(keywords)
            for keyword in keywords:

                hits = self.get_complex_ids_matching_protein_name(keyword, tsession = tsession)
                if hits:
                    if not quiet:
                        colortext.warning('Partial match on "{0}".'.format(keyword))
                    possible_match_ids.extend(hits)

        possible_match_idses = sorted(set(possible_match_ids))
        return [self.get_complex_details(id) for id in possible_match_ids]


    @ppi_data_entry
    def add_complex_structure_pair(self, complex_structure_definition_pair, keywords = None, force = False, previously_added = set(), trust_database_content = False, update_sections = set(), allow_missing_params_files = False, debug = False, minimum_sequence_identity = 95.0):
        '''Wrapper function for add_designed_pdb and add_complex.

           complex_structure_definition_pair should be a dict with the structure:
                dict(
                    Structure = <see the definition in import_api:add_designed_pdb>,
                    Complex = <see the definition in ppi_api:add_complex>,
                )

            To simplify the logic, we treat this function call as an atomic call i.e. it creates its own session and rolls back or commits.
        '''

        # Sanity checks
        assert(complex_structure_definition_pair['Complex']['structure_id'] == complex_structure_definition_pair['Structure']['db_id'])
        if 'chain_mapping' in complex_structure_definition_pair['Structure']:
            assert(sorted(complex_structure_definition_pair['Structure']['chain_mapping'].keys()) == sorted(complex_structure_definition_pair['Complex']['LChains'] + complex_structure_definition_pair['Complex']['RChains']))

        # Create a new session
        tsession = self.importer.get_session(new_session = True, utf = False)
        try:
            # Add the structure
            self.importer.add_designed_pdb(complex_structure_definition_pair['Structure'], previously_added = previously_added, trust_database_content = trust_database_content,
                                           update_sections = update_sections, allow_missing_params_files = allow_missing_params_files,
                                           minimum_sequence_identity = minimum_sequence_identity, tsession = tsession, debug = debug)
            if debug:
                tsession.rollback()
            else:
                tsession.commit()
            tsession.close()
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise

        tsession = self.importer.get_session(new_session = True, utf = True)
        try:
            # Add the complex definition and PDB definition
            api_response = self.add_complex(complex_structure_definition_pair['Complex'], keywords = keywords, force = force, debug = debug, tsession = tsession)
            if api_response['success']:
                str(api_response['PPIPDBSet']) # this forced lookup of partner_chains seems to be crucial when accessing it later (which should only be done for printing as the data cannot be guaranteed to be up-to-date)
                tsession.expunge_all() # note: we only need to expunge api_response['PPIPDBSet'].partner_chains (it is loaded lazily/deferred)
                if debug:
                    api_response = dict(success = False, error = 'Debug call - rolling back the transaction.')
                    tsession.rollback()
                else:

                    tsession.commit()
            else:
                tsession.rollback()
            tsession.close()
            return api_response
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise


    def lookup_pdb_set(self, tsession, passed_pdb_set, allow_partial_matches = True, complex_id = None):
        '''Takes a dict {'L' -> List(Tuple(PDB ID, Chain ID)), 'R' -> List(Tuple(PDB ID, Chain ID))} and returns all PDB
           sets (complex_id, set_number, reverse_match) which have either partial or exact matches depending on
           whether allow_partial_matches is True or False respectively. If reverse_match is True it means that the
           partner definitions are reversed (left partner = right partner,...).

           The matching is symmetric over the partner definitions i.e. if L1 matches R2 and R1 matches L2 then we consider this a match.
           If complex_id is specified then we restrict matches to that particular ID (PPComplex.ID). Otherwise, all definitions
           in the database are considered.

           If allow_partial_matches is True then we return hits if there is at least one common chain in each partner.
           Otherwise, we return hits if there are exact matches (modulo chain ordering)
           '''

        defined_sets = {}
        if complex_id != None:
            # Consider sets for a specific complex
            defined_sets[complex_id] = {}
            for r in tsession.query(dbmodel.PPIPDBPartnerChain).filter(dbmodel.PPIPDBPartnerChain.PPComplexID == complex_id):
                set_number = r.SetNumber
                defined_sets[complex_id][set_number] = defined_sets[complex_id].get(set_number, {'L' : [], 'R' : []})
                defined_sets[complex_id][set_number][r.Side].append((r.PDBFileID, r.Chain))
        else:
            # Consider all sets
            for r in tsession.query(dbmodel.PPIPDBPartnerChain):
                set_number = r.SetNumber
                c_id = r.PPComplexID
                defined_sets[c_id] = defined_sets.get(c_id, {})
                defined_sets[c_id][set_number] = defined_sets[c_id].get(set_number, {'L' : [], 'R' : []})
                defined_sets[c_id][set_number][r.Side].append((r.PDBFileID, r.Chain))

        set_number_hits = set()
        for c_id, set_definitions in sorted(defined_sets.iteritems()):
            for set_number, set_partners in sorted(set_definitions.iteritems()):
                # Check for matches against the stored PDB sets. Check for the symmetric definition as well
                if allow_partial_matches:
                    # Partial matching
                    if set(passed_pdb_set['L']).intersection(set_partners['L']) and set(passed_pdb_set['R']).intersection(set_partners['R']):
                        set_number_hits.add((c_id, set_number, False))
                    if set(passed_pdb_set['L']).intersection(set_partners['R']) and set(passed_pdb_set['R']).intersection(set_partners['L']):
                        set_number_hits.add((c_id, set_number, True))
                else:
                    # Exact matching
                    if (sorted(passed_pdb_set['L']) == sorted(set_partners['L'])) and (sorted(passed_pdb_set['R']) == sorted(set_partners['R'])):
                        set_number_hits.add((c_id, set_number, False))
                    if (sorted(passed_pdb_set['L']) == sorted(set_partners['R'])) and (sorted(passed_pdb_set['R']) == sorted(set_partners['L'])):
                        set_number_hits.add((c_id, set_number, True))

        if len(set([t[2] for t in set_number_hits])) > 1:
            raise colortext.Exception('WARNING: the complex definition has at least two PDB sets where the left and right partners are in the reverse direction. This indicates a redundancy in the database.')

        return sorted(set_number_hits)


    def lookup_complex_by_details(self, tsession, complex_details, allow_partial_matches = True):
        '''Takes a complex_details dict (as defined in add_complex) for a bound complex (i.e. a single PDB ID) and returns
           the corresponding complex(es) and PDB set details if the defined complex exists in the database.

           There are two paths. First, we check whether a complex exists with an exact match on all fields in the PPComplex
           table. This case is probably only likely in the case where the same complex definition is being added repeatedly
           e.g. if a data import script is being run over and over again. Next, we check whether a complex exists based on
           the PDB set i.e. whether a complex using the same PDB chains exists in the database.

           Note that this function will NOT detect cases where the same complex is being used as an existing complex in the
           database but where there are differences in the partner names and a different PDB file is being specified. Therefore,
           care must still be taken when adding complexes to the database to ensure that we do not store duplicate definitions.

           This function is mainly useful as a helper function for add_complex to avoid hitting fail branches when force == False
           in that function. It results in cleaner handling of attempts to re-add existing data.

           Note: We ignore the ChainIndex field in PPIPDBPartnerChain - i.e. we treat partner definitions as bags, not sequences

           Returns: a dict mapping:
                complex_id ->  Dict(reverse_match -> Boolean, # reverse_match is None, True, or False and indicates whether or not the matched complex names (L, R) are in the same order
                                    set_numbers -> List(dict(set_number -> set_number, reverse_match = Boolean))) # reverse_match here is True or False and indicates whether or not the matched PDB sets (L, R) are in the same order

        '''

        # todo: this part of the function currently only allows bound complexes as there is a single structure_id parameter
        # todo: this is the part of the code to change to allow the function to handle unbound complexes
        passed_pdb_set = dict(
            L = sorted([(complex_details['structure_id'], c) for c in complex_details['LChains']]),
            R = sorted([(complex_details['structure_id'], c) for c in complex_details['RChains']])
        )

        complex_id = None
        complex_reverse_match = None

        # Try for an exact match
        # This branch is only useful when the user is adding the same definition multiple times i.e. the same names for the complex.
        # This is mostly hit when import scripts are run multiple times.
        complex_record = get_or_create_in_transaction(tsession, dbmodel.PPComplex, complex_details, variable_columns = ['ID'], only_use_supplied_columns = True, read_only = True)
        if complex_record:
            results = [r for r in tsession.query(dbmodel.PPComplex).filter(and_(dbmodel.PPComplex.LName == complex_details['LName'], dbmodel.PPComplex.RName == complex_details['RName']))]
            results += [r for r in tsession.query(dbmodel.PPComplex).filter(and_(dbmodel.PPComplex.LShortName == complex_details['LShortName'], dbmodel.PPComplex.RShortName == complex_details['RShortName']))]
            results += [r for r in tsession.query(dbmodel.PPComplex).filter(and_(dbmodel.PPComplex.LHTMLName == complex_details['LHTMLName'], dbmodel.PPComplex.RHTMLName == complex_details['RHTMLName']))]
            complex_ids = sorted(set([r.ID for r in results]))
            if complex_ids:
                if not len(complex_ids) == 1:
                    raise colortext.Exception('WARNING: Multiple complex definitions (PPComplex.ID = {0}) share the same partner names. This indicates a redundancy in the database.'.format(', '.join(complex_ids)))
                complex_id = complex_ids[0]
                complex_record = tsession.query(dbmodel.PPComplex).filter(dbmodel.PPComplex.ID == complex_id).one()
                complex_reverse_match = False

            results = [r for r in tsession.query(dbmodel.PPComplex).filter(and_(dbmodel.PPComplex.LName == complex_details['RName'], dbmodel.PPComplex.RName == complex_details['LName']))]
            results += [r for r in tsession.query(dbmodel.PPComplex).filter(and_(dbmodel.PPComplex.LShortName == complex_details['RShortName'], dbmodel.PPComplex.RShortName == complex_details['LShortName']))]
            results += [r for r in tsession.query(dbmodel.PPComplex).filter(and_(dbmodel.PPComplex.LHTMLName == complex_details['RHTMLName'], dbmodel.PPComplex.LHTMLName == complex_details['LHTMLName']))]
            complex_ids = sorted(set([r.ID for r in results]))
            if complex_ids:
                if (complex_id != None) or (len(complex_ids) != 1):
                    raise colortext.Exception('WARNING: Multiple complex definitions (PPComplex.ID = {0}) share the same partner names. This indicates a redundancy in the database.'.format(', '.join(complex_ids)))
                complex_id = complex_ids[0]
                complex_record = tsession.query(dbmodel.PPComplex).filter(dbmodel.PPComplex.ID == complex_id).one()
                complex_reverse_match = True

        if complex_record:
            # We found an associated PPComplex record. Now we check to see whether an associated PPIPDBSet exists
            complex_id = complex_record.ID

            # todo: this part of the function allows unbound complexes and does not need to be updated
            set_number_hits = self.lookup_pdb_set(tsession, passed_pdb_set, allow_partial_matches = allow_partial_matches, complex_id = complex_id)

            # One exact hit for the complex definition with one or many PDB sets
            l = []
            for h in set_number_hits:
                assert(h[0] == complex_id)
                set_number = h[1]
                reverse_match = h[2]
                assert(complex_reverse_match == reverse_match)
                l.append(dict(set_number = set_number, reverse_match = reverse_match))
            return {complex_id : dict(reverse_match = complex_reverse_match, set_numbers = l)}
        else:
            # The complex did not exactly match a PPComplex record however there may simply be differences in the partner names.
            # We proceed by looking for a match based on the PDB chains by checking all PDB sets.
            set_number_hits = self.lookup_pdb_set(tsession, passed_pdb_set, allow_partial_matches = allow_partial_matches)
            results_by_complex = {}
            for h in set_number_hits:
                complex_id = h[0]
                set_number = h[1]
                reverse_match = h[2]
                results_by_complex[complex_id] = results_by_complex.get(complex_id, dict(reverse_match = None, set_numbers = []))
                results_by_complex[complex_id]['set_numbers'].append(dict(set_number = set_number, reverse_match = reverse_match))
            return results_by_complex
        return None


    @ppi_data_entry
    def add_complex(self, complex_details, keywords = [], force = False, debug = False, tsession = None):
        '''Add a complex to the database using a defined dict structure.

        :param complex_details: A dict fitting the defined structure (see below).
        :param keywords: A list of keywords used to search existing complexes for an existing match. Not necessary but
                         advised, particularly when adding a small number of complexes.
        :param force: If a potentially similar complex is found and force is False then then the function returns with a
                      message and without adding the complex. The ForceAddition setting in the Complex dict (see below)
                      will have the same effect as setting this variable.
        :param debug: If debug is set to True then the transaction used to insert the complex into the database will be
                      rolled back and a message stating that the insertion would have been successful is returned in the
                      return dict.
        :return: On successful import, the dict
                   {success = True, ComplexID -> Long, SetNumber -> Long, ReverseMatch -> Boolean}
                 corresponding to the database PPIPDBSet primary key is returned. ReverseMatch is True if the complex was
                 found in the database with the same partner ordering (Left = Left, Right = Right) and False otherwise.
                 If a similar complex is detected and force is False then a dict
                    {success = False, ComplexID -> Long, SetNumber -> Long, ReverseMatch -> Boolean, message -> String}
                 will be returned instead.
                 On error, a dict {success = False, error -> String} is returned.

        The database uses Unicode to encode the strings, allowing us to use e.g. Greek characters
        For this reason, please contain all structure definitions in a file encoded as Unicode. On Linux, you can add the
            # -*- coding: utf-8 -*-
        declaration at the top of the file (with no leading whitespace).

        One example of the dict structure is as follows:

            dict(
                # There are two cases - the complex exists in the database or we will be adding a new complex.
                # Note: Before adding new complexes, you should make sure that there is no existing complex in the
                #       database. This will help to reduce redundancy and provide us with better data.

                # These fields are required in both cases and specify the partners of the complex
                # Note: Please ensure that the LChains (resp. RChains) chains correspond to the protein/complex
                # identified by LName, LShortName, LHTMLName (resp. RName, RShortName, RHTMLName)
                structure_id = '1A2K_TP0',
                LChains = ['A'],
                RChains = ['C'],

                # Case 1: These fields should be used if there is an existing complex in the database.
                ComplexID = 202,

                # Case 2: These fields should only be used if there is no existing complex in the database.
                AdditionalKeywords = ['GSP1'], # Used to search for existing complexes. The PDB ID, LName, LShortName, etc. fields will automatically be used for the search so there is no need to specify those.
                LName = 'Ras-related nuclear protein', # the full protein name for the left partner. This is a Unicode field.
                LShortName = 'RAN', # the short-hand name commonly used
                LHTMLName = 'RAN', # a version of the short-hand name converted to HTML e.g. &alpha; used in place of an alpha character. This is an ASCII field.
                RName = 'Ran-specific GTPase-activating protein', # similar
                RShortName = 'RanGAP1', # similar
                RHTMLName = 'RanGAP1', # similar
                FunctionalClassID = 'OG', # One of A (Antibody-antigen), AB (Antigen/Bound Antibody), EI (Enzyme/inhibitor), ER (Enzyme containing complex),
                                          # ES (Enzyme containing complex), OG (G-proteins), OR (Receptors), or OX (Miscellaneous)
                PPDBMFunctionalClassID = 'O', # One of A (Antibody-antigen), AB (Antigen/Bound Antibody), E (Enzyme/Inhibitor or Enzyme/Substrate), or O (Miscellaneous)
                PPDBMDifficulty = None,   # specific to the protein-protein docking benchmark i.e. use None here
                IsWildType = True,        # if this is the wildtype sequence
                WildTypeComplexID = None, # if this is not wildtype sequence and the wildtype complex is in the database, please specify that complex ID here
                Notes = '...'             # any notes on the complex e.g. 'There is a related complex in the database (complex #119 at the time of writing) with all three unique chains from 1K5D (AB|C).'
                Warnings = None,          # any warnings about the complex in general. Note: Structural warnings belong in the description field of the Structure dict.

                # Optional fields for either case
                PDBComplexNotes = '...'   # any notes specific to the particular PDB structure rather than the complex
                DatabaseKeys = [                        # Used when adding complexes from databases to help map them back to that database
                    dict(
                        DatabaseName = "SKEMPI",
                        DatabaseKey = "1NCA_N_LH",
                    ),
                    ...
                ]

            )
        '''
        # todo: this function currently only adds bound complexes (which is the typical case). It is straightforward to generalize the structure above for unbound complexes e.g. by changing LChains and RChains to include structure ids

        existing_session = not(not(tsession))
        tsession = tsession or self.importer.get_session(new_session = True, utf = True)

        # Search for exact matches first, then partial matches
        pp_complex = None
        reverse_match = None
        for match_param in [False, True]:
            existing_complexes = self.lookup_complex_by_details(tsession, complex_details, allow_partial_matches = match_param)
            if existing_complexes:
                if len(existing_complexes) == 1:
                    existing_complex_id = existing_complexes.keys()[0]
                    pp_complex = tsession.query(dbmodel.PPComplex).filter(dbmodel.PPComplex.ID == existing_complex_id)
                    if 'ComplexID' in complex_details:
                        if complex_details['ComplexID'] != pp_complex.ID:
                            raise colortext.Exception('ComplexID {0} was passed but complex #{1} was found which seems to match the complex definition.'.format(complex_details['ComplexID'], pp_complex.ID))
                    reverse_match = existing_complexes[existing_complex_id]['reverse_match']
                    existing_pdb_sets = existing_complexes[existing_complex_id]['set_numbers']
                    if existing_pdb_sets:
                        if len(existing_pdb_sets) == 1:
                            existing_pdb_set = existing_pdb_sets[0]
                            msg = None
                            if match_param == True:
                                msg = 'A match was found on the partner/PDB set definition but the complex fields had different values e.g. different names of each partner.'
                                if not force:
                                    return dict(success = False, message = msg, ComplexID = existing_complex_id, SetNumber = existing_pdb_set['set_number'], ReverseMatch = existing_pdb_set['reverse_match'])
                                else:
                                    colortext.warning(msg)
                                    return dict(success = True, message = msg, ComplexID = existing_complex_id, SetNumber = existing_pdb_set['set_number'], ReverseMatch = existing_pdb_set['reverse_match'])
                            return dict(success = True, ComplexID = existing_complex_id, SetNumber = existing_pdb_set['set_number'], ReverseMatch = existing_pdb_set['reverse_match'])
                        else:
                            raise colortext.Exception('The complex definition exists in the database but multiple PDB sets / partner definitions match the passed parameters. Check this case manually.')
                    else:
                        # If force is not passed, raise an exception. Else, cascade into the new partner definition creation below.
                        if not force:
                            raise colortext.Exception('The complex definition exists in the database although no PDB sets / partner definitions corresponding EXACTLY to the partner definition were found. Check this case manually to see whether existing definitions would suit better than the passed definition (else, the force parameter can be passed to force creation of a new definition).')
                else:
                    raise colortext.Exception('Multiple complex definitions exists in the database which match the passed complex definition. Check this case manually.')

        # We have not found an exact match or (if force == True) a similar match has been found.
        # If force is False and a similar complex was found, we should have raise an exception above.
        try:
            assert('DatabaseKeys' not in complex_details) # todo: write this code

            # Check parameters
            passed_keys = sorted(complex_details.keys())
            expected_keys = ['structure_id', 'LChains', 'RChains']
            for k in expected_keys:
                assert(k in complex_details)
            structure_id, LChains, RChains = complex_details['structure_id'], complex_details['LChains'], complex_details['RChains']

            # Check that the structure is already in the database
            structure_record = None
            try:
                structure_record = tsession.query(dbmodel.PDBFile).filter(dbmodel.PDBFile.ID == structure_id).one()
            except:
                raise Exception('The structure "{0}" does not exist in the database.'.format(structure_id))

            # Add the PPComplex record
            if pp_complex:
                if reverse_match == True:
                    raise Exception('Write this case. We should add the passed chains in the opposite order (L = R, R = L) since the found complex has the opposite partner ordering.')
                else:
                    assert(force)
                    assert(reverse_match == False) # i.e. it is not equal to None
            else:
                pp_complex = None
                if 'ComplexID' in complex_details:
                    expected_keys.append('ComplexID')
                    if (('PDBComplexNotes' in complex_details) and len(complex_details) != 5) or (('PDBComplexNotes' not in complex_details) and (len(complex_details) != 4)):
                        raise Exception('As the ComplexID was specified, the only expected fields were "{0}" but "{1}" were passed.'.format('", "'.join(sorted(expected_keys)), '", "'.join(passed_keys)))
                    pp_complex = tsession.query(dbmodel.PPComplex).filter(dbmodel.PPComplex.ID == complex_details['ComplexID']).one()
                else:
                    keywords = keywords + [complex_details['LName'], complex_details['LShortName'], complex_details['LHTMLName'], complex_details['RName'], complex_details['RShortName'], complex_details['RHTMLName']]
                    if complex_details.get('AdditionalKeywords'):
                        keywords.extend(complex_details['AdditionalKeywords'])

                    possible_matches = self.find_complex([structure_id], keywords, tsession = tsession)
                    if possible_matches:
                        if not force:
                            return dict(success = False, debug = debug, error = 'Complexes exist in the database which may be related. Please check whether any of these complexes match your case.', possible_matches = possible_matches)
                        colortext.warning('Complexes exist in the database which may be related. Continuing to add a new complex regardless.')

                    pp_complex = get_or_create_in_transaction(tsession, dbmodel.PPComplex, dict(
                        LName = complex_details['LName'],
                        LShortName = complex_details['LShortName'],
                        LHTMLName = complex_details['LHTMLName'],
                        RName = complex_details['RName'],
                        RShortName = complex_details['RShortName'],
                        RHTMLName = complex_details['RHTMLName'],
                        FunctionalClassID = complex_details['FunctionalClassID'],
                        PPDBMFunctionalClassID = complex_details['PPDBMFunctionalClassID'],
                        PPDBMDifficulty = complex_details['PPDBMDifficulty'],
                        IsWildType = complex_details['IsWildType'],
                        WildTypeComplexID = complex_details['WildTypeComplexID'],
                        Notes = complex_details['Notes'],
                        Warnings = complex_details['Warnings'],
                    ), missing_columns = ['ID'])

            # Search for an existing PDB set. Read the current definitions, treating them as bags then sorting lexically
            pdb_sets = {}
            for pschain in tsession.query(dbmodel.PPIPDBPartnerChain).filter(dbmodel.PPIPDBPartnerChain.PPComplexID == pp_complex.ID):
                pdb_sets[pschain.SetNumber] = pdb_sets.get(pschain.SetNumber, {'L' : [], 'R' : []})
                pdb_sets[pschain.SetNumber][pschain.Side].append((pschain.PDBFileID, pschain.Chain))

            # Create a bag from the new definition then sort lexically
            new_pdb_set = dict(L = sorted([(structure_id, c) for c in LChains]),
                               R = sorted([(structure_id, c) for c in RChains]))

            # Check whether an exact match already exists
            matching_set, reverse_match = None, None
            for set_number, set_def in pdb_sets.iteritems():
                set_def['L'] = sorted(set_def['L'])
                set_def['R'] = sorted(set_def['R'])
                if set_def['L'] == new_pdb_set['L'] and set_def['R'] == new_pdb_set['R']:
                    matching_set, reverse_match = True, False
                elif set_def['L'] == new_pdb_set['R'] and set_def['R'] == new_pdb_set['L']:
                    matching_set, reverse_match = True, True
                if matching_set:
                    pdb_set = tsession.query(dbmodel.PPIPDBSet).filter(and_(dbmodel.PPIPDBSet.PPComplexID == pp_complex.ID, dbmodel.PPIPDBSet.SetNumber == set_number)).one()
                    return dict(success = True, ReverseMatch = reverse_match, ComplexID = pp_complex.ID, SetNumber = set_number) # this used to also return PPIPDBSet = pdb_set

            # No match. Create a new set by adding a PPIPDBSet record.
            if pdb_sets:
                new_set_number = max(pdb_sets.keys()) + 1
            else:
                new_set_number = 0
            assert(tsession.query(dbmodel.PPIPDBSet).filter(and_(dbmodel.PPIPDBSet.PPComplexID == pp_complex.ID, dbmodel.PPIPDBSet.SetNumber == new_set_number)).count() == 0) # Sanity check
            pdb_complex_notes = None
            if 'PDBComplexNotes' in complex_details:
                pdb_complex_notes = complex_details['PDBComplexNotes']
            pdb_set_object = get_or_create_in_transaction(tsession, dbmodel.PPIPDBSet,
                dict(
                    PPComplexID = pp_complex.ID,
                    SetNumber = new_set_number,
                    IsComplex = True, # todo: change when we allow unbound complexes
                    Notes = pdb_complex_notes,
                ))

            # Create the associated PPIPDBPartnerChain records
            for set_side, side_chains in sorted(new_pdb_set.iteritems()):
                chain_index = 0
                for pc in sorted(side_chains):
                    get_or_create_in_transaction(tsession, dbmodel.PPIPDBPartnerChain,
                        dict(
                            PPComplexID = pp_complex.ID,
                            SetNumber = new_set_number,
                            Side = set_side,
                            ChainIndex = chain_index,
                            PDBFileID = pc[0],
                            Chain = pc[1],
                            NMRModel = None, # todo
                        ), missing_columns = ['ID'])
                    chain_index += 1

            # Return the API response
            api_response = dict(success = True, ReverseMatch = False, PPIPDBSet = pdb_set_object, ComplexID = pp_complex.ID, SetNumber = new_set_number) # this used to also return PPIPDBSet = pdb_set_object
            if not(existing_session):
                if debug:
                    api_response = dict(success = False, debug = debug, error = 'Debug call - rolling back the transaction.')
                    tsession.rollback()
                    tsession.close()
                else:
                    tsession.commit()
                    tsession.close()
            return api_response
        except:
            colortext.error('Failure.')
            print(traceback.format_exc())
            tsession.rollback()
            tsession.close()
            raise


    @ppi_data_entry
    def add_user_dataset_case(self, tsession, user_dataset_case, user_dataset_name_to_id_map = {}):
        '''Add a user dataset case to the database using a defined dict structure.

        :param tsession: A transaction session. This must be created and passed into this function as user datasets should
                         be added in one transaction.
        :param user_dataset_case: A single case for the user dataset matching the structure defined below.
        :param user_dataset_name_to_id_map: Used to cache the mapping from user dataset names to their integer IDs
        :return: On success, the UserDataSetExperiment object is returned.


        user_dataset_case should be structured as in the following example:

            dict(

                # These records are used to create a PPMutagenesis record and the associated mutagenesis details

                Mutagenesis = dict(
                    RecognizableString = 'TinaGSP_32',
                    PPComplexID = -1,
                ),

                Mutations = [
                    # There is one dict per mutation
                    dict(
                        MutagenesisMutation = dict(
                            # PPMutagenesisID will be filled in when the PPMutagenesis record is created.
                            RecordKey = 'A D123E', # chain_id, wildtype_aa, residue_id.strip(), mutant_aa
                            ProteinID = None, # todo
                            ResidueIndex = None, # todo
                            WildTypeAA = 'D',
                            MutantAA = 'E',
                        ),
                        MutagenesisPDBMutation = dict(
                            # PPMutagenesisID and PPMutagenesisMutationID will be filled in when the PPMutagenesisMutation record is created.
                            # PPComplexID is taken from the PPMutagenesis section. WildTypeAA and MutantAA are taken from the PPMutagenesisMutation section.
                            SetNumber = -1,
                            PDBFileID = '1A2K_TP0',
                            Chain = 'A',
                            ResidueID = ' 123 ',
                        ),
                    ),
                ],

                # This field is used to create the UserPPDataSetExperiment record. All other fields can be derived from the above.
                # Note: We use the human-readable label here. The database ID is retrieved using e.g. ppi_api.get_defined_user_datasets()[<UserDataSetTextID>]['ID']
                UserDataSetTextID = 'RAN-GSP',
            )
        '''

        udc = user_dataset_case

        # Extract the PDB file and complex set number
        pdb_file_id = set([m['MutagenesisPDBMutation']['PDBFileID'] for m in udc['Mutations']])
        assert(len(pdb_file_id) == 1)
        pdb_file_id = pdb_file_id.pop()
        set_number = set([m['MutagenesisPDBMutation']['SetNumber'] for m in udc['Mutations']])
        assert(len(set_number) == 1)
        set_number = set_number.pop()

        is_wildtype = 1
        if udc['Mutations']:
            is_wildtype = 0

        # 1. Create the mutagenesis record
        pp_mutagenesis = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesis, dict(
                    PPComplexID = udc['Mutagenesis']['PPComplexID'],
                    SKEMPI_KEY = udc['Mutagenesis']['RecognizableString'],
                    WildType = is_wildtype,
                ), missing_columns = ['ID'])

        # 2. Create the PPMutagenesisMutation and PPMutagenesisPDBMutation records
        for m in udc['Mutations']:

            # 2a. Create the PPMutagenesisMutation record
            mmut = m['MutagenesisMutation']
            mmut['PPMutagenesisID'] = pp_mutagenesis.ID

            # Sanity check existing records
            existing_record = tsession.query(dbmodel.PPMutagenesisMutation).filter(and_(
                    dbmodel.PPMutagenesisMutation.PPMutagenesisID == mmut['PPMutagenesisID'], dbmodel.PPMutagenesisMutation.RecordKey == mmut['RecordKey']))
            if existing_record.count() > 0:
                existing_record = existing_record.one()
                assert(existing_record.MutantAA == mmut['MutantAA'])
                assert(existing_record.WildTypeAA == mmut['WildTypeAA'])

            # Add the new record
            pp_mutagenesis_mutation = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesisMutation, mmut, missing_columns = ['ID'])

            # 2b. Create the PPMutagenesisPDBMutation record
            pmut = m['MutagenesisPDBMutation']
            pmut['PPMutagenesisID'] = pp_mutagenesis.ID
            pmut['PPMutagenesisMutationID'] = pp_mutagenesis_mutation.ID
            pmut['PPComplexID'] = pp_mutagenesis.PPComplexID
            pmut['WildTypeAA'] = pp_mutagenesis_mutation.WildTypeAA
            pmut['MutantAA'] = pp_mutagenesis_mutation.MutantAA
            pmut['ResidueID'] = PDB.ResidueID2String(pmut['ResidueID']) # handle stripped strings

            # Sanity check existing records
            existing_record = tsession.query(dbmodel.PPMutagenesisPDBMutation).filter(and_(
                    dbmodel.PPMutagenesisPDBMutation.PPMutagenesisMutationID == pmut['PPMutagenesisMutationID'],
                    dbmodel.PPMutagenesisPDBMutation.PDBFileID == pdb_file_id,
                    dbmodel.PPMutagenesisPDBMutation.SetNumber == set_number,
                    dbmodel.PPMutagenesisPDBMutation.Chain == pmut['Chain'],
                    dbmodel.PPMutagenesisPDBMutation.ResidueID == pmut['ResidueID'],
            ))
            if existing_record.count() > 0:
                existing_record = existing_record.one()
                assert(existing_record.PPMutagenesisID == pmut['PPMutagenesisID'])
                assert(existing_record.PPComplexID == pmut['PPComplexID'])
                assert(existing_record.WildTypeAA == pmut['WildTypeAA'])
                assert(existing_record.MutantAA == pmut['MutantAA'])

            # Add the new record
            pp_mutagenesis_pdb_mutation = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesisPDBMutation, pmut, missing_columns = ['ID'])

        # 3. Create the UserPPDataSetExperiment record
        user_dataset_name = udc['UserDataSetTextID']
        if not user_dataset_name_to_id_map.get(user_dataset_name):
            user_dataset_name_to_id_map[user_dataset_name] = tsession.query(dbmodel.UserDataSet).filter(dbmodel.UserDataSet.TextID == user_dataset_name).one().ID
        user_dataset_id = user_dataset_name_to_id_map[user_dataset_name]

        new_record = True
        if tsession.query(dbmodel.UserPPDataSetExperiment).filter(and_(
                dbmodel.UserPPDataSetExperiment.UserDataSetID == user_dataset_id,
                dbmodel.UserPPDataSetExperiment.PPMutagenesisID == pp_mutagenesis.ID,
                dbmodel.UserPPDataSetExperiment.PDBFileID == pdb_file_id,
                dbmodel.UserPPDataSetExperiment.PPComplexID == pp_mutagenesis.PPComplexID,
                dbmodel.UserPPDataSetExperiment.SetNumber == set_number)).count() > 0:
            new_record = False
        user_dataset_experiment = get_or_create_in_transaction(tsession, dbmodel.UserPPDataSetExperiment, dict(
                UserDataSetID = user_dataset_id,
                PPMutagenesisID = pp_mutagenesis.ID,
                PDBFileID = pdb_file_id,
                PPComplexID = pp_mutagenesis.PPComplexID,
                SetNumber = set_number,
                IsComplex = True,
            ), missing_columns = ['ID'])
        if new_record:
            colortext.wgreen('.')
        else:
            colortext.wcyan('.')


    @general_data_entry
    def add_de_dataset(self, user_id, long_id, short_id, description, ddg_convention, dataset_creation_start_date = None, dataset_creation_end_date = None, publication_ids = [], existing_session = None):
        '''Convenience wrapper for add_dataset for DeltaE-only datasets.'''
        return self.add_dataset(user_id, long_id, short_id, description, False, False, True, ddg_convention, dataset_creation_start_date = dataset_creation_start_date, dataset_creation_end_date = dataset_creation_end_date, publication_ids = publication_ids, existing_session = existing_session)


    @ppi_data_entry
    def add_ssm_dataset(self, dataset_short_id, user_dataset_id, complex_id, set_number, mutations_dataframe, existing_session = None, debug = True):
        '''Import SSM data from an RCSB PDB file. Non-RCSB files are not currently handled. Some data (DataSet and UserDataSet)
           must be set up before calling this function.

        :param dataset_short_id: The short ID of the existing dataset in the database (DataSet.ShortID)
        :param user_dataset_id: The ID of the existing user dataset in the database (UserDataSet.ID)
        :param pp_complex_id: The complex ID used in the database (PPComplex.ID). This will be used to add the structure to the database.
        :param set_number: The set_number of the complex used in the database (PPIPDBSet.SetNumber). This is used to determine the choice of chains in predictions.
        :param mutations_dataframe: A pandas dataframe in the intermediate input format described below.
        :param debug: If True then the transaction is rolled back. This is set to True by default to reduce data-entry errors i.e. you should do a test-run of add_ssm_dataset first and then do a run with debug = False.
        :return: Dict {success : <True/False>, DataSetID : dataset_id, [errors : <list of error strings if failed>]}

        This function requires the complex, DataSet, and UserDataSet records to have been created. Those records can be added using
        the appropriate functions e.g.

            ppi_api = get_ppi_interface(read_file('pw'))

            # If the complex structure has not been added to the database:
            ppi_api.importer.add_pdb_from_rcsb(pdb_id, trust_database_content = True)

            # If the complex has not been added to the database:
            complex_ids = ppi_api.search_complexes_by_pdb_id(pdb_id)
            if complex_ids:
                colortext.warning('The PDB file {0} has associated complexes: {1}'.format(pdb_id, ', '.join(map(str, complex_ids))))
            api_response = ppi_api.add_complex(json.loads(read_file('my_complex.json')[path][to][complex_definition])) # The structure of the JSON file is described in the docstring for add_complex
            if not api_response['success']:
                raise Exception(api_response['error'])
            pp_complex_id, set_number = api_response['ComplexID'], api_response['SetNumber']

            # else if the complex already exists in the database:
            pp_complex_id, set_number = ..., ...

            # Add dataset publications
            publication_ids = [
                ppi_api.add_publication(...).ID, # currently not implemented
                ...
                ppi_api.add_publication(...).ID, # currently not implemented
            ]

            # Add the dataset and user dataset records
            dataset = ppi_api.add_de_dataset('oconchus', 'SSM_Psd95-CRIPT_Rama_10.1038/nature11500', 'Psd95-CRIPT', 'description...', ddg_convention, dataset_creation_start_date = datetime.date(...), dataset_creation_end_date = datetime.date(...), publication_ids = [...])
            user_dataset = ppi_api.add_de_user_dataset('oconchus', 'SSM-Psd95-CRIPT', '...')

            # Finally, import the SSM dataset
            add_ssm_dataset(dataset.ShortID, user_dataset.ID, pp_complex_id, set_number, mutations_dataframe)

        @todo: write the add_publication function (using the RIS parsing module in klab and the PubMed/DOI downloading modules).

        mutations_dataframe should be a complete (either a value or null at all positions in the m x n array) pandas
        dataframe with a standardized structure.
        This simplifies the data import. The dataframe should be indexed/row-indexed by residue type and column-indexed
        by a string chain ID + <underscore> + residue ID without spaces e.g. 'A_311' is residue ' 311 ' of chain A and 'A_312B' is residue ' 312B' of chain A.
        We include an underscore in the format to reduce confusion for cases where the PDB chain ID is an integer.
        For example, if the input file is a TSV formatted like:

            Pos/aa	A_311	A_312	...
            A	0.131	-0.42	...
            C	0.413	-0.022	...
            ...

        then a valid mutations_dataframe can be constructed via

            mutations_dataframe = pandas.read_csv(ssm_input_data_path, sep = '\t', header = 0, index_col = 0)

        '''

        tsession = existing_session or self.get_session(new_session = True, utf = False)

        # Sanity checks
        assert(complex_id != None and set_number != None)
        dataset_id = None
        try:
            dataset_id = tsession.query(dbmodel.DataSet).filter(dbmodel.DataSet.ShortID == dataset_short_id).one().ID
        except:
            raise Exception('No dataset with ShortID "{0}" exists in the database.'.format(dataset_short_id))
        try:
            tsession.query(dbmodel.UserDataSet).filter(dbmodel.UserDataSet.ID== user_dataset_id).one()
        except:
            raise Exception('No user dataset with TextID "{0}" exists in the database.'.format(user_dataset_id))

        # Retrieve the mapping from chain -> residue ID -> wildtype residue
        pdb_id, complex_chains = self.get_bound_pdb_set_details(complex_id, set_number)
        chain_wt_residue_by_pos = self.get_pdb_residues_by_pos(pdb_id, strip_res_ids = True)

        # Sanity checks on column indices
        chain_ids = set()
        for v in mutations_dataframe.columns.values:
            error_msg = 'The column index "{0}" does not have the expected format: <chain>_<residue id> e.g. "A_123".'.format(v)
            if v.find('_') == -1 or len(v.split('_')) != 2:
                raise colortext.Exception(error_msg)
            tokens = v.split('_')
            chain_id = tokens[0]
            residue_id = tokens[1]
            if len(chain_id) != 1 or (not(residue_id.strip().isdigit()) and not(residue_id.strip()[:-1].isdigit())):
                raise colortext.Exception(error_msg)
            chain_ids.add(chain_id)

        # Sanity checks on row indices
        mut_aas = sorted(mutations_dataframe.index)
        expected_mut_aas = set(residue_type_1to3_map.keys())
        expected_mut_aas.remove('X')
        assert(len(expected_mut_aas) == 20)
        if set(mut_aas).difference(expected_mut_aas):
            raise colortext.Exception('The row indices contain values which are non canonical residue types: "{0}".'.format('", "'.join(sorted(set(mut_aas).difference(expected_mut_aas)))))

        # Extract the data into a list of point mutations, iterating by column/position then row/AA
        # Add a single wildtype PPMutagenesis record (essentially a Complex with no corresponding mutation records)
        # For all single PDB mutations in the list
        #     if not wildtype
        #        add a PPMutagenesis record and corresponding mutation records
        #     add a PPIDataSetDE record to represent the original data (experimental data) in the database
        #     add a UserPPDataSetExperiment record to be used to create prediction runs
        #     add a UserPPAnalysisSetDE record to be used when analyzing prediction runs against the experimental data
        #
        # Note that there will be one UserPPAnalysisSetDE record for each mutant but only one record for wildtype even though
        # the wildtype sequence has exactly one corresponding DeltaE for each position. There will be exactly one UserPPAnalysisSetDE
        # record per mutant and one wildtype record for each position however all of the wildtype UserPPAnalysisSetDE records
        # will be associated to the sole wildtype UserPPAnalysisSetDE record.
        colortext.warning('Adding data for complex #{0}, dataset "{1}", user dataset #{2}.'.format(complex_id, dataset_id, user_dataset_id))
        record_number = 0
        mut_aas = list(mutations_dataframe.index)
        res_ids = list(mutations_dataframe.columns.values)
        try:
            # Add a PPMutagenesis record with no mutation records i.e. the wildtype/null 'mutagenesis'
            pp_wt_mutagenesis = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesis, dict(
                PPComplexID = complex_id,
                SKEMPI_KEY = 'SSM {0}| WildType'.format(pdb_id), # todo: this format is ambiguous if we start to store multiple SSM datasets with different choices of bound partners. We should ideally check all PPMutagenesisMutation/PPMutagenesisPDBMutation records on the complex for a match. At present (2016), it is unlikely that we will have many SSM datasets for consideration, never mind overlapping sets.
                WildType = 1,
            ), missing_columns = ['ID',])
            pp_wt_mutagenesis_id = pp_wt_mutagenesis.ID
            first_wt_record_number = None
            for chain_res_id in res_ids:

                tokens = chain_res_id.split('_')
                assert(len(tokens) == 2)
                chain_id = tokens[0]
                assert(len(chain_id) == 1)
                assert(chain_id in chain_wt_residue_by_pos)
                res_id = tokens[1]
                assert(res_id in chain_wt_residue_by_pos[chain_id])
                wt_aa = chain_wt_residue_by_pos[chain_id][res_id]
                for mut_aa in mut_aas:
                    record_number += 1
                    if record_number % 10 == 0:
                        colortext.wgreen('.')
                        sys.stdout.flush()

                    # Add the PPMutagenesis records for mutant cases
                    if mut_aa == wt_aa:
                        ppi_dataset_de_key = 'SSM {0}| WildType'.format(pdb_id)
                        if first_wt_record_number == None:
                            first_wt_record_number = record_number
                        analysis_set_record_number = first_wt_record_number
                        pp_mutagenesis_id = pp_wt_mutagenesis_id
                    else:
                        ppi_dataset_de_key = 'SSM {0}| {1} {2} {3} {4}'.format(pdb_id, chain_id, wt_aa, res_id, mut_aa) # SKEMPI_KEY is a bad name for a field!,
                        analysis_set_record_number = record_number

                        # Add a PPMutagenesis record with no mutation records i.e. the wildtype/null 'mutagenesis'
                        pp_mutagenesis = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesis, dict(
                            PPComplexID = complex_id,
                            SKEMPI_KEY = 'SSM {0}| {1} {2} {3} {4}'.format(pdb_id, chain_id, wt_aa, res_id, mut_aa), # SKEMPI_KEY is a bad name for a field!,
                            WildType = 0,
                        ), missing_columns = ['ID'])
                        pp_mutagenesis_id = pp_mutagenesis.ID
                        #pprint.pprint(pp_mutagenesis.__dict__)

                        pp_mutagenesis_mutation = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesisMutation, dict(
                            PPMutagenesisID = pp_mutagenesis_id,
                            RecordKey = '{0} {1}{2}{3}'.format(chain_id, wt_aa, res_id, mut_aa),
                            ProteinID = None,
                            ResidueIndex = None,
                            WildTypeAA = wt_aa,
                            MutantAA = mut_aa,
                        ), missing_columns = ['ID',])
                        pp_mutagenesis_mutation_id = pp_mutagenesis_mutation.ID
                        #pprint.pprint(pp_mutagenesis_mutation.__dict__)

                        pp_mutagenesis_pdb_mutation = get_or_create_in_transaction(tsession, dbmodel.PPMutagenesisPDBMutation, dict(
                            PPMutagenesisID = pp_mutagenesis_id,
                            PPMutagenesisMutationID = pp_mutagenesis_mutation_id,
                            PPComplexID = complex_id,
                            SetNumber = set_number,
                            PDBFileID = pdb_id,
                            Chain = chain_id,
                            WildTypeAA = wt_aa,
                            ResidueID = PDB.ResidueID2String(res_id),
                            MutantAA = mut_aa,
                        ), missing_columns = ['ID',])
                        pp_mutagenesis_pdb_mutation_id = pp_mutagenesis_pdb_mutation.ID
                        #pprint.pprint(pp_mutagenesis_pdb_mutation.__dict__)

                    # Add a DeltaE measurement record (PPIDataSetDE)
                    ppi_dataset_de = get_or_create_in_transaction(tsession, dbmodel.PPIDataSetDE, dict(
                        SecondaryID = ppi_dataset_de_key, # optional field
                        DataSetID = dataset_id,
                        Section = 'Supplementary Information II',
                        RecordNumber = record_number,
                        DE = mutations_dataframe[chain_res_id][mut_aa],
                        DEUnit = 'DeltaE (see DataSet.Description)',
                        PublishedError = None,
                        NumberOfMeasurements = None,
                        PPMutagenesisID = pp_mutagenesis_id,
                        PPComplexID = complex_id,
                        SetNumber = set_number,
                        PublishedPDBFileID = pdb_id,
                        PossibleError = False,
                        Remarks = None,
                        IsABadEntry = 0,
                        AddedBy = 'oconchus',
                        AddedDate = datetime.datetime.now(),
                        LastModifiedBy = 'oconchus',
                        LastModifiedDate = datetime.datetime.now(),
                    ), missing_columns = ['ID',], variable_columns = ['AddedDate', 'LastModifiedDate'])
                    ppi_dataset_de_id = ppi_dataset_de.ID

                    # Add a record (UserPPDataSetExperiment) to be included in the associated prediction run
                    user_pp_dataset_experiment = get_or_create_in_transaction(tsession, dbmodel.UserPPDataSetExperiment, dict(
                        UserDataSetID = user_dataset_id,
                        PPMutagenesisID = pp_mutagenesis_id,
                        PDBFileID = pdb_id,
                        PPComplexID = complex_id,
                        SetNumber = set_number,
                        IsComplex = 1
                    ), missing_columns = ['ID',])
                    user_pp_dataset_experiment_id = user_pp_dataset_experiment.ID

                    # dd a record (UserPPAnalysisSetDE) to be used in the analysis, linking the UserPPDataSetExperiment with the DeltaE (PPIDataSetDE) record
                    user_pp_analysis_set_de = get_or_create_in_transaction(tsession, dbmodel.UserPPAnalysisSetDE, dict(
                        Subset = 'Psd95-Cript',
                        Section = 'McLaughlin2012',
                        RecordNumber = analysis_set_record_number,
                        UserPPDataSetExperimentID = user_pp_dataset_experiment_id,
                        PPIDataSetDEID = ppi_dataset_de_id,
                        PPMutagenesisID = pp_mutagenesis_id,
                    ), missing_columns = ['ID',])
                    user_pp_analysis_set_de_id = user_pp_analysis_set_de.ID

            if debug:
                colortext.warning('\nDEBUG MODE IS SET. THE CODE RAN SUCCESSFULLY BUT THE DATASET WILL NOT BE ADDED. RE-RUN THIS FUNCTION WITH debug = False.')
                tsession.rollback()
                tsession.close()
            else:
                tsession.commit()
                tsession.close()
        except Exception, e:
            tsession.rollback()
            tsession.close()
            colortext.warning(traceback.format_exc())
            raise colortext.Exception(str(e))
