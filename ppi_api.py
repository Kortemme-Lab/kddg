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
import traceback
import gzip
import shutil
import sqlite3
import cPickle as pickle
import datetime

import numpy
from sqlalchemy import and_

from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation
from klab.fs.fsio import read_file, write_temp_file
from klab.rosetta.input_files import Mutfile, Resfile
from klab.benchmarking.analysis.ddg_binding_affinity_analysis import DBBenchmarkRun as BindingAffinityBenchmarkRun
from klab.bio.alignment import ScaffoldModelChainMapper
from klab.bio.clustalo import PDBChainSequenceAligner
from klab.db.sqlalchemy_interface import row_to_dict, get_or_create_in_transaction

import db_schema as dbmodel
from api_layers import *
from db_api import ddG, PartialDataException, SanityCheckException


DeclarativeBase = dbmodel.DeclarativeBase


def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, port = 3306):
    '''This is the function that should be used to get a BindingAffinityDDGInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(BindingAffinityDDGInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port)


def get_interface_with_config_file(host_config_name = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None, get_interface_factory = get_interface, passed_port = None):
    # Uses ~/.my.cnf to get authentication information
    ### Example .my.cnf (host_config_name will equal guybrush2):
    ### [clientguybrush2]
    ### user=myname
    ### password=notmyrealpass
    ### host=guybrush.ucsf.edu
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


    def __init__(self, passwd = None, username = 'kortemmelab', hostname = None, rosetta_scripts_path = None, rosetta_database_path = None, port = 3306):
        super(BindingAffinityDDGInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port)
        self.prediction_data_path = self.DDG_db.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionPPIDataPath"')[0]['Value']


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
    def get_chains_for_mutatagenesis(self, mutagenesis_id, pdb_file_id, pdb_set_number, complex_id = None):
        '''Returns a dictionary mapping 'L' to the list of left chains and 'R' to the list of right chains.
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.
        '''

        pp_mutagenesis = self.DDG_db.execute_select("SELECT * FROM PPMutagenesis WHERE ID=%s", parameters = (mutagenesis_id,))

        # Sanity checks
        assert(len(pp_mutagenesis) == 1)
        if complex_id:
            assert(pp_mutagenesis[0]['PPComplexID'] == complex_id)
            pdb_set = self.DDG_db.execute_select("SELECT * FROM PPIPDBSet WHERE PPComplexID=%s AND SetNumber=%s", parameters = (complex_id, pdb_set_number))
            assert(len(pdb_set) == 1 and pdb_set[0]['IsComplex'] == 1) # complex structure check
        else:
            complex_id = pp_mutagenesis[0]['PPComplexID']

        complex_chains = dict(L = [], R = [])
        crecords = self.DDG_db.execute_select("SELECT * FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND SetNumber=%s ORDER BY ChainIndex", parameters = (complex_id, pdb_set_number))
        for c in crecords:
            assert(c['PDBFileID'] == pdb_file_id) # complex structure check
            complex_chains[c['Side']].append(c['Chain'])
        assert(complex_chains['L'] and complex_chains['R'])
        assert(len(set(complex_chains['L']).intersection(set(complex_chains['R']))) == 0) # in one unbound case, the same chain appears twice on one side (2CLR_DE|1CD8_AA, may be an error since this was published as 1CD8_AB but 1CD8 has no chain B) but it seems reasonable to assume that a chain should only appear on one side
        return complex_chains


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
        return json.dumps(self.get_prediction_set_case_details(prediction_set_id, retrieve_references = retrieve_references))

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
        self.DDG_db._get_connection()
        con = self.DDG_db.connection
        with con:
            cur = con.cursor()
            for prediction_id in prediction_ids:
                qry = 'UPDATE PredictionPPI SET DevelopmentProtocolID=%s WHERE ID=%s'
                cur.execute(qry, (dev_protocol_id, prediction_id))

                if rosetta_script:
                    # Add the Rosetta script to the database, not using cursor
                    file_content_id = self._add_prediction_file(prediction_id, rosetta_script, os.path.basename(rosetta_script_file), 'RosettaScript', 'RosettaScript', 'Input', rm_trailing_line_whitespace = True, forced_mime_type = 'text/xml', file_content_id = file_content_id )


    @job_creator
    def add_job(self, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = None):
        '''This function inserts a prediction into the database.
            The parameters define:
                - the prediction set id used to group this prediction with other predictions for analysis;
                - the protocol to be used to run the prediction;
                - the set of mutations and PDB complex associated with the mutagenesis experiment;
                - whether HETATM lines are to be kept or not.
                - additional Rosetta flags e.g. "-ignore_zero_occupancy false" used to determine the mapping from PDB to Rosetta numbering. These flags should correspond to those used in the protocol otherwise errors could occur.
            We strip the PDB based on the chains defined by the complex and keep_hetatm_lines and store the PDB in the database.
            Next, the mapping from Rosetta numbering to PDB numbering is determined and stored in the database.
            Then, the appropriate input files e.g. resfiles or mutfiles are generated and stored in the database.
            Finally, we add the prediction record and associate it with the generated files.'''
        return self._add_job(prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)


    @job_creator
    def add_job_by_user_dataset_record(self, prediction_set_id, user_dataset_name, user_dataset_experiment_id, protocol_id, extra_rosetta_command_flags = None, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = None):
        '''Add a prediction job based on a user dataset record. This is typically called during add_prediction_run rather than directly by the user.
           user_dataset_name is implied by user_dataset_experiment_id but we include it for sanity checking errors in data-entry.

           The extra_rosetta_command_flags variable is used to add additional flags e.g. "-ignore_zero_occupancy false". These should be added if they are used in the protocol.'''

        try:
            user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
        except:
            raise colortext.Exception('The user dataset "%s" does not exist for this API.' % user_dataset_name)

        ude = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(user_dataset_experiment_id, user_dataset_id))
        if not len(ude) == 1:
            raise colortext.Exception('User dataset experiment %d does not exist for/correspond to this user dataset.' % user_dataset_experiment_id)
        ude = ude[0]
        #colortext.message(pprint.pformat(ude))
        return self._add_job(prediction_set_id, protocol_id, ude['PPMutagenesisID'], ude['PPComplexID'], ude['PDBFileID'], ude['SetNumber'], extra_rosetta_command_flags = extra_rosetta_command_flags, user_dataset_experiment_id = user_dataset_experiment_id, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)


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
    def add_prediction_run(self, prediction_set_id, user_dataset_name, extra_rosetta_command_flags = None, protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, input_files = {}, quiet = False, test_only = False, only_single_mutations = False, short_run = False, test_run_first = True, show_full_errors = False):
        '''Adds all jobs corresponding to a user dataset e.g. add_prediction_run("my first run", "AllBindingAffinityData", tagged_subset = "ZEMu").
           If keep_hetatm_lines is False then all HETATM records for the PDB prediction chains will be removed. Otherwise, they are kept.
           input_files is a global parameter for the run which is generally empty. Any files added here will be associated to all predictions in the run.

           The extra_rosetta_command_flags parameter e.g. "-ignore_zero_occupancy false" is used to determine the mapping
           from PDB to Rosetta numbering. These flags should correspond to those used in the protocol otherwise errors could occur.

           Returns False if no predictions were added to the run else return True if all predictions (and there were some) were added to the run.'''

        # Check preconditions
        assert(not(input_files)) # todo: do something with input_files when we use that here - call self._add_file_content, associate the filenames with the FileContent IDs, and pass that dict to add_job which will create PredictionPPIFile records
        assert(only_single_mutations == False) # todo: support this later? it may make more sense to just define new UserDataSets
        self._add_prediction_run_preconditions(prediction_set_id, user_dataset_name, tagged_subset)

        # Get the list of user dataset experiment records
        user_dataset_experiments = self.get_user_dataset_experiment_ids(user_dataset_name, tagged_subset = tagged_subset)
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
        pdb_residues_to_rosetta_cache = {}

        # Test all predictions before creating records
        if test_only or test_run_first:
            if not quiet:
                colortext.message('Testing %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
            # Progress counter setup
            count, records_per_dot = 0, 50
            showprogress = not(quiet) and len(user_dataset_experiments) > 300
            if showprogress: print("|" + ("*" * (int(len(user_dataset_experiments)/records_per_dot)-2)) + "|")
            for ude in user_dataset_experiments:
                # If the mutagenesis already exists in the prediction set, do not test it again
                if protocol_id:
                    existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
                else:
                    existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID IS NULL", parameters=(prediction_set_id, ude['ID']))
                if len(existing_results) > 0: continue

                # Test the prediction setup
                prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = True, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)
                # Progress counter
                count += 1
                if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
                if short_run and count > 4: break
            if not quiet: print('')

        if test_only:
            return

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
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
                    prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = False, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)
                except Exception, e:
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
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
            if showprogress and count % records_per_dot == 0: colortext.write(".", "green", flush = True)
            if short_run and count > 4: break

        if failed_jobs:
            colortext.error('Some jobs failed to run:\n%s' % pprint.pformat(failed_jobs))
        if not quiet: print('')
        return True


    def _create_pdb_residues_to_rosetta_cache_mp(self, pdb_residues_to_rosetta_cache, pdb_file_id, pdb_chains_to_keep, extra_rosetta_command_flags, keep_hetatm_lines):
        # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
        raise Exception('Shane should finish this')
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
        assert(not(input_files)) # todo: do something with input_files when we use that here - call self._add_file_content, associate the filenames with the FileContent IDs, and pass that dict to add_job which will create PredictionPPIFile records
        assert(only_single_mutations == False) # todo: support this later? it may make more sense to just define new UserDataSets
        self._add_prediction_run_preconditions(prediction_set_id, user_dataset_name, tagged_subset)

        # Get the list of user dataset experiment records
        user_dataset_experiments = self.get_user_dataset_experiment_ids(user_dataset_name, tagged_subset = tagged_subset)
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
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
                    prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = False, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)
                except Exception, e:
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
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


    def _add_job(self, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, user_dataset_experiment_id = None, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = {}):
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
        pp_complex = self.DDG_db.execute_select("SELECT * FROM PPComplex WHERE ID=%s", parameters = (pp_complex_id,))

        # Determine the list of PDB chains that will be kept
        pdb_chains = self.get_chains_for_mutatagenesis(pp_mutagenesis_id, pdb_file_id, pp_complex_pdb_set_number, complex_id = pp_complex_id)
        pdb_chains_to_keep = set(pdb_chains['L'] + pdb_chains['R'])
        if self.rosetta_database_path:
            cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, self.rosetta_database_path, extra_rosetta_command_flags)
        else:
            cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, extra_rosetta_command_flags)

        if cache_maps and pdb_residues_to_rosetta_cache.get(cache_key):
            stripped_p = pdb_residues_to_rosetta_cache[cache_key]['stripped_p']
        else:
            # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
            pdb_file = self.DDG_db.execute_select("SELECT * FROM PDBFile WHERE ID=%s", parameters = (pdb_file_id,))
            p = PDB(pdb_file[0]['Content'])
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
            colortext.error('%s: %s, chains %s' % (str(e), str(stripped_p.pdb_id), str(pdb_chains_to_keep)))

        # Assert that there are no empty sequences
        assert(sorted(stripped_p.atom_sequences.keys()) == sorted(pdb_chains_to_keep))
        for chain_id, sequence in stripped_p.atom_sequences.iteritems():
            assert(len(sequence) > 0)

        # Get the PDB mutations and check that they make sense in the context of the stripped PDB file
        # Note: the schema assumes that at most one set of mutations can be specified per PDB file per complex per mutagenesis. We may want to relax that in future by adding the SetNumber to the PPMutagenesisPDBMutation table
        complex_mutations = self.DDG_db.execute_select('SELECT * FROM PPMutagenesisMutation WHERE PPMutagenesisID=%s', parameters=(pp_mutagenesis_id,))
        pdb_complex_mutations = self.DDG_db.execute_select('SELECT * FROM PPMutagenesisPDBMutation WHERE PPMutagenesisID=%s AND PPComplexID=%s AND PDBFileID=%s', parameters=(pp_mutagenesis_id, pp_complex_id, pdb_file_id))
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

        # Make mutfile and resfile
        rosetta_mutations = stripped_p.map_pdb_residues_to_rosetta_residues(mutations)
        mf = Mutfile.from_mutagenesis(rosetta_mutations)
        mutfile_name = 'mutations.mutfile'
        rf = Resfile.from_mutageneses(mutations)
        resfile_name = 'mutations.resfile'

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
        parameter_sub_dict = {
            '%%input_pdb%%' : pdb_filename,
            '%%chainstomove%%' : pdb_chains_to_move_str,
            '%%pathtoresfile%%' : resfile_name,
            '%%pathtomutfile%%' : mutfile_name,
        }

        if test_only:
            return

        # All functions within the next with block should use the same database cursor.
        # The commands then function as parts of a transaction which is rolled back if errors occur within the block
        # or else is committed.
        self.DDG_db._get_connection()
        con = self.DDG_db.connection
        with con:
            cur = con.cursor()

            if protocol_id:
                qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID=%%s' % self._get_prediction_table()
                cur.execute(qry, (prediction_set_id, user_dataset_experiment_id, protocol_id))
            else:
                qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID IS NULL' % self._get_prediction_table()
                cur.execute(qry, (prediction_set_id, user_dataset_experiment_id))
            existing_record = cur.fetchall()
            assert(len(existing_record) == 0)

            prediction_record = dict(
                PredictionSet = prediction_set_id,
                PPMutagenesisID = pp_mutagenesis_id,
                UserPPDataSetExperimentID = user_dataset_experiment_id,
                ProtocolID = protocol_id,
                JSONParameters = json.dumps(parameter_sub_dict),
                DevelopmentProtocolID = None,
                ExtraParameters = extra_rosetta_command_flags,
                Status = 'queued',
                Cost = total_num_residues,
                KeptHETATMLines = keep_hetatm_lines,
            )
            sql, params, record_exists = self.DDG_db.create_insert_dict_string(self._get_prediction_table(), prediction_record, ['PredictionSet', 'UserPPDataSetExperimentID', 'ProtocolID'])
            cur.execute(sql, params)
            prediction_id = cur.lastrowid

            # Add the stripped PDB file
            self._add_prediction_file(prediction_id, '\n'.join(stripped_p.lines), pdb_filename, 'PDB', 'StrippedPDB', 'Input', db_cursor = cur, rm_trailing_line_whitespace = True, forced_mime_type = 'chemical/x-pdb')

            # Add the mutfile
            self._add_prediction_file(prediction_id, str(mf), mutfile_name, 'Mutfile', 'Mutfile', 'Input', db_cursor = cur, rm_trailing_line_whitespace = True)
            # Add the resfile
            self._add_prediction_file(prediction_id, str(rf), resfile_name, 'Resfile', 'Resfile', 'Input', db_cursor = cur, rm_trailing_line_whitespace = True)

            # Add the residue mappings
            self._add_prediction_file(prediction_id, rosetta_to_atom_residue_map, 'rosetta2pdb.resmap.json', 'RosettaPDBMapping', 'Rosetta residue->PDB residue map', 'Input', db_cursor = cur, forced_mime_type = "application/json")
            self._add_prediction_file(prediction_id, atom_to_rosetta_residue_map, 'pdb2rosetta.resmap.json', 'RosettaPDBMapping', 'PDB residue->Rosetta residue map', 'Input', db_cursor = cur, forced_mime_type = "application/json")

        if protocol_id:
            qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID=%%s' % self._get_prediction_table()
            existing_record = self.DDG_db.execute_select(qry, parameters=(prediction_set_id, user_dataset_experiment_id, protocol_id))
        else:
            qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID IS NULL' % self._get_prediction_table()
            existing_record = self.DDG_db.execute_select(qry, parameters=(prediction_set_id, user_dataset_experiment_id))
        assert(len(existing_record) == 1)
        prediction_id = existing_record[0]['ID']
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
    def determine_best_pair(self, prediction_id, score_method_id, expectn = None, returnn = 1):
        '''This returns the best wildtype/mutant pair for a prediction given a scoring method. NOTE: Consider generalising this to the n best pairs.'''
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        mutant_complexes = []
        wildtype_complexes = []
        for structure_id, scores in scores.iteritems():
            if scores.get('MutantComplex'):
                mutant_complexes.append((scores['MutantComplex']['total'], structure_id))
            if scores.get('WildTypeComplex'):
                wildtype_complexes.append((scores['WildTypeComplex']['total'], structure_id))
        wildtype_complexes = sorted(wildtype_complexes)
        mutant_complexes = sorted(mutant_complexes)
        if wildtype_complexes and mutant_complexes:
            return wildtype_complexes[0][1], mutant_complexes[0][1]
        return None, None


    def get_pubs_job_archive(self, *args, **kw):
        prediction_type = kw.get('PredictionType')
        prediction_id = kw.get('PredictionID')
        if prediction_id and prediction_type:
            if prediction_type == 'Binding affinity':
                source_path = os.path.join('/home/kyleb/gits/ddg/job_output/151106-kyleb_dipubs/{0}-ddg/'.format(str(prediction_id)))
                return self.zip_job_files(source_path, 'job_data_{0}_ba.zip'.format(prediction_id))
            elif prediction_type == 'Monomeric stability':
                source_filepath = os.path.join('/kortemmelab', 'shared', 'DDG', 'jobs', '{0}.zip'.format(str(prediction_id)))
                try:
                    assert(os.path.exists(source_filepath))
                    archive = read_file(source_filepath)
                    tg.response.headers['Content-type'] = 'application/zip'
                    tg.response.headers['Content-Disposition'] = 'attachment;filename=job_data_{0}_ms.zip'.format(prediction_id)
                    return archive
                except Exception, e:
                    return {'success' : False, 'message' : traceback.format_exc()}
        else:
            return {'success' : False}

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


    def _get_prediction_data(self, prediction_id, score_method_id, main_ddg_analysis_type, top_x = 3, expectn = None, extract_data_for_case_if_missing = False, root_directory = None, dataframe_type = "Binding affinity", prediction_data = {}):
        try:
            top_x_ddg = self.get_top_x_ddg(prediction_id, score_method_id, top_x = top_x, expectn = expectn)
        except Exception, e:
            colortext.pcyan(str(e))
            colortext.warning(traceback.format_exc())
            if extract_data_for_case_if_missing:
                self.extract_data_for_case(prediction_id, root_directory = root_directory, force = True, score_method_id = score_method_id)
            try:
                top_x_ddg = self.get_top_x_ddg(prediction_id, score_method_id, top_x = top_x, expectn = expectn)
            except PartialDataException, e:
                raise
            except Exception, e:
                raise
        top_x_ddg_stability = self.get_top_x_ddg_stability(prediction_id, score_method_id, top_x = top_x, expectn = expectn)

        prediction_data[main_ddg_analysis_type] = top_x_ddg
        prediction_data['DDGStability_Top%d' % top_x] = top_x_ddg_stability
        return prediction_data


    @analysis_api
    def get_top_x_ddg(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''Returns the TopX value for the prediction. Typically, this is the mean value of the top X predictions for a
           case computed using the associated Score records in the database.'''

        # scores is a mapping from nstruct -> ScoreType -> score record where ScoreType is one of 'DDG', 'WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex'
        # if we do the calculation in Python, pull scores out to the top level first
        # otherwise, we can add a stored procedure to determine the TopX
        # if we go the Python route, we can implement different variations on TopX (including a stored procedure) and pass the function pointers as an argument to the main analysis function

        # Make sure that we have as many cases as we expect
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        # todo: move get_top_x_ddg_total_score into this function
        return self.get_top_x_ddg_total_score(scores, top_x)


    def get_top_x_ddg_total_score(self, scores, top_x):
        if scores == None:
            return None

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


    @analysis_api
    def get_top_x_ddg_stability(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''Returns the TopX value for the prediction only considering the complex scores. This computation may work as a
           measure of a stability DDG value.'''
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
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
            take_lowest = 3,
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
    def analyze(self, prediction_set_ids, score_method_id,
            experimental_data_exists = True,
            analysis_set_ids = [],
            prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_credits = {}, prediction_set_colors = {}, prediction_set_alphas = {},
            use_published_data = False,
            allow_failures = False,
            use_existing_benchmark_data = True, recreate_graphs = False,
            include_derived_mutations = False,
            expectn = 50,
            use_single_reported_value = False,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            output_directory = None,
            generate_plots = True,
            generate_matplotlib_plots = False,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            debug = False,
            restrict_to = set(),
            remove_cases = set(),
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
           burial_cutoff defines what should be considered buried (DSSPExposure field). Values around 1.0 are fully exposed, values of 0.0 are fully buried. For technical reasons, the DSSP value can exceed 1.0 but usually not by much.
           stability_classication_experimental_cutoff AKA x_cutoff. This defines the neutral mutation range for experimental values in kcal/mol i.e. values between -1.0 and 1.0 kcal/mol are considered neutral by default.
           stability_classication_predicted_cutoff AKA y_cutoff. This defines the neutral mutation range for predicted values in energy units.

           * Reporting arguments *

           output_directory : The directory in which to save plots and reports.
           generate_plots   : if plots are not needed, setting this to False can shorten the analysis time.
           report_analysis  : Whether or not to print analysis to stdout.
           silent = False   : Whether or not anything should be printed to stdout (True is useful for webserver interaction).
        '''
        assert(take_lowest > 0 and (int(take_lowest) == take_lowest))
        assert(0 <= burial_cutoff <= 2.0)
        assert(stability_classication_experimental_cutoff > 0)
        assert(stability_classication_predicted_cutoff > 0)
        assert(expectn > 0 and (int(expectn) == expectn))

        for prediction_set_id in prediction_set_ids:
            benchmark_run = self.get_analysis_dataframe(prediction_set_id,
                experimental_data_exists = experimental_data_exists,
                prediction_set_series_name = prediction_set_series_names.get(prediction_set_id),
                prediction_set_description = prediction_set_descriptions.get(prediction_set_id),
                prediction_set_credit = prediction_set_credits.get(prediction_set_id),
                prediction_set_color = prediction_set_colors.get(prediction_set_id),
                prediction_set_alpha = prediction_set_alphas.get(prediction_set_id),
                use_existing_benchmark_data = use_existing_benchmark_data,
                include_derived_mutations = include_derived_mutations,
                use_single_reported_value = use_single_reported_value,
                take_lowest = take_lowest,
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
            analysis_sets_to_run = sorted(analysis_sets_to_run)
            if experimental_data_exists:
                #todo: hack. this currently seems to expect all datapoints to be present. handle the case when we are missing data e.g. prediction set "ZEMu run 1"
                analysis_sets_to_run = ['ZEMu'] # ['BeAtMuSiC', 'SKEMPI', 'ZEMu']

            for analysis_set_id in analysis_sets_to_run:
                colortext.message(analysis_set_id)

                benchmark_run.calculate_metrics(analysis_set = analysis_set_id, analysis_directory = output_directory)
                benchmark_run.write_dataframe_to_csv(os.path.join(output_directory, 'data.csv'))
                benchmark_run.plot(analysis_set = analysis_set_id, analysis_directory = output_directory, matplotlib_plots = generate_matplotlib_plots)

                self.output_score_method_information(
                    score_method_id, output_directory,
                    analysis_set_id = analysis_set_id,
                    take_lowest = take_lowest,
                    expectn = expectn,
                )

                return benchmark_run
                # recreate_graphs
                # analysis_directory = output_directory

                # colors, alpha, and default series name and descriptions are taken from PredictionSet records
                # The order (if p1 before p2 then p1 will be on the X-axis in comparative plots) in comparative analysis plots is determined by the order in PredictionSets
                # assert PredictionSet for PredictionSet in PredictionSets is in the database

                # calls get_analysis_dataframe(options) over all PredictionSets
                # if output_directory is set, save files
                # think about how to handle this in-memory. Maybe return a dict like:
                    #"run_analyis" -> benchmark_name -> {analysis_type -> object}
                    #"comparative_analysis" -> (benchmark_name_1, benchmark_name_2) -> {analysis_type -> object}
                # comparative analysis
                #   only compare dataframes with the exact same points
                #   allow cutoffs, take_lowest to differ but report if they do so



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
    def get_complex_ids_matching_protein_name(self, partial_name):
        '''Returns a list of PPComplex IDs where at least one of the partner names matches partial_name.'''
        partial_name = '%{0}%'.format(partial_name)
        results = self.DDG_db.execute_select('''
            SELECT ID FROM PPComplex
            WHERE
            LName LIKE %s
            OR LShortName LIKE %s
            OR LHTMLName LIKE %s
            OR RName LIKE %s
            OR RShortName LIKE %s
            OR RHTMLName LIKE %s ORDER BY ID''', parameters=(partial_name, partial_name, partial_name, partial_name, partial_name, partial_name))
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
    def find_complex(self, pdbs, keywords = []):
        for pdb_id in pdb_ids:
            existing_records = self.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
            if existing_records:
                colortext.warning('The PDB file {0} exists in the database.'.format(pdb_id))
            complex_ids = self.search_complexes_by_pdb_id(pdb_id)
            if complex_ids:
                colortext.warning('The PDB file {0} has associated complexes: {1}'.format(pdb_id, ', '.join(map(str, complex_ids))))
                assert(len(complex_ids) == 1)
                complex_id = complex_ids[0]
                colortext.warning('Complex #{0}'.format(complex_id))
                pprint.pprint(self.get_complex_details(complex_id))

            assert(type(keywords) == list)
            keywords = set(keywords)

            ids = []
            for keyword in keywords:
                ids.extend(self.get_complex_ids_matching_protein_name(keyword))
            ids = sorted(set(ids))

            for id in ids:
                d = self.get_complex_details(id)
                colortext.warning(id)
                print('{0}, {1}, {2}'.format(d['LName'].encode('utf-8').strip(), d['LShortName'].encode('utf-8').strip(), d['LHTMLName'].encode('utf-8').strip()))
                print('{0}, {1}, {2}'.format(d['RName'].encode('utf-8').strip(), d['RShortName'].encode('utf-8').strip(), d['RHTMLName'].encode('utf-8').strip()))


    @ppi_data_entry
    def add_complex(self, complex_definition, keywords = [], force = False, allow_missing_params_files = False, debug = True, minimum_sequence_identity = 95.0):
        '''Add a complex to the database using a defined dict structure.

        :param complex_definition: A dict fitting the defined structure (see below).
        :param keywords: A list of keywords used to search existing complexes for an existing match. Not necessary but
                         advised, particularly when adding a small number of complexes.
        :param force: If a potentially similar complex is found and force is False then then the function returns with a
                      message and without adding the complex. The ForceAddition setting in the Complex dict (see below)
                      will have the same effect as setting this variable.
        :param allow_missing_params_files: If the PDB file contains ligands with no associated params files in complex_definition
                                           then the function will return with a message. If this option is set to True
                                           then the file will be added but the ligand is not guaranteed to be kept
                                           by the protocol (this depends on the version of Rosetta - versions from February
                                           2016 onwards should have better ligand support).
        :param minimum_sequence_identity: For non-RCSB files, we require a chain mapping from chains to the corresponding
                                          RCSB chains. Clustal is run to ensure that the sequence identity is at least this
                                          value.
        :param debug: If debug is set to True then the transaction used to insert the complex into the database will be
                      rolled back and a message stating that the insertion would have been successful is returned in the
                      return dict.
        :return: On successful import, the dict {success = True, ComplexID -> Long, SetNumber -> Long} corresponding to
                 the database PPIPDBSet primary key is returned. If a similar complex is detected and force is False then
                 a dict {success = False, debug = debug, message -> String} will be returned instead. On error, a dict
                 {success = False, error -> String} is returned.

        The database uses Unicode to encode the strings, allowing us to use e.g. Greek characters
        For this reason, please contain all structure definitions in a file encoded as Unicode. On Linux, you can add the
            # -*- coding: utf-8 -*-
        declaration at the top of the file (with no leading whitespace).

        One example of the dict structure is as follows:

            dict(
                Structure = dict(

                    # Fields required for all complexes. We require non-RCSB structures to be based on RCSB structures at least in this import interface.
                    rcsb_id = '1K5D',

                    # Exactly one of these fields must be specified for non-RCSB complexes.
                    pdb_object = p,
                    file_path = 'pdbs/1K5D2.pdb',

                    # Fields required for non-RCSB complexes.
                    db_id = '1K5D_TP0', # must be between 5-10 characters
                    techniques = "PDB_REDO",
                    file_source = "Roger Wilco",
                    user_id = "sierra",

                    # Required for non-RCSB complexes and optional for RCSB complexes
                    description = 'GSP1 complex (RanGAP1) from Tina Perica. This file was taken from PDB_REDO. Removed residues 180-213 from chain A (RAN/GSP1) (the ones wrapping around YRB1).',

                    # Exactly one of these fields must be specified for non-RCSB complexes.
                    # If identical_chains is passed then it must be set to True and all protein, DNA, and RNA chains in
                    # the input file must have the same chain ID as in the RCSB file.
                    # If chain_mapping is instead passed then the keys must cover all protein, DNA, and RNA chains in the input file.
                    # Since ligands are sometimes reassigned to a new chain e.g. 'X', we do not require a mapping for them in
                    # chain_mapping and we do not consider them in the sanity check for the identical_chains case.
                    # Note: This logic does not currently handle cases where multiple chains in the RCSB file were merged into one
                    #       chain in the input file i.e. chain A in the input file corresponds to chains L and H in the RCSB
                    #       structure.
                    # Clustal is run over the chain mapping to ensure a sequence identity of minimum_sequence_identity% (95% by default).
                    identical_chains = True
                        # or
                    chain_mapping = dict(
                        A = 'C', # To help other uses read the input file, it is useful to mention any other choices available i.e. "choice of A, D, G, J." followed by the chain name e.g. "RAN"
                        C = 'A', # e.g. "choice of C, F, I, L. Ran GTPase activating protein 1"
                    ),

                    # Fields required for non-RCSB complexes if there are ligands/ions in the input structure.
                    # First, all ligands and ions will be found in the input structure. Each (non-water) ligand and ion code must be
                    # covered in either the ligand_instance_mapping, the ligand_code_mapping, or the unchanged_ligand_codes mapping
                    # In practice, you will probably one want to use ligand_code_mapping, unchanged_ligand_codes, and unchanged_ion_codes
                    # as it is rare for ion codes to be changed.
                    # LigandMap is defined in klab.bio.ligand.
                    # The three types of mapping are:

                    # 1. detailed mappings - one mapping per instance. You will probably not want to bother using this format as it can become unwieldy with large numbers of ligands
                    ligand_instance_mapping = LigandMap.from_tuples_dict({ # Input PDB's HET code, residue ID -> RCSB HET code, RCSB residue ID
                        ('G13', 'X   1 ') : ('GNP', 'A1250 '),
                    }),
                    ion_instance_mapping = LigandMap.from_tuples_dict({ # Input PDB's HET code, residue ID -> RCSB HET code, RCSB residue ID
                        ('UNX', 'A1249 ') : ('MG', 'A1250 '),
                    }),

                    # 2. ligand type mappings - one mapping per ligand code. This is the easiest to specify when there are changes to ligand codes
                    ligand_code_mapping = {
                        'G13' : 'GNP',
                    )
                    ion_code_mapping = {
                        'UNX' : 'MG',
                    )

                    # 3. white lists - lists of ligand and ion codes which have not changed. Most of your cases will probably fall into this category
                    unchanged_ligand_codes = ['MSE'],
                    unchanged_ion_codes = ['FE2'],

                    # Fields currently advised when there are ligands in the structure in order for the protocols to handle
                    # ligands properly. This may be a non-issue in the future; at the time of writing, February 5th 2016,
                    # there was a large push by the Rosetta XRW to address the problem of handling arbitrary ligands so this
                    # problem may be mainly solved.
                    params_file_paths = {
                        'G13' : 'temp/pdbs/1K5D2.params'
                    },
                ),
                Complex = dict(
                    # There are two cases - the complex exists in the database or we will be adding a new complex.
                    # Note: Before adding new complexes, you should make sure that there is no existing complex in the
                    #       database. This will help to reduce redundancy and provide us with better data.

                    # These fields are required in both cases and specify the partners of the complex
                    # Note: Please ensure that the LChains (resp. RChains) chains correspond to the protein/complex
                    # identified by LName, LShortName, LHTMLName (resp. RName, RShortName, RHTMLName)
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
                )
            )
        '''

        # complex_definition, keywords = [], force = False, allow_missing_params_files = False, debug = True

        structural_details = complex_definition['Structure']
        complex_details = complex_definition['Complex']

        tsession = self.importer.get_session(new_session = True)
        rcsb_id = structural_details['rcsb_id']
        params_file_paths = structural_details.get('params_file_paths', {})
        for k, v in params_file_paths.iteritems():
            params_file_paths[k] = os.path.abspath(v)

        # Step 1: Add the associated PDB files
        if 'file_path' or 'pdb_object' in structural_details:


            from klab.bio.ligand import LigandMap


            pdb2pdb_chain_maps = []

            # Required fields
            database_id = structural_details['db_id']
            techniques = structural_details['techniques']
            assert(5 <= len(database_id) <= 10)
            if 'description' not in structural_details:
                 raise colortext.Exception('A description is required for non-RCSB files. This should include any details on the structure preparation.')
            if 'file_source' not in structural_details:
                 raise colortext.Exception('''A file_source is required for non-RCSB files. This should describe where the file came from e.g. "PDB REDO", "Rosetta", or the creator's name e.g. "Jon Snow".''')
            if 'user_id' not in structural_details:
                 raise colortext.Exception('A user_id is required for non-RCSB files. This should correspond to a record in the User table.')

            # Read the input PDB
            designed_pdb_object = None
            if structural_details.get('pdb_object'):
                if structural_details.get('file_path'):
                    raise colortext.Exception('Only one of pdb_object and file_path should be specified, not both.')
                designed_pdb_object = structural_details.get('pdb_object')
            else:
                if not structural_details.get('file_path'):
                    raise colortext.Exception('Exactly one of pdb_object or file_path must be specified.')
                pdb_file_path = os.path.abspath(structural_details['file_path'])
                if not os.path.exists(pdb_file_path):
                    raise colortext.Exception('Could not locate the file "{0}".'.format(pdb_file_path))
                designed_pdb_object = PDB(read_file(pdb_file_path))

            # Read the chain IDs
            all_chain_ids, main_chain_ids = [], []
            for k, v in sorted(designed_pdb_object.chain_types.iteritems()):
                all_chain_ids.append(k)
                if v != 'Unknown' and v != 'Solution' and v != 'Ligand':
                    main_chain_ids.append(k)

            # Check the chain mapping
            chain_mapping, chain_mapping_keys = {}, set()
            if 'identical_chains' in structural_details:
                # Create a partial chain mapping (ignoring ligand chains which may have a corresponding match)
                assert('chain_mapping' not in structural_details)
                assert(structural_details['identical_chains'] == True)
                for c in main_chain_ids:
                    chain_mapping[c] = c
                    chain_mapping_keys.add(c)
            else:
                # Check that the chain mapping domain is a complete mapping over the main chain IDs ( main_chain_ids <= chain_mapping.keys() <= all_chain_ids where "<=" is non-strict subset of)
                assert('chain_mapping' in structural_details)
                chain_mapping = structural_details['chain_mapping']
                chain_mapping_keys = set(chain_mapping.keys())
                assert(len(chain_mapping_keys.intersection(set(main_chain_ids))) == len(main_chain_ids)) # main_chain_ids is a subset of chain_mapping_keys
                assert(len(set(all_chain_ids).intersection(chain_mapping_keys)) == len(chain_mapping_keys)) # chain_mapping_keys is a subset of all_chain_ids
            unmapped_chains = sorted(set(all_chain_ids).difference(chain_mapping_keys))

            # Add the RCSB structure to the database and then use that object to run Clustal
            self.importer.add_pdb_from_rcsb(rcsb_id, trust_database_content = True, ligand_params_file_paths = {})
            rcsb_object = self.importer.get_pdb_object(rcsb_id, tsession = None)
            design_sequences = None
            if designed_pdb_object.seqres_sequences:
                design_sequences = designed_pdb_object.seqres_sequences
            else:
                design_sequences = designed_pdb_object.atom_sequences
            pcsa = PDBChainSequenceAligner()
            for chain_id, sequence in sorted(rcsb_object.seqres_sequences.iteritems()):
                pcsa.add(rcsb_id, chain_id, str(sequence))
            for chain_id, sequence in sorted(design_sequences.iteritems()):
                pcsa.add(database_id, chain_id, str(sequence))
            output, best_matches = pcsa.align()
            colortext.warning(pprint.pformat(best_matches))
            for dc, rc in sorted(chain_mapping.iteritems()):
                # todo: this mapping is untested on ligand chains. We could do: for all c in unmapped_chains, see if the 3-letter sequences match and, if so, and add those
                k1 = '{0}_{1}'.format(database_id, dc)
                k2 = '{0}_{1}'.format(rcsb_id, rc)
                print(k1, k2)
                assert(k1 in best_matches and k2 in best_matches[k1])
                sequence_identity = best_matches[k1][k2]
                pdb2pdb_chain_maps.append(dict(
                    PDBFileID1 = database_id,
                    Chain1 = dc,
                    PDBFileID2 = rcsb_id,
                    Chain2 = rc,
                    SequenceIdentity = sequence_identity,
                ))
                if sequence_identity < minimum_sequence_identity:
                    raise colortext.Exception('The chain mapping is defined as {0}, chain {1} -> {2}, chain {3} but the sequence identity for these chains is only {4}%.'.format(database_id, dc, rcsb_id, rc, sequence_identity))

            # Make sure that no params files were specified that do not fit the PDB file
            ligands_not_present = set(params_file_paths.keys()).difference(designed_pdb_object.get_ligand_codes())
            if len(ligands_not_present) > 0:
                raise colortext.Exception('Params files were specified for ligands which do not exist in the PDB file: "{0}".'.format('", "'.join(ligands_not_present)))

            # Sanity-check the ligand mapping for consistency
            ligand_code_map = {}
            if 'ligand_instance_mapping' in structural_details:
                pprint.pprint(structural_details['ligand_instance_mapping'].code_map)
                for k, v in structural_details['ligand_instance_mapping'].code_map.iteritems():
                    assert(ligand_code_map.get(k) == None or ligand_code_map[k] == v)
                    ligand_code_map[k] = v
            if 'ligand_code_mapping' in structural_details:
                pprint.pprint(structural_details['ligand_code_mapping'])
                for k, v in structural_details['ligand_code_mapping'].iteritems():
                    assert(ligand_code_map.get(k) == None or ligand_code_map[k] == v)
                    ligand_code_map[k] = v
            if 'unchanged_ligand_codes' in structural_details:
                for k in structural_details['unchanged_ligand_codes']:
                    assert(ligand_code_map.get(k) == None or ligand_code_map[k] == k)
                    ligand_code_map[k] = k
            pprint.pprint(ligand_code_map)

            # Sanity-check the ligand mapping for completeness
            ligand_mapping = LigandMap.from_code_map(ligand_code_map)
            print(sorted(ligand_code_map.keys()), sorted(designed_pdb_object.get_ligand_codes()))
            assert(sorted(ligand_code_map.keys()) == sorted(designed_pdb_object.get_ligand_codes()))

            # todo: currently unhandled
            assert('ion_instance_mapping' not in structural_details and 'ion_code_mapping' not in structural_details)

            # todo: update this function call
            importer.add_designed_pdb(designed_pdb_object, database_id, rcsb_id,
                                  structural_details['file_source'], structural_details['description'] , structural_details['user_id'],
                                  chain_mapping = chain_mapping, ligand_mapping = ligand_mapping,
                                  ligand_params_file_paths = structural_details.get('params_file_paths', {}), techniques=techniques)

            for pdb2pdb_chain_map in pdb2pdb_chain_maps:
                pdb2pdb_chain_map = get_or_create_in_transaction(tsession, dbmodel.PDB2PDBChainMap, pdb2pdb_chain_map)

        else:
            # RCSB file
            raise colortext.Exception('This should work but has not been tested')
            # Read the RCSB PDB
            rcsb_pdb_object = self.importer.get_pdb_object(rcsb_id, tsession = None)

            # Params files
            print(rcsb_pdb_object.get_ligand_codes())
            assert(len(set(params_file_paths.keys()).difference(rcsb_pdb_object.get_ligand_codes())) == 0)

            self.importer.add_pdb_from_rcsb(rcsb_id, trust_database_content = True, ligand_params_file_paths = params_file_paths)


        colortext.message('Importing {0} as {1}'.format(tina_pdb_id, tina_db_id))

        '''

                    # Fields required for non-RCSB complexes if there are ligands/ions in the input structure.
                    # First, all ligands and ions will be found in the input structure. Each (non-water) ligand and ion code must be
                    # covered in either the ligand_instance_mapping, the ligand_code_mapping, or the unchanged_ligand_codes mapping
                    # In practice, you will probably one want to use ligand_code_mapping, unchanged_ligand_codes, and unchanged_ion_codes
                    # as it is rare for ion codes to be changed.
                    # LigandMap is defined in klab.bio.ligand.
                    # The three types of mapping are:

                    # 1. detailed mappings - one mapping per instance. You will probably not want to bother using this format as it can become unwieldy with large numbers of ligands
                    ligand_instance_mapping = LigandMap.from_tuples_dict({ # Input PDB's HET code, residue ID -> RCSB HET code, RCSB residue ID
                        ('G13', 'X   1 ') : ('GNP', 'A1250 '),
                    }),
                    ion_instance_mapping = LigandMap.from_tuples_dict({ # Input PDB's HET code, residue ID -> RCSB HET code, RCSB residue ID
                        ('UNX', 'A1249 ') : ('MG', 'A1250 '),
                    }),

                    # 2. ligand type mappings - one mapping per ligand code. This is the easiest to specify when there are changes to ligand codes
                    ligand_code_mapping = {
                        'G13' : 'GNP',
                    )
                    ion_code_mapping = {
                        'UNX' : 'MG',
                    )

                    # 3. white lists - lists of ligand and ion codes which have not changed. Most of your cases will probably fall into this category
                    unchanged_ligand_codes = ['MSE'],
                    unchanged_ion_codes = ['FE2'],

                    # Fields currently advised when there are ligands in the structure in order for the protocols to handle
                    # ligands properly. This may be a non-issue in the future; at the time of writing, February 5th 2016,
                    # there was a large push by the Rosetta XRW to address the problem of handling arbitrary ligands so this
                    # problem may be mainly solved.
                    params_files = {
                        'G13' : 'temp/pdbs/1K5D2.params'
                    },
                ),
        '''




        assert(sorted(details['Structure']['chain_mapping'].keys()) == sorted(details['Complex']['LChains'] + details['Complex']['RChains']))

        assert(details['Structure']['unchanged_ligand_codes'] or details['Structure']['ligand_mapping'])

        for k in details['Structure']['params_files'].keys():
            print(k)
            assert(k in details['Structure']['ligand_mapping'.code_map])


        assert ComplexID or ComplexDetails - LChains and RChains


        tina_pdb_id = tina_pdb_id.upper()
        rcsb_pdb_id = tina_pdb_id_to_rcsb_pdb_id[tina_pdb_id]
        assert(structural_details['rcsb_id'] == rcsb_pdb_id)

        tina_pdb_object = tina_pdb_objects[tina_pdb_id]
        rcsb_pdb_object = rcsb_pdb_objects[rcsb_pdb_id]
        tina_db_id = structural_details['db_id']

        assert((tina_db_id != tina_pdb_id) and (tina_db_id != rcsb_pdb_id) and (len(tina_db_id) > 7) and (tina_db_id[4:7] == '_TP'))

        #    Step 1: Add complexes
        complex_id = None
        if 'complex_id' in complex_details:
           complex_id = complex_details['complex_id']
        else:


            Complex = dict(
                ComplexID = 202,
                LChains = ['A'],
                RChains = ['B'],
            )

        ppi_api.add_complex

        mut_complex_3H7P = dict(
            LName = 'Ubiquitin (yeast) K63R',
            LShortName = 'Ubiquitin K63R',
            LHTMLName = 'Ubiquitin (yeast) K63R',
            RName = 'Ubiquitin (yeast)',
            RShortName = 'Ubiquitin',
            RHTMLName = 'Ubiquitin (yeast)',
            FunctionalClassID = 'OX',
            PPDBMFunctionalClassID = 'O',
            PPDBMDifficulty = None,
            IsWildType = False,
            WildTypeComplexID = wt_complex_id,
            Notes = None,
            Warnings = None,
        )
        DDGdb.insertDictIfNew('PPComplex', mut_complex_3H7P, ['LName', 'RName'])







        #complex_definition, keywords = [], force = False, allow_missing_params_files = False)


        self.find_complex(pdbs, keywords)

        #pdb_ids

        #complex = complex_definition['Complex']
        #Keywords = [complex['LName'], complex['LShortName'], complex['LHTMLName'], complex['RName'], complex['RShortName'], complex['RHTMLName']]


        mut_complex_3H7P = dict(
            LName = 'Ubiquitin (yeast) K63R',
            LShortName = 'Ubiquitin K63R',
            LHTMLName = 'Ubiquitin (yeast) K63R',
            RName = 'Ubiquitin (yeast)',
            RShortName = 'Ubiquitin',
            RHTMLName = 'Ubiquitin (yeast)',
            FunctionalClassID = 'OX',
            PPDBMFunctionalClassID = 'O',
            PPDBMDifficulty = None,
            IsWildType = False,
            WildTypeComplexID = wt_complex_id,
            Notes = None,
            Warnings = None,
        )
        DDGdb.insertDictIfNew('PPComplex', mut_complex_3H7P, ['LName', 'RName'])

        '''

    for tina_pdb_id, details in sorted(complex_definitions.iteritems()):

        structural_details = details['Structure']
        complex_details = details['Complex']

        if tina_pdb_id in ['1A2K']:
            continue

        if tina_pdb_id != '1I2M':
            assert(details['Structure']['chain_mapping'])

        assert(sorted(details['Structure']['chain_mapping'].keys()) == sorted(details['Complex']['LChains'] + details['Complex']['RChains']))

        assert(details['Structure']['unchanged_ligand_codes'] or details['Structure']['ligand_mapping'])

        for k in details['Structure']['params_files'].keys():
            print(k)
            assert(k in details['Structure']['ligand_mapping'.code_map])


        assert ComplexID or ComplexDetails - LChains and RChains


        tina_pdb_id = tina_pdb_id.upper()
        rcsb_pdb_id = tina_pdb_id_to_rcsb_pdb_id[tina_pdb_id]
        assert(structural_details['rcsb_id'] == rcsb_pdb_id)

        tina_pdb_object = tina_pdb_objects[tina_pdb_id]
        rcsb_pdb_object = rcsb_pdb_objects[rcsb_pdb_id]
        tina_db_id = structural_details['db_id']

        assert((tina_db_id != tina_pdb_id) and (tina_db_id != rcsb_pdb_id) and (len(tina_db_id) > 7) and (tina_db_id[4:7] == '_TP'))

        #    Step 1: Add complexes
        complex_id = None
        if 'complex_id' in complex_details:
           complex_id = complex_details['complex_id']
        else:


        Complex = dict(
            ComplexID = 202,
            LChains = ['A'],
            RChains = ['B'],
        )

        ppi_api.add_complex

        mut_complex_3H7P = dict(
            LName = 'Ubiquitin (yeast) K63R',
            LShortName = 'Ubiquitin K63R',
            LHTMLName = 'Ubiquitin (yeast) K63R',
            RName = 'Ubiquitin (yeast)',
            RShortName = 'Ubiquitin',
            RHTMLName = 'Ubiquitin (yeast)',
            FunctionalClassID = 'OX',
            PPDBMFunctionalClassID = 'O',
            PPDBMDifficulty = None,
            IsWildType = False,
            WildTypeComplexID = wt_complex_id,
            Notes = None,
            Warnings = None,
        )
        DDGdb.insertDictIfNew('PPComplex', mut_complex_3H7P, ['LName', 'RName'])'''
