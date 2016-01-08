#!/usr/bin/python2.4
# encoding: utf-8
"""
import_api.py
High-level functions for importing data into the DDG database.

Created by Shane O'Connor 2015.
Copyright (c) 2015 Shane O'Connor. All rights reserved.
"""

import sys
import pprint
from io import BytesIO
import os
import copy
import json
import zipfile
import traceback
import gzip
import shutil
import sqlite3
import cPickle as pickle
from types import NoneType

import numpy
import pandas

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import create_engine
from sqlalchemy import inspect as sqlalchemy_inspect

if __name__ == '__main__':
    sys.path.insert(0, '../../klab')

from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation
from klab.fs.fsio import read_file, write_temp_file, open_temp_file, write_file
from klab.bio.pfam import Pfam
from klab.bio.dssp import MonomerDSSP, ComplexDSSP, MissingAtomException
from klab.bio.ligand import Ligand, PDBLigand

from db_schema import test_schema_against_database_instance
from db_schema import PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue, LigandDescriptor, LigandIdentifier, LigandSynonym, PDBLigand
from db_schema import Ligand as DBLigand
from api_layers import *
from db_api import ddG, PartialDataException, SanityCheckException
import ddgdbapi


class DataImportInterface(object):
    '''This is the data import API class which should be used when adding basic data (PDB files, complex definitions, etc.)
       to the database.

            from ddglib.interface_api import DataImportInterface
            importer = DataImportInterface(read_file('ddgdb.pw'))
            e.g.
            importer.add_pdb_from_rcsb('1A2K')

       Objects of this class and derived subclasses has three main members:

          self.DDG_db - a database interface used to interact directly with the database via MySQL commands
          self.DDG_db_utf - the same interface but with UTF support. This should be used when dealing with UTF fields e.g. publication data
          self.prediction_data_path - this is the location on the file server where output form jobs of the derived class type (e.g. binding affinity jobs) should be stored.
    '''


    ##################
    #                #
    #  Constructors  #
    #                #
    ##################


    def __init__(self, passwd, connect_string, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):

        # Set up MySQLdb connections
        passwd = passwd.strip()
        self.DDG_db = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname)
        self.DDG_db_utf = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname, use_utf = True)
        test_schema_against_database_instance(self.DDG_db)

        # Set up SQLAlchemy connections
        self.connect_string = connect_string
        self.engine, self.session = None, None
        self.get_engine()
        self.get_session()

        self.rosetta_scripts_path = rosetta_scripts_path
        self.rosetta_database_path = rosetta_database_path

        # Parse PDB chain -> Pfam mapping
        self.pfam_api = Pfam()


    @classmethod
    def get_interface_with_config_file(cls, database = 'ddg', host_config_name = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None, my_cnf_path = None):
        # Uses ~/.my.cnf to get authentication information
        ### Example .my.cnf (host_config_name will equal guybrush2):
        ### [clientguybrush2]
        ### user=kyleb
        ### password=notmyrealpass
        ### host=guybrush.ucsf.edu
        if not my_cnf_path:
            my_cnf_path = os.path.expanduser(os.path.join('~', '.my.cnf'))
        if not os.path.isfile(os.path.expanduser(my_cnf_path)):
            raise Exception("A .my.cnf file must exist at: " + my_cnf_path)

        # These four variables must be set in a section of .my.cnf named host_config_name
        user = None
        password = None
        host = None
        connection_string = None
        connection_string_key = 'sqlalchemy.{0}.url'.format(database)
        with open(my_cnf_path, 'r') as f:
            parsing_config_section = False
            for line in f:
                if line.strip() == '[client%s]' % host_config_name:
                    parsing_config_section = True
                elif line.strip() == '':
                    parsing_config_section = False
                elif parsing_config_section:
                    if '=' in line:
                        key, val = line.strip().split('=')
                        key, val = key.strip(), val.strip()
                        if key == 'user':
                            user = val
                        elif key == 'password':
                            password = val
                        elif key == 'host':
                            host = val
                        elif key == connection_string_key:
                            connection_string = val
                    else:
                        parsing_config_section = False

        if not user or not password or not host or not connection_string:
            raise Exception("Couldn't find host(%s), username(%s), password, or connection string in section %s in %s" % (host, user, host_config_name, my_cnf_path) )

        return cls(password, connection_string, username = user, hostname = host, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)


    def __del__(self):
        pass         #self.DDG_db.close()         #self.ddGDataDB.close()



    #############################
    #                           #
    #  SQLAlchemy Engine setup  #
    #                           #
    #############################



    def get_engine(self):
        if not self.engine:
            self.engine = create_engine(self.connect_string)
        return self.engine


    def get_connection(self):
        # e.g. connection = importer.get_connection(); connection.execute("SELECT * FROM User")
        self.get_engine()
        return self.engine.connect()


    def get_session(self, new_session = False, autoflush = True, autocommit = False):
        self.get_engine()
        if new_session or not(self.session):
            maker_ddgdatabase = sessionmaker(autoflush = autoflush, autocommit = autocommit)
            s = scoped_session(maker_ddgdatabase)
            DeclarativeBaseDDG = declarative_base()
            metadata_ddg = DeclarativeBaseDDG.metadata
            s.configure(bind=self.engine)
            metadata_ddg.bind = self.engine
            if new_session:
                return s
            else:
                self.session = s
        return self.session


    def renew(self):
        self.session = self.get_session(new_session = True)


    def get_or_create_in_transaction(self, tsession, model, values, missing_columns = []):
        '''
        Uses the SQLAlchemy model to retrieve an existing record based on the supplied field values or, if there is no
        existing record, to create a new database record.
        :param tsession: An SQLAlchemy transactioned session
        :param model: The name of the SQLAlchemy class representing the table
        :param values: A dict of values which will be used to populate the fields of the model
        :param missing_columns: Elements of missing_columns are expected to be fields in the model but are left blank regardless of whether they exist in values. This is useful for auto_increment fields.
        :return:
        '''

        values = copy.deepcopy(values) # todo: this does not seem to be necessary since we do not seem to be writing

        fieldnames = [c.name for c in list(sqlalchemy_inspect(model).columns)]
        for c in missing_columns:
            fieldnames.remove(c)

        pruned_values = {}
        for k in set(values.keys()).intersection(set(fieldnames)):
            v = values[k]
            pruned_values[k] = v
        assert(sorted(pruned_values.keys()) == sorted(fieldnames))

        instance = tsession.query(model).filter_by(**pruned_values)
        if instance.count() > 1:
            raise Exception('Multiple records were found with the search criteria.')
        instance = instance.first()

        if instance:
            return instance
        else:
            instance = model(**pruned_values)
            tsession.add(instance)
            tsession.flush()
            return instance


    #################################
    #                               #
    #  Ligand entry - public API    #
    #                               #
    #  Missing tables:              #
    #      LigandPrice              #
    #      LigandReference          #
    #                               #
    #################################


    def add_ligand_by_pdb_code(self, pdb_code):
        '''This function adds a ligand to the database using the ligand's PDB code. The insertion is handled by a transaction.
           Touched tables:
               Ligand
               LigandDescriptor
               LigandIdentifier
               LigandSynonym
        '''
        colortext.message('Adding ligand {0}'.format(pdb_code))
        l = Ligand.retrieve_data_from_rcsb(pdb_code, cached_dir = '/tmp')
        colortext.ppurple(l)

        tsession = self.get_session(new_session = True) # As this may be called by another function, we want to keep the ligand entry separate from other transactions.
        try:
            # Create the main ligand record
            if l.InChI == None:
                # Error handling for the unknown ligands
                assert(l.PDBCode == 'UNL' or l.PDBCode == 'UNK' or l.PDBCode == 'UNX')
                l.InChI = l.PDBCode
                l.InChIKey = l.PDBCode
            db_ligand = self.get_or_create_in_transaction(tsession, DBLigand, l.__dict__, missing_columns = ['ID'])

            # Create the ligand descriptor records
            descriptor_fieldnames = [c.name for c in list(sqlalchemy_inspect(LigandDescriptor).columns)]

            for descriptor in l.descriptors:
                descriptor = copy.deepcopy(descriptor)
                descriptor['LigandID'] = db_ligand.ID
                db_ligand_descriptor = self.get_or_create_in_transaction(tsession, LigandDescriptor, descriptor, missing_columns = ['ID'])

            for identifier in l.identifiers:
                identifier = copy.deepcopy(identifier)
                identifier['LigandID'] = db_ligand.ID
                db_ligand_identifier = self.get_or_create_in_transaction(tsession, LigandIdentifier, identifier, missing_columns = ['ID'])

            for synonym in l.synonyms:
                db_ligand_synonym = self.get_or_create_in_transaction(tsession, LigandSynonym, dict(LigandID = db_ligand.ID, Synonym = synonym.strip()))

            tsession.commit()
            print('Success.\n')
        except:
            colortext.error('Failure.')
            tsession.rollback()
            raise


    #################################
    #                               #
    #  PDB file entry - public API  #
    #                               #
    #  Missing tables:              #
    #      FileContent              #
    #      LigandPrice              #
    #      LigandReference          #
    #      PDBChain                 #
    #      PDBFile                  #
    #      PDBMolecule              #
    #      PDBMoleculeChain         #
    #      PDBResidue               #
    #                               #
    #      Protein                  #
    #      ProteinDatabaseIdentifier#
    #      ProteinName              #
    #      ProteinOrganism          #
    #      ProteinResidue           #
    #      ProteinSegment           #
    #      Publication              #
    #                               #
    #################################


    def add_pdb_from_rcsb(self, pdb_id):
        '''NOTE: This API is used to create and analysis predictions or retrieve information from the database.
                 This function adds new raw data to the database and does not seem to belong here. It should be moved into
                 an admin API instead.
           This function adds imports a PDB into the database, creating the associated molecule, chain and residue etc. records.'''

        assert(len(pdb_id) == 4)
        file_source = 'RCSB'

        contains_membrane_protein = None
        protein = None
        UniProtAC = None
        UniProtID = None
        testonly = False
        force = False
        techniques = None
        derived_from = None
        notes = None
        allow_missing_molecules = False


        # for all ligand codes in the PDB (if RCSB PDB):
        #    call importer.add_ligand_by_pdb_code(l)

        # for all ligand codes in the PDB (if not RCSB PDB), require a mapping to PDB ligand IDs and then:
        #    call importer.add_ligand_by_pdb_code(l)
        #    create a PDBLigand record with the Ligand.ID associated with the PDB ligand ID

        #-------
        p = PDB.retrieve(pdb_id)

        pdbtmo = PDBTM(read_file('/kortemmelab/shared/mirror/PDBTM/pdbtmall.xml'))
        pdb_id_map = pdbtmo.get_pdb_id_map()
        uc_pdb_id_map = {}
        for k, v in pdb_id_map.iteritems():
            uc_pdb_id_map[k.upper()] = v
        if pdb_id.upper() in uc_pdb_id_map:
            print('{0} is a transmembrane protein.'.format(pdb_id))
        return

        user_ids = set([str(p.upper()) for p in json.loads(kw['pdb_ids'])])
        common_pdb_ids = sorted(pdbtm_ids.intersection(user_ids))

        comment = '%d of the PDB IDs are marked as membrane proteins in the PDBTM' % len(pdbtm_ids.intersection(user_ids))
        #-------

        assert(file_source)
        if filepath:
            if not os.path.exists(filepath):
                raise Exception("The file %s does not exist." % filepath)
            filename = os.path.split(filepath)[-1]
            rootname, extension = os.path.splitext(filename)
            if not extension.lower() == ".pdb":
                raise Exception("Aborting: The file does not have a .pdb extension.")
        if pdb_id:
            rootname = pdb_id

        try:
            dbp = PDBFileImporter(self.DDG_db, rootname, contains_membrane_protein = contains_membrane_protein, protein = protein, file_source = file_source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID, testonly = testonly, techniques = techniques, derived_from = derived_from, notes = notes)
            #Structure.getPDBContents(self.DDG_db)
            results = self.DDG_db.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters = (rootname,))

            if results:
                #ddgdbapi.getUniProtMapping(pdb_id, storeInDatabase = True)
                #raise Exception("There is already a structure in the database with the ID %s." % rootname)
                if force:
                    dbp.commit(testonly = testonly, allow_missing_molecules = allow_missing_molecules)
                return None
            dbp.commit(testonly = testonly)
            return rootname
        except Exception, e:
            colortext.error(str(e))
            colortext.error(traceback.format_exc())
            raise Exception("An exception occurred committing %s to the database." % filepath)


    def _add_pdb_ligands(self, pdb_object, pdb_tsession):
        '''This function associates the ligands of a PDB file (which may be arbitrarily named) with ligands entered in
           the database using the ligand's PDB code. The insertion is handled by a transaction which should be set up
           by the caller.
           Touched tables:
               PDBLigand
        '''
        pass


    def _update_pdb_residues(self, check_chains = True):
        '''This function should not usually be run but it is handy to keep around.'''

        from ppi_api import get_interface as get_ppi_interface
        ppi_api = get_ppi_interface(read_file('../pw'))

        records = self.DDG_db.execute_select('SELECT DISTINCT PDBFileID FROM PDBResidue WHERE MonomericExposure IS NULL OR MonomericDSSP IS NULL OR ComplexExposure IS NULL OR ComplexDSSP IS NULL')
        records = self.DDG_db.execute_select('SELECT DISTINCT ID AS PDBFileID FROM PDBFile')
        pdb_ids = [r['PDBFileID'] for r in records]
        num_pdb_ids = len(pdb_ids)
        counter = 0

        pdb_objects = {}

        missing_chains = True
        for pdb_id in pdb_ids:

            counter += 1

            #426/1474: 1GPW PDB chains: A,B,C,D,E,F	 DB chains: A,B,C,D,E,F,P
            #if pdb_id != '1GTX':
            #    continue
            if counter < 434:
                continue

            colortext.message('{0}/{1}: {2}'.format(counter, num_pdb_ids, pdb_id))
            pdb_object = PDB(self.DDG_db.execute_select("SELECT Content FROM PDBFile WHERE ID=%s", parameters = (pdb_id,))[0]['Content'])
            pdb_objects[pdb_id] = pdb_object

            #pprint.pprint(pdb_object.chain_types)

            chains = sorted(set(pdb_object.seqres_sequences.keys() + pdb_object.atom_sequences.keys() + pdb_object.chain_types.keys()))
            if check_chains:
                db_chains = [r['Chain'] for r in self.DDG_db.execute_select("SELECT Chain FROM PDBChain WHERE PDBFileID=%s ORDER BY Chain", parameters = (pdb_id,))]
                if chains != db_chains:
                    colortext.warning('PDB chains: {0}\t DB chains: {1}'.format(','.join(chains), ','.join(db_chains)))
                    colortext.error('Missing chains.')
                    if pdb_id == '1GTX':
                        for c in pdb_object.chain_types.keys():
                            pdbc = dict(
                                PDBFileID = pdb_id,
                                Chain = c,
                                MoleculeType = pdb_object.chain_types[c],
                            )
                            self.DDG_db.insertDictIfNew('PDBChain', pdbc, ['PDBFileID', 'Chain'])
                    break
                    for c in chains:
                        # todo:
                        pdbc = dict(
                            PDBFileID = pdb_id,
                            Chain = c,
                            MoleculeType = pdb_object.chain_types[c],
                        )
                        self.ddGdb.insertDictIfNew('PDBChain', pdbc, ['PDBFileID', 'Chain'])

            for chain_id in pdb_object.chain_types.keys():
                self.DDG_db.execute("UPDATE PDBChain SET MoleculeType=%s WHERE PDBFileID=%s AND Chain=%s", parameters = (pdb_object.chain_types[chain_id], pdb_id, chain_id))

            if False:
                # Add PDBMolecule and PDBMoleculeChain records
                try:
                    molecules = pdb_object.get_molecules_and_source()
                except MissingRecordsException:
                    molecules = []
                for molecule in molecules:
                    chains = molecule['Chains']
                    molecule['PDBFileID'] = pdb_id
                    molecule['Organism'] = molecule['OrganismScientificName'] or molecule['OrganismCommonName']
                    md = {}
                    for k in ['PDBFileID', 'MoleculeID', 'Name', 'Organism', 'Fragment', 'Synonym', 'Engineered', 'EC', 'Mutation', 'OtherDetails']:
                        md[k] = molecule[k]
                    self.ddGdb.insertDictIfNew('PDBMolecule', md, ['PDBFileID', 'MoleculeID'])
                    for c in chains:
                        try:
                            mcd = dict(
                                PDBFileID = pdb_id,
                                MoleculeID = md['MoleculeID'],
                                Chain = c
                            )
                            self.ddGdb.insertDictIfNew('PDBMoleculeChain', mcd, ['PDBFileID', 'MoleculeID', 'Chain'])
                        except:
                            if allow_missing_molecules: pass
                            else: raise



            # Run DSSP
            dssp_complex_d, dssp_monomer_d = None, None
            try:
                # This fails for some PDB e.g. if they only have CA atoms
                dssp_complex_d = ComplexDSSP(pdb_object)
            except MissingAtomException, e:
                print('DSSP (complex) failed for this case.')
            try:
                # This fails for some PDB e.g. if they only have CA atoms
                dssp_monomer_d = MonomerDSSP(pdb_object)
            except MissingAtomException, e:
                print('DSSP (monomer) failed for this case.')

            # Make sure that the residue records exist for all results of DSSP
            if dssp_monomer_d:
                for chain_id in pdb_object.atom_sequences.keys():
                    if pdb_object.chain_types[chain_id] == 'Protein':
                        for chain_id, mapping in dssp_monomer_d:
                            for residue_id, residue_details in sorted(mapping.iteritems()):
                                residue_record = self.DDG_db.execute_select('SELECT ID FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                                assert(len(residue_record) == 1)

                # Add the monomeric DSSP results
                for chain_id in pdb_object.atom_sequences.keys():
                    if pdb_object.chain_types[chain_id] == 'Protein':
                        colortext.warning('\tRunning monomeric DSSP on chain %s' % chain_id)
                        for chain_id, mapping in dssp_monomer_d:
                            for residue_id, residue_details in sorted(mapping.iteritems()):
                                residue_record = self.DDG_db.execute_select('SELECT ID, MonomericExposure, MonomericDSSP FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                                assert(len(residue_record) == 1)
                                PDBResidueID = residue_record[0]['ID']
                                if residue_record[0]['MonomericDSSP'] == None:
                                    self.DDG_db.execute('UPDATE PDBResidue SET MonomericDSSP=%s WHERE ID=%s', parameters=(residue_details['ss'], PDBResidueID))
                                else:
                                    #print(pdb_id, chain_id, residue_id, residue_record[0]['MonomericDSSP'], residue_details['ss'])
                                    assert(residue_record[0]['MonomericDSSP'] == residue_details['ss'])
                                if residue_record[0]['MonomericExposure'] == None:
                                    self.DDG_db.execute('UPDATE PDBResidue SET MonomericExposure=%s WHERE ID=%s', parameters=(residue_details['exposure'], PDBResidueID))
                                else:
                                    #print(pdb_id, chain_id, residue_id, residue_record[0]['MonomericExposure'], residue_details['exposure'])
                                    #assert(residue_record[0]['MonomericExposure'] == residue_details['exposure'])
                                    assert(abs(residue_record[0]['MonomericExposure'] - residue_details['exposure']) < 0.0001)

            # Make sure that the residue records exist for all results of DSSP
            if dssp_complex_d:
                for chain_id in pdb_object.atom_sequences.keys():
                    if pdb_object.chain_types[chain_id] == 'Protein':
                        for chain_id, mapping in dssp_complex_d:
                            for residue_id, residue_details in sorted(mapping.iteritems()):
                                residue_record = self.DDG_db.execute_select('SELECT ID FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                                assert(len(residue_record) == 1)

                # Add the complex DSSP results
                for chain_id in pdb_object.atom_sequences.keys():
                    if pdb_object.chain_types[chain_id] == 'Protein':
                        colortext.warning('\tRunning complex DSSP on chain %s' % chain_id)
                        for chain_id, mapping in dssp_complex_d:
                            for residue_id, residue_details in sorted(mapping.iteritems()):
                                residue_record = self.DDG_db.execute_select('SELECT ID, ComplexExposure, ComplexDSSP FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                                assert(len(residue_record) == 1)
                                PDBResidueID = residue_record[0]['ID']
                                if residue_record[0]['ComplexDSSP'] == None:
                                    self.DDG_db.execute('UPDATE PDBResidue SET ComplexDSSP=%s WHERE ID=%s', parameters=(residue_details['ss'], PDBResidueID))
                                else:
                                    assert(residue_record[0]['ComplexDSSP'] == residue_details['ss'])
                                if residue_record[0]['ComplexExposure'] == None:
                                    self.DDG_db.execute('UPDATE PDBResidue SET ComplexExposure=%s WHERE ID=%s', parameters=(residue_details['exposure'], PDBResidueID))
                                else:
                                    #assert(residue_record[0]['ComplexExposure'] == residue_details['exposure'])
                                    assert(abs(residue_record[0]['ComplexExposure'] - residue_details['exposure']) < 0.0001)



            # Extract and store the coordinates
            for chain_id in pdb_object.atom_sequences.keys():
                chain_dataframe = pdb_object.extract_xyz_matrix_from_chain(chain_id)
                if isinstance(chain_dataframe, NoneType):
                    raise Exception('The coordinates dataframe could not be created for {0}, chain {1}'.format(pdb_id, chain_id))
                ufname, cfname = None, None
                if isinstance(ppi_api.get_pdb_chain_coordinates(pdb_id, chain_id), NoneType):
                    try:
                        f, ufname = open_temp_file('/tmp', suffix = '.hdf5')
                        f.close()
                        f, cfname = open_temp_file('/tmp', suffix = '.hdf5.gz')
                        f.close()
                        store = pandas.HDFStore(ufname)
                        store['dataframe'] = chain_dataframe
                        store.close()
                        content = read_file(ufname, binary = True)
                        with gzip.open(cfname, 'wb') as f:
                            f.write(content)
                        f = open(cfname)
                        zipped_contents = f.read()
                        f.close()
                        ppi_api.DDG_db.execute('UPDATE PDBChain SET Coordinates=%s WHERE PDBFileID=%s AND Chain=%s', parameters=(zipped_contents, pdb_id, chain_id))
                        os.remove(ufname)
                        os.remove(cfname)
                    except Exception, e:
                        print('Failed to add coordinates for {0}, chain {1}'.format(pdb_id, chain_id))
                        if ufname: os.remove(ufname)
                        if cfname: os.remove(cfname)
                        print(str(e))
                        print(traceback.format_exc())
                else:
                    print(pdb_id + chain_id + ' has coordinates')

            continue

            if True:
                # Add PDBResidue records
                for c, seq in pdb_object.atom_sequences.iteritems():
                    count = 1
                    for s in seq:
                        res_id, r = s
                        assert(len(r.ResidueID) == 5)
                        assert(c == r.Chain)

                        residue_record = self.DDG_db.execute_select('SELECT * FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, c, r.ResidueID))
                        if len(residue_record) != 1:
                            print('Missing {0}, chain {1}, {2}'.format(pdb_id, c, r.ResidueID))



                        return
                        continue
                        db_res = dict(
                            PDBFileID = pdb_id,
                            Chain = c,
                            ResidueID = r.ResidueID,
                            ResidueAA = r.ResidueAA,
                            ResidueType = r.residue_type,
                            IndexWithinChain = count,
                            CoordinatesExist = True,
                            RecognizedByRosetta = None,
                            BFactorMean = None,
                            BFactorDeviation = None,
                            SecondaryStructurePosition = None,
                            AccessibleSurfaceArea = None,
                            MonomericExposure = None,
                            MonomericDSSP = None,
                            ComplexExposure = None,
                            ComplexDSSP = None,
                        )
                        self.ddGdb.insertDictIfNew('PDBResidue', db_res, ['PDBFileID', 'Chain', 'ResidueID'])
                        count += 1

            if False:
                # Update PDBResidue records
                for c, seq in pdb_object.atom_sequences.iteritems():
                    count = 1
                    for s in seq:
                        res_id, r = s
                        assert(len(r.ResidueID) == 5)
                        assert(c == r.Chain)
                        db_res = dict(
                            PDBFileID = pdb_id,
                            Chain = c,
                            ResidueID = r.ResidueID,
                            ResidueAA = r.ResidueAA,
                            ResidueType = r.residue_type,
                            IndexWithinChain = count,
                            CoordinatesExist = True,
                            RecognizedByRosetta = None,
                            BFactorMean = None,
                            BFactorDeviation = None,
                            SecondaryStructurePosition = None,
                            AccessibleSurfaceArea = None,
                            MonomericExposure = None,
                            MonomericDSSP = None,
                            ComplexExposure = None,
                            ComplexDSSP = None,
                        )
                        self.ddGdb.insertDictIfNew('PDBResidue', db_res, ['PDBFileID', 'Chain', 'ResidueID'])
                        count += 1




        return



    def update_pdb_chains(self):
        '''
        '''

        pdb_files = self.DDG_db.execute_select('SELECT ID, FileSource FROM PDBFile ORDER BY ID')
        num_files = len(pdb_files)

        x = 0
        ligand_codes = set()
        colortext.message('Updating {0} PDB files.'.format(num_files))
        rescount = {}
        for pdb_file in pdb_files:
            x += 1
            #if True and x < 67:
            #    continue

            pdb_id = pdb_file['ID']

            #if pdb_id < '1DAN':
            #    continue

            #print('')
            #print(pdb_id)
            rescount[pdb_id] = self.DDG_db.execute_select('SELECT COUNT(ID) AS NumResidues FROM PDBResidue WHERE PDBFileID=%s AND ComplexDSSP IS NOT NULL', parameters = (pdb_id,))[0]['NumResidues']
            continue
            file_source = pdb_file['FileSource']
            if len(pdb_id) == 4:
                assert(file_source == 'RCSB')
            else:
                assert(file_source != 'RCSB')
            #print('{0}/{1}'.format(x, num_files))
            sys.stdout.write('{0}/{1}\r'.format(x, num_files)); sys.stdout.flush()

            pdb_id = pdb_file['ID']
            p = PDB(self.DDG_db.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters = (pdb_id,))[0]['Content'])

            print(p.hetatm_formulae)
            ligand_codes = ligand_codes.union(p.get_ligand_codes())
            pprint.pprint(p.ligands)
            print(ligand_codes)
            continue

            print('here')
            break

            chains = set(p.seqres_sequences.keys()).union(set(p.atom_sequences.keys()))
            for c in chains:
                self.DDG_db.execute_select('UPDATE PDBChain SET MoleculeType=%s WHERE PDBFileID=%s AND Chain=%s', parameters = (p.chain_types[c], pdb_id, c))
            #break
            print('.')
            #print(p.seqres_sequences)
            #print(p.atom_sequences)

        pprint.pprint(rescount)
        print(len(rescount))
        for k, v in rescount.iteritems():
            if v < 50:
                print(k, v)

        print(sorted(ligand_codes))
        #colortext.warning(self.pfam_api.get_pfam_accession_numbers_from_pdb_chain('1TVA', 'A'))


    ##################
    # Adding PDB files
    ##################





def test():
    importer = DataImportInterface.get_interface_with_config_file()
    session = importer.session

    if True:
        ligand_codes = ['0Z6', '0G6', '0EF', '0QE', '13P', '15P', '1AC', '1PE', '1PG', '2AB', '2GP', '3AA', '3GP', '3PG', '544', '5IU', '9CR', 'A80', 'ABA', 'AC9', 'ACA', 'ACE', 'ACH', 'ACT', 'ACY', 'ADC', 'ADP', 'ADU', 'AE3', 'AF3', 'AGF', 'AGL', 'AGP', 'AI2', 'ALA', 'ALC', 'ALF', 'AMB', 'AMP', 'ANL', 'ANP', 'ANS', 'AOA', 'AP5', 'APB', 'AR7', 'ARA', 'ARB', 'ASP', 'ATP', 'AXP', 'B3I', 'B3P', 'B98', 'BCT', 'BEN', 'BEO', 'BGC', 'BHD', 'BID', 'BLA', 'BM6', 'BMA', 'BME', 'BNG', 'BNZ', 'BOG', 'BR', 'BRL', 'BTW', 'C1O', 'C2O', 'C8E', 'CA', 'CAC', 'CAM', 'CCS', 'CD', 'CEM', 'CF0', 'CFM', 'CGU', 'CIT', 'CL', 'CLF', 'CME', 'CMP', 'CO', 'CO3', 'CO9', 'COT', 'CSO', 'CSS', 'CSW', 'CSX', 'CTP', 'CU', 'CU1', 'CUA', 'CYN', 'CYS', 'DAF', 'DDE', 'DDF', 'DGA', 'DGN', 'DIO', 'DKA', 'DLY', 'DMS', 'DPG', 'DSN', 'DVA', 'DVR', 'E64', 'E6C', 'EBP', 'EDO', 'EMO', 'EPE', 'EQP', 'EQU', 'ETX', 'F09', 'F25', 'FAD', 'FCY', 'FE', 'FES', 'FFO', 'FK5', 'FKP', 'FLA', 'FLC', 'FME', 'FMN', 'FMT', 'FOK', 'FOS', 'FUC', 'FUL', 'G6D', 'GAI', 'GAL', 'GCP', 'GDP', 'GL0', 'GLA', 'GLC', 'GLU', 'GLY', 'GLZ', 'GNH', 'GNP', 'GOL', 'GSH', 'GSP', 'GTP', 'GTS', 'GTT', 'GVE', 'HAS', 'HCA', 'HEA', 'HEC', 'HED', 'HEM', 'HG', 'HIC', 'HMC', 'HPR', 'HYP', 'IAS', 'IBR', 'IET', 'IGP', 'IMD', 'IOD', 'IPA', 'IPH', 'IUM', 'K', 'LAR', 'LAT', 'LBT', 'LDA', 'LNO', 'LPX', 'LYZ', 'M3L', 'MAK', 'MAL', 'MAN', 'MES', 'MG', 'MHO', 'MLA', 'MN', 'MN3', 'MNA', 'MNQ', 'MOH', 'MPC', 'MPD', 'MPT', 'MRD', 'MSE', 'MTX', 'MYR', 'NA', 'NAD', 'NAG', 'NAI', 'NAP', 'NDG', 'NDP', 'NES', 'NH2', 'NH4', 'NHE', 'NI', 'NO3', 'NPS', 'O', 'OCT', 'OLA', 'OXL', 'OXP', 'P34', 'P5P', 'P6G', 'PBM', 'PCA', 'PG4', 'PGA', 'PGE', 'PHQ', 'PLG', 'PLM', 'PLP', 'PMP', 'PMS', 'PO4', 'PON', 'PRO', 'PTR', 'RAP', 'RB', 'RDA', 'REA', 'RET', 'RIP', 'RS2', 'RTL', 'SCH', 'SE4', 'SEP', 'SF4', 'SIA', 'SIN', 'SM', 'SMC', 'SN1', 'SO2', 'SO3', 'SO4', 'SOC', 'SPA', 'SR', 'STI', 'STU', 'SVA', 'TAC', 'TAD', 'TAR', 'TBA', 'TEM', 'THP', 'THU', 'TLA', 'TMP', 'TPO', 'TRN', 'TRP', 'TRQ', 'TRS', 'TXP', 'TYS', 'UNL', 'UNX', 'URA', 'V7O', 'VAL', 'XYP', 'Y1', 'YBT', 'ZAP', 'ZN', 'ZNH']
        num_ligands = len(ligand_codes)
        c = 0
        for l in ligand_codes:
            c += 1
            #if c < 280:
            #    continue

            print('{0}/{1}'.format(c, num_ligands))
            importer.add_ligand_by_pdb_code(l)

    return
    importer._update_pdb_residues()
    return

    pdbfile = session.query(PDBFile).filter(PDBFile.ID == '1A2K')[0]

    for m in session.query(PDBMoleculeChain).filter(PDBMoleculeChain.PDBFileID == '1A2K'):
        print('')
        print(m)
        print(m.pdb_molecule)
        print(m.pdb_chain)

    print('')

    for m in session.query(PDBMolecule).filter(PDBMolecule.PDBFileID == '1A2K'):
        print(m)
        for mc in m.chains:
            for r in mc.pdb_chain.residues:
                print(r)

    print('')

test()

#from ppi_api import get_interface as get_ppi_interface
#ppi_api = get_ppi_interface(read_file('../pw'))


sys.exit(0)

class PDBFileImporter(object):
    '''This class is responsible for adding all records related to a PDBFile record i.e.:
          PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue
    '''

    # At the time of writing, these PDB IDs had no UniProt entries
    NoUniProtIDs = set(['1GTX', '1UOX', '1WSY', '1YYJ', '2IMM', '2MBP'])

    # At the time of writing, these PDB IDs had no JRNL lines
    NoPublicationData = ['2FX5']


    @staticmethod
    def retrieve(self, rcsb_pdb_id):
        '''Retrieves a PDB file from the RCSB and adds the associated metadata.'''
        pass


    @staticmethod
    def from_filepath(self, pdb_id, filepath, file_source):
        '''Retrieves a PDB file from the RCSB and adds the associated metadata.'''
        assert(file_source != 'RCSB')


    def __init__(self,
                 pdb_id, content = None, protein = None, contains_membrane_protein = None, file_source = None, filepath = None,
                 UniProtAC = None, UniProtID = None, testonly = False, techniques = None, derived_from = None, notes = None):
        '''UniProtACs have forms like 'P62937' whereas UniProtIDs have forms like 'PPIA_HUMAN.'''

        if derived_from != None:
            assert(isinstance(derived_from, str) and len(derived_from) == 4)

        self.dict = dict(
            ID = pdb_id,
            FileSource = file_source,
            Content = content,
            FASTA = None,
            Resolution = None,
            Techniques = techniques,
            BFactors = None,
            Publication =  None,
            Transmembrane = contains_membrane_protein,
            Notes = notes,
            DerivedFrom = derived_from,
        )

        self.testonly = testonly
        self.filepath = filepath
        self.UniProtAC = UniProtAC
        self.UniProtID = UniProtID
        self.ACtoID_mapping = None
        self.PDBtoAC_mapping = None
        self.chains = None
        self.pdb_object = None


    def get_contents(self):
        if self.dict['Content']:
            return self.dict['Content']

        ddGdb = self.ddGdb
        d = self.dict
        pdb_id = d['ID']
        if len(pdb_id) != 4:
            print(pdb_id)
        assert(len(pdb_id) <= 10)

        if self.filepath:
            filename = self.filepath
        else:
            filename = os.path.join("../pdbs", pdb_id + ".pdb")
        contents = None
        chains = {}

        if not os.path.exists(filename):
            sys.stdout.write("The file for %s is missing. Retrieving it now from RCSB: " % (pdb_id))
            sys.stdout.flush()
            try:
                contents = rcsb.retrieve_pdb(pdb_id)
                self.dict['FileSource'] = 'RCSB'
                write_file(os.path.join("../pdbs", pdb_id + ".pdb"), contents)
            except:
                print(traceback.format_exc())
                raise Exception("Error retrieving %s." % filename)
        else:
            contents = read_file(filename)

        self.dict['Content'] = contents
        return contents


    def parse_pdb_object(self):

        pdb_object = PDB(self.get_contents())
        self.pdb_object = pdb_object

        # Resolution
        self.dict['Resolution'] = pdb_object.get_resolution()
        if not self.dict['Resolution']:
            colortext.error("Could not determine resolution for %s." % self.filepath)
            #raise Exception("Could not determine resolution for %s." % filename)
        if self.dict['Resolution'] == "N/A":
            self.dict['Resolution'] = None

        # Techniques
        self.dict['Techniques'] = pdb_object.get_techniques() or self.dict['Techniques']

        # FASTA
        self.dict['FASTA'] = pdb_object.create_fasta()

        # Checks
        self.chains = pdb_object.atom_chain_order
        if not self.chains:
            raise Exception('No ATOM chains were found in the PDB file.')

        #UniqueIDs[pdb_id] = True
        pdb_id = self.dict['ID']
        ref_pdb_id = self.dict['DerivedFrom'] or self.dict['ID']
        if ref_pdb_id not in self.NoUniProtIDs:
            read_UniProt_map(self.ddGdb)
            if not PDBToUniProt.get(ref_pdb_id):
                if not (self.UniProtAC and self.UniProtID):
                    ACtoID_mapping, PDBtoAC_mapping = None, None
                    try:
                        getUniProtMapping(ref_pdb_id, storeInDatabase = False, ddGdb = self.ddGdb)
                        if not (PDBtoAC_mapping and ACtoID_mapping):
                            raise Exception("Could not find a UniProt mapping for %s in %s." % (ref_pdb_id, uniprotmapping))
                    except:
                        colortext.error("Could not find a UniProt mapping for %s in %s." % (ref_pdb_id, uniprotmapping))
                    self.ACtoID_mapping = ACtoID_mapping
                    self.PDBtoAC_mapping = PDBtoAC_mapping

        if pdb_id not in self.NoPublicationData:
            self.get_publication()

        foundRes = pdb_object.CheckForPresenceOf(["CSE", "MSE"])
        if foundRes:
            colortext.error("The PDB %s contains residues which could affect computation (%s)." % (pdb_id, join(foundRes, ", ")))
            if "CSE" in foundRes:
                colortext.error("The PDB %s contains CSE. Check." % pdb_id)
            if "MSE" in foundRes:
                colortext.error("The PDB %s contains MSE. Check." % pdb_id)

        self.dict['BFactors'] = pickle.dumps(pdb_object.ComputeBFactors())
        return pdb_object


    def get_publication(self):
        '''Extracts the PDB source information.'''
        ddGdb = self.ddGdb
        d = self.dict

        PUBTYPES = ['ISSN', 'ESSN']

        p = PDB(d['Content'].split("\n"))
        j = p.get_journal()
        if j:
            pdbID = d['ID'].strip()

            # We identify the sources for a PDB identifier with that identifier
            PublicationID = "PDB:%s" % pdbID
            sourceExists = self.ddGdb.locked_execute("SELECT ID FROM Publication WHERE ID=%s", parameters=(PublicationID,))
            if not sourceExists:
                if not self.testonly:
                    self.ddGdb.insertDict('Publication', {'ID' : PublicationID})

            d['Publication'] = PublicationID

            locations = self.ddGdb.locked_execute("SELECT * FROM PublicationIdentifier WHERE SourceID=%s", parameters=(PublicationID,))
            publocations = [location for location in locations if location['Type'] in PUBTYPES]
            doilocations = [location for location in locations if location['Type'] == "DOI"]
            assert(len(publocations) <= 1)
            assert(len(doilocations) <= 1)
            if j["published"]:
                skip = False
                if publocations:
                    location = publocations[0]
                    if j["REFN"]["type"] == location['Type']:
                        if j["REFN"]["ID"] != location['ID']:
                            colortext.warning("REFN: Check that the PublicationIdentifier data ('%s') matches the PDB REFN data ('%s')." % (str(location), j["REFN"]))
                else:
                    if j.get('REFN'):
                        assert(j["REFN"]["type"] in PUBTYPES)
                        source_location_dict = dict(
                            SourceID	= PublicationID,
                            ID			= j["REFN"]["ID"],
                            Type		= j["REFN"]["type"],
                        )
                        if not self.testonly:
                            self.ddGdb.insertDict('PublicationIdentifier', source_location_dict)
            if j["DOI"]:
                if doilocations:
                    location = doilocations[0]
                    if j["DOI"] != location['ID']:
                        colortext.warning("DOI: Check that the PublicationIdentifier data ('%s') matches the PDB DOI data ('%s')." % (str(doilocations), j["DOI"]))
                else:
                    source_location_dict = dict (
                        SourceID	= PublicationID,
                        ID			= j["DOI"],
                        Type		= "DOI",
                    )
                    if not self.testonly:
                        self.ddGdb.insertDict('PublicationIdentifier', source_location_dict)


    def parse_FASTA(self):
        pdbID = self.dict[self.ddGdb.FlatFieldNames.PDB_ID]
        results = self.ddGdb.execute_select("SELECT FASTA FROM Structure WHERE PDB_ID = %s", parameters = (pdbID,))
        if results:
            assert(len(results) == 1)
            return rcsb.parseFASTAs(results[0]["FASTA"])
        else:
            return rcsb.parseFASTAs(pdbID)

    def find(self):
        results = self.ddGdb.locked_execute('''SELECT PDB_ID FROM Structure WHERE PDB_ID=%s''', parameters = (self.PDB_ID, ))
        if results:
            return self.PDB_ID, self.PDB_ID
        return None, None

    def test(self):
        # #todo: Implement this later
        pass

    def commit(self, testonly = False, allow_missing_molecules = False):
        '''Returns the database record ID if an insert occurs but will typically return None if the PDB is already in the database.
           allow_missing_molecules should be set if molecules have been removed from the file e.g. Rosetta-generated structure
           but the headers still contain this information.
        '''
        d = self.dict
        pdb_id = d['ID']
        ddGdb = self.ddGdb
        testonly = testonly or self.testonly

        pdb_object = self.parse_pdb_object()

        if self.UniProtAC and self.UniProtID:
            # todo: Append to uniprotmapping.csv file
            results = self.ddGdb.locked_execute("SELECT * FROM UniProtKB WHERE UniProtKB_AC=%s", parameters=(self.UniProtAC,))
            if results:
                if results[0]['UniProtKB_ID'] != self.UniProtID:
                    raise Exception("Existing UniProt mapping (%s->%s) does not agree with the passed-in parameters (%s->%s)." % (results[0]['UniProtKB_AC'],results[0]['UniProtKB_ID'],self.UniProtAC,self.UniProtID))
            else:
                UniProtMapping = dict(
                    UniProtKB_AC = self.UniProtAC,
                    UniProtKB_ID = self.UniProtID,
                )
                if not testonly:
                    self.ddGdb.insertDict('UniProtKB', UniProtMapping)

        # Either add a new PDBFile record or update an existing one. todo: not all fields are currently updated
        results = self.ddGdb.locked_execute("SELECT * FROM PDBFile WHERE ID=%s", parameters = (d['ID'],))
        if results:
            assert(len(results) == 1)
            result = results[0]
            pdbID = result['ID']
            for k, v in d.iteritems():
                if k != 'PDBFileID':
                    if k == 'Techniques' and result[k] == "":
                        SQL = "UPDATE PDBFile SET %s" % k
                        SQL += "=%s WHERE PDB_ID=%s"
                        if not testonly:
                            results = self.ddGdb.locked_execute(SQL, parameters = (v, pdbID))
                    if d[k] and not(result[k]):
                        SQL = "UPDATE PDBFile SET %s" % k
                        SQL += "=%s WHERE PDB_ID=%s"
                        if not testonly:
                            results = self.ddGdb.locked_execute(SQL, parameters = (v, pdbID))
        else:
            if not testonly:
                if self.dict['FileSource'] == None:
                    raise Exception('A file source must be specified.')

                self.ddGdb.insertDict('PDBFile', self.dict)
                self.databaseID = self.ddGdb.getLastRowID()

        # Add PDBChain records
        for c in self.chains:
            pdbc = dict(
                PDBFileID = pdb_id,
                Chain = c,
            )
            self.ddGdb.insertDictIfNew('PDBChain', pdbc, ['PDBFileID', 'Chain'])

        # Add PDBMolecule and PDBMoleculeChain records
        try:
            molecules = pdb_object.get_molecules_and_source()
        except MissingRecordsException:
            molecules = []
        for molecule in molecules:
            chains = molecule['Chains']
            molecule['PDBFileID'] = pdb_id
            molecule['Organism'] = molecule['OrganismScientificName'] or molecule['OrganismCommonName']
            md = {}
            for k in ['PDBFileID', 'MoleculeID', 'Name', 'Organism', 'Fragment', 'Synonym', 'Engineered', 'EC', 'Mutation', 'OtherDetails']:
                md[k] = molecule[k]
            self.ddGdb.insertDictIfNew('PDBMolecule', md, ['PDBFileID', 'MoleculeID'])
            for c in chains:
                try:
                    mcd = dict(
                        PDBFileID = pdb_id,
                        MoleculeID = md['MoleculeID'],
                        Chain = c
                    )
                    self.ddGdb.insertDictIfNew('PDBMoleculeChain', mcd, ['PDBFileID', 'MoleculeID', 'Chain'])
                except:
                    if allow_missing_molecules: pass
                    else: raise

        # Add PDBResidue records
        for c, seq in pdb_object.atom_sequences.iteritems():
            count = 1
            for s in seq:
                res_id, r = s
                assert(len(r.ResidueID) == 5)
                assert(c == r.Chain)
                db_res = dict(
                    PDBFileID = pdb_id,
                    Chain = c,
                    ResidueID = r.ResidueID,
                    ResidueAA = r.ResidueAA,
                    ResidueType = r.residue_type,
                    IndexWithinChain = count,
                    CoordinatesExist = True,
                    RecognizedByRosetta = None,
                    BFactorMean = None,
                    BFactorDeviation = None,
                    SecondaryStructurePosition = None,
                    AccessibleSurfaceArea = None,
                    MonomericExposure = None,
                    MonomericDSSP = None,
                    ComplexExposure = None,
                    ComplexDSSP = None,
                )
                self.ddGdb.insertDictIfNew('PDBResidue', db_res, ['PDBFileID', 'Chain', 'ResidueID'])
                count += 1

        # Add UniProt mapping
        if self.UniProtAC and self.UniProtID:
            results = self.ddGdb.locked_execute("SELECT * FROM UniProtKBMapping WHERE UniProtKB_AC=%s", parameters=(self.UniProtAC,))
            if results:
                if results[0]['PDBFileID'] != d['ID']:
                    raise Exception("Existing UniProt mapping (%s->%s) does not agree with the passed-in parameters (%s->%s)." % (results[0]['UniProtKB_AC'],results[0]['PDBFileID'],self.UniProtAC,d['ID']))
            else:
                UniProtPDBMapping = dict(
                    UniProtKB_AC = self.UniProtAC,
                    PDBFileID = d[FieldNames_.PDB_ID],
                )
                if not testonly:
                    self.ddGdb.insertDict('UniProtKBMapping', UniProtPDBMapping)

        # Store the UniProt mapping in the database
        ref_pdb_id = self.dict['DerivedFrom'] or self.dict['ID']
        if ref_pdb_id not in self.NoUniProtIDs:
            if not (self.ACtoID_mapping and self.PDBtoAC_mapping):
                try:
                    self.ACtoID_mapping, self.PDBtoAC_mapping = getUniProtMapping(ref_pdb_id, storeInDatabase = testonly)
                except:
                    if False and self.dict['Techniques'] != 'Rosetta model':
                        raise
            if False and self.dict['Techniques'] != 'Rosetta model':
                assert(self.ACtoID_mapping and self.PDBtoAC_mapping)
            if not testonly:
                if self.ACtoID_mapping and self.PDBtoAC_mapping:
                    commitUniProtMapping(self.ddGdb, self.ACtoID_mapping, self.PDBtoAC_mapping)

        # Run DSSP
        dssp_complex_d, dssp_monomer_d = None, None
        try:
            # This fails for some PDB e.g. if they only have CA atoms
            dssp_complex_d = ComplexDSSP(self.pdb_object)
        except MissingAtomException, e:
            print('DSSP (complex) failed for this case.')
        try:
            # This fails for some PDB e.g. if they only have CA atoms
            dssp_monomer_d = MonomerDSSP(self.pdb_object)
        except MissingAtomException, e:
            print('DSSP (monomer) failed for this case.')

        # Make sure that the residue records exist for all results of DSSP
        if dssp_monomer_d:
            for chain_id in self.pdb_object.atom_sequences.keys():
                if self.pdb_object.chain_types[chain_id] == 'Protein':
                    for chain_id, mapping in dssp_monomer_d:
                        for residue_id, residue_details in sorted(mapping.iteritems()):
                            residue_record = ddGdb.execute_select('SELECT ID FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                            assert(len(residue_record) == 1)

            # Add the monomeric DSSP results
            for chain_id in self.pdb_object.atom_sequences.keys():
                if self.pdb_object.chain_types[chain_id] == 'Protein':
                    colortext.warning('\tRunning DSSP on chain %s' % chain_id)
                    for chain_id, mapping in dssp_monomer_d:
                        for residue_id, residue_details in sorted(mapping.iteritems()):
                            residue_record = ddGdb.execute_select('SELECT ID, MonomericExposure, MonomericDSSP FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                            assert(len(residue_record) == 1)
                            PDBResidueID = residue_record[0]['ID']
                            if residue_record[0]['MonomericDSSP'] == None:
                                ddGdb.execute('UPDATE PDBResidue SET MonomericDSSP=%s WHERE ID=%s', parameters=(residue_details['ss'], PDBResidueID))
                            if residue_record[0]['MonomericExposure'] == None:
                                ddGdb.execute('UPDATE PDBResidue SET MonomericExposure=%s WHERE ID=%s', parameters=(residue_details['exposure'], PDBResidueID))
        # Make sure that the residue records exist for all results of DSSP
        if dssp_complex_d:
            for chain_id in self.pdb_object.atom_sequences.keys():
                if self.pdb_object.chain_types[chain_id] == 'Protein':
                    for chain_id, mapping in dssp_complex_d:
                        for residue_id, residue_details in sorted(mapping.iteritems()):
                            residue_record = ddGdb.execute_select('SELECT ID FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                            assert(len(residue_record) == 1)

            # Add the complex DSSP results
            for chain_id in self.pdb_object.atom_sequences.keys():
                if self.pdb_object.chain_types[chain_id] == 'Protein':
                    colortext.warning('\tRunning DSSP on chain %s' % chain_id)
                    for chain_id, mapping in dssp_complex_d:
                        for residue_id, residue_details in sorted(mapping.iteritems()):
                            residue_record = ddGdb.execute_select('SELECT ID, ComplexExposure, ComplexDSSP FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(pdb_id, chain_id, residue_id))
                            assert(len(residue_record) == 1)
                            PDBResidueID = residue_record[0]['ID']
                            if residue_record[0]['ComplexDSSP'] == None:
                                ddGdb.execute('UPDATE PDBResidue SET ComplexDSSP=%s WHERE ID=%s', parameters=(residue_details['ss'], PDBResidueID))
                            if residue_record[0]['ComplexExposure'] == None:
                                ddGdb.execute('UPDATE PDBResidue SET ComplexExposure=%s WHERE ID=%s', parameters=(residue_details['exposure'], PDBResidueID))

        if not testonly:
            return self.databaseID
        else:
            return None

        #raise Exception('Before using this again, add the functionality from ddgadmin/updatedb/compute_all_dssp.py to add the molecules and DSSP values.')

        return self.databaseID

    def remove(self):
        # We do not usually want to remove PDBs from the database
        pass

    def __repr__(self):
        ddGdb = self.ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        d = self.dict
        str = []
        str.append("PDBFileID: %s" % d['ID'])
        str.append("Protein: %s" % d['FileSource'])
        return join(str, "\n")