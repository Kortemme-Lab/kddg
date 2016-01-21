#!/usr/bin/python2.4
# encoding: utf-8
"""
import_api.py
High-level functions for importing data into the DDG database.

Created by Shane O'Connor 2015.
Copyright (c) 2015 Shane O'Connor. All rights reserved.
"""

# todo: add_pdb_from_rcsb: FileContent (see ticket 1489)
# todo: get_pdb_details (see below)
# todo: PDB Chains (_add_pdb_chains)
#       add self.pfam_api.? call to get mapping from RCSB PDB files to the UniProt sequences. Add this to a separate function, _add_pdb_uniprot_mapping.
#       add self.pfam_api.get_pfam_accession_numbers_from_pdb_chain(database_pdb_id, c)) calls. Use this in _add_pdb_uniprot_mapping.
#       add self.scope_api.get_chain_details(database_pdb_id, c))) calls
#       add coordinates
# todo: _add_pdb_uniprot_mapping. implement UniProt mapping. The old approach was flawed. Use the new approach (similar to how I used PPComplex as an abstraction layer)

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
import time

import numpy
import pandas

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import create_engine, and_
from sqlalchemy import inspect as sqlalchemy_inspect
from sqlalchemy.exc import TimeoutError as SQLAlchemyTimeoutError
from MySQLdb import OperationalError as MySQLOperationalError

if __name__ == '__main__':
    sys.path.insert(0, '../../klab')

from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation
from klab.fs.fsio import read_file, write_temp_file, open_temp_file, write_file
from klab.bio.pfam import Pfam
from klab.bio.dssp import MonomerDSSP, ComplexDSSP, MissingAtomException
from klab.bio.ligand import Ligand, PDBLigand
from klab.bio.pdbtm import PDBTM
from klab.db.sqlalchemy_interface import get_single_record_from_query, get_or_create_in_transaction
from klab.bio import rcsb

from db_schema import test_schema_against_database_instance
from db_schema import PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue, LigandDescriptor, LigandIdentifier, LigandSynonym, LigandPrice, LigandReference, PDBLigand, PDBIon, FileContent
from db_schema import Ligand as DBLigand
from db_schema import Ion as DBIon
from db_schema import Publication, PublicationAuthor, PublicationIdentifier
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


    def __init__(self, passwd, connect_string, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, cache_dir = None, echo_sql = False):
        '''
        :param passwd:
        :param connect_string:
        :param username:
        :param hostname:
        :param rosetta_scripts_path:
        :param rosetta_database_path:
        :param cache_dir: Used to cache downloaded files e.g. PDB files from the RCSB servers. Particularly useful during testing to avoid spamming their servers with requests.
        :param echo_sql: If echo_sql is set then all executed SQL commands are printed to stdout (by default) which may be useful for debugging.
        :return:
        '''

        # Set up MySQLdb connections
        passwd = passwd.strip()
        self.DDG_db = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname)
        self.DDG_db_utf = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname, use_utf = True)

        self.cache_dir = cache_dir
        self.echo_sql = echo_sql

        test_schema_against_database_instance(self.DDG_db)

        if self.cache_dir:
            self.initialize_cache_directory()

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
    def get_interface_with_config_file(cls, database = 'ddg', host_config_name = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None, my_cnf_path = None, cache_dir = None, echo_sql = False):
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

        return cls(password, connection_string, username = user, hostname = host, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, cache_dir = cache_dir, echo_sql = echo_sql)


    def __del__(self):
        pass         #self.DDG_db.close()         #self.ddGDataDB.close()



    #############################
    #                           #
    #  SQLAlchemy Engine setup  #
    #                           #
    #############################



    def get_engine(self):
        if not self.engine:
            self.engine = create_engine(self.connect_string, echo = self.echo_sql)
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


    #################################
    #                               #
    #  Data removal API             #
    #                               #
    #################################


    def remove_ligands(self):
        '''This function should not generally be called. It was added while testing the ligand addition code.
           Removals are protected by sessions to prevent partial deletions and ligands will only be removed with this
           function if there are no corresponding records in other tables e.g. PDBLigand.

           I initially added ligands based on IDs extracted from PDBs however this list included ions like FE2 so I wanted
           to scrub all existing records and only add ligands with >1 atoms to the Ligand table. Ions are now added to the
           Ion table.
        '''

        tsession = self.get_session()
        ligand_ids = [l.ID for l in tsession.query(DBLigand)]
        for ligand_id in ligand_ids:
            tsession = self.get_session(new_session = True) # do not allow partial deletions
            try:
                colortext.message('Removing ligand {0}.'.format(ligand_id))
                for ltbl in [LigandDescriptor, LigandIdentifier, LigandPrice, LigandReference, LigandSynonym]:
                    tsession.query(ltbl).filter(ltbl.LigandID == ligand_id).delete()
                tsession.query(DBLigand).filter(DBLigand.ID == ligand_id).delete()
                tsession.commit()
                tsession.close()
                print('Success.\n')
            except Exception, e:
                colortext.error('Failure.')
                print(str(e))
                print(traceback.format_exc() + '\n')
                tsession.rollback()
                tsession.close()
            print('')


    #################################
    #                               #
    #  Data update API              #
    #                               #
    #################################


    def update_pdbs(self, pdb_ids = [], update_sections = set(), start_at = None, restrict_to_file_source = None):
        '''Updates all or selected data for all or selected PDB files in the database, enumerating in alphabetical order.
           If start_at is specified, the update begins at the specified PDB identifier.
        '''
        if not pdb_ids:
            pdb_ids = [r for r in self.DDG_db.execute_select('SELECT ID, FileSource FROM PDBFile ORDER BY ID')]

        counter = 0
        num_pdb_ids = len(pdb_ids)
        hit_starting_pdb = False
        for r in pdb_ids:
            counter += 1
            pdb_id = r['ID']
            if (not restrict_to_file_source) or (r['FileSource'] == restrict_to_file_source):
                if (not start_at) or (r['ID'].upper() == start_at):
                    hit_starting_pdb = True
                if hit_starting_pdb:
                    colortext.message('Updating data for {0} ({1}/{2}).'.format(pdb_id, counter, num_pdb_ids))
                    tsession = self.get_session(new_session = True)
                    self.add_pdb_data(tsession, pdb_id, update_sections = update_sections)
                    tsession.commit()
                    tsession.close()
        if not hit_starting_pdb:
            raise Exception('We never hit the starting PDB "{0}".'.format(start_at))

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
           The value of the ID field of the Ligand record is returned.
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
            db_ligand = get_or_create_in_transaction(tsession, DBLigand, l.__dict__, missing_columns = ['ID'])

            # Create the ligand descriptor records
            descriptor_fieldnames = [c.name for c in list(sqlalchemy_inspect(LigandDescriptor).columns)]

            for descriptor in l.descriptors:
                descriptor = copy.deepcopy(descriptor)
                descriptor['LigandID'] = db_ligand.ID
                db_ligand_descriptor = get_or_create_in_transaction(tsession, LigandDescriptor, descriptor, missing_columns = ['ID'])

            for identifier in l.identifiers:
                identifier = copy.deepcopy(identifier)
                identifier['LigandID'] = db_ligand.ID
                db_ligand_identifier = get_or_create_in_transaction(tsession, LigandIdentifier, identifier, missing_columns = ['ID'])

            for synonym in l.synonyms:
                db_ligand_synonym = get_or_create_in_transaction(tsession, LigandSynonym, dict(LigandID = db_ligand.ID, Synonym = synonym.strip()))

            db_ligand_id = db_ligand.ID
            tsession.commit()
            tsession.close()
            print('Success.\n')
            return db_ligand_id
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise


    #################################
    #                               #
    #  PDB data retrieval API       #
    #                               #
    #################################


    def get_pdb_details(self, pdb_ids, cached_pdb_details = None):
        '''Returns the details stored in the database about the PDB files associated with pdb_ids e.g. chains, resolution,
           technique used to determine the structure etc.'''
        pdbs = {}
        cached_pdb_ids = []
        if cached_pdb_details:
            cached_pdb_ids = set(cached_pdb_details.keys())
        for pdb_id in pdb_ids:
            if pdb_id in cached_pdb_ids:
                pdbs[pdb_id] = cached_pdb_details[pdb_id]
            else:
                record = self.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))[0]
                p = PDB(record['Content'])
                pdb_chain_lengths = {}
                for chain_id, s in p.atom_sequences.iteritems():
                    pdb_chain_lengths[chain_id] = len(s)
                # todo: get the list of protein chains and PDB residues from the database and assert that they are the same
                #       as what were extracted from the PDB file.
                #       maybe change 'chains' below to 'protein_chains'
                pdbs[pdb_id] = dict(
                    chains = pdb_chain_lengths,
                    TM = record['Transmembrane'],
                    Technique = record['Techniques'],
                    XRay = record['Techniques'].find('X-RAY') != -1,
                    Resolution = record['Resolution'],
                )
        return pdbs


    #################################
    #                               #
    #  PDB file entry - public API  #
    #                               #
    #  Missing tables:              #
    #      FileContent              #
    #                               #
    #      Protein                  #
    #      ProteinDatabaseIdentifier#
    #      ProteinName              #
    #      ProteinOrganism          #
    #      ProteinResidue           #
    #      ProteinSegment           #
    #                               #
    #################################


    def initialize_cache_directory(self):
        pdb_dir = os.path.join(self.cache_dir, 'pdbs')
        try:
            os.makedirs(pdb_dir)
        except: pass
        assert(os.path.exists(pdb_dir))


    def _retrieve_pdb_contents(self, pdb_id, fresh = False):
        if (not fresh) and (self.cache_dir):
            cached_filepath = os.path.join(self.cache_dir, 'pdbs', '{0}.pdb'.format(pdb_id))
            if os.path.exists(cached_filepath):
                print('Retrieving locally cached file for {0}.'.format(pdb_id))
                return read_file(cached_filepath)
            elif os.path.exists(cached_filepath + '.gz'):
                print('Retrieving locally cached file for {0}.'.format(pdb_id))
                return read_file(cached_filepath + '.gz')

        print('Retrieving {0} from RCSB.'.format(pdb_id))
        contents = rcsb.retrieve_pdb(pdb_id)
        if self.cache_dir:
            write_file(os.path.join(self.cache_dir, 'pdbs', '{0}.pdb'.format(pdb_id)), contents)
        return contents


    def add_pdb_from_rcsb(self, pdb_id, previously_added = set(), update_sections = set()):
        '''NOTE: This API is used to create and analysis predictions or retrieve information from the database.
                 This function adds new raw data to the database and does not seem to belong here. It should be moved into
                 an admin API instead.
           This function adds imports a PDB into the database, creating the associated molecule, chain and residue etc. records.
           If previously_added contains pdb_id then we return. Otherwise, we step through the full PDB file import. This
           is added purely for optimistic efficiency e.g. if were to add 100 designed files based off the same RCSB PDB,
           we would not want to run the full import code for the RCSB PDB file more than once.

           Touched tables:
               PDBFile
               todo: FileContent (see ticket 1489)
        '''

        if pdb_id in previously_added:
            return pdb_id
        assert(len(pdb_id) == 4)
        pdb_id = pdb_id.upper()

        tsession = self.get_session(new_session = True)
        try:
            db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == pdb_id))
            existing_record = db_record != None
            if existing_record:
                print('Retrieving {0} from database.'.format(pdb_id))
                pdb_object = PDB(db_record.Content)
                assert(db_record.FileSource == 'RCSB')
            else:
                db_record = PDBFile()
                contents = self._retrieve_pdb_contents(pdb_id)
                pdb_object = PDB(contents)
                db_record.ID = pdb_id
                db_record.FileSource = 'RCSB'
                db_record.Content = contents
                update_sections = set() # add all related data

            # Checks
            self.chains = pdb_object.atom_chain_order
            if not self.chains:
                raise Exception('No ATOM chains were found in the PDB file.')
            foundRes = pdb_object.CheckForPresenceOf(["CSE", "MSE"])
            if foundRes:
                colortext.error("The PDB %s contains residues which could affect computation (%s)." % (pdb_id, join(foundRes, ", ")))
                if "CSE" in foundRes:
                    colortext.error("The PDB %s contains CSE. Check." % pdb_id)
                if "MSE" in foundRes:
                    colortext.error("The PDB %s contains MSE. Check." % pdb_id)

            # FASTA
            if not(existing_record) or (not db_record.FASTA):
                db_record.FASTA = pdb_object.create_fasta()

            # Resolution
            if not(existing_record) or (not db_record.Resolution):
                resolution = pdb_object.get_resolution()
                if not resolution:
                    colortext.error("Could not determine resolution for %s." % self.filepath)
                if resolution == "N/A":
                    resolution = None
                db_record.Resolution = resolution

            # Techniques
            if not(existing_record) or (not db_record.Techniques):
                db_record.Techniques = pdb_object.get_techniques()

            # Transmembrane
            if not(existing_record) or (db_record.Transmembrane == None):
                pdbtmo = PDBTM(read_file('/kortemmelab/shared/mirror/PDBTM/pdbtmall.xml'))
                pdb_id_map = pdbtmo.get_pdb_id_map()
                uc_pdb_id_map = {}
                for k, v in pdb_id_map.iteritems():
                    uc_pdb_id_map[k.upper()] = v
                if pdb_id in uc_pdb_id_map:
                    colortext.warning('{0} is a transmembrane protein.'.format(pdb_id))
                db_record.Transmembrane = pdb_id in pdb_id_map


            # Friday - add bfactors back
            overall_bfactors = pdb_object.get_B_factors().get('Overall')
            db_record.BFactorMean = overall_bfactors['mean']
            db_record.BFactorDeviation = overall_bfactors['stddev']

            # Add a new PDBFile record
            if not(existing_record):
                db_record.UserID = None
                db_record.Notes = None
                db_record.DerivedFrom = None
                tsession.add(db_record)
                tsession.flush()
                print(db_record)

            # Publication
            if not(existing_record) or (not db_record.Publication):
                self._add_pdb_publication(tsession, db_record.ID, pdb_object = pdb_object)

            # add all other data
            self.add_pdb_data(tsession, pdb_id, update_sections = update_sections)

            previously_added.add(pdb_id)

            print('Success.\n')
            tsession.commit()
            tsession.close()
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise


    def get_pdb_object(self, database_pdb_id, tsession = None):
        '''Create a PDB object from content in the database.'''
        tsession = tsession or self.get_session()
        db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id))
        assert(db_record)
        return PDB(db_record.Content)


    def add_pdb_data(self, tsession, database_pdb_id, update_sections = set(), ligand_mapping = {}, chain_mapping = {}):
        '''
           database_pdb_id is the RCSB ID for RCSB files and a custom ID for other (designed) structures.
           If transaction_session is None (e.g. if this was called directly outside of a transaction), create a transaction
           session for the remaining inner calls. If update_sections is non-empty, just call those specific inner functions.
           Note: If the caller creates a transaction then it is responsible for committing/rolling back the transaction.
           ligand_mapping = {}, chain_mapping = {}
        :param tsession:
        :param database_pdb_id:
        :param update_sections:
        :param ligand_mapping: A mapping from ligand IDs (e.g. "FPP") to RCSB IDs. This is only necessary for non-RCSB files as these files may have modified IDs.
        :param chain_mapping: A mapping from chain IDs (e.g. "A") to the chain in the original RCSB file. This is only necessary for non-RCSB files as these files may have modified chain IDs (e.g. changing the ligand chain to "X").
        :return:
        '''

        if not tsession:
            tsession = self.get_session(new_session = True)

        # Retrieve the PDB object
        pdb_object = self.get_pdb_object(database_pdb_id, tsession = tsession)

        if not(update_sections) or ('Chains' in update_sections):
            colortext.warning('*** Chains ***')
            self._add_pdb_chains(tsession, database_pdb_id, pdb_object, chain_mapping)
        if not(update_sections) or ('Molecules' in update_sections):
            colortext.warning('*** Molecules ***')
            self._add_pdb_molecules(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Residues' in update_sections):
            colortext.warning('*** Residues ***')
            self._add_pdb_residues(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Ligands' in update_sections):
            colortext.warning('*** Ligands ***')
            self._add_pdb_rcsb_ligands(tsession, database_pdb_id, pdb_object, ligand_mapping)
        if not(update_sections) or ('Ions' in update_sections):
            colortext.warning('*** Ions ***')
            self._add_pdb_rcsb_ions(tsession, database_pdb_id, pdb_object)
        #if not(update_sections) or ('UniProt' in update_sections):
        #    colortext.warning('*** UniProt ***')
        #    self._add_pdb_uniprot_mapping(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Publication' in update_sections):
            colortext.warning('*** Publication ***')
            self._add_pdb_publication(tsession, database_pdb_id, pdb_object)


    def _add_pdb_chains(self, tsession, database_pdb_id, pdb_object = None, chain_mapping = {}):
        '''
           Touched tables:
               PDBChain
        '''
        from db_schema import PDBFile, PDBChain
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)

        db_chains = {}
        for r in tsession.query(PDBChain).filter(PDBChain.PDBFileID == database_pdb_id).order_by(PDBChain.Chain):
            db_chains[r.Chain] = r

        db_chains_ids = sorted(db_chains.keys())
        chain_ids = sorted(set(pdb_object.seqres_sequences.keys() + pdb_object.atom_sequences.keys() + pdb_object.chain_types.keys()))

        if chain_ids != db_chains_ids:
            #colortext.warning('PDB chains: {0}\t DB chains: {1}'.format(','.join(chain_ids), ','.join(db_chains_ids)))
            #colortext.error('Missing chains.')
            new_chain_ids = sorted(set(chain_ids).difference(db_chains_ids))
            for c in new_chain_ids:
                db_chain = get_or_create_in_transaction(tsession, PDBChain, dict(
                    PDBFileID = database_pdb_id,
                    Chain = c,
                    MoleculeType = pdb_object.chain_types[c]
                ), missing_columns = ['WildtypeProteinID', 'FullProteinID', 'SegmentProteinID', 'WildtypeAlignedProteinID', 'AcquiredProteinID', 'Coordinates'])
            db_chains = {}
            for r in tsession.query(PDBChain).filter(PDBChain.PDBFileID == database_pdb_id).order_by(PDBChain.Chain):
                db_chains[r.Chain] = r
            db_chains_ids = sorted(db_chains.keys())
            assert(chain_ids == db_chains_ids)

        for chain_id in pdb_object.chain_types.keys():
            if pdb_object.chain_types[chain_id] != db_chains[chain_id].MoleculeType:
                db_chain = tsession.query(PDBChain).filter(and_(PDBChain.PDBFileID == database_pdb_id, PDBChain.Chain == chain_id)).one() # we expect exactly one record
                db_chain.MoleculeType = pdb_object.chain_types[chain_id]
                tsession.flush()

        return

        # todo: Extract and store the coordinates
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


    def _add_pdb_molecules(self, tsession, database_pdb_id, pdb_object = None, allow_missing_molecules = True):
        '''
           Add PDBMolecule and PDBMoleculeChain records
           Touched tables:
               PDBMolecule
               PDBMoleculeChain
        '''
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)
        try:
            molecules = pdb_object.get_molecules_and_source()
        except MissingRecordsException:
            molecules = []
        for molecule in molecules:
            chains = molecule['Chains']
            molecule['PDBFileID'] = database_pdb_id
            molecule['Organism'] = molecule['OrganismScientificName'] or molecule['OrganismCommonName']
            md = {}
            for k in ['PDBFileID', 'MoleculeID', 'Name', 'Organism', 'Fragment', 'Synonym', 'Engineered', 'EC', 'Mutation', 'OtherDetails']:
                md[k] = molecule[k]

            #db_molecule = get_or_create_in_transaction(tsession, PDBMolecule, md, updatable_columns = ['Synonym', 'OtherDetails']) # this was used to fix a problem where the db entries had truncated values. I have since resized the field and updated the records.
            db_molecule = get_or_create_in_transaction(tsession, PDBMolecule, md)
            for c in chains:
                try:
                    db_molecule_chain = get_or_create_in_transaction(tsession, PDBMoleculeChain, dict(
                        PDBFileID = database_pdb_id,
                        MoleculeID = md['MoleculeID'],
                        Chain = c
                    ))
                except:
                    if allow_missing_molecules: pass
                    else: raise


    def _add_pdb_residues(self, tsession, database_pdb_id, pdb_object = None):
        '''
           Touched tables:
               PDBResidue
        '''
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)
        residue_bfactors = pdb_object.get_B_factors().get('PerResidue')

        # Run DSSP over the entire structure
        dssp_complex_d, dssp_monomer_d = None, None
        try:
            # This fails for some PDB e.g. if they only have CA atoms
            dssp_complex_d = ComplexDSSP(pdb_object, read_only = True)
        except MissingAtomException, e:
            print('DSSP (complex) failed for this case: {0}.'.format(database_pdb_id))
            for db_chain in tsession.query(PDBChain).filter(PDBChain.PDBFileID == database_pdb_id):
                print(db_chain.Chain, db_chain.MoleculeType)
                assert(db_chain.MoleculeType != 'Protein') # we should always pass on protein chains

        # Run DSSP over the individual chains
        try:
            # This fails for some PDB e.g. if they only have CA atoms
            dssp_monomer_d = MonomerDSSP(pdb_object, read_only = True)
        except MissingAtomException, e:
            print('DSSP (monomer) failed for this case: {0}.'.format(database_pdb_id))
            for db_chain in tsession.query(PDBChain).filter(PDBChain.PDBFileID == database_pdb_id):
                print(db_chain.Chain, db_chain.MoleculeType)
                assert(db_chain.MoleculeType != 'Protein') # we should always pass on protein chains

        # Generate a list of residues with coordinates in the PDB file
        parsed_residues = set()
        for c, seq in pdb_object.atom_sequences.iteritems():
            for s in seq:
                res_id, r = s
                parsed_residues.add(c + r.ResidueID)

        # Sanity checks: make sure that the residue records exist for all results of DSSP
        monomeric_records, complex_records = {}, {}
        if dssp_monomer_d:
            for chain_id, mapping in dssp_monomer_d:
                if pdb_object.chain_types[chain_id] == 'Protein' or pdb_object.chain_types[chain_id] == 'Protein skeleton':
                   for residue_id, residue_details in sorted(mapping.iteritems()):
                        chain_residue_id = chain_id + residue_id
                        assert(chain_residue_id in parsed_residues)
                        monomeric_records[chain_residue_id] = residue_details
        if dssp_complex_d:
            for chain_id, mapping in dssp_complex_d:
                if pdb_object.chain_types[chain_id] == 'Protein' or pdb_object.chain_types[chain_id] == 'Protein skeleton':
                    for residue_id, residue_details in sorted(mapping.iteritems()):
                        chain_residue_id = chain_id + residue_id
                        assert(chain_residue_id in parsed_residues)
                        complex_records[chain_residue_id] = residue_details

        # Read existing data from the database
        existing_residues = {}
        for r in tsession.query(PDBResidue).filter(PDBResidue.PDBFileID == database_pdb_id):
            chain_residue_id = r.Chain + r.ResidueID
            assert(chain_residue_id not in existing_residues)
            existing_residues[chain_residue_id] = r

        # Add PDBResidue records
        # dssp_monomer_d and dssp_complex_d are maps: chain -> residue_id -> DSSP record
        for c, seq in sorted(pdb_object.atom_sequences.iteritems()):
            count = 1
            for s in seq:
                res_id, r = s
                assert(len(r.ResidueID) == 5)
                assert(c == r.Chain)

                chain_residue_id = c + r.ResidueID
                dssp_res_complex_ss, dssp_res_complex_exposure, dssp_res_monomer_ss, dssp_res_monomer_exposure = ' ', None, ' ', None
                monomeric_record = monomeric_records.get(chain_residue_id)
                if monomeric_record:
                    dssp_res_monomer_ss = monomeric_record['ss']
                    dssp_res_monomer_exposure = monomeric_record['exposure']
                complex_record = complex_records.get(chain_residue_id)
                if complex_record:
                    dssp_res_complex_ss = complex_record['ss']
                    dssp_res_complex_exposure = complex_record['exposure']

                average_bfactors = residue_bfactors.get(chain_residue_id, {})
                existing_residue_record = existing_residues.get(chain_residue_id)
                if not existing_residue_record:
                    db_residue = get_or_create_in_transaction(tsession, PDBResidue, dict(
                        PDBFileID = database_pdb_id,
                        Chain = c,
                        ResidueID = r.ResidueID,
                        ResidueAA = r.ResidueAA,
                        ResidueType = r.residue_type,
                        IndexWithinChain = count,
                        CoordinatesExist = True,
                        RecognizedByRosetta = None,
                        BFactorMean = average_bfactors.get('mean'),
                        BFactorDeviation = average_bfactors.get('stddev'),
                        MonomericExposure = dssp_res_monomer_exposure,
                        MonomericDSSP = dssp_res_monomer_ss,
                        ComplexExposure = dssp_res_complex_exposure,
                        ComplexDSSP = dssp_res_complex_ss,
                    ))
                else:
                    # Sanity check: make sure that the current data matches the database
                    #print('EXISTING RESIDUE')
                    #pprint.pprint(existing_residue_record.__dict__)
                    if existing_residue_record.BFactorMean != None:
                        assert(abs(float(existing_residue_record.BFactorMean) - average_bfactors.get('mean')) < 0.001)
                    if existing_residue_record.BFactorDeviation != None:
                        assert(abs(float(existing_residue_record.BFactorDeviation) - average_bfactors.get('stddev')) < 0.001)
                    if existing_residue_record.MonomericExposure != None:
                        assert(abs(float(existing_residue_record.MonomericExposure) - dssp_res_monomer_exposure) < 0.001)
                    if existing_residue_record.MonomericDSSP != None:
                        assert(r.ResidueAA == 'X' or (existing_residue_record.MonomericDSSP == dssp_res_monomer_ss))
                    if existing_residue_record.ComplexExposure != None:
                        assert(abs(float(existing_residue_record.ComplexExposure) - dssp_res_complex_exposure) < 0.001)
                    if existing_residue_record.ComplexDSSP != None:
                        assert(r.ResidueAA == 'X' or (existing_residue_record.ComplexDSSP == dssp_res_complex_ss))

                    # Update data (add new data if is was previously missing)
                    existing_residue_record.BFactorMean = average_bfactors.get('mean')
                    existing_residue_record.BFactorDeviation = average_bfactors.get('stddev')
                    existing_residue_record.MonomericExposure = dssp_res_monomer_exposure
                    existing_residue_record.MonomericDSSP = dssp_res_monomer_ss
                    existing_residue_record.ComplexExposure = dssp_res_complex_exposure
                    existing_residue_record.ComplexDSSP = dssp_res_complex_ss
                    tsession.flush()
                #self.ddGdb.insertDictIfNew('PDBResidue', db_res, ['PDBFileID', 'Chain', 'ResidueID'])
                count += 1
        #print(count)


    def _add_pdb_rcsb_ligands(self, tsession, database_pdb_id, pdb_object = None, ligand_mapping = {}, ligand_params_files = {}):
        '''This function associates the ligands of a PDB file (which may be arbitrarily named) with ligands entered in
           the database using the ligand's PDB code. The insertion is handled by a transaction which should be set up
           by the caller.
           Touched tables:
               PDBLigand
               Other Ligand tables by proxy (via add_ligand_by_pdb_code)
        '''

        assert(not(ligand_params_files)) # write this functionality when needed. First, add a new file using FileContent. Next, associate this file with the PDBLigand record.

        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)
        db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id))

        db_ligand_ids = {}
        if db_record.FileSource == 'RCSB':
            # This structure came straight from the RCSB. We trust the ligand codes to be correct
            assert(not(ligand_mapping))

            # Add the ligand description using data from the RCSB if they do not exist.
            for ligand_code in pdb_object.get_ligand_codes():
                db_ligand_ids[ligand_code] = self.add_ligand_by_pdb_code(ligand_code)
        else:
            raise Exception('implement this - try S9G10_best')
            # This structure is not from the RCSB and may use non-standard ligand codes.
            # We therefore require a mapping from all ligand codes in the PDB to RCSB ligand codes.

            ligand_codes = pdb_object.get_ligand_codes()

            # Check all codes have a mapping
            for ligand_code in ligand_codes:
                assert(ligand_code in ligand_mapping)

            # Check whether any codes exist in the mapping which have corresponding Ion records in the database.
            # Since non-RCSB files may be missing headers and we cannot assume an heterogen is an ion purely due to
            # the number of ATOM records totaling one record (e.g. missing coordinates), we need to make this check.
            # Note: this check assumes that the sets of ligand and ion PDB codes are mutually exclusive.
            for ligand_code in ligand_codes:
                existing_db_record = tsession.query(DBIon).filter(DBIon.PDBCode == ligand_code)
                if existing_db_record.count() > 0:
                    raise Exception('Handle this case and add a PDBIon record below instead.')

            # Add the ligands using the RCSB code
            for ligand_code in ligand_codes:
                db_ligand_ids[ligand_code] = self.add_ligand_by_pdb_code(ligand_mapping[ligand_code])

        # Record all instances of ligands in the PDB file (add PDBLigand records).
        # PDBLigandCode is the code used by the PDB file regardless of whether the structure came from the RCSB i.e. it
        # may not be the standard code. The standard code can be found by looking up the associated Ligand record.
        for chain_id, chain_ligands in sorted(pdb_object.ligands.iteritems()):
            for het_seq_id, lig in sorted(chain_ligands.iteritems()):
                try:
                    assert(lig.PDBCode in db_ligand_ids)
                    db_ion = get_or_create_in_transaction(tsession, PDBLigand, dict(
                        PDBFileID = database_pdb_id,
                        Chain = chain_id,
                        SeqID = het_seq_id,
                        PDBLigandCode = lig.PDBCode,
                        LigandID = db_ligand_ids[lig.PDBCode],
                        ParamsFileContentID = None,
                    ))
                except Exception, e:
                    colortext.error(str(e))
                    colortext.error(traceback.format_exc())
                    raise Exception('An exception occurred committing ligand "{0}" from {1} to the database.'.format(lig.PDBCode, database_pdb_id))


    def _add_pdb_rcsb_ions(self, tsession, database_pdb_id, pdb_object = None):
        '''This function associates the ions of a PDB file with ions entered in the database using PDB codes from RCSB
           PDB files.
           For simplicity, we make the assumption that ion codes in all PDB files are not modified from the original PDB file.
           We support this assumption with a couple of checks:
             - we check that the elemental code is the same for the ion's atom and for the database record;
             - we only allow the addition of Ion records from RCSB PDB files. Since we require that the RCSB PDB file be
               added prior to adding derived/designed structures, any ions should have a corresponding record in the database
               unless: i) the ion was manually or otherwise added to the structure; or ii) the ion code was indeed changed.

           The insertion is handled by a transaction which should be set up
           by the caller.
           Touched tables:
               PDBIon
               Ion
        '''

        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)
        db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id))

        # Create a set of Ion records from the PDB file. ions maps PDB codes to dicts containing Ion table fields
        ions = {}
        for c, cions in pdb_object.ions.iteritems():
            for seq_id, pdb_ion_object in cions.iteritems():
                pprint.pprint(pdb_ion_object.__dict__)
                pprint.pprint(pdb_ion_object.get_db_records(database_pdb_id))
                if not ions.get(pdb_ion_object.PDBCode):
                    ions[pdb_ion_object.PDBCode] = copy.deepcopy(pdb_ion_object.get_db_records(database_pdb_id)['Ion'])
                else:
                    # Make sure that all ions in the PDB file have the same formula, description, etc.
                    subsequent_instance = copy.deepcopy(pdb_ion_object.get_db_records(database_pdb_id)['Ion'])
                    for k, v in ions[pdb_ion_object.PDBCode].iteritems():
                        assert(v == subsequent_instance[k])
        colortext.warning(pprint.pformat(ions))
        # Make sure that the ions in the PDB file have the same formula, description, etc. as currently in the database
        existing_ion_codes = set()
        for pdb_code, d in ions.iteritems():
            colortext.pcyan(pdb_code)
            existing_db_record = tsession.query(DBIon).filter(DBIon.PDBCode == pdb_code)
            colortext.pcyan(d)
            if existing_db_record.count() > 0:
                assert(existing_db_record.count() == 1)
                existing_db_record = existing_db_record.one()
                colortext.pcyan(existing_db_record.__dict__)
                pprint.pprint(ions[pdb_ion_object.PDBCode])
                if ions[pdb_ion_object.PDBCode]['PDBCode'] == existing_db_record.PDBCode:
                    assert(ions[pdb_ion_object.PDBCode]['Description'] == existing_db_record.Description) # This can differ e.g. CL in 127L is 3(CL 1-) since there are 3 ions but in PDB files with 2 ions, this can be 2(CL 1-). We can assert this if we do extra parsing.
                    assert(ions[pdb_ion_object.PDBCode]['Formula'] == existing_db_record.Formula) # This can differ e.g. CL in 127L is 3(CL 1-) since there are 3 ions but in PDB files with 2 ions, this can be 2(CL 1-). We can assert this if we do extra parsing.
                existing_ion_codes.add(pdb_code)
            elif db_record.FileSource != 'RCSB':
                raise Exception('We cannot add ion "{0}" from this file as there is no Ion record in the database. This ion should have been added by an underlying RCSB file.'.format(pdb_code))

        # Create the main Ion records, only creating records for ions in RCSB files.
        if db_record.FileSource == 'RCSB':
            for pdb_code, ion_record in ions.iteritems():
                if pdb_code not in existing_ion_codes:
                    # Do not add existing records
                    colortext.message('Adding ion {0}'.format(pdb_code))
                    pprint.pprint(ion_record)
                    db_ion = get_or_create_in_transaction(tsession, DBIon, ion_record, missing_columns = ['ID'])

        # Get the mapping from PDB code to Ion objects
        db_ions = {}
        for pdb_code in ions.keys():
            existing_db_record = tsession.query(DBIon).filter(DBIon.PDBCode == pdb_code)
            assert(existing_db_record.count() == 1)
            db_ions[pdb_code] = existing_db_record.one()

        # Record all instances of ions in the PDB file (add PDBIon records).
        for c, cions in pdb_object.ions.iteritems():
            for seq_id, pdb_ion_object in cions.iteritems():

                assert(pdb_ion_object.get_db_records(None)['Ion']['PDBCode'] == db_ions[pdb_ion_object.PDBCode].PDBCode)
                #assert(pdb_ion_object.get_db_records(None)['Ion']['Formula'] == db_ions[pdb_ion_object.PDBCode].Formula) # not always true e.g. see CL comment above. We can assert this if we do extra parsing.
                assert(pdb_ion_object.get_db_records(None)['Ion']['Description'] == db_ions[pdb_ion_object.PDBCode].Description)
                pdb_ion_record = pdb_ion_object.get_db_records(database_pdb_id, ion_id = db_ions[pdb_ion_object.PDBCode].ID)['PDBIon']
                try:
                    db_ion = get_or_create_in_transaction(tsession, PDBIon, pdb_ion_record)
                except Exception, e:
                    colortext.error(str(e))
                    colortext.error(traceback.format_exc())
                    raise Exception('An exception occurred committing ion "{0}" from {1} to the database.'.format(pdb_ion_object.get_db_records(None)['Ion']['PDBCode'], database_pdb_id))


    def _add_pdb_uniprot_mapping(self, tsession, database_pdb_id, pdb_object = None):
        '''UniProtACs have forms like 'P62937' whereas UniProtIDs have forms like 'PPIA_HUMAN.'''
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)
        return
        #protein = None
        #UniProtAC = None
        #UniProtID = None
        # todo: add UniProt mapping here. The old approach was flawed.
        #ref_pdb_id = self.dict['DerivedFrom'] or self.dict['ID']
        #if ref_pdb_id not in self.NoUniProtIDs:
        #    read_UniProt_map(self.ddGdb)
        #    if not PDBToUniProt.get(ref_pdb_id):
        #        if not (self.UniProtAC and self.UniProtID):
        #            ACtoID_mapping, PDBtoAC_mapping = None, None
        #            try:
        #                getUniProtMapping(ref_pdb_id, storeInDatabase = False, ddGdb = self.ddGdb)
        #                if not (PDBtoAC_mapping and ACtoID_mapping):
        #                    raise Exception("Could not find a UniProt mapping for %s in %s." % (ref_pdb_id, uniprotmapping))
        #            except:
        #                colortext.error("Could not find a UniProt mapping for %s in %s." % (ref_pdb_id, uniprotmapping))
        #            self.ACtoID_mapping = ACtoID_mapping
        #            self.PDBtoAC_mapping = PDBtoAC_mapping
        #
        # Add UniProt mapping
        #if self.UniProtAC and self.UniProtID:
        #    results = self.ddGdb.locked_execute("SELECT * FROM UniProtKBMapping WHERE UniProtKB_AC=%s", parameters=(self.UniProtAC,))
        #    if results:
        #        if results[0]['PDBFileID'] != d['ID']:
        #            raise Exception("Existing UniProt mapping (%s->%s) does not agree with the passed-in parameters (%s->%s)." % (results[0]['UniProtKB_AC'],results[0]['PDBFileID'],self.UniProtAC,d['ID']))
        #    else:
        #        UniProtPDBMapping = dict(
        #            UniProtKB_AC = self.UniProtAC,
        #            PDBFileID = d[FieldNames_.PDB_ID],
        #        )
        # Store the UniProt mapping in the database
        #ref_pdb_id = self.dict['DerivedFrom'] or self.dict['ID']
        #if ref_pdb_id not in self.NoUniProtIDs:
        #    if not (self.ACtoID_mapping and self.PDBtoAC_mapping):
        #        try:
        #            self.ACtoID_mapping, self.PDBtoAC_mapping = getUniProtMapping(ref_pdb_id, storeInDatabase = testonly)
        #        except:
        #            if False and self.dict['Techniques'] != 'Rosetta model':
        #                raise
        #    if False and self.dict['Techniques'] != 'Rosetta model':
        #        assert(self.ACtoID_mapping and self.PDBtoAC_mapping)
        #    if not testonly:
        #        if self.ACtoID_mapping and self.PDBtoAC_mapping:
        #            commitUniProtMapping(self.ddGdb, self.ACtoID_mapping, self.PDBtoAC_mapping)


    def _add_pdb_publication(self, tsession, database_pdb_id, pdb_object = None):
        '''Extracts the PDB source information.
           Touched tables:
               Publication
               PublicationIdentifier
        '''

        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)
        pdb_record = tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id).one()
        PUBTYPES = ['ISSN', 'ESSN']

        j = pdb_object.get_journal()
        if not j:
            return
        database_pdb_id = database_pdb_id.strip()

        # We identify the sources for a PDB identifier with that identifier
        PublicationID = "PDB:%s" % database_pdb_id
        publication_record = get_or_create_in_transaction(tsession, Publication, dict(ID = PublicationID), only_use_supplied_columns = True)
        pdb_record.Publication = publication_record.ID
        tsession.flush()

        locations = tsession.query(PublicationIdentifier).filter(PublicationIdentifier.SourceID == publication_record.ID)
        pub_locations = [location for location in locations if location.Type in PUBTYPES]
        doi_locations = [location for location in locations if location.Type == 'DOI']
        assert(len(pub_locations) <= 1)
        assert(len(doi_locations) <= 1)
        if j["published"]:
            if pub_locations:
                location = pub_locations[0]
                if j["REFN"]["type"] == location.Type:
                    if j["REFN"]["ID"] != location.ID:
                        colortext.warning("REFN: Check that the PublicationIdentifier data ('{0}, {1}') matches the PDB REFN data ({2}).".format(location.ID, location.Type, j["REFN"]))
            elif j.get('REFN'):
                assert(j["REFN"]["type"] in PUBTYPES)
                db_pub_id = get_or_create_in_transaction(tsession, PublicationIdentifier, dict(
                    SourceID    = PublicationID,
                    ID          = j["REFN"]["ID"],
                    Type        = j["REFN"]["type"],
                ))
        if j["DOI"]:
            if doi_locations:
                location = doi_locations[0]
                if j["DOI"] != location.ID:
                    colortext.warning("DOI: Check that the PublicationIdentifier data  ('{0}, {1}') matches the PDB DOI data ({2}).".format(location.ID, location.Type, j["DOI"]))
            else:
                db_pub_id = get_or_create_in_transaction(tsession, PublicationIdentifier, dict(
                    SourceID    = PublicationID,
                    ID          = j["DOI"],
                    Type        = "DOI",
                ))


    def add_designed_pdb_file(self, designed_pdb_filepath, design_pdb_id, original_pdb_id, description, username, chain_mapping = {}, ligand_mapping = {}, previously_added = set()):
        '''Wrapper for add_designed_pdb.'''
        return self.add_designed_pdb(self, PDB.from_filepath(designed_pdb_filepath), design_pdb_id, original_pdb_id, description, username, chain_mapping = chain_mapping, ligand_mapping = ligand_mapping, previously_added = previously_added)


    def add_designed_pdb(self, designed_pdb_object, design_pdb_id, original_pdb_id, description, username, chain_mapping = {}, ligand_mapping = {}, previously_added = set()):
        assert(5 <= len(design_pdb_id) <= 10)

        # I used to have a trust_database_content boolean parameter but am not sure whether or not we should use it
        colortext.message('Adding designed PDB file {0} based off {1}.'.format(design_pdb_id, original_pdb_id))

        colortext.pcyan('Adding the original PDB file using a separate transaction.')
        self.add_pdb_from_rcsb(original_pdb_id, previously_added = previously_added)

        raise Exception('implement this')

        #dbp = PDBFileImporter(self.DDG_db, rootname, contains_membrane_protein = contains_membrane_protein, protein = protein, file_source = file_source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID, testonly = testonly, techniques = techniques, derived_from = derived_from, notes = notes)

        #---
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
        #---


        p = designed_pdb_object
        pdb_ligand_codes = p.get_ligand_codes()
        if True or pdb_ligand_codes and not ligand_mapping:
            d = {}
            for pdb_ligand_code in pdb_ligand_codes: d[pdb_ligand_codes] = pdb_ligand_codes
            sample_mapping = pprint.pformat(d)
            raise Exception('The PDB file contains ligands {0} but no mapping is given from the ligand codes to RCSB ligand codes. We require this mapping as it is not unusual for ligand codes to be renamed e.g. to " X ". An example mapping would be:\nligand_mapping = {1}'.format(pdb_ligand_codes, sample_mapping))

        self.get_pdb_details(design_pdb_id)
        self.get_pdb_details(original_pdb_id)

        # regardless of chain_mapping being specified, run clustal and assert that the sequence for chain c in design_pdb_id and chain_mapping.get(c, c) in original_pdb_id matches >90%

        # start data entry
        # open a transaction, call inner methods, then close the transaction
        # Particularly, call  add_pdb_data which then calls other functions. The point of separating all of these sub-functions into one big inner function is so we have an API to update the information for specific PDB files.





def test():

    echo_sql = False
    importer = DataImportInterface.get_interface_with_config_file(cache_dir = '/kortemmelab/data/oconchus/ddgcache', echo_sql = echo_sql)
    session = importer.session

    #importer.update_pdbs(update_sections = set(['Chains', 'Molecules']), start_at = '1W99', restrict_to_file_source = 'RCSB')
    #importer.update_pdbs(update_sections = set(['Publication']), start_at = None, restrict_to_file_source = 'RCSB')
    importer.update_pdbs(update_sections = set(['Residues', 'Publication']), start_at = None, restrict_to_file_source = 'RCSB')

test()

test_cases = '2FX5' #  at the time of writing, these PDB IDs had no JRNL lines
test_cases = '1GTX', '1UOX', '1WSY', '1YYJ', '2IMM', '2MBP' # at the time of writing, these PDB IDs had no UniProt entries


