#!/usr/bin/python2.4
# encoding: utf-8
"""
import_api.py
High-level functions for importing data into the DDG database.


Example usage:
    # Create an import API instance
    importer = DataImportInterface.get_interface_with_config_file(cache_dir = '/kortemmelab/data/oconchus/ddgcache', echo_sql = False)

    # Access the SQLAlchemy session directly
    session = importer.session

    # Access the MySQLdb interface layer directly
    DDG_db = importer.DDG_db # or importerDDG_db_utf

    # Access an RCSB PDB file to the database
    importer.add_pdb_from_rcsb('1A2K')

    # Access ligand details to the database (note: this will be called by add_pdb_from_rcsb)
    importer.add_ligand_by_pdb_code('GTP')

    # Update certain properties of RCSB files in the database
    importer.update_pdbs(update_sections = set(['Residues', 'Publication']), start_at = None, restrict_to_file_source = 'RCSB')


@todo list:
  - ticket 1489: add_pdb_from_rcsb: FileContent (see ticket 1489)
  - get_pdb_details (see below)
  - PDB Chains (_add_pdb_chains)
    ticket 1463: SCOP/SCOPe classifications
    add self.pfam_api.? call to get mapping from RCSB PDB files to the UniProt sequences. Add this to a separate function, _add_pdb_uniprot_mapping.
    add self.pfam_api.get_pfam_accession_numbers_from_pdb_chain(database_pdb_id, c)) calls. Use this in _add_pdb_uniprot_mapping.
    add self.scope_api.get_chain_details(database_pdb_id, c))) calls
    ticket 1472: add coordinates
  - ticket 1488, 1473: _add_pdb_uniprot_mapping. implement UniProt mapping. The old approach was flawed. Use the new approach (similar to how I used PPComplex as an abstraction layer)
  - ticket 1493: add b-factors and DSSP for modified residues
  - ticket 1491: electron densities

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

from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation
from klab.fs.fsio import read_file, write_temp_file, open_temp_file, write_file
from klab.bio.pfam import Pfam
from klab.bio.dssp import MonomerDSSP, ComplexDSSP, MissingAtomException
from klab.bio.ligand import Ligand, PDBLigand, LigandMap
from klab.bio.pdbtm import PDBTM
from klab.bio.clustalo import PDBChainSequenceAligner
from klab.db.sqlalchemy_interface import get_single_record_from_query, get_or_create_in_transaction
from klab.bio import rcsb
from klab.general.strutil import remove_trailing_line_whitespace
from klab.hash.md5 import get_hexdigest

from db_schema import test_schema_against_database_instance
from db_schema import PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue, LigandDescriptor, LigandIdentifier, LigandSynonym, LigandPrice, LigandReference, PDBLigand, PDBLigandFile, PDBIon, FileContent, PDB2PDBChainMap
from db_schema import Ligand as DBLigand
from db_schema import Ion as DBIon
from db_schema import User as DBUser
from db_schema import Publication, PublicationAuthor, PublicationIdentifier, DeclarativeBase
from api_layers import *
import ddgdbapi

try:
    import magic
except ImportError:
    colortext.error('Failed to import magic package. This failure will prevent you from being able to import new file content into the database.')
    pass


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


    def __init__(self, passwd, connect_string, connect_string_utf, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, cache_dir = None, echo_sql = False, port = 3306):
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
        self.DDG_db = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname, port = port)
        self.DDG_db_utf = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname, use_utf = True, port = port)

        self.cache_dir = cache_dir
        self.echo_sql = echo_sql

        test_schema_against_database_instance(self.DDG_db)

        if self.cache_dir:
            self.initialize_cache_directory()
        else:
            colortext.warning('Warning: No cache directory has been specified in your configuration file e.g.\n  cache_dir = /kortemmelab/data/username/ddgcache\nThis may result in files being retrieved from the RCSB servers multiple times.')

        # Set up SQLAlchemy connections
        self.connect_string = connect_string
        self.connect_string_utf = connect_string_utf
        self.engine, self.session = None, None
        self.engine_utf, self.session_utf = None, None
        self.get_engine(utf = False)
        self.get_engine(utf = True)
        self.get_session(utf = False)
        self.get_session(utf = True)

        self.rosetta_scripts_path = rosetta_scripts_path
        self.rosetta_database_path = rosetta_database_path

        # Parse PDB chain -> Pfam mapping
        self.pfam_api = Pfam()


    @classmethod
    def get_interface_with_config_file(cls, database = 'ddg', host_config_name = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None, my_cnf_path = None, cache_dir = None, echo_sql = False, port = 3306):
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
        connection_string_key_utf = 'sqlalchemy.{0}.url.utf'.format(database)
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
                        elif key == 'cache_dir':
                            cache_dir = val
                        elif key == 'host':
                            host = val
                        elif key == 'port':
                            port = int(val)
                        elif key == connection_string_key:
                            connection_string = val
                        elif key == connection_string_key_utf:
                            connect_string_utf = val
                    else:
                        parsing_config_section = False

        if not user or not password or not host or not connection_string or not connect_string_utf:
            raise Exception("Couldn't find host(%s), username(%s), password, or connection string in section %s in %s" % (host, user, host_config_name, my_cnf_path) )

        return cls(password, connection_string, connect_string_utf, username = user, hostname = host, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, cache_dir = cache_dir, echo_sql = echo_sql, port = port)


    def __del__(self):
        pass         #self.DDG_db.close()         #self.ddGDataDB.close()



    #############################
    #                           #
    #  SQLAlchemy Engine setup  #
    #                           #
    #############################



    def get_engine(self, utf = False):
        if utf:
            if not self.engine_utf:
                self.engine_utf = create_engine(self.connect_string_utf, echo = self.echo_sql)
            return self.engine_utf
        else:
            if not self.engine:
                self.engine = create_engine(self.connect_string, echo = self.echo_sql)
            return self.engine


    def get_connection(self, utf = False):
        # e.g. connection = importer.get_connection(); connection.execute("SELECT * FROM User")
        engine = self.get_engine(utf = utf)
        return engine.connect()


    def get_session(self, new_session = False, autoflush = True, autocommit = False, utf = False):
        engine = self.get_engine(utf = utf)
        if new_session or ((not(utf) and not(self.session)) or (utf and not(self.session_utf))):
            maker_ddgdatabase = sessionmaker(autoflush = autoflush, autocommit = autocommit)
            s = scoped_session(maker_ddgdatabase)
            DeclarativeBaseDDG = declarative_base()
            metadata_ddg = DeclarativeBaseDDG.metadata
            s.configure(bind=engine)
            metadata_ddg.bind = engine
            if new_session:
                return s
            else:
                if utf:
                    self.session_utf = s
                else:
                    self.session = s
        if utf:
            return self.session_utf
        else:
            return self.session


    def renew(self, utf = False):
        self.session = self.get_session(new_session = True, utf = utf)


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


    def update_pdbs(self, pdb_ids = [], update_sections = set(), start_at = None, restrict_to_file_source = None, pdb_ligand_params_files = {}):
        '''Updates all or selected data for all or selected PDB files in the database, enumerating in alphabetical order.
           If start_at is specified, the update begins at the specified PDB identifier.
           pdb_ligand_params_files should be a mapping from pdb_id -> ligand_code -> params_file_path
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
                    ligand_params_file_paths = pdb_ligand_params_files.get(pdb_id, {})
                    self.add_pdb_data(tsession, pdb_id, update_sections = update_sections, ligand_params_file_paths = ligand_params_file_paths)
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


    ###########################################################################################
    ## File management layer
    ##
    ## This part of the API is responsible for file content abstraction
    ###########################################################################################


    def get_file_id(self, content, tsession = None, hexdigest = None):
        '''Searches the database to see whether the FileContent already exists. The search uses the digest and filesize as
           heuristics to speed up the search. If a file has the same hex digest and file size then we do a straight comparison
           of the contents.
           If the FileContent exists, the value of the ID field is returned else None is returned.
           '''
        tsession = tsession or self.get_session()
        existing_filecontent_id = None
        hexdigest = hexdigest or get_hexdigest(content)
        filesize = len(content)
        for r in tsession.query(FileContent).filter(and_(FileContent.MD5HexDigest == hexdigest, FileContent.Filesize == filesize)):
            if r.Content == content:
                assert(existing_filecontent_id == None) # content uniqueness check
                existing_filecontent_id = r.ID
        return existing_filecontent_id


    @informational_file
    def get_file_id_using_old_interface(self, content, db_cursor = None, hexdigest = None):
        '''Searches the database to see whether the FileContent already exists. The search uses the digest and filesize as
           heuristics to speed up the search. If a file has the same hex digest and file size then we do a straight comparison
           of the contents.
           If the FileContent exists, the value of the ID field is returned else None is returned.
           '''
        # @todo: This function is deprecated. Use get_file_id instead.

        existing_filecontent_id = None
        hexdigest = hexdigest or get_hexdigest(content)
        filesize = len(content)
        if db_cursor:
            db_cursor.execute('SELECT * FROM FileContent WHERE MD5HexDigest=%s AND Filesize=%s', (hexdigest, filesize))
            results = db_cursor.fetchall()
        else:
            results = self.DDG_db.execute_select('SELECT * FROM FileContent WHERE MD5HexDigest=%s AND Filesize=%s', parameters=(hexdigest, filesize))
        for r in results:
            if r['Content'] == content:
                assert(existing_filecontent_id == None) # content uniqueness check
                existing_filecontent_id = r['ID']
        return existing_filecontent_id


    def _add_file_content(self, file_content, tsession = None, rm_trailing_line_whitespace = False, forced_mime_type = None):
        '''Takes file file_content (and an option to remove trailing whitespace from lines e.g. to normalize PDB files), adds
           a new record if necessary, and returns the associated FileContent.ID value.'''

        tsession = tsession or self.get_session()
        if rm_trailing_line_whitespace:
            file_content = remove_trailing_line_whitespace(file_content)

        # Check to see whether the file has been uploaded before
        hexdigest = get_hexdigest(file_content)
        existing_filecontent_id = self.get_file_id(file_content, tsession = tsession, hexdigest = hexdigest)

        # Create the FileContent record if the file is a new file
        if existing_filecontent_id == None:

            # Determing the MIME type
            mime_type = forced_mime_type
            if not mime_type:
                # Note: in case the wrong mime-types are being returned, try saving to file first and then calling magic.from_file.
                # See commit c62883b58649bd813bf022f7d1193abb06f1676d for the code. This used to be necessary for some odd reason.
                mime_type = magic.from_buffer(file_content, mime = True)

            # Create the database record
            file_content_record = get_or_create_in_transaction(tsession, FileContent, dict(
                Content = file_content,
                MIMEType = mime_type,
                Filesize = len(file_content),
                MD5HexDigest = hexdigest
            ), missing_columns = ['ID'])
            existing_filecontent_id = file_content_record.ID

        assert(existing_filecontent_id != None)
        return existing_filecontent_id


    def _add_file_content_using_old_interface(self, file_content, db_cursor = None, rm_trailing_line_whitespace = False, forced_mime_type = None):
        '''Takes file content (and an option to remove trailing whitespace from lines e.g. to normalize PDB files), adds
           a new record if necessary, and returns the associated FileContent.ID value.'''

        # todo: This function is deprecated. Use _add_file_content instead.

        if rm_trailing_line_whitespace:
            file_content = remove_trailing_line_whitespace(file_content)

        # Check to see whether the file has been uploaded before
        hexdigest = get_hexdigest(file_content)
        existing_filecontent_id = self.get_file_id_using_old_interface(file_content, db_cursor = db_cursor, hexdigest = hexdigest)

        # Create the FileContent record if the file is a new file
        if existing_filecontent_id == None:

            # Determing the MIME type
            mime_type = forced_mime_type
            if not mime_type:
                # Note: in case the wrong mime-types are being returned, try saving to file first and then calling magic.from_file.
                # See commit c62883b58649bd813bf022f7d1193abb06f1676d for the code. This used to be necessary for some odd reason.
                mime_type = magic.from_buffer(file_content, mime = True)

            # Create the database record
            d = dict(
                Content = file_content,
                MIMEType = mime_type,
                Filesize = len(file_content),
                MD5HexDigest = hexdigest
            )
            if db_cursor:
                sql, params, record_exists = self.DDG_db.create_insert_dict_string('FileContent', d, ['Content'])
                db_cursor.execute(sql, params)
            else:
                self.DDG_db.insertDictIfNew('FileContent', d, ['Content'], locked = False)
            existing_filecontent_id = self.get_file_id(file_content, db_cursor = db_cursor, hexdigest = hexdigest)

        assert(existing_filecontent_id != None)
        return existing_filecontent_id


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
                p = PDB(record['Content'], parse_ligands = True)
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


    def get_rcsb_record(self, pdbfile_db_record, tsession = None):
        '''pdbfile_db_record should be a db_schema.py::PDBFile object.
           Winds up the 'derived from' tree to find the RCSB file that this file originated from.
           Throws an exception if there are no such files.
           This is useful for a number of reasons:
             - Looking up the resolution of the original structure, the determination technique, and its b-factors
               We do not copy this information into derived structures as it is generally meaningless (e.g. for PDB_REDO
               structures or structures minimized or repacked with some force-field).
             - Determining the name of the molecules in the derived PDB file.
           etc.
        '''
        if not tsession:
            tsession = self.get_session()
        try:
            c = 0
            while (pdbfile_db_record.DerivedFrom) and (pdbfile_db_record.FileSource != 'RCSB') and (c < 40): # the last expression should be unnecessary but just in case...
                pdbfile_db_record = tsession.query(PDBFile).filter(PDBFile.ID == pdbfile_db_record.DerivedFrom).one()
                c += 1
            assert(pdbfile_db_record.FileSource == 'RCSB')
            return pdbfile_db_record
        except Exception, e:
            raise Exception('Failed to retrieve an RCSB record corresponding to "{0}".'.format(pdbfile_db_record.ID))


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
        if not os.path.exists(pdb_dir):
            raise colortext.Exception('The cache directory "{0}" could not be created.'.format(self.cache_dir))


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


    def add_pdb_from_rcsb(self, pdb_id, previously_added = set(), update_sections = set(), trust_database_content = False, ligand_params_file_paths = {}, debug = False):
        '''NOTE: This API is used to create and analysis predictions or retrieve information from the database.
                 This function adds new raw data to the database and does not seem to belong here. It should be moved into
                 an admin API instead.
           This function adds imports a PDB into the database, creating the associated molecule, chain and residue etc. records.
           If previously_added contains pdb_id then we return. Otherwise, we step through the full PDB file import. This
           is added purely for optimistic efficiency e.g. if were to add 100 designed files based off the same RCSB PDB,
           we would not want to run the full import code for the RCSB PDB file more than once.

           If trust_database_content is True then we return if we find a PDBFile record without delving into the related tables.
           This is useful in certain circumstances e.g. when adding derived PDB files using add_designed_pdb.

           Touched tables:
               PDBFile
               todo: FileContent (see ticket 1489)
        '''

        assert(not ligand_params_file_paths) # todo: handle this case when we need to. See the implementation for add_designed_pdb

        if pdb_id in previously_added:
            return pdb_id
        assert(len(pdb_id) == 4)
        pdb_id = pdb_id.upper()

        # RCSB files should be a straightforward case so we will use a new session and commit all changes
        tsession = self.get_session(new_session = True)
        try:
            pdb_object = None
            db_record_object = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == pdb_id))
            is_new_record = db_record_object == None
            if not is_new_record:
                print('Retrieving {0} from database.'.format(pdb_id))
                pdb_object = PDB(db_record_object.Content, parse_ligands = True)
                assert(db_record_object.FileSource == 'RCSB')
                if trust_database_content:
                    print('Trusting the existing data and early-outing.')
                    return pdb_id
            else:
                db_record_object = PDBFile()
                contents = self._retrieve_pdb_contents(pdb_id)
                pdb_object = PDB(contents, parse_ligands = True)
                db_record_object.ID = pdb_id
                db_record_object.FileSource = 'RCSB'
                db_record_object.Content = contents
                update_sections = set() # add all related data

            # Fill in the FASTA, Resolution, Techniques, Transmembrane, and b-factor fields of a PDBFile record
            self._add_pdb_file_information(db_record_object, pdb_object, pdb_id, is_new_record, is_rcsb_pdb = True)

            # Add a new PDBFile record
            if is_new_record:
                db_record_object.UserID = None
                db_record_object.Notes = None
                db_record_object.DerivedFrom = None
                tsession.add(db_record_object)
            tsession.flush()

            # Publication
            if is_new_record or (not db_record_object.Publication):
                self._add_pdb_publication(tsession, db_record_object.ID, pdb_object = pdb_object)

            # add all other data
            self.add_pdb_data(tsession, pdb_id, update_sections = update_sections, ligand_params_file_paths = ligand_params_file_paths)

            previously_added.add(pdb_id)

            print('Success.\n')
            if debug:
                print('Debug call - rolling back the transaction.\n')
                tsession.rollback()
            else:
                tsession.commit()
            tsession.close()
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise


    def _add_pdb_file_information(self, db_record_object, pdb_object, pdb_id, is_new_record, is_rcsb_pdb = False):
        '''Fills in the FASTA, Resolution, Techniques, Transmembrane, and b-factor fields of a PDBFile record.'''

        # Checks
        pdb_atom_chains = pdb_object.atom_chain_order
        if not pdb_atom_chains:
            raise Exception('No ATOM chains were found in the PDB file.')
        foundRes = pdb_object.CheckForPresenceOf(["CSE", "MSE"])
        if foundRes:
            colortext.error("The PDB %s contains residues which could affect computation (%s)." % (pdb_id, ', '.join(foundRes, )))
            if "CSE" in foundRes:
                colortext.error("The PDB %s contains CSE. Check." % pdb_id)
            if "MSE" in foundRes:
                colortext.error("The PDB %s contains MSE. Check." % pdb_id)

        # FASTA
        if is_new_record or (not db_record_object.FASTA):
            db_record_object.FASTA = pdb_object.create_fasta()

        # Resolution
        if is_new_record or (not db_record_object.Resolution):
            resolution = pdb_object.get_resolution()
            if not resolution:
                colortext.error("Could not determine resolution for {0}.".format(pdb_id))
            if resolution == "N/A":
                resolution = None
            db_record_object.Resolution = resolution

        # Techniques
        if is_rcsb_pdb:
            if is_new_record or (not db_record_object.Techniques):
                db_record_object.Techniques = pdb_object.get_techniques()

        # Transmembrane
        if is_rcsb_pdb:
            if is_new_record or (db_record_object.Transmembrane == None):
                pdbtmo = PDBTM(read_file('/kortemmelab/shared/mirror/PDBTM/pdbtmall.xml')) # todo: set this up as a parameter in the configuration file
                pdb_id_map = pdbtmo.get_pdb_id_map()
                uc_pdb_id_map = {}
                for k, v in pdb_id_map.iteritems():
                    uc_pdb_id_map[k.upper()] = v
                if pdb_id in uc_pdb_id_map:
                    colortext.warning('{0} is a transmembrane protein.'.format(pdb_id))
                db_record_object.Transmembrane = pdb_id in pdb_id_map

        # B-factors
        overall_bfactors = pdb_object.get_B_factors().get('Overall')
        db_record_object.BFactorMean = overall_bfactors['mean']
        db_record_object.BFactorDeviation = overall_bfactors['stddev']


    def is_session_utf(self, tsession):
        return str(tsession.bind.url).lower().find('utf') != -1


    def get_pdb_object(self, database_pdb_id, tsession = None):
        '''Create a PDB object from content in the database.'''
        if tsession:
            assert(not(self.is_session_utf(tsession)))
        tsession = tsession or self.get_session()
        assert(not(self.is_session_utf(tsession)))
        db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id))
        assert(db_record)
        return PDB(db_record.Content, parse_ligands = True)


    def add_pdb_data(self, tsession, database_pdb_id, update_sections = set(), ligand_mapping = {}, chain_mapping = {}, ligand_params_file_paths = {}):
        '''
            The point of separating the data entry into these sub-functions and calling them from this function is so we
            have an API to update the information for specific PDB files.

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
            self._add_pdb_chains(tsession, database_pdb_id, pdb_object, chain_mapping = chain_mapping)
        if not(update_sections) or ('Molecules' in update_sections):
            colortext.warning('*** Molecules ***')
            self._add_pdb_molecules(tsession, database_pdb_id, pdb_object, chain_mapping = chain_mapping)
        if not(update_sections) or ('Residues' in update_sections):
            colortext.warning('*** Residues ***')
            self._add_pdb_residues(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Ligands' in update_sections):
            colortext.warning('*** Ligands ***')
            self._add_pdb_rcsb_ligands(tsession, database_pdb_id, pdb_object, ligand_mapping, ligand_params_file_paths = ligand_params_file_paths)
        if not(update_sections) or ('Ions' in update_sections):
            colortext.warning('*** Ions ***')
            self._add_pdb_rcsb_ions(tsession, database_pdb_id, pdb_object)
        #if not(update_sections) or ('UniProt' in update_sections):
        #    colortext.warning('*** UniProt ***')
        #    self._add_pdb_uniprot_mapping(tsession, database_pdb_id, pdb_object)


    def _add_pdb_chains(self, tsession, database_pdb_id, pdb_object = None, chain_mapping = {}):
        '''
           Touched tables:
               PDBChain
        '''

        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)

        db_chains = {}
        for r in tsession.query(PDBChain).filter(PDBChain.PDBFileID == database_pdb_id).order_by(PDBChain.Chain):
            db_chains[r.Chain] = r

        db_chains_ids = sorted(db_chains.keys())
        chain_ids = sorted(set(pdb_object.seqres_sequences.keys() + pdb_object.atom_sequences.keys() + pdb_object.chain_types.keys()))

        # Sanity checks for derived structures
        self._check_derived_record_against_rcsb_record(tsession, database_pdb_id, pdb_object, chain_mapping)

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
            for r in tsession.query(PDBChain).filter(PDBChain.PDBFileID == database_pdb_id):
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


    def _add_pdb_molecules(self, tsession, database_pdb_id, pdb_object = None, allow_missing_molecules = False, chain_mapping = {}):
        '''
           Add PDBMolecule and PDBMoleculeChain records
           Touched tables:
               PDBMolecule
               PDBMoleculeChain
        '''

        assert(allow_missing_molecules == False) # todo: do we ever use allow_missing_molecules? We can inspect that case when it presents itself

        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id, tsession = tsession)

        # Sanity checks for derived structures
        self._check_derived_record_against_rcsb_record(tsession, database_pdb_id, pdb_object, chain_mapping)

        if not(chain_mapping):
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

                # Add the PDBMolecule record
                db_molecule = get_or_create_in_transaction(tsession, PDBMolecule, md)

                # Add the PDBMoleculeChain records
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
        else:
            # Copy the molecule information from the original RCSB structure
            # First, get the DB records for the derived structure and the original RCSB structure
            db_record = tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id).one()
            assert(db_record.FileSource != 'RCSB')
            rcsb_record = self.get_rcsb_record(db_record, tsession = tsession)

            # Get the list of RCSB chains
            rcsb_chains = chain_mapping.values()

            # Get the list of PDB molecules associated with chains that are in the derived PDB (accounting for chain renaming)
            rcsb_molecule_chains = {} # a dict mapping RCSB chain IDs to the associated PDBMoleculeChain record
            rcsb_molecule_ids = set()
            for r in tsession.query(PDBMoleculeChain).filter(PDBMoleculeChain.PDBFileID == rcsb_record.ID):
                if r.Chain in rcsb_chains:
                    rcsb_molecule_chains[r.Chain] = r
                    rcsb_molecule_ids.add(r.MoleculeID)

            # Add the PDBMolecule records
            for r in tsession.query(PDBMolecule).filter(PDBMolecule.PDBFileID == rcsb_record.ID):
                if r.MoleculeID in rcsb_molecule_ids:
                    db_molecule = get_or_create_in_transaction(tsession, PDBMolecule, dict(
                        PDBFileID = database_pdb_id,
                        MoleculeID = r.MoleculeID,
                        Name = r.Name,
                        Organism = r.Organism,
                        Fragment = r.Fragment,
                        Synonym = r.Synonym,
                        Engineered = r.Engineered,
                        EC = r.EC,
                        Mutation = r.Mutation,
                        OtherDetails = r.OtherDetails,
                    ))

            # Add the PDBMoleculeChain records
            for derived_chain_id, rcsb_chain_id in sorted(chain_mapping.iteritems()):
                associated_molecule_id = rcsb_molecule_chains[rcsb_chain_id].MoleculeID
                try:
                    db_molecule_chain = get_or_create_in_transaction(tsession, PDBMoleculeChain, dict(
                        PDBFileID = database_pdb_id,
                        MoleculeID = associated_molecule_id,
                        Chain = derived_chain_id
                    ))
                except:
                    if allow_missing_molecules: pass
                    else: raise


    def _add_pdb_residues(self, tsession, database_pdb_id, pdb_object = None):
        '''
           The code here is the same for both RCSB and non-RCSB structures.

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
                        if residue_details['3LC'] != 'UNK': # todo: should we handle these residues?
                            chain_residue_id = chain_id + residue_id
                            assert(chain_residue_id in parsed_residues)
                            monomeric_records[chain_residue_id] = residue_details
        if dssp_complex_d:
            for chain_id, mapping in dssp_complex_d:
                if pdb_object.chain_types[chain_id] == 'Protein' or pdb_object.chain_types[chain_id] == 'Protein skeleton':
                    for residue_id, residue_details in sorted(mapping.iteritems()):
                        if residue_details['3LC'] != 'UNK': # todo: should we handle these residues?
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
                    ), missing_columns = ['ID'])
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


    def _add_pdb_rcsb_ligands(self, tsession, database_pdb_id, pdb_object = None, ligand_mapping = {}, ligand_params_file_paths = {}):
        '''This function associates the ligands of a PDB file (which may be arbitrarily named) with ligands entered in
           the database using the ligand's PDB code. The insertion is handled by a transaction which should be set up
           by the caller.
           Touched tables:
               PDBLigand
               Other Ligand tables by proxy (via add_ligand_by_pdb_code)
        '''

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
            # This structure is not from the RCSB and may use non-standard ligand codes.
            # We therefore require a mapping from all ligand codes in the PDB to RCSB ligand codes.
            ligand_codes = pdb_object.get_ligand_codes()
            if ligand_codes:
                assert(ligand_mapping)

            # Check all codes have a mapping and that the codomain values already exist in the database (the underlying
            # RCSB file and its ligands should already have been added)
            for ligand_code in ligand_codes:
                assert(ligand_mapping.map_code(ligand_code))
                db_ligand_record = tsession.query(DBLigand).filter(DBLigand.PDBCode == ligand_mapping.map_code(ligand_code)).one()
                db_ligand_ids[ligand_code] = db_ligand_record.ID

            # Check whether any codes exist in the mapping which have corresponding Ion records in the database.
            # Since non-RCSB files may be missing headers and we cannot assume an heterogen is an ion purely due to
            # the number of ATOM records totaling one record (e.g. missing coordinates), we need to make this check.
            # Note: this check assumes that the sets of ligand and ion PDB codes are mutually exclusive.
            for ligand_code in ligand_codes:
                existing_db_record = tsession.query(DBIon).filter(DBIon.PDBCode == ligand_code)
                if existing_db_record.count() > 0:
                    raise Exception('Handle this case and add a PDBIon record below instead.')

        # Record all instances of ligands in the PDB file (add PDBLigand records).
        # PDBLigandCode is the code used by the PDB file regardless of whether the structure came from the RCSB i.e. it
        # may not be the standard code. The standard code can be found by looking up the associated Ligand record.
        for chain_id, chain_ligands in sorted(pdb_object.ligands.iteritems()):
            for het_seq_id, lig in sorted(chain_ligands.iteritems()):
                try:
                    assert(lig.PDBCode in db_ligand_ids)
                    pdb_ligand = get_or_create_in_transaction(tsession, PDBLigand, dict(
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

        # Params files
        if ligand_params_file_paths:
            if not(0 < max(map(len, ligand_params_file_paths.keys())) <= 3):
                bad_keys = sorted([k for k in ligand_params_file_paths.keys() if len(k) > 3])
                raise colortext.Exception('The ligand codes "{0}" are invalid - all codes must be between 1 and 3 characters e.g. "CIT".'.format('", "'.join(bad_keys)))
        bad_keys = sorted(set(ligand_params_file_paths.keys()).difference(pdb_object.get_ligand_codes()))
        if bad_keys:
            raise colortext.Exception('The ligand codes "{0}" were specified but were not found in the PDB file.'.format('", "'.join(bad_keys)))

        ligand_params_file_content = {}
        if ligand_params_file_paths:

            # Read all params files
            for ligand_code, params_filepath in ligand_params_file_paths.iteritems():
                ligand_params_file_content[ligand_code] = read_file(params_filepath)

            for ligand_code, params_file_content in ligand_params_file_content.iteritems():
                # First, add a new file using FileContent.
                file_content_id = self._add_file_content(params_file_content, tsession = tsession, rm_trailing_line_whitespace = True, forced_mime_type = 'text/plain')

                # Next, associate this file with the PDBLigand record.
                pdb_ligand_file = get_or_create_in_transaction(tsession, PDBLigandFile, dict(
                    PDBFileID = database_pdb_id,
                    PDBLigandCode = lig.PDBCode,
                    ParamsFileContentID = file_content_id,
                ))


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
                if not ions.get(pdb_ion_object.PDBCode):
                    ions[pdb_ion_object.PDBCode] = copy.deepcopy(pdb_ion_object.get_db_records(database_pdb_id)['Ion'])
                else:
                    # Make sure that all ions in the PDB file have the same formula, description, etc.
                    subsequent_instance = copy.deepcopy(pdb_ion_object.get_db_records(database_pdb_id)['Ion'])
                    for k, v in ions[pdb_ion_object.PDBCode].iteritems():
                        assert(v == subsequent_instance[k])

        # Make sure that the ions in the PDB file have the same formula, description, etc. as currently in the database
        existing_ion_codes = set()
        for pdb_code, d in ions.iteritems():
            existing_db_record = tsession.query(DBIon).filter(DBIon.PDBCode == pdb_code)
            if existing_db_record.count() > 0:
                assert(existing_db_record.count() == 1)
                existing_db_record = existing_db_record.one()
                if ions[pdb_ion_object.PDBCode]['PDBCode'] == existing_db_record.PDBCode:
                    if db_record.FileSource == 'RCSB':
                        assert(ions[pdb_ion_object.PDBCode]['Description'] == existing_db_record.Description) # This can differ e.g. CL in 127L is 3(CL 1-) since there are 3 ions but in PDB files with 2 ions, this can be 2(CL 1-). We can assert this if we do extra parsing.
                        assert(ions[pdb_ion_object.PDBCode]['Formula'] == existing_db_record.Formula) # This can differ e.g. CL in 127L is 3(CL 1-) since there are 3 ions but in PDB files with 2 ions, this can be 2(CL 1-). We can assert this if we do extra parsing.
                    else:
                        if ions[pdb_ion_object.PDBCode]['Description'] != existing_db_record.Description:
                            colortext.warning('The description for {0} ("{1}") does not match the database record. This may occur if the PDB headers are missing. However, it also may indicate that the code "{0}" for the ion was manually altered which is not handled by our pipeline and could cause errors.'.format(pdb_code, ions[pdb_ion_object.PDBCode]['Description']))
                        if ions[pdb_ion_object.PDBCode]['Formula'] != existing_db_record.Formula:
                            colortext.warning('The formula for {0} ("{1}") does not match the database record. This may occur if the PDB headers are missing. However, it also may indicate that the code "{0}" for the ion was manually altered which is not handled by our pipeline and could cause errors.'.format(pdb_code, ions[pdb_ion_object.PDBCode]['Formula']))
                existing_ion_codes.add(pdb_code)


        # Create the main Ion records, only creating records for ions in RCSB files.
        if db_record.FileSource == 'RCSB':
            for pdb_code, ion_record in ions.iteritems():
                if pdb_code not in existing_ion_codes:
                    # Do not add existing records
                    colortext.message('Adding ion {0}'.format(pdb_code))
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
                if db_record.FileSource == 'RCSB':
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


    def add_designed_pdb(self, structural_details,
                           allow_missing_params_files = False,
                           minimum_sequence_identity = 95.0,
                           previously_added = set(),
                           trust_database_content = False,
                           update_sections = set(),
                           tsession = None,
                           debug = True):
        '''

        :param structural_details: A dict fitting the defined structure (see below).
        :param allow_missing_params_files: If the PDB file contains ligands with no associated params files in structural_details
                                           then the function will return with a message. If this option is set to True
                                           then the file will be added but the ligand is not guaranteed to be kept
                                           by the protocol (this depends on the version of Rosetta - versions from February
                                           2016 onwards should have better ligand support).
        :param minimum_sequence_identity:  For non-RCSB files, we require a chain mapping from chains to the corresponding
                                           RCSB chains. Clustal is run to ensure that the sequence identity is at least this
                                           value.
        :param previously_added:
        :param trust_database_content:
        :param update_sections:
        :param tsession:
        :param debug: If debug is set to True then the transaction used to insert the structure into the database will be
                      rolled back and a message stating that the insertion would have been successful is returned in the
                      return dict.
        :return:

        One example of the dict structure is as follows:

            dict(
                # Fields required for all structures. We require non-RCSB structures to be based on RCSB structures at least in this import interface.
                rcsb_id = '1K5D',

                # Exactly one of these fields must be specified for non-RCSB structures.
                pdb_object = p,
                file_path = 'pdbs/1K5D2.pdb',

                # Fields required for non-RCSB structures.
                db_id = '1K5D_TP0', # must be between 5-10 characters
                techniques = "PDB_REDO",
                file_source = "Roger Wilco",
                user_id = "sierra",

                # Required for non-RCSB structures and optional for RCSB structures
                description = 'GSP1 complex (RanGAP1) from Tina Perica. This file was taken from PDB_REDO. Removed residues 180-213 from chain A (RAN/GSP1) (the ones wrapping around YRB1).',

                # Optional fields for both non-RCSB structures and optional for RCSB structures
                transmembrane = None, # None or True or False
                publication = None, # this should have a corresponding entry (Publication.ID) in the database
                resolution = None, # should be either None or a floating-point number

                # Exactly one of these fields must be specified for non-RCSB structures.
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
                ligand_params_file_paths = {
                    'G13' : 'temp/pdbs/1K5D2.params'
                },
            )
            '''


        ################################
        # Checks and balances
        ################################

        rcsb_id = structural_details['rcsb_id']
        ligand_params_file_paths = structural_details.get('ligand_params_file_paths', {})
        assert(isinstance(ligand_params_file_paths, dict))
        for k, v in ligand_params_file_paths.iteritems():
            ligand_params_file_paths[k] = os.path.abspath(v)

        # Type checks
        assert((isinstance(rcsb_id, str) or isinstance(rcsb_id, unicode)) and (4 == len(rcsb_id.strip())))
        rcsb_id = str(rcsb_id.strip())

        # Adding an RCSB structure - cascade into add_pdb_from_rcsb
        if not('file_path' in structural_details or 'pdb_object' in structural_details):

            # Read the RCSB PDB
            rcsb_pdb_object = self.get_pdb_object(rcsb_id)

            # Params files
            assert(len(set(ligand_params_file_paths.keys()).difference(rcsb_pdb_object.get_ligand_codes())) == 0)

            return self.add_pdb_from_rcsb(rcsb_id, previously_added = previously_added, trust_database_content = trust_database_content,
                                            update_sections = update_sections, ligand_params_file_paths = ligand_params_file_paths, debug = debug)

        # The remainder of this function adds a designed structure to the database
        pdb2pdb_chain_maps = []

        # Required fields
        design_pdb_id = structural_details['db_id']
        techniques = structural_details['techniques']
        assert((isinstance(design_pdb_id, str) or isinstance(design_pdb_id, unicode)) and (5 <= len(design_pdb_id.strip()) <= 10))
        design_pdb_id = str(design_pdb_id.strip())
        techniques = (techniques or '').strip() or None
        if (techniques == None) or not(isinstance(techniques, str) or isinstance(techniques, unicode)):
            raise colortext.Exception('The technique for generating the PDB file must be specified e.g. "Rosetta model" or "PDB_REDO structure" or "Manual edit".')
        if 'description' not in structural_details:
             raise colortext.Exception('A description is required for non-RCSB files. This should include any details on the structure preparation.')
        if 'file_source' not in structural_details:
             raise colortext.Exception('''A file_source is required for non-RCSB files. This should describe where the file came from e.g. "PDB REDO", "Rosetta", or the creator's name e.g. "Jon Snow".''')
        if 'user_id' not in structural_details:
             raise colortext.Exception('A user_id is required for non-RCSB files. This should correspond to a record in the User table.')
        file_source, description, user_id = structural_details['file_source'], structural_details['description'] , structural_details['user_id']
        assert((isinstance(file_source, str) or isinstance(file_source, unicode)) and (file_source.strip()))
        assert((isinstance(description, str) or isinstance(description, unicode)) and (description.strip()))
        assert((isinstance(user_id, str) or isinstance(user_id, unicode)) and (user_id.strip()))
        file_source = file_source.strip()
        description = description.strip()
        user_id = str(user_id.strip())

        # User checks. Make sure the user has a record in the database
        try:
            user_record = self.get_session().query(DBUser).filter(DBUser.ID == user_id).one()
        except:
            colortext.error('Could not find user "{0}" in the database.'.format(user_id))
            raise
        colortext.warning('User: {1} ({0})'.format(user_record.ID, ' '.join([n for n in [user_record.FirstName, user_record.MiddleName, user_record.Surname] if n])))

        # Optional fields
        resolution = structural_details.get('resolution')
        assert(resolution == None or isinstance(resolution, float))
        transmembrane = structural_details.get('transmembrane')
        assert(transmembrane == None or isinstance(transmembrane, bool))

        # Publication checks
        # todo: if publication, assert that the publication record exists
        publication = structural_details.get('publication')
        assert(publication == None) #todo: handle

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
            assert(isinstance(chain_mapping, dict) and chain_mapping)
            chain_mapping_keys = set(chain_mapping.keys())
            assert(len(chain_mapping_keys.intersection(set(main_chain_ids))) == len(main_chain_ids)) # main_chain_ids is a subset of chain_mapping_keys
            assert(len(set(all_chain_ids).intersection(chain_mapping_keys)) == len(chain_mapping_keys)) # chain_mapping_keys is a subset of all_chain_ids
        unmapped_chains = sorted(set(all_chain_ids).difference(chain_mapping_keys))

        # Add the RCSB structure to the database and then use that object to run Clustal
        colortext.pcyan('Adding the original PDB file using a separate transaction.')
        self.add_pdb_from_rcsb(rcsb_id, previously_added = previously_added, trust_database_content = True, ligand_params_file_paths = {})

        # Now that we have added the RCSB structure, create a new session which will be aware of that structure
        tsession = tsession or self.get_session(new_session = True)
        try:
            attempts = 0
            rcsb_object = None
            while (not rcsb_object) and (attempts < 10):
                # Hacky but there seems to be a race condition here between closing the previous transaction in add_pdb_from_rcsb and creating the new transaction (tsession) above
                try:
                    rcsb_object = self.get_pdb_object(rcsb_id, tsession = tsession)
                except Exception, e:
                    colortext.warning(str(e))
                    colortext.warning(traceback.format_exc())
                    attempts += 1
                    time.sleep(1)
            if not rcsb_object:
                raise colortext.Exception('Race condition detected. Try the import again.')

            # Run Clustal to check whether the chain mapping is correct
            rcsb_db_record_object = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == rcsb_id))
            design_sequences = None

            if designed_pdb_object.seqres_sequences:
                design_sequences = designed_pdb_object.seqres_sequences
            else:
                design_sequences = designed_pdb_object.atom_sequences
            pcsa = PDBChainSequenceAligner()
            for chain_id, sequence in sorted(rcsb_object.seqres_sequences.iteritems()):
                pcsa.add(rcsb_id, chain_id, str(sequence))
            for chain_id, sequence in sorted(design_sequences.iteritems()):
                pcsa.add(design_pdb_id, chain_id, str(sequence))
            output, best_matches = pcsa.align(ignore_bad_chains = True)
            for dc, rc in sorted(chain_mapping.iteritems()):
                # todo: this mapping is untested on ligand chains. We could do: for all c in unmapped_chains, see if the 3-letter sequences match and, if so, and add those
                k1 = '{0}_{1}'.format(design_pdb_id, dc)
                k2 = '{0}_{1}'.format(rcsb_id, rc)
                assert(k1 in best_matches and k2 in best_matches[k1])
                assert(k2 in best_matches and k1 in best_matches[k2])

                sequence_identity = best_matches[k1][k2]
                pdb2pdb_chain_maps.append(dict(
                    PDBFileID1 = design_pdb_id,
                    Chain1 = dc,
                    PDBFileID2 = rcsb_id,
                    Chain2 = rc,
                    SequenceIdentity = sequence_identity,
                ))
                if sequence_identity < minimum_sequence_identity:
                    raise colortext.Exception('The chain mapping is defined as {0}, chain {1} -> {2}, chain {3} but the sequence identity for these chains is only {4}%.'.format(design_pdb_id, dc, rcsb_id, rc, sequence_identity))

                sequence_identity = best_matches[k2][k1]
                pdb2pdb_chain_maps.append(dict(
                    PDBFileID1 = rcsb_id,
                    Chain1 = rc,
                    PDBFileID2 = design_pdb_id,
                    Chain2 = dc,
                    SequenceIdentity = sequence_identity,
                ))
                if sequence_identity < minimum_sequence_identity:
                    raise colortext.Exception('The chain mapping is defined as {0}, chain {1} -> {2}, chain {3} but the sequence identity for these chains is only {4}%.'.format(rcsb_id, rc, design_pdb_id, dc, sequence_identity))

            # Make sure that no params files were specified that do not fit the PDB file
            ligands_not_present = set(ligand_params_file_paths.keys()).difference(designed_pdb_object.get_ligand_codes())
            if len(ligands_not_present) > 0:
                raise colortext.Exception('Params files were specified for ligands which do not exist in the PDB file: "{0}".'.format('", "'.join(ligands_not_present)))

            # Sanity-check the ligand mapping for consistency
            ligand_code_map = {}
            if 'ligand_instance_mapping' in structural_details:
                for k, v in structural_details['ligand_instance_mapping'].code_map.iteritems():
                    assert(ligand_code_map.get(k) == None or ligand_code_map[k] == v)
                    ligand_code_map[k] = v
            if 'ligand_code_mapping' in structural_details:
                for k, v in structural_details['ligand_code_mapping'].iteritems():
                    assert(ligand_code_map.get(k) == None or ligand_code_map[k] == v)
                    ligand_code_map[k] = v
            if 'unchanged_ligand_codes' in structural_details:
                for k in structural_details['unchanged_ligand_codes']:
                    assert(ligand_code_map.get(k) == None or ligand_code_map[k] == k)
                    ligand_code_map[k] = k

            # Sanity-check the ligand mapping for completeness
            ligand_mapping = LigandMap.from_code_map(ligand_code_map)
            assert(isinstance(ligand_mapping, LigandMap) and ligand_mapping)
            if sorted(ligand_code_map.keys()) != sorted(designed_pdb_object.get_ligand_codes()):
                raise colortext.Exception('Incomplete mapping or unexpected entries: The ligand mapping contains mappings for "{0}" but the PDB file contains ligand codes "{1}".'.format('", "'.join(sorted(ligand_code_map.keys())), '", "'.join(sorted(designed_pdb_object.get_ligand_codes()))))

            # todo: currently unhandled
            # todo: add ion mapping support - this seems less important as it is probably less likely that users will rename ion codes
            assert('ion_instance_mapping' not in structural_details and 'ion_code_mapping' not in structural_details)
            # Check to make sure that the set of ions is a subset of those in the original PDB
            # todo: Ideally, we should pass an ion mapping like for the ligand mapping. However, this may be annoying for users
            #       to have to specify. Instead, we could check to make sure that the mapping ion_code -> atom_type matches in
            #       both structures which would be a more specific check than the assertion below. This would disallow renaming
            #       of ions but this seems a reasonable trade-off.
            assert(len(set(designed_pdb_object.get_ion_codes()).difference(set(rcsb_object.get_ion_codes()))) == 0)


            ################################
            # Data entry
            ################################


            # Add the original PDB file to the database
            colortext.message('Adding designed PDB file {0} based off {1}.'.format(design_pdb_id, rcsb_id))
            db_record_object = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == design_pdb_id))
            is_new_record = db_record_object == None
            if not is_new_record:
                print('Retrieving designed PDB {0} from database.'.format(design_pdb_id))
                db_pdb_object = PDB(db_record_object.Content, parse_ligands = True)
                assert(designed_pdb_object.lines == db_pdb_object.lines)
                assert(db_record_object.FileSource != 'RCSB')
                if trust_database_content:
                    print('Trusting the existing data and early-outing.')
                    return design_pdb_id
            else:
                db_record_object = PDBFile(**dict(
                    ID = design_pdb_id,
                    FileSource = file_source,
                    Content = str(designed_pdb_object),
                    Techniques = techniques,
                    UserID = user_record.ID,
                    Notes = description,
                    DerivedFrom = rcsb_id,
                ))
                update_sections = set() # add all related data

            # Fill in the FASTA, Resolution, and b-factor fields of a PDBFile record.
            # We do not use the information from the RCSB object - derived PDB files may have useful new data e.g. b-factors
            # from PDB_REDO and may have different sequences.
            self._add_pdb_file_information(db_record_object, designed_pdb_object, design_pdb_id, is_new_record, is_rcsb_pdb = False)

            # We copy the transmembrane classification (assuming that the designed protein keeps the same characteristic)
            if transmembrane == None:
                db_record_object.Transmembrane = rcsb_db_record_object.Transmembrane
            else:
                db_record_object.Transmembrane = transmembrane

            # Add a new PDBFile record
            if is_new_record:
                tsession.add(db_record_object)
            tsession.flush()

            assert(not publication) # todo: add publication entry here

            # add all other data
            self.add_pdb_data(tsession, design_pdb_id, update_sections = set(), ligand_mapping = ligand_mapping, chain_mapping = chain_mapping, ligand_params_file_paths = ligand_params_file_paths)

            # add the designed PDB -> RCSB PDB chain mapping
            print('Creating the chain mapping')
            for pdb2pdb_chain_map in pdb2pdb_chain_maps:
                pdb2pdb_chain_map = get_or_create_in_transaction(tsession, PDB2PDBChainMap, pdb2pdb_chain_map, missing_columns = ['ID'])

            previously_added.add(design_pdb_id)

            print('Success.\n')
            if debug:
                print('Debug call - rolling back the transaction.\n')
                tsession.rollback()
            else:
                tsession.commit()
            tsession.close()
        except:
            colortext.error('Failure.')
            tsession.rollback()
            tsession.close()
            raise


    def _check_derived_record_against_rcsb_record(self, tsession, database_pdb_id, pdb_object, chain_mapping):
        '''Sanity checks for derived structures compared to their RCSB ancestor.'''
        rcsb_chains = None
        rcsb_record = None
        chain_ids = sorted(set(pdb_object.seqres_sequences.keys() + pdb_object.atom_sequences.keys() + pdb_object.chain_types.keys()))
        if chain_mapping:
            db_record = tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id).one()
            assert(db_record.FileSource != 'RCSB')

            rcsb_chains = {}
            rcsb_record = self.get_rcsb_record(db_record, tsession = tsession)
            for r in tsession.query(PDBChain).filter(PDBChain.PDBFileID == rcsb_record.ID):
                rcsb_chains[r.Chain] = r

            for chain_id in chain_ids:
                if chain_id in chain_mapping:
                    assert(pdb_object.chain_types[chain_id] == rcsb_chains[chain_mapping[chain_id]].MoleculeType)
                else:
                    # We cannot assert(chain_ids == sorted(chain_mapping.keys())) as this can fail e.g. if a user splits
                    # chain C (protein + ligand) into chain A (protein) and chain X (ligand). Instead, we use this weaker assertion.
                    assert(pdb_object.chain_types[chain_id] != 'Protein')






def _test():

    # Create an import API instance
    importer = DataImportInterface.get_interface_with_config_file(cache_dir = '/kortemmelab/data/oconchus/ddgcache', echo_sql = False)

    # Access the SQLAlchemy session directly
    session = importer.session

    # Access the MySQLdb interface layer directly
    DDG_db = importer.DDG_db # or importerDDG_db_utf

    # Update certain properties of RCSB files in the database
    importer.update_pdbs(update_sections = set(['Residues', 'Publication']), start_at = None, restrict_to_file_source = 'RCSB')
    #importer.update_pdbs(update_sections = set(['Ligands', 'Ions']), start_at = None, restrict_to_file_source = 'RCSB')


