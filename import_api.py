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
from sqlalchemy import create_engine, and_
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
from klab.bio.pdbtm import PDBTM
from klab.db.sqlalchemy import get_single_record_from_query, get_or_create_in_transaction

from db_schema import test_schema_against_database_instance
from db_schema import PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue, LigandDescriptor, LigandIdentifier, LigandSynonym, PDBLigand
from db_schema import Ligand as DBLigand
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

            tsession.commit()
            print('Success.\n')
        except:
            colortext.error('Failure.')
            tsession.rollback()
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
            PDBFile()
            db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == pdb_id))
            existing_record = db_record != None
            if existing_record:
                print('Retrieving {0} from database.'.format(pdb_id))
                assert(db_record.count() == 1)
                db_record = db_record[0]
                pdb_object = PDB(db_record.Content)
                assert(db_record.FileSource == 'RCSB')
            else:
                db_record = PDBFile()
                print('Retrieving {0} from RCSB.'.format(pdb_id))
                contents = rcsb.retrieve_pdb(pdb_id)
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

            # B-factors
            if not(existing_record) or (not db_record.BFactors):
                db_record.BFactors = pickle.dumps(pdb_object.ComputeBFactors())

            # Publication
            if not(existing_record) or (not db_record.Publication):
                self._add_pdb_publication(tsession, db_record.ID)

            # Transmembrane
            if not(existing_record) or (db_record.Transmembrane == None):
                pdbtmo = PDBTM(read_file('/kortemmelab/shared/mirror/PDBTM/pdbtmall.xml'))
                pdb_id_map = pdbtmo.get_pdb_id_map()
                uc_pdb_id_map = {}
                for k, v in pdb_id_map.iteritems():
                    uc_pdb_id_map[k.upper()] = v
                if pdb_id in uc_pdb_id_map:
                    print('{0} is a transmembrane protein.'.format(pdb_id))
                else:
                    print('{0} is not a transmembrane protein.'.format(pdb_id))
                db_record.Transmembrane = pdb_id in pdb_id_map

            # Add a new PDBFile record
            if not(existing_record):
                db_record.UserID = None
                db_record.Notes = None
                db_record.DerivedFrom = None
                tsession.add(db_record)
                tsession.flush()

            # add all other data
            self.add_pdb_data(tsession, pdb_id, update_sections = update_sections)

            previously_added.add(pdb_id)
            #tsession.commit()
            print('Success.\n')
        except:
            colortext.error('Failure.')
            tsession.rollback()
            raise


    def get_pdb_object(self, database_pdb_id, tsession = None):
        '''Create a PDB object from content in the database.'''
        tsession = tsession or self.get_session()
        db_record = get_single_record_from_query(tsession.query(PDBFile).filter(PDBFile.ID == database_pdb_id))
        assert(db_record)
        return PDB(db_record.Content)


    def update_pdbs(self, pdb_ids = [], update_sections = {}):
        '''Updates all or selected data for all or selected PDB files in the database.'''
        if not pdb_ids:
            pdb_ids = [r['ID'] for r in self.DDG_db.execute_select('SELECT ID FROM PDBFile')]
        for pdb_id in pdb_ids:
            #if pdb_id not in ['1QM4', '1SEE', '1WSY', '2y2W9N', 'S9G10_best', 'uby_1UBQ', 'uby_CUE', 'uby_OTU', 'uby_RPN13', 'uby_SH3', 'y1AAR']:
            #if pdb_id.upper() < '1GPW':#. not in ['1QM4']:
            #    continue
            #if pdb_id.upper() != '2VLO':
            #    continue
            colortext.message('Updating data for {0}.'.format(pdb_id))
            tsession = self.get_session(new_session = True)
            self.add_pdb_data(tsession, pdb_id, update_sections = update_sections)
            tsession.commit()


    def add_pdb_data(self, tsession, database_pdb_id, update_sections = {}):
        '''database_pdb_id is the RCSB ID for RCSB files and a custom ID for other (designed) structures.
           If transaction_session is None (e.g. if this was called directly outside of a transaction), create a transaction
           session for the remaining inner calls. If update_sections is non-empty, just call those specific inner functions.
        '''
        if not tsession:
            tsession = self.get_session(new_session = True)

        # Create a PDB object
        pdb_object = self.get_pdb_object(database_pdb_id)

        if not(update_sections) or ('Chains' in update_sections):
            self._add_pdb_chains(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Molecules' in update_sections):
            self._add_pdb_molecules(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Residues' in update_sections):
            self._add_pdb_residues(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Ligands' in update_sections):
            self._add_pdb_rcsb_ligands(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('UniProt' in update_sections):
            self._add_pdb_uniprot_mapping(tsession, database_pdb_id, pdb_object)
        if not(update_sections) or ('Publication' in update_sections):
            self._add_pdb_publication(tsession, database_pdb_id, pdb_object)


    def _add_pdb_chains(self, tsession, database_pdb_id, pdb_object = None):
        '''
           Touched tables:
               PDBChain
        '''
        from db_schema import PDBFile, PDBChain
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id)

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
                instance = PDBChain(**dict(
                    PDBFileID = database_pdb_id,
                    Chain = c,
                    MoleculeType = pdb_object.chain_types[c]
                ))
                tsession.add(instance)
                tsession.flush()
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


    def _add_pdb_molecules(self, tsession, database_pdb_id, pdb_object = None, allow_missing_molecules = True):
        '''
           Add PDBMolecule and PDBMoleculeChain records
           Touched tables:
               PDBMolecule
               PDBMoleculeChain
        '''
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id)
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
            instance = PDBMolecule(md)
            tsession.add(instance)
            tsession.flush()
            for c in chains:
                try:
                    instance = PDBMoleculeChain(dict(
                        PDBFileID = database_pdb_id,
                        MoleculeID = md['MoleculeID'],
                        Chain = c
                    ))
                    tsession.add(instance)
                    tsession.flush()
                except:
                    if allow_missing_molecules: pass
                    else: raise


    def _add_pdb_residues(self, tsession, database_pdb_id, pdb_object = None):
        '''
           Touched tables:
               PDBResidue
        '''
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id)
        raise Exception('Implement')


        #allow_missing_molecules = False


    def _add_pdb_rcsb_ligands(self, tsession, database_pdb_id, pdb_object = None):
        '''This function associates the ligands of a PDB file (which may be arbitrarily named) with ligands entered in
           the database using the ligand's PDB code. The insertion is handled by a transaction which should be set up
           by the caller.
           Touched tables:
               PDBLigand
        '''
        # for all ligand codes in the PDB (if RCSB PDB):
        #    call importer.add_ligand_by_pdb_code(l)

        # for all ligand codes in the PDB (if not RCSB PDB), require a mapping to PDB ligand IDs and then:
        #    call importer.add_ligand_by_pdb_code(l)
        #    create a PDBLigand record with the Ligand.ID associated with the PDB ligand ID
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id)
        pass

        try: pass
        except Exception, e:
            colortext.error(str(e))
            colortext.error(traceback.format_exc())
            raise Exception("An exception occurred committing %s to the database." % filepath)


    def _add_pdb_uniprot_mapping(self, tsession, database_pdb_id, pdb_object = None):
        '''UniProtACs have forms like 'P62937' whereas UniProtIDs have forms like 'PPIA_HUMAN.'''
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id)
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

        raise Exception('test here')
        pdb_object = pdb_object or self.get_pdb_object(database_pdb_id)
        pdb_record = tsession.query(PDBFile).filter(PDBFile.ID == pdb_id)
        PUBTYPES = ['ISSN', 'ESSN']

        pdb_object = PDB(pdb_record.Content)
        j = pdb_object.get_journal()
        if not j:
            return
        pdb_id = pdb_id.strip()
        print(j["REFN"])
        raise Exception('test here')

        # We identify the sources for a PDB identifier with that identifier
        PublicationID = "PDB:%s" % pdbID
        publication_record = get_or_create_in_transaction(tsession, Publication, dict(ID = PublicationID))
        pdb_record.Publication = publication_record.ID

        locations = tsession.query(PublicationIdentifier).filter(PublicationIdentifier.SourceID == publication_record.ID)
        pub_locations = [location for location in locations if location.Type in PUBTYPES]
        doi_locations = [location for location in locations if location.Type == 'DOI']
        assert(len(pub_locations) <= 1)
        assert(len(doi_locations) <= 1)
        if j["published"]:
            if pub_locations:
                location = pub_locations[0]
                if j["REFN"]["type"] == location['Type']:
                    if j["REFN"]["ID"] != location['ID']:
                        colortext.warning("REFN: Check that the PublicationIdentifier data ('%s') matches the PDB REFN data ('%s')." % (str(location), j["REFN"]))
            elif j.get('REFN'):
                assert(j["REFN"]["type"] in PUBTYPES)
                instance = PublicationIdentifier(dict(
                    SourceID    = PublicationID,
                    ID          = j["REFN"]["ID"],
                    Type        = j["REFN"]["type"],
                ))
                tsession.add(instance)
                tsession.flush()
        if j["DOI"]:
            if doi_locations:
                location = doi_locations[0]
                if j["DOI"] != location['ID']:
                    colortext.warning("DOI: Check that the PublicationIdentifier data ('%s') matches the PDB DOI data ('%s')." % (str(doi_locations), j["DOI"]))
            else:
                instance = PublicationIdentifier(dict(
                    SourceID    = PublicationID,
                    ID          = j["DOI"],
                    Type        = "DOI",
                ))
                tsession.add(instance)
                tsession.flush()


    def add_designed_pdb_file(self, designed_pdb_filepath, design_pdb_id, original_pdb_id, description, username, chain_mapping = {}, ligand_mapping = {}):
        '''Wrapper for add_designed_pdb.'''
        return self.add_designed_pdb(self, PDB.from_filepath(designed_pdb_filepath), design_pdb_id, original_pdb_id, description, username, chain_mapping = chain_mapping, ligand_mapping = ligand_mapping)


    def add_designed_pdb(self, designed_pdb_object, design_pdb_id, original_pdb_id, description, username, chain_mapping = {}, ligand_mapping = {}, trust_database_content = True):
        assert(5 <= len(design_pdb_id) <= 10)

        colortext.message('Adding designed PDB file {0} based off {1}.'.format(design_pdb_id, original_pdb_id))

        colortext.pcyan('Adding the original PDB file using a separate transaction.')
        self.add_pdb_from_rcsb(original_pdb_id, trust_database_content = trust_database_content)
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

    if False:
        ligand_codes = ['0Z6', '0G6', '0EF', '0QE', '13P', '15P', '1AC', '1PE', '1PG', '2AB', '2GP', '3AA', '3GP', '3PG', '544', '5IU', '9CR', 'A80', 'ABA', 'AC9', 'ACA', 'ACE', 'ACH', 'ACT', 'ACY', 'ADC', 'ADP', 'ADU', 'AE3', 'AF3', 'AGF', 'AGL', 'AGP', 'AI2', 'ALA', 'ALC', 'ALF', 'AMB', 'AMP', 'ANL', 'ANP', 'ANS', 'AOA', 'AP5', 'APB', 'AR7', 'ARA', 'ARB', 'ASP', 'ATP', 'AXP', 'B3I', 'B3P', 'B98', 'BCT', 'BEN', 'BEO', 'BGC', 'BHD', 'BID', 'BLA', 'BM6', 'BMA', 'BME', 'BNG', 'BNZ', 'BOG', 'BR', 'BRL', 'BTW', 'C1O', 'C2O', 'C8E', 'CA', 'CAC', 'CAM', 'CCS', 'CD', 'CEM', 'CF0', 'CFM', 'CGU', 'CIT', 'CL', 'CLF', 'CME', 'CMP', 'CO', 'CO3', 'CO9', 'COT', 'CSO', 'CSS', 'CSW', 'CSX', 'CTP', 'CU', 'CU1', 'CUA', 'CYN', 'CYS', 'DAF', 'DDE', 'DDF', 'DGA', 'DGN', 'DIO', 'DKA', 'DLY', 'DMS', 'DPG', 'DSN', 'DVA', 'DVR', 'E64', 'E6C', 'EBP', 'EDO', 'EMO', 'EPE', 'EQP', 'EQU', 'ETX', 'F09', 'F25', 'FAD', 'FCY', 'FE', 'FES', 'FFO', 'FK5', 'FKP', 'FLA', 'FLC', 'FME', 'FMN', 'FMT', 'FOK', 'FOS', 'FUC', 'FUL', 'G6D', 'GAI', 'GAL', 'GCP', 'GDP', 'GL0', 'GLA', 'GLC', 'GLU', 'GLY', 'GLZ', 'GNH', 'GNP', 'GOL', 'GSH', 'GSP', 'GTP', 'GTS', 'GTT', 'GVE', 'HAS', 'HCA', 'HEA', 'HEC', 'HED', 'HEM', 'HG', 'HIC', 'HMC', 'HPR', 'HYP', 'IAS', 'IBR', 'IET', 'IGP', 'IMD', 'IOD', 'IPA', 'IPH', 'IUM', 'K', 'LAR', 'LAT', 'LBT', 'LDA', 'LNO', 'LPX', 'LYZ', 'M3L', 'MAK', 'MAL', 'MAN', 'MES', 'MG', 'MHO', 'MLA', 'MN', 'MN3', 'MNA', 'MNQ', 'MOH', 'MPC', 'MPD', 'MPT', 'MRD', 'MSE', 'MTX', 'MYR', 'NA', 'NAD', 'NAG', 'NAI', 'NAP', 'NDG', 'NDP', 'NES', 'NH2', 'NH4', 'NHE', 'NI', 'NO3', 'NPS', 'O', 'OCT', 'OLA', 'OXL', 'OXP', 'P34', 'P5P', 'P6G', 'PBM', 'PCA', 'PG4', 'PGA', 'PGE', 'PHQ', 'PLG', 'PLM', 'PLP', 'PMP', 'PMS', 'PO4', 'PON', 'PRO', 'PTR', 'RAP', 'RB', 'RDA', 'REA', 'RET', 'RIP', 'RS2', 'RTL', 'SCH', 'SE4', 'SEP', 'SF4', 'SIA', 'SIN', 'SM', 'SMC', 'SN1', 'SO2', 'SO3', 'SO4', 'SOC', 'SPA', 'SR', 'STI', 'STU', 'SVA', 'TAC', 'TAD', 'TAR', 'TBA', 'TEM', 'THP', 'THU', 'TLA', 'TMP', 'TPO', 'TRN', 'TRP', 'TRQ', 'TRS', 'TXP', 'TYS', 'UNL', 'UNX', 'URA', 'V7O', 'VAL', 'XYP', 'Y1', 'YBT', 'ZAP', 'ZN', 'ZNH']
        num_ligands = len(ligand_codes)
        c = 0
        for l in ligand_codes:
            c += 1
            #if c < 280:
            #    continue
            print('{0}/{1}'.format(c, num_ligands))
            importer.add_ligand_by_pdb_code(l)

    #importer.update_pdbs(pdb_ids = ['1a22'], update_sections = ['Chains'])
    importer.update_pdbs(update_sections = ['Chains'])
    return

    importer.add_pdb_from_rcsb('1bn3')

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

test_cases = '2FX5' #  at the time of writing, these PDB IDs had no JRNL lines
test_cases = '1GTX', '1UOX', '1WSY', '1YYJ', '2IMM', '2MBP' # at the time of writing, these PDB IDs had no UniProt entries

sys.exit(0)

