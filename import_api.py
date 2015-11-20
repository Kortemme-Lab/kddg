#!/usr/bin/python2.4
# encoding: utf-8
"""
import_api.py
High-level functions for importing data into the DDG database.

Created by Shane O'Connor 2015.
Copyright (c) 2015 Shane O'Connor. All rights reserved.
"""

import sys

if __name__ == '__main__':
    sys.path.insert(0, '../../klab')

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import create_engine

from klab.fs.fsio import read_file

from db_schema import PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue



###########################
#                         #
#  SQLAlchemy Engine      #
#                         #
###########################

engine = None
session = None

def get_engine(connect_string):
    global engine
    if not engine:
        engine = create_engine(connect_string)
    return engine


def get_connection(connect_string):
    # e.g. connection = get_connection(); connection.execute("SELECT * FROM User")
    global engine
    engine = get_engine(connect_string)
    return engine.connect()


def get_session(connect_string, new_session = False):
    global engine
    global session
    print(connect_string)
    engine = get_engine(connect_string)
    if new_session or not(session):
        maker_ddgdatabase = sessionmaker(autoflush = True, autocommit = False)
        s = scoped_session(maker_ddgdatabase)
        DeclarativeBaseDDG = declarative_base()
        metadata_ddg = DeclarativeBaseDDG.metadata
        s.configure(bind=engine)
        metadata_ddg.bind = engine
        if new_session:
            return s
        else:
            session = s
            return s
    return session


def test():
    session = get_session(read_file('db_connection_str'))
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

class PDBFile(object):
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


    def __init__(self, ddGdb, pdb_id, content = None, protein = None, contains_membrane_protein = None, file_source = None, filepath = None, UniProtAC = None, UniProtID = None, testonly = False, techniques = None, derived_from = None, notes = None):
        '''UniProtACs have forms like 'P62937' whereas UniProtIDs have forms like 'PPIA_HUMAN.'''
        super(PDBStructure, self).__init__(ddGdb)

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