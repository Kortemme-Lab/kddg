#!/usr/bin/python2.4
# encoding: utf-8
"""
db_schema.py
SQLAlchemy representation

Created by Shane O'Connor 2015.
Copyright (c) 2015 Shane O'Connor. All rights reserved.
"""

import sys
import os
import inspect
import pprint

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy import inspect as sqlalchemy_inspect
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, LONGBLOB
from sqlalchemy.types import DateTime, Enum, Integer, TIMESTAMP, Text, Unicode, String

if __name__ == '__main__':
    sys.path.insert(0, '../../klab')

from klab.db.sqlalchemy import MySQLSchemaConverter
from klab.fs.fsio import read_file
from klab import colortext


DeclarativeBase = declarative_base()


#############################
#                           #
#  Foundational constructs  #
#                           #
#############################



class AminoAcid(DeclarativeBase):
    __tablename__ = 'AminoAcid'

    Code = Column(Unicode(1), nullable=False, primary_key=True)
    LongCode = Column(Unicode(3), nullable=False)
    Name = Column(Unicode(32), nullable=False)
    Polarity = Column(Enum('polar','non-polar','charged'), nullable=True)
    Aromaticity = Column(Enum('aromatic','aliphatic','neither'), nullable=True)
    Hydrophobicity_pH7 = Column(Enum('hydrophobic','hydrophilic'), nullable=True)
    SideChainAcidity = Column(Enum('acidic','basic','neutral'), nullable=True)
    pKa = Column(DOUBLE, nullable=True)
    AverageMass = Column(DOUBLE, nullable=True)
    Volume = Column(DOUBLE, nullable=True)
    Size = Column(Enum('small','large'), nullable=True)
    Tiny = Column(TINYINT(1), nullable=True, default=0)


######################################
#                                    #
#  Ligands and associated records    #
#                                    #
######################################


class Ligand(DeclarativeBase):
    __tablename__ = 'Ligand'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBCode = Column(String(3), nullable=False)
    LigandCode = Column(String(512), nullable=False)
    Formula = Column(String(256), nullable=False)
    MolecularWeight = Column(DOUBLE, nullable=False)
    LigandType = Column(String(256), nullable=False)
    Solubility = Column(String(256), nullable=True)
    CellPermeability = Column(Enum('Yes','No','Yes (probably)', nullable=True))
    AssaysToDetermineConcentrationInCells = Column(String(256), nullable=True)
    ProductionInCells = Column(Enum('Yes','No'), nullable=True)
    ProductionInCellsNotes = Column(String(64), nullable=True)
    Diagram = Column(LONGBLOB, nullable=True)
    SimilarCompoundsDiagram = Column(LONGBLOB, nullable=True)
    InChI = Column(Text, nullable=False)
    InChIKey = Column(String(27), nullable=False)


class LigandDescriptor(DeclarativeBase):
    __tablename__ = 'LigandDescriptor'

    ID = Column(Integer, nullable=False, primary_key=True)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    Descriptor = Column(Text, nullable=False)
    DescriptorType = Column(String(128), nullable=False)
    Program = Column(String(64), nullable=False)
    Version = Column(String(16), nullable=False)


class LigandIdentifier(DeclarativeBase):
    __tablename__ = 'LigandIdentifier'

    ID = Column(Integer, nullable=False, primary_key=True)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    Identifier = Column(String(1024), nullable=False)
    IDType = Column(String(128), nullable=False)
    Program = Column(String(64), nullable=False)
    Version = Column(String(16), nullable=False)


class LigandPrice(DeclarativeBase):
    __tablename__ = 'LigandPrice'

    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False, primary_key=True)
    PriceDate = Column(DateTime, nullable=False, primary_key=True)
    USDPricePerGram = Column(DOUBLE, nullable=False)
    PriceNote = Column(String(256), nullable=True)


class LigandReference(DeclarativeBase):
    __tablename__ = 'LigandReference'

    ID = Column(Integer, nullable=False, primary_key=True)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    PublicationID = Column(String(64), ForeignKey('Publication.ID'), nullable=True)
    Type = Column(Enum('Reference','Assay'), nullable=False, default=u'Reference')
    Notes = Column(Text, nullable=True)


class LigandSynonym(DeclarativeBase):
    __tablename__ = 'LigandSynonym'

    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False, primary_key=True)
    Synonym = Column(String(256), nullable=False, primary_key=True)


######################################
#                                    #
#  PDB files and associated records  #
#                                    #
######################################


class PDBFile(DeclarativeBase):
    __tablename__ = 'PDBFile'

    ID = Column(Unicode(10), nullable=False, primary_key=True, default=u'')
    FileSource = Column(Unicode(64), nullable=False, default=u'RCSB')
    Content = Column(Text, nullable=False)
    FASTA = Column(Text, nullable=False)
    Resolution = Column(DOUBLE, nullable=True)
    Techniques = Column(Unicode(256), nullable=False)
    BFactors = Column(Text, nullable=False)
    Publication = Column(Unicode(64), nullable=True)
    Transmembrane = Column(TINYINT(1), nullable=True)
    UserID = Column(Unicode(64), nullable=True)
    Notes = Column(Text, nullable=True)
    DerivedFrom = Column(Unicode(4), nullable=True)

    # Relationships
    # todo: add publication relationship publication = relationship('Publication', foreign_keys=[Publication])

    def __repr__(self):
        return 'PDBFile: {0}. Source: {1}. Resolution {2}. {3}'.format(self.ID, self.FileSource, self.Resolution, self.Notes or '')


class PDBChain(DeclarativeBase):
    __tablename__ = 'PDBChain'

    PDBFileID = Column(Unicode(10), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    Chain = Column(Unicode(1), nullable=False, primary_key=True)
    MoleculeType = Column(Enum('Protein','DNA','RNA','Ligand'), nullable=True)
    WildtypeProteinID = Column(Unicode(18), nullable=True)
    FullProteinID = Column(Unicode(18), nullable=True)
    SegmentProteinID = Column(Unicode(18), nullable=True)
    WildtypeAlignedProteinID = Column(Unicode(18), nullable=True)
    AcquiredProteinID = Column(Unicode(18), nullable=True)
    Coordinates = Column(LONGBLOB, nullable=True)

    # Parent relationships
    pdb_file = relationship('PDBFile', primaryjoin="PDBChain.PDBFileID==PDBFile.ID")

    # Children relationships
    residues = relationship('PDBResidue', primaryjoin="and_(PDBResidue.PDBFileID==PDBChain.PDBFileID, PDBResidue.Chain==PDBChain.Chain)")

    def __repr__(self):
        return 'PDBChain: {0}, {1}'.format(self.PDBFileID, self.Chain)


class PDBMolecule(DeclarativeBase):
    __tablename__ = 'PDBMolecule'

    PDBFileID = Column(Unicode(8), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    MoleculeID = Column(Integer, nullable=False, primary_key=True)
    Name = Column(Unicode(256), nullable=False)
    Organism = Column(Unicode(256), nullable=True)
    Fragment = Column(Unicode(256), nullable=True)
    Synonym = Column(Unicode(256), nullable=True)
    Engineered = Column(TINYINT(1), nullable=True)
    EC = Column(Unicode(32), nullable=True)
    Mutation = Column(TINYINT(1), nullable=True)
    OtherDetails = Column(Unicode(256), nullable=True)

    # Parent relationships
    pdb_file = relationship('PDBFile', primaryjoin="PDBMolecule.PDBFileID==PDBFile.ID")

    # Children relationships
    chains = relationship("PDBMoleculeChain", primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBMolecule.PDBFileID, PDBMoleculeChain.MoleculeID==PDBMolecule.MoleculeID)")

    def __repr__(self):
        return 'PDBMolecule ({0}-{1}). Name: {2}. Organism: {3}'.format(self.PDBFileID, self.MoleculeID, '/'.join([s for s in [self.Name or '', self.Synonym or ''] if s]), self.Organism or 'N/A')


class PDBMoleculeChain(DeclarativeBase):
    __tablename__ = 'PDBMoleculeChain'

    PDBFileID = Column(Unicode(8), ForeignKey('PDBMolecule.PDBFileID'), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    MoleculeID = Column(Integer, ForeignKey('PDBMolecule.MoleculeID'), nullable=False, primary_key=True)
    Chain = Column(Unicode(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)

    # Parent relationships
    pdb_molecule = relationship('PDBMolecule', primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBMolecule.PDBFileID, PDBMoleculeChain.MoleculeID==PDBMolecule.MoleculeID)")
    pdb_chain = relationship('PDBChain', primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBChain.PDBFileID, PDBMoleculeChain.Chain==PDBChain.Chain)")

    def __repr__(self):
        return 'PDBMoleculeChain ({0}-{1}-{2}).'.format(self.PDBFileID, self.MoleculeID, self.Chain)


class PDBResidue(DeclarativeBase):
    __tablename__ = 'PDBResidue'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBFileID = Column(Unicode(10), ForeignKey('PDBChain.PDBFileID'), nullable=False)
    Chain = Column(Unicode(1), ForeignKey('PDBChain.Chain'), nullable=False)
    ResidueID = Column(Unicode(5), nullable=False)
    ResidueAA = Column(Unicode(1), ForeignKey('AminoAcid.Code'), nullable=False)
    ResidueType = Column(Enum('Protein','DNA','RNA','Ligand'), nullable=False, default=u'Protein')
    IndexWithinChain = Column(Integer, nullable=False)
    CoordinatesExist = Column(TINYINT(1), nullable=False)
    RecognizedByRosetta = Column(TINYINT(1), nullable=True)
    BFactorMean = Column(DOUBLE, nullable=True)
    BFactorDeviation = Column(DOUBLE, nullable=True)
    SecondaryStructurePosition = Column(Enum('Coil','Helix','Sheet','Turn','3_10-Helix'), nullable=True)
    AccessibleSurfaceArea = Column(DOUBLE, nullable=True)
    MonomericExposure = Column(DOUBLE, nullable=True)
    MonomericDSSP = Column(Unicode(1), nullable=True)
    ComplexExposure = Column(DOUBLE, nullable=True)
    ComplexDSSP = Column(Unicode(1), nullable=True)

    # Parent relationships
    pdb_chain = relationship('PDBChain', primaryjoin="and_(PDBResidue.PDBFileID==PDBChain.PDBFileID, PDBResidue.Chain==PDBChain.Chain)")

    # Child relationships
    residue = relationship('AminoAcid', primaryjoin="PDBResidue.ResidueAA==AminoAcid.Code")

    def __repr__(self):
        return 'PDBResidue {0}. {1} {2} {3} ({4}). Exposure: {5}. DSSP: {6}.'.format(self.residue.LongCode, self.PDBFileID, self.Chain, (self.ResidueAA + self.ResidueID.strip()).ljust(5), self.ResidueType, self.ComplexExposure, self.ComplexDSSP)


class PDBLigand(DeclarativeBase):
    __tablename__ = 'PDBLigand'

    PDBFileID = Column(Unicode(10), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    Chain = Column(Unicode(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)
    SeqID = Column(Unicode(5), nullable=False, primary_key=True)
    PDBLigandCode = Column(Unicode(3), nullable=False)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    ParamsFileContentID = Column(Integer, ForeignKey('FileContent.ID'), nullable=False)


###########################
#                         #
#  Bookkeeping functions  #
#                         #
###########################




def generate_sqlalchemy_definition(tablenames = []):
    '''This function generates the SQLAlchemy class definitions from the database. The generation does not parse the
       entire definition - it omits unique keys, foreign key constraints etc. but it saves a lot of manual work setting
       up the boilerplate field definitions. When the database schema changes, call this function to update the
       SQLAlchemy class definitions. You may want/need to reuse any existing relationships defined between tables.'''
    sc = MySQLSchemaConverter('kortemmelab', 'kortemmelab.ucsf.edu', 'ddG', read_file(os.path.join('..', 'pw')).strip(), 3306, "/var/lib/mysql/mysql.sock")
    #sc.get_sqlalchemy_schema(['PDBFile', 'PDBChain', 'PDBMolecule', 'PDBMoleculeChain', 'PDBResidue'])
    sc.get_sqlalchemy_schema(tablenames)


def test_schema_against_database_instance(DDG_db):
    '''Make sure that our SQLAlchemy definitions match the database. This should be run by the API prior to connection
       as it lets the admin know that they need to update the schema here (use generate_sqlalchemy_definition to update
       the schema).'''
    database_to_class_mapping = {}
    clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
    clsmembers = [c[1] for c in clsmembers if issubclass(c[1], DeclarativeBase) and c[1] != DeclarativeBase]
    for c in clsmembers:
        database_to_class_mapping[c.__tablename__] = c

    inconsistencies = []
    for tblname, pcls in sorted(database_to_class_mapping.iteritems()):
        represented_columns = set([c.name for c in list(sqlalchemy_inspect(pcls).columns)])
        tbl_columns = set([c['Field'] for c in DDG_db.execute_select('SHOW COLUMNS FROM {0}'.format(tblname))])
        if sorted(represented_columns) != sorted(tbl_columns):
            inconsistencies.append(pcls.__name__)
            colortext.error('The SQLAlchemy class {0} does not match the database schema.'.format(pcls.__name__))
            if represented_columns.difference(tbl_columns):
                colortext.warning('  The SQLAlchemy class has extra columns: {0}'.format(', '.join(sorted(represented_columns.difference(tbl_columns)))))
            if tbl_columns.difference(represented_columns):
                colortext.pcyan('  The MySQL schema definition has extra columns: {0}'.format(', '.join(sorted(tbl_columns.difference(represented_columns)))))
    if inconsistencies:
        generate_sqlalchemy_definition(inconsistencies)
        raise colortext.Exception('The following SQLAlchemy classes do not match the database schema: {0}.'.format(', '.join(inconsistencies)))



if __name__ == '__main__':
    generate_sqlalchemy_definition(['AminoAcid'])
    sys.exit(0)
    #generate_sqlalchemy_definition(['PDBFile', 'PDBChain', 'PDBMolecule', 'PDBMoleculeChain', 'PDBResidue'])
    from ppi_api import get_interface as get_ppi_interface
    ppi_api = get_ppi_interface(read_file(os.path.join('..', 'pw')).strip())
    test_schema_against_database_instance(ppi_api.DDG_db)