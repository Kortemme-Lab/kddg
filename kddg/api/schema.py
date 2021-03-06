#!/usr/bin/python2.4

# encoding: utf-8
"""
schema.py
SQLAlchemy representation

Created by Shane O'Connor 2015.
Copyright (c) 2015 Shane O'Connor. All rights reserved.
"""

import sys
import os
import inspect
import pprint
import StringIO
import gzip
import pandas
import traceback
import getpass

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy import inspect as sqlalchemy_inspect
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, LONGBLOB
from sqlalchemy.types import DateTime, Enum, Integer, TIMESTAMP, Text, Unicode, String
from sqlalchemy.orm import deferred

from klab.db.sqlalchemy_interface import MySQLSchemaConverter
from klab.fs.fsio import read_file
from klab import colortext
from klab.fs.fsio import read_file, write_file, write_temp_file

from kddg.api import settings

sys_settings = settings.load()

DeclarativeBase = declarative_base()


#############################
#                           #
#  Foundational constructs  #
#                           #
#############################



class AminoAcid(DeclarativeBase):
    __tablename__ = 'AminoAcid'

    Code = Column(String(1), nullable=False, primary_key=True)
    LongCode = Column(String(3), nullable=False)
    Name = Column(String(32), nullable=False)
    Polarity = Column(Enum('polar','non-polar','charged'), nullable=True)
    Aromaticity = Column(Enum('aromatic','aliphatic','neither'), nullable=True)
    Hydrophobicity_pH7 = Column(Enum('hydrophobic','hydrophilic'), nullable=True)
    SideChainAcidity = Column(Enum('acidic','basic','neutral'), nullable=True)
    pKa = Column(DOUBLE, nullable=True)
    AverageMass = Column(DOUBLE, nullable=True)
    Volume = Column(DOUBLE, nullable=True)
    Size = Column(Enum('small','large'), nullable=True)
    Tiny = Column(TINYINT(1), nullable=True, default=0)


class FileContent(DeclarativeBase):
    __tablename__ = 'FileContent'

    ID = Column(Integer, nullable=False, primary_key=True)
    Content = deferred(Column(LONGBLOB, nullable=False))
    MIMEType = Column(String(64), nullable=False)
    Filesize = Column(Integer, nullable=False)
    MD5HexDigest = Column(String(32), nullable=False)


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
    CellPermeability = Column(Enum('Yes','No','Yes (probably)'), nullable=True)
    AssaysToDetermineConcentrationInCells = Column(String(256), nullable=True)
    ProductionInCells = Column(Enum('Yes','No'), nullable=True)
    ProductionInCellsNotes = Column(String(64), nullable=True)
    Diagram = deferred(Column(LONGBLOB, nullable=True))
    SimilarCompoundsDiagram = deferred(Column(LONGBLOB, nullable=True))
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
#  Ions and associated records       #
#                                    #
######################################


class Ion(DeclarativeBase):
    __tablename__ = 'Ion'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBCode = Column(String(3), nullable=False)
    Formula = Column(String(256), nullable=False)
    Description = Column(String(256), nullable=True)


######################################
#                                    #
#  PDB files and associated records  #
#                                    #
######################################


class PDBFile(DeclarativeBase):
    __tablename__ = 'PDBFile'

    ID = Column(String(10), nullable=False, primary_key=True, default=u'')
    FileSource = Column(String(64), nullable=False, default=u'RCSB')
    Content = deferred(Column(Text, nullable=False))
    FASTA = deferred(Column(Text, nullable=False))
    Resolution = Column(DOUBLE, nullable=True)
    Techniques = Column(String(256), nullable=False)
    BFactorMean = Column(DOUBLE, nullable=True)
    BFactorDeviation = Column(DOUBLE, nullable=True)
    Publication = Column(String(64), ForeignKey('Publication.ID'), nullable=True)
    Transmembrane = Column(TINYINT(1), nullable=True)
    UserID = Column(String(64), nullable=True)
    Notes = Column(Text, nullable=True)
    DerivedFrom = Column(String(4), nullable=True)

    # Relationships
    publication = relationship('Publication', viewonly=True, primaryjoin="PDBFile.Publication==Publication.ID")

    def __repr__(self):
        notes = ''
        resolution = 'unknown resolution'
        if self.Resolution:
            resolution = str(self.Resolution) + 'A'
        if self.Notes and self.Notes.strip():
            notes = self.Notes.strip()
            if notes[-1] != '.':
                notes += '.'
        return 'PDBFile: {0}, {1} ({2}). Source: {3}. B-Factors: {4} (mean), {5} (stddev). Transmembrane protein: {6}. {7}'.format(self.ID, resolution, self.Techniques, self.FileSource, self.BFactorMean, self.BFactorDeviation, self.Transmembrane, notes)


class PDBChain(DeclarativeBase):
    __tablename__ = 'PDBChain'

    PDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    Chain = Column(String(1), nullable=False, primary_key=True)
    MoleculeType = Column(Enum('Protein', 'DNA', 'RNA', 'Ligand', 'Protein skeleton', 'Heterogen', 'Solution', 'Unknown'), nullable=True)
    WildtypeProteinID = Column(String(18), nullable=True)
    FullProteinID = Column(String(18), nullable=True)
    SegmentProteinID = Column(String(18), nullable=True)
    WildtypeAlignedProteinID = Column(String(18), nullable=True)
    AcquiredProteinID = Column(String(18), nullable=True)
    Coordinates = deferred(Column(LONGBLOB, nullable=True))

    # Parent relationships
    pdb_file = relationship('PDBFile', viewonly=True, primaryjoin="PDBChain.PDBFileID==PDBFile.ID")

    # Children relationships
    residues = relationship('PDBResidue', viewonly=True, primaryjoin="and_(PDBResidue.PDBFileID==PDBChain.PDBFileID, PDBResidue.Chain==PDBChain.Chain)")

    def __init__(self, **kwargs):
        super(PDBChain, self).__init__(**kwargs)
        # do custom initialization here

    def __repr__(self):
        return 'PDBChain: {0}, {1} ({2})'.format(self.PDBFileID, self.Chain, self.MoleculeType)


class PDBMolecule(DeclarativeBase):
    __tablename__ = 'PDBMolecule'

    PDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    MoleculeID = Column(Integer, nullable=False, primary_key=True)
    Name = Column(String(256), nullable=False)
    Organism = Column(String(256), nullable=True)
    Fragment = Column(String(256), nullable=True)
    Synonym = Column(String(256), nullable=True)
    Engineered = Column(TINYINT(1), nullable=True)
    EC = Column(String(32), nullable=True)
    Mutation = Column(TINYINT(1), nullable=True)
    OtherDetails = Column(String(256), nullable=True)

    # Parent relationships
    pdb_file = relationship('PDBFile', viewonly=True, primaryjoin="PDBMolecule.PDBFileID==PDBFile.ID")

    # Children relationships
    chains = relationship("PDBMoleculeChain", viewonly=True, primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBMolecule.PDBFileID, PDBMoleculeChain.MoleculeID==PDBMolecule.MoleculeID)")

    def __repr__(self):
        return 'PDBMolecule ({0}-{1}). Name: {2}. Organism: {3}'.format(self.PDBFileID, self.MoleculeID, '/'.join([s for s in [self.Name or '', self.Synonym or ''] if s]), self.Organism or 'N/A')


class PDBMoleculeChain(DeclarativeBase):
    __tablename__ = 'PDBMoleculeChain'

    # example of how to specify multiple foreign key relationships on one field (relationships defined below)

    PDBFileID = Column(String(10), ForeignKey('PDBMolecule.PDBFileID'), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    MoleculeID = Column(Integer, ForeignKey('PDBMolecule.MoleculeID'), nullable=False, primary_key=True)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)

    # Parent relationships
    pdb_molecule = relationship('PDBMolecule', viewonly=True, primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBMolecule.PDBFileID, PDBMoleculeChain.MoleculeID==PDBMolecule.MoleculeID)")
    pdb_chain = relationship('PDBChain', viewonly=True, primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBChain.PDBFileID, PDBMoleculeChain.Chain==PDBChain.Chain)")

    def __repr__(self):
        return 'PDBMoleculeChain ({0}-{1}-{2}).'.format(self.PDBFileID, self.MoleculeID, self.Chain)


class PDBResidue(DeclarativeBase):
    __tablename__ = 'PDBResidue'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBFileID = Column(String(10), ForeignKey('PDBChain.PDBFileID'), nullable=False)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False)
    ResidueID = Column(String(5), nullable=False)
    ResidueAA = Column(String(1), ForeignKey('AminoAcid.Code'), nullable=False)
    ResidueType = Column(Enum('Protein', 'DNA', 'RNA'), nullable=False, default=u'Protein')
    IndexWithinChain = Column(Integer, nullable=False)
    CoordinatesExist = Column(TINYINT(1), nullable=False)
    RecognizedByRosetta = Column(TINYINT(1), nullable=True)
    BFactorMean = Column(DOUBLE, nullable=True)
    BFactorDeviation = Column(DOUBLE, nullable=True)
    MonomericExposure = Column(DOUBLE, nullable=True)
    MonomericDSSP = Column(String(1), nullable=True)
    ComplexExposure = Column(DOUBLE, nullable=True)
    ComplexDSSP = Column(String(1), nullable=True)

    # Parent relationships
    pdb_chain = relationship('PDBChain', viewonly=True, primaryjoin="and_(PDBResidue.PDBFileID==PDBChain.PDBFileID, PDBResidue.Chain==PDBChain.Chain)")

    # Child relationships
    residue = relationship('AminoAcid', viewonly=True, primaryjoin="PDBResidue.ResidueAA==AminoAcid.Code")

    def __repr__(self):
        return 'PDBResidue {0}. {1} {2} {3} ({4}). Exposure: {5}. DSSP: {6}.'.format(self.residue.LongCode, self.PDBFileID, self.Chain, (self.ResidueAA + self.ResidueID.strip()).ljust(5), self.ResidueType, self.ComplexExposure, self.ComplexDSSP)


class PDBLigand(DeclarativeBase):
    __tablename__ = 'PDBLigand'

    PDBFileID = Column(String(10), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)
    SeqID = Column(String(5), nullable=False, primary_key=True)
    PDBLigandCode = Column(String(3), nullable=False)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)


class PDBLigandFile(DeclarativeBase):
    __tablename__ = 'PDBLigandFile'

    PDBFileID = Column(String(10), ForeignKey('PDBLigand.PDBFileID'), nullable=False, primary_key=True)
    PDBLigandCode = Column(String(3), ForeignKey('PDBLigand.PDBLigandCode'), nullable=False, primary_key=True)
    ParamsFileContentID = Column(Integer, ForeignKey('FileContent.ID'), nullable=False)


class PDBIon(DeclarativeBase):
    __tablename__ = 'PDBIon'

    PDBFileID = Column(String(10), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)
    SeqID = Column(String(5), nullable=False, primary_key=True)
    PDBIonCode = Column(String(3), nullable=False)
    IonID = Column(Integer, ForeignKey('Ion.ID'), nullable=False)
    ParamsFileContentID = Column(Integer, ForeignKey('FileContent.ID'), nullable=True)
    Element = Column(String(2), nullable=False)


class PDB2PDBChainMap(DeclarativeBase):
    __tablename__ = 'PDB2PDBChainMap'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBFileID1 = Column(String(10), nullable=False)
    Chain1 = Column(String(1), nullable=False)
    PDBFileID2 = Column(String(10), nullable=False)
    Chain2 = Column(String(1), nullable=False)
    SequenceIdentity = Column(DOUBLE, nullable=False)


class Project(DeclarativeBase):
    __tablename__ = 'Project'

    ID = Column(Unicode(64, collation='utf8_unicode_ci'), nullable=False, primary_key=True)
    Description = Column(Text(collation='utf8_unicode_ci'), nullable=True)


class ProjectPDBFile(DeclarativeBase):
    __tablename__ = 'ProjectPDBFile'

    PDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    ProjectID = Column(Unicode(64, collation='utf8_unicode_ci'), ForeignKey('Project.ID'), nullable=False, primary_key=True)
    Notes = Column(Text)


#########################################
#                                       #
#  Publications and associated records  #
#                                       #
#########################################


class Publication(DeclarativeBase):
    __tablename__ = 'Publication'

    ID = Column(String(64), nullable=False, primary_key=True)
    DGUnit = Column(Enum('kJ/mol', 'kcal/mol', 'cal/mol', 'fitness'), nullable=True)
    DDGConvention = Column(Enum('Rosetta','ProTherm','Unknown','Not applicable'), nullable=True)
    Notes = Column(Text, nullable=True)
    DGNotes = Column(Unicode(1024), nullable=True)
    DGUnitUsedInProTherm = Column(Enum('kcal/mol','kJ/mol'), nullable=True)
    DDGProThermSignNotes = Column(String(1024), nullable=True)
    DDGValuesNeedToBeChecked = Column(TINYINT(1), nullable=False, default=0)
    RIS = deferred(Column(Text, nullable=True))
    Title = Column(Unicode(256), nullable=True)
    Publication = Column(String(256), nullable=True)
    Volume = Column(String(8), nullable=True)
    Issue = Column(String(8), nullable=True)
    StartPage = Column(String(16), nullable=True)
    EndPage = Column(String(16), nullable=True)
    PublicationYear = Column(Integer, nullable=True)
    PublicationDate = Column(DateTime, nullable=True)
    DOI = Column(String(64), nullable=True)
    URL = Column(String(128), nullable=True)

    authors = relationship('PublicationAuthor', viewonly=True, primaryjoin="PublicationAuthor.PublicationID==Publication.ID", order_by="PublicationAuthor.AuthorOrder")

    def get_authors(self):
        '''Warning: This function should be called via a UTF-friendly session/cursor as the strings are stored as unicode in the database.'''
        athrs = []
        if self.authors:
            for a in self.authors:
                initials = u''.join([n[0] for n in (u'{0} {1}'.format(a.FirstName or u'', a.MiddleNames or u'')).strip().split()])
                if initials:
                    athrs.append(u'{0} {1}'.format(a.Surname, initials))
                else:
                    athrs.append(a.Surname)
        return ', '.join(athrs)


class PublicationAuthor(DeclarativeBase):
    __tablename__ = 'PublicationAuthor'

    PublicationID = Column(String(64), ForeignKey('Publication.ID'), nullable=False, primary_key=True)
    AuthorOrder = Column(Integer, nullable=False, primary_key=True)
    FirstName = Column(Unicode(64), nullable=False)
    MiddleNames = Column(Unicode(64), nullable=True)
    Surname = Column(Unicode(64), nullable=True)


class PublicationIdentifier(DeclarativeBase):
    __tablename__ = 'PublicationIdentifier'

    SourceID = Column(String(64), nullable=False, primary_key=True)
    ID = Column(String(256), nullable=False, primary_key=True)
    Type = Column(Enum('URL','DOI','ISSN','ESSN','PMID','MANUAL'), nullable=False)


class PublicationDDGValueLocation(DeclarativeBase):
    __tablename__ = 'PublicationDDGValueLocation'

    SourceID = Column(String(64), nullable=False, primary_key=True)
    Location = Column(String(256), nullable=False, primary_key=True)
    Notes = Column(String(512), nullable=True)


#########################################
#                                       #
#  Users                                #
#                                       #
#########################################


class User(DeclarativeBase):
    __tablename__ = 'User'

    ID = Column(String(64), nullable=False, primary_key=True)
    FirstName = Column(Unicode(64), nullable=False)
    MiddleName = Column(Unicode(64), nullable=True)
    Surname = Column(Unicode(64), nullable=True)
    Email = Column(String(80), nullable=True)


#############################################
#                                           #
#  Protein-protein complex classifications  #
#                                           #
#############################################


class PPDBMFunctionalClass(DeclarativeBase):
    __tablename__ = 'PPDBMFunctionalClass'

    ID = Column(String(2), nullable=False, primary_key=True)
    Description = Column(Unicode(128), nullable=False)


class PPIFunctionalClass(DeclarativeBase):
    __tablename__ = 'PPIFunctionalClass'

    ID = Column(String(2), nullable=False, primary_key=True)
    Description = Column(Unicode(128), nullable=False)


#########################################
#                                       #
#  Protein-protein complex definitions  #
#                                       #
#########################################


class PPComplex(DeclarativeBase):
    __tablename__ = 'PPComplex'

    ID = Column(Integer, nullable=False, primary_key=True)
    LName = Column(Unicode(256), nullable=False)
    LShortName = Column(Unicode(127), nullable=False)
    LHTMLName = Column(String(127), nullable=False)
    RName = Column(Unicode(255), nullable=False)
    RShortName = Column(Unicode(127), nullable=False)
    RHTMLName = Column(String(127), nullable=False)
    FunctionalClassID = Column(String(2), ForeignKey('PPIFunctionalClass.ID'), nullable=True)
    PPDBMFunctionalClassID = Column(String(2), ForeignKey('PPDBMFunctionalClass.ID'), nullable=True)
    PPDBMDifficulty = Column(Enum('Difficult','Medium','Rigid-body'), nullable=True)
    IsWildType = Column(TINYINT(1), nullable=False)
    WildTypeComplexID = Column(Integer, ForeignKey('PPComplex.ID'), nullable=True)
    Notes = Column(Unicode(1024), nullable=True)
    Warnings = Column(Unicode(1024), nullable=True)


    def __repr__(self):
        functional_class = None
        if self.FunctionalClassID:
            functional_class = self.FunctionalClassID
        elif self.FunctionalClassID:
            functional_class = self.PPDBMFunctionalClassID
        if functional_class:
            functional_class = '({0})'.format(functional_class)
        return ('{0} | {1} {2}'.format(self.LName, self.RName, functional_class or '')).strip()


class PPIPDBSet(DeclarativeBase):
    __tablename__ = 'PPIPDBSet'

    PPComplexID = Column(Integer, ForeignKey('PPComplex.ID'), nullable=False, primary_key=True)
    SetNumber = Column(Integer, nullable=False, primary_key=True)
    IsComplex = Column(TINYINT(1), nullable=False)
    Notes = Column(String(1024), nullable=True)

    # Relationships
    partner_chains = relationship('PPIPDBPartnerChain', viewonly=True, primaryjoin="and_(PPIPDBSet.PPComplexID==PPIPDBPartnerChain.PPComplexID, PPIPDBSet.SetNumber==PPIPDBPartnerChain.SetNumber)", order_by="PPIPDBPartnerChain.Side, PPIPDBPartnerChain.ChainIndex")

    def __repr__(self):
        pdb_ids = set([pc.PDBFileID for pc in self.partner_chains])
        if len(pdb_ids) == 1:
            d = dict(L = '', R = '')
            for pc in self.partner_chains:
                d[pc.Side] += pc.Chain
            return '#{0} {1}|{2}'.format(self.SetNumber, d['L'], d['R'])
        else:
            d = dict(L = [], R = [])
            for pc in self.partner_chains:
                d[pc.Side].append('{0} {1}'.format(pc.PDBFileID, pc.Chain))
            return '#{0} {1}|{2}'.format(self.SetNumber, ','.join(d['L']), ','.join(d['R']))


class PPIPDBPartnerChain(DeclarativeBase):
    __tablename__ = 'PPIPDBPartnerChain'

    ID = Column(Integer, nullable=False, primary_key=True)
    PPComplexID = Column(Integer, ForeignKey('PPIPDBSet.PPComplexID'), nullable=False)
    SetNumber = Column(Integer, ForeignKey('PPIPDBSet.SetNumber'), nullable=False)
    Side = Column(Enum('L','R'), nullable=False)
    ChainIndex = Column(Integer, nullable=False)
    PDBFileID = Column(String(10), ForeignKey('PDBChain.PDBFileID'), nullable=False)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False)
    NMRModel = Column(Integer, nullable=True)


class PPIConformationalChange(DeclarativeBase):
    __tablename__ = 'PPIConformationalChange'

    PPComplexID = Column(Integer, ForeignKey('PPIPDBSet.PPComplexID'), nullable=False, primary_key=True)
    ComplexSetNumber = Column(Integer, ForeignKey('PPIPDBSet.SetNumber'), nullable=False, primary_key=True)
    UnboundSetNumber = Column(Integer, ForeignKey('PPIPDBSet.SetNumber'), nullable=False, primary_key=True)
    DASA = Column(DOUBLE, nullable=True)
    I_RMSD = Column(DOUBLE, nullable=True)


##############################################
#                                            #
#  Protein-protein complex dataset mappings  #
#                                            #
##############################################


class PPIDataSetCrossmap(DeclarativeBase):
    __tablename__ = 'PPIDataSetCrossmap'

    ID = Column(Integer, nullable=False, primary_key=True)
    FromDataSetID = Column(String(128), ForeignKey('PPIDataSetDDG.DataSetID'), nullable=False)
    FromSection = Column(String(64), ForeignKey('PPIDataSetDDG.Section'), nullable=False)
    FromRecordNumber = Column(Integer, ForeignKey('PPIDataSetDDG.RecordNumber'), nullable=False)
    ToDataSetID = Column(String(128), ForeignKey('PPIDataSetDDG.DataSetID'), nullable=False)
    ToSection = Column(String(64), ForeignKey('PPIDataSetDDG.Section'), nullable=False)
    ToRecordNumber = Column(Integer, ForeignKey('PPIDataSetDDG.RecordNumber'), nullable=False)


class PPIDatabaseComplex(DeclarativeBase):
    __tablename__ = 'PPIDatabaseComplex'

    DatabaseName = Column(Enum('Kortemme & Baker', 'Kastritis et al.', 'Protein Protein Docking Benchmark v4.0', 'SKEMPI', 'ZEMu', 'CC/PBSA', 'Ben Stranges'), nullable=False, primary_key=True)
    DatabaseKey = Column(String(32), nullable=False, primary_key=True)
    PPComplexID = Column(Integer, ForeignKey('PPComplex.ID'), nullable=False)


#########################################
#                                       #
#  Protein-protein complex mutageneses  #
#                                       #
#########################################


class PPMutagenesis(DeclarativeBase):
    __tablename__ = 'PPMutagenesis'

    ID = Column(Integer, nullable=False, primary_key=True)
    PPComplexID = Column(Integer, ForeignKey('PPComplex.ID'), nullable=False)
    SKEMPI_KEY = Column(String(256), nullable=True)
    WildType = Column(TINYINT(1), nullable=False, default=0)

    # Relationships
    complex = relationship('PPComplex', viewonly=True, primaryjoin="PPMutagenesis.PPComplexID==PPComplex.ID")
    pdb_mutations = relationship('PPMutagenesisPDBMutation', viewonly=True, primaryjoin="PPMutagenesis.ID==PPMutagenesisPDBMutation.PPMutagenesisID")


class PPMutagenesisMutation(DeclarativeBase):
    __tablename__ = 'PPMutagenesisMutation'

    ID = Column(Integer, nullable=False, primary_key=True)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesis.ID'), nullable=False)
    RecordKey = Column(String(255), nullable=False)
    ProteinID = Column(String(18), nullable=True)
    ResidueIndex = Column(Integer, nullable=True)
    WildTypeAA = Column(String(1), ForeignKey('AminoAcid.Code'), nullable=False)
    MutantAA = Column(String(1), ForeignKey('AminoAcid.Code'), nullable=False)


class PPMutagenesisPDBMutation(DeclarativeBase):
    __tablename__ = 'PPMutagenesisPDBMutation'

    ID = Column(Integer, nullable=False, primary_key=True)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesisMutation.PPMutagenesisID'), ForeignKey('PPMutagenesis.ID'), nullable=False)
    PPMutagenesisMutationID = Column(Integer, ForeignKey('PPMutagenesisMutation.ID'), nullable=False)
    PPComplexID = Column(Integer, ForeignKey('PPIPDBPartnerChain.PPComplexID'), ForeignKey('PPMutagenesis.PPComplexID'), nullable=False)
    SetNumber = Column(Integer, ForeignKey('PPIPDBPartnerChain.SetNumber'), nullable=False)
    PDBFileID = Column(String(10), ForeignKey('PDBResidue.PDBFileID'), ForeignKey('PPIPDBPartnerChain.PDBFileID'), nullable=False)
    Chain = Column(String(1), ForeignKey('PDBResidue.Chain'), ForeignKey('PPIPDBPartnerChain.Chain'), nullable=False)
    WildTypeAA = Column(String(1), ForeignKey('AminoAcid.Code'), ForeignKey('PDBResidue.ResidueAA'), ForeignKey('PPMutagenesisMutation.WildTypeAA'), nullable=False)
    ResidueID = Column(String(5), ForeignKey('PDBResidue.ResidueID'), nullable=False)
    MutantAA = Column(String(1), ForeignKey('AminoAcid.Code'), ForeignKey('PPMutagenesisMutation.MutantAA'), nullable=False)

    # Parent relationships
    pdb_residue = relationship('PDBResidue', viewonly=True, primaryjoin="and_(PDBResidue.PDBFileID==PPMutagenesisPDBMutation.PDBFileID, PDBResidue.Chain==PPMutagenesisPDBMutation.Chain, PDBResidue.ResidueID==PPMutagenesisPDBMutation.ResidueID, PDBResidue.ResidueAA==PPMutagenesisPDBMutation.WildTypeAA)")
    wt_aa = relationship('AminoAcid', viewonly=True, primaryjoin="and_(AminoAcid.Code==PPMutagenesisPDBMutation.WildTypeAA)")
    mut_aa = relationship('AminoAcid', viewonly=True, primaryjoin="and_(AminoAcid.Code==PPMutagenesisPDBMutation.WildTypeAA)")
    pp_mutagenesis_mutation = relationship('PPMutagenesisMutation', viewonly=True, primaryjoin="and_(PPMutagenesisMutation.ID==PPMutagenesisPDBMutation.PPMutagenesisMutationID, PPMutagenesisMutation.PPMutagenesisID==PPMutagenesisPDBMutation.PPMutagenesisID, PPMutagenesisMutation.WildTypeAA==PPMutagenesisPDBMutation.WildTypeAA, PPMutagenesisMutation.MutantAA==PPMutagenesisPDBMutation.MutantAA)")
    ppi_pdb_partner_chain = relationship('PPIPDBPartnerChain', viewonly=True, primaryjoin="and_(PPIPDBPartnerChain.PPComplexID==PPMutagenesisPDBMutation.PPComplexID, PPIPDBPartnerChain.SetNumber==PPMutagenesisPDBMutation.SetNumber, PPIPDBPartnerChain.PDBFileID==PPMutagenesisPDBMutation.PDBFileID, PPIPDBPartnerChain.Chain==PPMutagenesisPDBMutation.Chain)")
    pdb_residue = relationship('PPMutagenesis', viewonly=True, primaryjoin="and_(PPMutagenesis.ID==PPMutagenesisPDBMutation.PPMutagenesisID, PPMutagenesis.PPComplexID==PPMutagenesisPDBMutation.PPComplexID)")


#######################################################
#                                                     #
#  Protein-protein complex experimental measurements  #
#                                                     #
#######################################################

'''
PPIDataSetDDG
PPIDataSetDDGSource
PPIExperimentalMeasurements'''


class PPIDDG(DeclarativeBase):
    __tablename__ = 'PPIDDG'

    ID = Column(Integer, nullable=False, primary_key=True)
    SecondaryID = Column(String(32), nullable=False)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesis.ID'), nullable=False)
    Publication = Column(String(64), ForeignKey('Publication.ID'), nullable=True)
    LocationOfValueInPublication = Column(String(96), nullable=True)
    DuplicateOf = Column(Integer, nullable=True)
    DDG = Column(DOUBLE, nullable=False)
    Temperature = Column(DOUBLE, nullable=True)
    pH = Column(DOUBLE, nullable=True)
    PublishedDDG = Column(DOUBLE, nullable=True)
    PublishedUnit = Column(Enum('kcal/mol','kJ/mol','cal/mol'), nullable=True)
    PublishedError = Column(String(16), nullable=True)
    DDGWasCheckedAgainstPublication = Column(TINYINT(1), nullable=False, default=0)
    Kd_wt = Column(DOUBLE, nullable=True)
    Kd_mut = Column(DOUBLE, nullable=True)
    Kon_wt = Column(DOUBLE, nullable=True)
    Kon_mut = Column(DOUBLE, nullable=True)
    Koff_wt = Column(DOUBLE, nullable=True)
    Koff_mut = Column(DOUBLE, nullable=True)
    DH_wt = Column(DOUBLE, nullable=True)
    DH_mut = Column(DOUBLE, nullable=True)
    DS_wt = Column(DOUBLE, nullable=True)
    DS_mut = Column(DOUBLE, nullable=True)
    NumberOfMeasurements = Column(Integer, nullable=True)
    Remarks = Column(Text, nullable=True)
    IsABadEntry = Column(TINYINT(1), nullable=False, default=0)
    AddedBy = Column(String(64), nullable=False)
    AddedDate = Column(DateTime, nullable=False)
    LastModifiedBy = Column(String(64), nullable=False)
    LastModifiedDate = Column(DateTime, nullable=False)

    publication = relationship('Publication', viewonly=True, primaryjoin="PPIDDG.Publication==Publication.ID")


class UserPPAnalysisSet(DeclarativeBase):
    __tablename__ = 'UserPPAnalysisSet'

    ID = Column(Integer, nullable=False, primary_key=True)
    Subset = Column(String(128), nullable=True)
    Section = Column(String(64), nullable=True)
    RecordNumber = Column(Integer, nullable=False)
    UserPPDataSetExperimentID = Column(Integer, ForeignKey('UserPPDataSetExperiment.ID'), nullable=False)
    PositiveDependentPPIDDGID = Column(Integer, ForeignKey('PPIDDG.ID'), nullable=True)
    NegativeDependentPPIDDGID = Column(Integer, ForeignKey('PPIDDG.ID'), nullable=True)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesis.ID'), nullable=False)

    positive_ddg = relationship('PPIDDG', viewonly=True, primaryjoin="PPIDDG.ID==UserPPAnalysisSet.PositiveDependentPPIDDGID")
    negative_ddg = relationship('PPIDDG', viewonly=True, primaryjoin="PPIDDG.ID==UserPPAnalysisSet.PositiveDependentPPIDDGID")


################################################################
#                                                              #
#  Protein-protein complex experimental measurements - DeltaE  #
#                                                              #
################################################################


class PPIDataSetDE(DeclarativeBase):
    __tablename__ = 'PPIDataSetDE'

    ID = Column(Integer, nullable=False, primary_key=True)
    SecondaryID = Column(String(32), nullable=False)
    DataSetID = Column(String(128), ForeignKey('DataSet.ID'), nullable=False)
    Section = Column(String(64), nullable=False)
    RecordNumber = Column(Integer, nullable=False)
    DE = Column(DOUBLE, nullable=False)
    DEUnit = Column(String(32), nullable=True)
    PublishedError = Column(String(16), nullable=True)
    NumberOfMeasurements = Column(Integer, nullable=True)

    # Identifies the mutagenesis used by the dataset. PPComplexID below is used as additional referential integrity glue
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesisPDBMutation.PPMutagenesisID'), ForeignKey('PPMutagenesis.ID'), nullable=False)

    # Identifies the complex definition used by the dataset
    PPComplexID = Column(Integer, ForeignKey('PPMutagenesisPDBMutation.PPComplexID'), ForeignKey('PPMutagenesis.PPComplexID'), ForeignKey('PPIPDBSet.PPComplexID'), ForeignKey('PPComplex.ID'), nullable=False)
    SetNumber = Column(Integer, ForeignKey('PPMutagenesisPDBMutation.SetNumber'), ForeignKey('PPIPDBSet.SetNumber'), nullable=False)

    # The PDB ID used in the publication. This is not used for any joins except to the PDBFile table i.e. we store this data for the purposes of archival or comparison but do not use it
    PublishedPDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=True)

    PossibleError = Column(TINYINT(1), nullable=False)
    Remarks = Column(Text, nullable=True)
    IsABadEntry = Column(TINYINT(1), nullable=False, default=0)
    AddedBy = Column(String(64), ForeignKey('User.ID'), nullable=False)
    AddedDate = Column(DateTime, nullable=False)
    LastModifiedBy = Column(String(64), ForeignKey('User.ID'), nullable=False)
    LastModifiedDate = Column(DateTime, nullable=False)

    mutagenesis = relationship('PPMutagenesis', viewonly=True, primaryjoin="and_(PPIDataSetDE.PPMutagenesisID==PPMutagenesis.ID, PPIDataSetDE.PPComplexID==PPMutagenesis.PPComplexID)")
    ppi_pdb_set = relationship('PPIPDBSet', viewonly=True, primaryjoin="and_(PPIPDBSet.PPComplexID==PPIDataSetDE.PPComplexID, PPIPDBSet.SetNumber==PPIDataSetDE.SetNumber)")
    pdb_mutations = relationship('PPMutagenesisPDBMutation', viewonly=True, primaryjoin="and_(PPIDataSetDE.PPMutagenesisID==PPMutagenesisPDBMutation.PPMutagenesisID, PPIDataSetDE.PPComplexID==PPMutagenesisPDBMutation.PPComplexID, PPIDataSetDE.SetNumber==PPMutagenesisPDBMutation.SetNumber)")
    complex = relationship('PPComplex', viewonly=True, primaryjoin="PPIDataSetDE.PPComplexID==PPComplex.ID")


class UserPPAnalysisSetDE(DeclarativeBase):
    __tablename__ = 'UserPPAnalysisSetDE'

    ID = Column(Integer, nullable=False, primary_key=True)
    Subset = Column(String(128), nullable=True)
    Section = Column(String(64), nullable=True)
    RecordNumber = Column(Integer, nullable=False)
    UserPPDataSetExperimentID = Column(Integer, ForeignKey('UserPPDataSetExperiment.ID'), nullable=False)
    PPIDataSetDEID = Column(Integer, ForeignKey('PPIDataSetDE.ID'), nullable=False)
    PPMutagenesisID = Column(Integer, ForeignKey('UserPPDataSetExperiment.PPMutagenesisID'), ForeignKey('PPIDataSetDE.PPMutagenesisID'), nullable=False)

    user_pp_dataset_experiment = relationship('UserPPDataSetExperiment', viewonly=True, primaryjoin="and_(UserPPAnalysisSetDE.UserPPDataSetExperimentID==UserPPDataSetExperiment.ID, UserPPAnalysisSetDE.PPMutagenesisID==UserPPDataSetExperiment.PPMutagenesisID)")
    ppi_dataset_de = relationship('PPIDataSetDE', viewonly=True, primaryjoin="and_(UserPPAnalysisSetDE.PPIDataSetDEID==PPIDataSetDE.ID, UserPPAnalysisSetDE.PPMutagenesisID==PPIDataSetDE.PPMutagenesisID)")


#######################################################
#                                                     #
#  Datasets                                           #
#                                                     #
#######################################################


class DataSet(DeclarativeBase):
    __tablename__ = 'DataSet'

    ID = Column(String(128), nullable=False, primary_key=True)
    ShortID = Column(String(32), nullable=False)
    UserID = Column(String(64), ForeignKey('User.ID'), nullable=True)
    Description = Column(String(512), nullable=True)
    DatasetType = Column(Enum('Protein stability', 'Binding affinity', 'Protein stability and binding affinity'), nullable=False)
    ContainsStabilityDDG = Column(TINYINT(1), nullable=False, default=0)
    ContainsBindingAffinityDDG = Column(TINYINT(1), nullable=False, default=0)
    ContainsBindingAffinityDE = Column(TINYINT(1), nullable=False, default=0)
    CreationDateStart = Column(DateTime, nullable=True)
    CreationDateEnd = Column(DateTime, nullable=True)
    DDGConvention = Column(Enum('Rosetta','ProTherm'), nullable=False)


class DataSetCrossmap(DeclarativeBase):
    __tablename__ = 'DataSetCrossmap'

    ID = Column(Integer, nullable=False, primary_key=True)
    FromDataSetID = Column(String(128), ForeignKey('DataSetDDG.DataSetID'), nullable=False)
    FromSection = Column(String(64), ForeignKey('DataSetDDG.Section'), nullable=False)
    FromRecordNumber = Column(Integer, ForeignKey('DataSetDDG.RecordNumber'), nullable=False)
    ToDataSetID = Column(String(128), ForeignKey('DataSetDDG.DataSetID'), nullable=False)
    ToSection = Column(String(64), ForeignKey('DataSetDDG.Section'), nullable=False)
    ToRecordNumber = Column(Integer, ForeignKey('DataSetDDG.RecordNumber'), nullable=False)


class DataSetDDG(DeclarativeBase):
    __tablename__ = 'DataSetDDG'

    ID = Column(Integer, nullable=False, primary_key=True)
    DataSetID = Column(String(128), ForeignKey('DataSet.ID'), nullable=False)
    Section = Column(String(64), nullable=False)
    RecordNumber = Column(Integer, nullable=False)
    AggregateType = Column(Enum('SingleValue','MeanValue'), nullable=False)
    PublishedValue = Column(DOUBLE, nullable=False)
    MutationIsReversed = Column(TINYINT(1), nullable=False)
    PDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=True)
    PublishedPDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=True)
    PossibleError = Column(TINYINT(1), nullable=False)
    Remark = Column(Text, nullable=True)
    CorrectionRemark = Column(Text, nullable=True)


class DataSetDDGSource(DeclarativeBase):
    __tablename__ = 'DataSetDDGSource'

    DataSetDDGID = Column(Integer, ForeignKey('DataSetDDG.ID'), nullable=False, primary_key=True)
    ExperimentAssayID = Column(Integer, ForeignKey('ExperimentAssayDDG.ExperimentAssayID'), nullable=False, primary_key=True)
    Type = Column(Enum('Unknown','DDG','DDG_H2O'), ForeignKey('ExperimentAssayDDG.Type'), nullable=False, primary_key=True)


class DataSetReference(DeclarativeBase):
    __tablename__ = 'DataSetReference'

    DataSetID = Column(String(128), ForeignKey('DataSet.ID'), nullable=False, primary_key=True)
    Publication = Column(String(64), ForeignKey('Publication.ID'), nullable=False, primary_key=True, default=u'')


#######################################################
#                                                     #
#  User datasets                                      #
#                                                     #
#######################################################


class UserDataSet(DeclarativeBase):
    __tablename__ = 'UserDataSet'

    ID = Column(Integer, nullable=False, primary_key=True)
    TextID = Column(String(32), nullable=False)
    UserID = Column(String(64), nullable=True)
    Description = Column(String(512), nullable=True)
    DatasetType = Column(Enum('Protein stability','Binding affinity'), nullable=False)
    AnalyzeDDG = Column(TINYINT(1), nullable=False, default=1) # should be set to True if analysis includes DDG values
    AnalyzeDE = Column(TINYINT(1), nullable=False, default=0) # should be set to True if analysis includes DE (Delta energy e.g. SSM) values
    FirstCreated = Column(DateTime, nullable=True)
    LastModified = Column(DateTime, nullable=True)


class UserPPDataSetExperiment(DeclarativeBase):
    __tablename__ = 'UserPPDataSetExperiment'

    ID = Column(Integer, nullable=False, primary_key=True)
    UserDataSetID = Column(Integer, ForeignKey('UserDataSet.ID'), nullable=False)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesis.ID'), nullable=False)
    PDBFileID = Column(String(10), ForeignKey('PPIPDBPartnerChain.PDBFileID'), nullable=False)
    PPComplexID = Column(Integer, ForeignKey('PPIPDBPartnerChain.PPComplexID'), ForeignKey('PPIPDBSet.PPComplexID'), ForeignKey('PPComplex.ID'), nullable=False)
    SetNumber = Column(Integer, ForeignKey('PPIPDBPartnerChain.SetNumber'), ForeignKey('PPIPDBSet.SetNumber'), nullable=False)
    IsComplex = Column(TINYINT(1), ForeignKey('PPIPDBSet.IsComplex'), nullable=False)

    # Parent relationships
    ppi_pdb_partner_chain = relationship('PPIPDBPartnerChain', viewonly=True, primaryjoin="and_(PPIPDBPartnerChain.PDBFileID==UserPPDataSetExperiment.PDBFileID, PPIPDBPartnerChain.PPComplexID==UserPPDataSetExperiment.PPComplexID, PPIPDBPartnerChain.SetNumber==UserPPDataSetExperiment.SetNumber)")
    ppi_pdb_set = relationship('PPIPDBSet', viewonly=True, primaryjoin="and_(PPIPDBSet.PPComplexID==UserPPDataSetExperiment.PPComplexID, PPIPDBSet.SetNumber==UserPPDataSetExperiment.SetNumber, PPIPDBSet.IsComplex==UserPPDataSetExperiment.IsComplex)")
    complex = relationship('PPComplex', viewonly=True, primaryjoin="UserPPDataSetExperiment.PPComplexID==PPComplex.ID")
    user_dataset = relationship('UserDataSet', viewonly=True, primaryjoin="UserPPDataSetExperiment.UserDataSetID==UserDataSet.ID")
    mutagenesis = relationship('PPMutagenesis', viewonly=True, primaryjoin="UserPPDataSetExperiment.PPMutagenesisID==PPMutagenesis.ID")

    def __repr__(self):
        mutations = []
        for m in self.mutagenesis.pdb_mutations:
            if m.PPComplexID == self.PPComplexID and m.SetNumber == self.SetNumber and m.PDBFileID == self.PDBFileID:
                mutations.append('{0} {1}{2}{3}'.format(m.Chain, m.WildTypeAA, m.ResidueID.strip(), m.MutantAA))
        mutations = ', '.join(mutations)
        return 'UserPPDataSetExperiment #{0} ({1}). Complex = {2}, Partners = {3}. Mutations: {4}'.format(self.ID, self.user_dataset.TextID, self.complex, self.ppi_pdb_set, mutations or 'N/A')


class UserPPDataSetExperimentTag(DeclarativeBase):
    __tablename__ = 'UserPPDataSetExperimentTag'

    UserPPDataSetExperimentID = Column(Integer, nullable=False, primary_key=True)
    Tag = Column(String(64), nullable=False, primary_key=True)


#######################################################
#                                                     #
#  Protocol tables                                    #
#                                                     #
#######################################################


class Protocol(DeclarativeBase):
    __tablename__ = 'Protocol'

    ID = Column(String(256), nullable=False, primary_key=True)
    Description = Column(Text, nullable=False)
    ClassName = Column(String(256), nullable=True)
    Publication = Column(String(64), ForeignKey('Publication.ID'), nullable=True)

#todo: add missing protocol-related tables


#######################################################
#                                                     #
#  Scoring method tables                              #
#                                                     #
#######################################################


class ScoreMethod(DeclarativeBase):
    __tablename__ = 'ScoreMethod'

    ID = Column(Integer, nullable=False, primary_key=True)
    MethodName = Column(Unicode(64, collation='utf8_unicode_ci'), nullable=False)
    MethodType = Column(Unicode(64, collation='utf8_unicode_ci'), nullable=False)
    Parameters = Column(Unicode(64, collation='utf8_unicode_ci'), nullable=True)
    Authors = Column(Unicode(255, collation='utf8_unicode_ci'), nullable=False)
    Notes = Column(Unicode(512, collation='utf8_unicode_ci'), nullable=True)

    def __repr__(self):
        sp = self.Parameters or ''
        return 'Score method #{0}: {1}, {2}{3} ({4})'.format(self.ID, self.MethodName, self.MethodType, sp, self.Authors)


#######################################################
#                                                     #
#  Predictions tables                                 #
#                                                     #
#######################################################


class PredictionSet(DeclarativeBase):
    __tablename__ = 'PredictionSet'

    ID = Column(String(48), nullable=False, primary_key=True)
    Status = Column(Enum('halted','active'), nullable=False)
    Priority = Column(Integer, nullable=False, default=5)
    ProteinStability = Column(TINYINT(4), nullable=False, default=1)
    BindingAffinity = Column(TINYINT(4), nullable=False, default=0)
    BatchSize = Column(Integer, nullable=False, default=40)
    SeriesName = Column(String(128), nullable=True)
    SeriesColor = Column(String(6), nullable=False, default=u'ff0000')
    SeriesAlpha = Column(DOUBLE, nullable=False, default=1)
    Description = Column(String(256), nullable=False)
    CanBeDeleted = Column(TINYINT(1), nullable=False, default=0)
    EntryDate = Column(TIMESTAMP, nullable=False)

    # relationships
    ppi_predictions = relationship('PredictionPPI', viewonly=True, primaryjoin="PredictionSet.ID==PredictionPPI.PredictionSet")


class Prediction(DeclarativeBase):
    __tablename__ = 'Prediction'

    ID = Column(Integer, nullable=False, primary_key=True)
    ExperimentID = Column(Integer, nullable=False)
    UserDataSetExperimentID = Column(Integer, nullable=True)
    PredictionSet = Column(Unicode(48), nullable=False)
    ProtocolID = Column(Unicode(256), nullable=False)
    Cost = Column(DOUBLE, nullable=False, default=0)
    KeptHETATMLines = Column(TINYINT(1), nullable=True)
#StrippedPDB = Column(Text, nullable=True)
#ResidueMapping = Column(Text, nullable=True)
#InputFiles = Column(Text, nullable=True)
#Description = Column(Text, nullable=True)
#ScoreVersion = Column(DOUBLE, nullable=False, default=0.23)
#ddG = Column(Text, nullable=True)
#Scores = Column(Text, nullable=True)
    NumberOfMeasurements = Column(Integer, nullable=False, default=1)
    EntryDate = Column(TIMESTAMP, nullable=False)
    cryptID = Column(Unicode(50), nullable=True)
    StartDate = Column(DateTime, nullable=True)
    EndDate = Column(DateTime, nullable=True)
    Status = Column(Enum('queued','active','done','failed','postponed'), nullable=True)
    Errors = Column(Text, nullable=True)
    AdminCommand = Column(Unicode(20), nullable=True)
    ExtraParameters = Column(Text, nullable=False)
    maxvmem = Column(DOUBLE, nullable=False)
    DDGTime = Column(DOUBLE, nullable=False)
#StoreOutput = Column(TINYINT(1), nullable=False, default=0)


prediction_ppi_clone_null_fields = ['EntryDate', 'StartDate', 'EndDate', 'Errors', 'AdminCommand', 'maxvmem', 'DDGTime']


class PredictionPPI(DeclarativeBase):
    __tablename__ = 'PredictionPPI'

    ID = Column(Integer, nullable=False, primary_key=True)
    PredictionSet = Column(String(48), ForeignKey('PredictionSet.ID'), nullable=False)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesis.ID'), nullable=False)
    UserPPDataSetExperimentID = Column(Integer, ForeignKey('UserPPDataSetExperiment.ID'), nullable=False)
    ProtocolID = Column(String(256), ForeignKey('Protocol.ID'), nullable=True)
    JSONParameters = Column(Text, nullable=True) # todo: this should be deferred if possible but last time we set it as deferred it caused a bug (it was quicker to just remove the deferred call than debug)
    EntryDate = Column(TIMESTAMP, nullable=False)
    StartDate = Column(DateTime, nullable=True)
    EndDate = Column(DateTime, nullable=True)
    Status = Column(Enum('queued','active','done','failed','postponed'), nullable=True)
    Errors = deferred(Column(Text, nullable=True))
    AdminCommand = Column(String(20), nullable=True)
    Cost = Column(DOUBLE, nullable=False, default=0)
    maxvmem = Column(DOUBLE, nullable=True)
    DDGTime = Column(DOUBLE, nullable=True)
    ExtraParameters = deferred(Column(Text, nullable=True))
    KeptHETATMLines = Column(TINYINT(1), nullable=True)
    NumberOfMeasurements = Column(Integer, nullable=False, default=1)
    DevelopmentProtocolID = Column(Integer, nullable=True)

    # Relationships
    files = relationship('PredictionPPIFile', viewonly=True, primaryjoin="PredictionPPI.ID==PredictionPPIFile.PredictionPPIID")
    mutagenesis = relationship('PPMutagenesis', viewonly=True, primaryjoin="PredictionPPI.PPMutagenesisID==PPMutagenesis.ID")
    user_dataset_experiment = relationship('UserPPDataSetExperiment', viewonly=True, primaryjoin="PredictionPPI.UserPPDataSetExperimentID==UserPPDataSetExperiment.ID")


    def __repr__(self):
        try:
            return 'Prediction #{0}, {1} ({2}): Mutagenesis #{3}.\n{4}.'.format(self.ID, self.PredictionSet, self.Status, self.mutagenesis.ID, self.user_dataset_experiment)
        except:
            raise Exception('The __repr__ function failed on this Prediction. Was it created using the .clone() method? This function has a known issue.')


    def clone(self, prediction_set):
        '''Returns a new fresh PredictionSet object to be inserted into the database.

           Warning: This code returns a PredictionSet object but the relationships cannot be called on this object e.g.
              print new_prediction.mutagenesis.ID
           will fail. I am guessing that the relationships are set up on object instantiation.

           todo: I think the correct way to do this is to use the expunge and transient functions on detached objects.
                 See http://docs.sqlalchemy.org/en/rel_1_1/orm/session_api.html#sqlalchemy.orm.session.make_transient
                 and http://stackoverflow.com/questions/20112850/sqlalchemy-clone-table-row-with-relations?lq=1
        '''
        assert(prediction_set != self.PredictionSet)
        fieldnames = [c.name for c in list(sqlalchemy_inspect(PredictionPPI).columns)]
        new_prediction = PredictionPPI()
        for c in fieldnames:
            if c in prediction_ppi_clone_null_fields:
                setattr(new_prediction, c, None)
            else:
                setattr(new_prediction, c, getattr(self, c))
        new_prediction.ID = None
        new_prediction.PredictionSet = prediction_set
        new_prediction.Status = 'queued'
        return new_prediction


class PredictionPPIFile(DeclarativeBase):
    __tablename__ = 'PredictionPPIFile'

    ID = Column(Integer, nullable=False, primary_key=True)
    PredictionPPIID = Column(Integer, ForeignKey('PredictionPPI.ID'), nullable=False)
    FileContentID = Column(Integer, ForeignKey('FileContent.ID'), nullable=False)
    Filename = Column(String(255), nullable=False)
    Filetype = Column(Enum('Image','MOL','Mutfile','Other','Params','PDB','PDF','Resfile','Text','RosettaPDBMapping'), nullable=False)
    FileRole = Column(String(64), nullable=False)
    Stage = Column(Enum('Input','Output','Analysis'), nullable=True)

    # Relationships
    content = relationship('FileContent', viewonly=True, primaryjoin="PredictionPPIFile.FileContentID==FileContent.ID")

    def clone(self, prediction_id):
        '''Returns a new fresh PredictionPPIFile object to be inserted into the database.

           Warning: This code returns a PredictionPPIFile object but the relationships cannot be called on this object e.g.
              print new_prediction_file.content
           will fail. I am guessing that the relationships are set up on object instantiation. See comment on PredictionPPI.clone() above.
        '''

        assert(prediction_id != self.PredictionPPIID)
        fieldnames = [c.name for c in list(sqlalchemy_inspect(PredictionPPIFile).columns)]
        new_prediction_file = PredictionPPIFile()
        for c in fieldnames:
            setattr(new_prediction_file, c, getattr(self, c))
        new_prediction_file.ID = None
        new_prediction_file.PredictionPPIID = prediction_id
        return new_prediction_file


class PredictionPPIStructureScore(DeclarativeBase):
    __tablename__ = 'PredictionPPIStructureScore'

    ID = Column(Integer, nullable=False, primary_key=True)
    PredictionPPIID = Column(Integer, ForeignKey('PredictionPPI.ID'), nullable=False)
    ScoreMethodID = Column(Integer, ForeignKey('ScoreMethod.ID'), nullable=False)
    ScoreType = Column(Enum('DDG','WildTypeLPartner','WildTypeRPartner','WildTypeComplex','MutantLPartner','MutantRPartner','MutantComplex'), nullable=False)
    StructureID = Column(Integer, nullable=True)
    DDG = Column(DOUBLE, nullable=True)
    total = Column(DOUBLE, nullable=True)
    dslf_fa13 = Column(DOUBLE, nullable=True)
    dslf_ca_dih = Column(DOUBLE, nullable=True)
    dslf_cs_ang = Column(DOUBLE, nullable=True)
    dslf_ss_dih = Column(DOUBLE, nullable=True)
    dslf_ss_dst = Column(DOUBLE, nullable=True)
    fa_pair = Column(DOUBLE, nullable=True)
    fa_atr = Column(DOUBLE, nullable=True)
    fa_dun = Column(DOUBLE, nullable=True)
    fa_elec = Column(DOUBLE, nullable=True)
    fa_intra_rep = Column(DOUBLE, nullable=True)
    fa_rep = Column(DOUBLE, nullable=True)
    fa_sol = Column(DOUBLE, nullable=True)
    hbond_bb_sc = Column(DOUBLE, nullable=True)
    hbond_lr_bb = Column(DOUBLE, nullable=True)
    hbond_sc = Column(DOUBLE, nullable=True)
    hbond_sr_bb = Column(DOUBLE, nullable=True)
    omega = Column(DOUBLE, nullable=True)
    p_aa_pp = Column(DOUBLE, nullable=True)
    pro_close = Column(DOUBLE, nullable=True)
    rama = Column(DOUBLE, nullable=True)
    ref = Column(DOUBLE, nullable=True)
    yhh_planarity = Column(DOUBLE, nullable=True)

    # beta_nov_16 columns
    fa_dun_dev = Column(DOUBLE, nullable=True)
    fa_dun_rot = Column(DOUBLE, nullable=True)
    fa_dun_semi = Column(DOUBLE, nullable=True)
    fa_intra_atr_xover4 = Column(DOUBLE, nullable=True)
    fa_intra_elec = Column(DOUBLE, nullable=True)
    fa_intra_rep_xover4 = Column(DOUBLE, nullable=True)
    fa_intra_sol_xover4 = Column(DOUBLE, nullable=True)
    hxl_tors = Column(DOUBLE, nullable=True)
    lk_ball = Column(DOUBLE, nullable=True)
    lk_ball_bridge = Column(DOUBLE, nullable=True)
    lk_ball_bridge_uncpl = Column(DOUBLE, nullable=True)
    lk_ball_iso = Column(DOUBLE, nullable=True)
    rama_prepro = Column(DOUBLE, nullable=True)
    cart_bonded = Column(DOUBLE, nullable=True)


#######################################################
#                                                     #
#  Analysis tables                                    #
#                                                     #
#######################################################


class AnalysisDataFrame(DeclarativeBase):
    __tablename__ = 'AnalysisDataFrame'

    ID = Column(Integer, nullable=False, primary_key=True)
    PredictionSet = Column(String(48), ForeignKey('PredictionSet.ID'), nullable=False)
    DataFrameType = Column(Enum('Stability','Binding affinity'), nullable=False, default=u'Binding affinity')
    ContainsExperimentalData = Column(TINYINT(1), nullable=False, default=1)
    ScoreMethodID = Column(Integer, ForeignKey('ScoreMethod.ID'), nullable=False)
    UseSingleReportedValue = Column(TINYINT(1), nullable=True, default=0)
    TopX = Column(Integer, nullable=False, default=3)
    BurialCutoff = Column(DOUBLE, nullable=False, default=0.25)
    StabilityClassicationExperimentalCutoff = Column(DOUBLE, nullable=False, default=1)
    StabilityClassicationPredictedCutoff = Column(DOUBLE, nullable=False, default=1)
    IncludesDerivedMutations = Column(TINYINT(1), nullable=False, default=0)
    DDGAnalysisType = Column(String(16), nullable=True, default=u'DDG_Top3')
    SeriesName = Column(String(128), nullable=False)
    SeriesColor = Column(String(6), nullable=False)
    SeriesAlpha = Column(DOUBLE, nullable=False, default=1)
    Description = Column(String(256), nullable=False)
    Credit = Column(String(128), nullable=False)
    DDGAnalysisTypeDescription = Column(String(256), nullable=False)
    PandasHDFStore = deferred(Column(LONGBLOB, nullable=False))
    CreationDate = Column(TIMESTAMP, nullable=False)

    # Relationships
    score_method = relationship('ScoreMethod', viewonly=True, primaryjoin="ScoreMethod.ID ==AnalysisDataFrame.ScoreMethodID")
    prediction_set = relationship('PredictionSet', viewonly=True, primaryjoin="PredictionSet.ID ==AnalysisDataFrame.PredictionSet")


    def __repr__(self):
        s = []
        s.append('Analysis dataframe: {0}, {1} ({2})'.format(self.PredictionSet, self.DataFrameType, self.CreationDate))
        s.append(str(self.score_method))
        s.append('TopX: {0}'.format(self.TopX))
        s.append('Burial cutoff: {0}'.format(self.BurialCutoff))
        s.append('Experimental cutoff: {0}'.format(self.StabilityClassicationExperimentalCutoff))
        s.append('Prediction cutoff: {0}'.format(self.StabilityClassicationPredictedCutoff))
        return '\n'.join(s)


    def get_dataframe_info(self):
        mem_zip = StringIO.StringIO()
        mem_zip.write(self.PandasHDFStore)
        mem_zip.seek(0)
        hdf_store_blob = gzip.GzipFile(fileobj = mem_zip, mode='rb').read()

        try:
            # read_hdf does not currently (as of pandas.__version__ == 0.17.0) accept stream objects so we write to file
            #mem_unzipped = StringIO.StringIO()
            #mem_unzipped.write(hdf_store_blob)
            #mem_unzipped.seek(0)
            analysis_pandas_input_filepath = write_temp_file('/tmp', hdf_store_blob, ftype = 'wb')

            df = pandas.read_hdf(analysis_pandas_input_filepath, 'dataframe')
            store = pandas.HDFStore(analysis_pandas_input_filepath)

            # Defensive programming in case the format changes
            scalar_adjustments, ddg_analysis_type, ddg_analysis_type_description, analysis_sets = None, None, None, None
            try: scalar_adjustments = store['scalar_adjustments'].to_dict()
            except: pass
            try: ddg_analysis_type = store['ddg_analysis_type'].to_dict()['ddg_analysis_type']
            except: pass
            try: ddg_analysis_type_description = store['ddg_analysis_type_description'].to_dict()['ddg_analysis_type_description']
            except: pass
            if scalar_adjustments:
                analysis_sets = scalar_adjustments.keys()

            d = dict(
                AnalysisDataFrameID = self.ID,
                dataframe = df,
                scalar_adjustments = scalar_adjustments,
                analysis_type = ddg_analysis_type,
                analysis_type_description = ddg_analysis_type_description,
                analysis_sets = analysis_sets,
                score_method_id = self.score_method.ID,
                prediction_set = self.prediction_set.ID,
                top_x = self.TopX,
            )
            os.remove(analysis_pandas_input_filepath)
            return d
        except Exception, e:
            if os.path.exists(analysis_pandas_input_filepath):
                os.remove(analysis_pandas_input_filepath)
            raise Exception('An exception occurred reading the dataframe: {0}\n.{1}'.format(str(e), traceback.format_exc()))


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
    sc = MySQLSchemaConverter(sys_settings.database.username, sys_settings.database.hostname, sys_settings.database.database, sys_settings.database.password, sys_settings.database.port, sys_settings.database.socket)
    #sc.get_sqlalchemy_schema(['PDBFile', 'PDBChain', 'PDBMolecule', 'PDBMoleculeChain', 'PDBResidue'])
    sc.get_sqlalchemy_schema(tablenames)


def test_schema_against_database_instance(DDG_db):
    '''Make sure that our SQLAlchemy definitions match the database. This should be run by the API prior to connection
       as it lets the admin know that they need to update the schema here (use generate_sqlalchemy_definition to update
       the schema).'''
    if getpass.getuser() == 'kyleb':
        colortext.rastaprint('Hello Kyle!')
        return
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
