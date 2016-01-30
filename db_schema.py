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
from sqlalchemy.orm import deferred

#if __name__ == '__main__':
#    sys.path.insert(0, '../../klab')

from klab.db.sqlalchemy_interface import MySQLSchemaConverter
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
    residues = relationship('Publication', viewonly=True, primaryjoin="PDBFile.Publication==Publication.ID")

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

    PDBFileID = Column(String(10), nullable=False, primary_key=True)
    Chain = Column(String(1), nullable=False, primary_key=True)
    SeqID = Column(String(5), nullable=False, primary_key=True)
    PDBIonCode = Column(String(3), nullable=False)
    IonID = Column(Integer, nullable=False)
    ParamsFileContentID = Column(Integer, nullable=True)
    Element = Column(String(2), nullable=False)


#########################################
#                                       #
#  Publications and associated records  #
#                                       #
#########################################


class Publication(DeclarativeBase):
    __tablename__ = 'Publication'

    ID = Column(String(64), nullable=False, primary_key=True)
    DGUnit = Column(Enum('kJ/mol','kcal/mol','cal/mol'), nullable=True)
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


class PublicationAuthor(DeclarativeBase):
    __tablename__ = 'PublicationAuthor'

    PublicationID = Column(String(64), nullable=False, primary_key=True)
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
        d = dict(L = '', R = '')
        pdb_ids = set([pc.PDBFileID for pc in self.partner_chains])
        if len(pdb_ids) == 1:
            for pc in self.partner_chains:
                d[pc.Side] += pc.Chain
            return '{0}|{1}'.format(d['L'], d['R'])
        else:
            for pc in self.partner_chains:
                d[pc.Side] += '{0} {1}'.format(pc.PDBFileID, pc.Chain)
            return '{0}|{1}'.format(','.join(d['L']), ','.join(d['R']))


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

    DatabaseName = Column(Enum('Kortemme & Baker','Kastritis et al.','Protein Protein Docking Benchmark v4.0','SKEMPI','ZEMu','CC/PBSA','Ben Stranges'), nullable=False, primary_key=True)
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

'''PPIDDG
PPIDataSetDDG
PPIDataSetDDGSource
PPIExperimentalMeasurements'''

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
    DatasetType = Column(Enum('Protein stability','Binding affinity'), nullable=True)
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


class PredictionPPI(DeclarativeBase):
    __tablename__ = 'PredictionPPI'

    ID = Column(Integer, nullable=False, primary_key=True)
    PredictionSet = Column(String(48), ForeignKey('PredictionSet.ID'), nullable=False)
    PPMutagenesisID = Column(Integer, ForeignKey('PPMutagenesis.ID'), nullable=False)
    UserPPDataSetExperimentID = Column(Integer, ForeignKey('UserPPDataSetExperiment.ID'), nullable=False)
    ProtocolID = Column(String(256), ForeignKey('Protocol.ID'), nullable=True)
    JSONParameters = Column(Text, nullable=True)
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
        return 'Prediction #{0}, {1} ({2}): Mutagenesis #{3}.\n{4}.'.format(self.ID, self.PredictionSet, self.Status, self.mutagenesis.ID, self.user_dataset_experiment)


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
#qry = 'SELECT {0}File.*, FileContent.Content, FileContent.MIMEType, FileContent.Filesize, FileContent.MD5HexDigest FROM {0}File INNER JOIN FileContent ON FileContentID=FileContent.ID WHERE {0}ID=%s'.format(*params)


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
    generate_sqlalchemy_definition(['UserDataSet'])

    #generate_sqlalchemy_definition(['AminoAcid'])
    sys.exit(0)
    from ppi_api import get_interface as get_ppi_interface
    ppi_api = get_ppi_interface(read_file(os.path.join('..', 'pw')).strip())
    test_schema_against_database_instance(ppi_api.DDG_db)
