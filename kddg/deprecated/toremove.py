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

#if __name__ == '__main__':
#    sys.path.insert(0, '../../klab')

from klab import colortext
from klab.bio.pdb import PDB
from klab.bio.basics import ChainMutation
from klab.fs.fsio import read_file, write_temp_file, open_temp_file, write_file
from klab.bio.pfam import Pfam
from klab.bio.dssp import MonomerDSSP, ComplexDSSP, MissingAtomException
from klab.bio.ligand import Ligand, PDBLigand
from klab.bio.pdbtm import PDBTM
from klab.db.sqlalchemy_interface import get_single_record_from_query, get_or_create_in_transaction

from db_schema import test_schema_against_database_instance
from db_schema import PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue, LigandDescriptor, LigandIdentifier, LigandSynonym, PDBLigand
from db_schema import Ligand as DBLigand
#from db_schema import Publication, PublicationAuthor, PublicationIdentifier
from api_layers import *
from db_api import ddG, PartialDataException, SanityCheckException
import ddgdbapi

rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease'
rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database'
p = PDB(read_file('/kortemmelab/data/kyleb/ddg_numbering_for_shane/24548-data/1CBW_FGHI.pdb'))
#p.construct_pdb_to_rosetta_residue_map(rosetta_scripts_path, rosetta_database_path)
p.construct_pdb_to_rosetta_residue_map(rosetta_scripts_path, rosetta_database_path, extra_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res')
pprint.pprint(p.get_atom_sequence_to_rosetta_map())
pprint.pprint(p.rosetta_sequences)

from ppi_api import get_interface as get_ppi_interface
ppi_api = get_ppi_interface(read_file('../misc/ddgdb.pw'),
                                rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
                                rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database')
content = ppi_api.DDG_db.execute_select('SELECT Content FROM PDBFile WHERE ID="1CBW"')[0]['Content']
print(content)
write_file('/tmp/ddginterface/1CBW_FGHI_db.pdb', content)