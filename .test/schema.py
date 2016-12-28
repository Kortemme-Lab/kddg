#!/usr/bin/python2
"""
.test/schema.py
Test script to generate skeleton SQLAlchemy code for schema.py. This is useful for when the schema is changed.

Created by Shane O'Connor 2016.
"""

import sys

from klab.db.sqlalchemy_interface import MySQLSchemaConverter
from kddg.api.data import DataImportInterface
from kddg.api import schema as dbmodel
from kddg.api import settings


sys_settings = settings.load()


data_api = None


def get_fresh_session():
    global data_api
    if not data_api:
        data_api = DataImportInterface.get_interface_with_config_file(cache_dir = sys_settings.cache.cache_dir, echo_sql = False)
        dbmodel.test_schema_against_database_instance(data_api.DDG_db)
    return data_api.get_session(new_session = True)


if __name__ == '__main__':
    tablenames = sys.argv[1:]
    if tablenames:
        dbmodel.generate_sqlalchemy_definition(tablenames)
    else:
        print('\nCall this script with a space-separated list of table names from the DDG database e.g.\n\tpython {0} AminoAcid PDBFile\n'.format(__file__))
        sys.exit(1)
