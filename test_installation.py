#!/usr/bin/python2.4
# encoding: utf-8
"""
test_installation.py
A simple test to check that the database is correctly configured. New installations will require the user to configure the
kddg/api/settings.json file which is created on the first attempt at connecting to a database.

Created by Shane O'Connor 2016.
Copyright (c) 2016 Shane O'Connor. All rights reserved.
"""

import pprint

print('\nImporting PPI module.')
from kddg.api.ppi import get_interface_with_config_file

print('\nAttempting to set up a database connection.')
ppi_api = get_interface_with_config_file()

print('\nRetrieving data from the database.')
pprint.pprint(ppi_api.get_amino_acid_details())

print('\nSuccess.\n')
