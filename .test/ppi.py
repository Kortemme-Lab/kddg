#!/usr/bin/python2
"""
.test/ppi.py
Simple tests of the ppi interface.

Created by Shane O'Connor 2016.
"""


import pprint

from kddg.api.ppi import get_interface_with_config_file


ppi_api = get_interface_with_config_file()

pprint.pprint(ppi_api.get_amino_acid_details())

ppi_api.help()