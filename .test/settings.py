#!/usr/bin/python2
"""
.test/settings.py
A simple test script for the settings module.

Created by Shane O'Connor 2016.
"""


import kddg.api.settings as settings
sys_settings = settings.load()

print(sys_settings.database)
print(sys_settings.database.username)
