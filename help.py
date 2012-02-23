#!/usr/bin/python2.4
# encoding: utf-8
"""
help.py
Help functionality for the ddG database.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import sys
from string import join
from common.ddgproject import FieldNames, StdCursor, ddGDatabase
from common import colortext
import ddglib#.filter import *
import inspect

dbfields = FieldNames()

def _get_ResultSetFilter_data():
	s_module = "ddglib.ddgfilters"
	#clsmembers = inspect.getmembers(sys.modules[s_module], lambda member: member.__module__ == s_module and inspect.isclass)
	#clsmembers = inspect.getmembers(sys.modules[s_module], lambda member: member.inspect.isclass(member))
	
	m_filters = []
	m_resultsets = []
	d_filters = {}
	
	s_module = "ddglib.ddgfilters"
	for m in inspect.getmembers(sys.modules[s_module]):
		o = m[1]
		if inspect.isclass(o) and o.__module__ == s_module:
			classnm = m[0]
			if classnm.find("Filter") != -1:
				d = {"name" : classnm, "class" : o}
				m_filters.append(d)
				d_filters[classnm] = d 
			elif classnm.find("ResultSet") != -1:
				e_Filter = "%sFilter" % classnm[:classnm.find("ResultSet")]
				e_Filter = o.allowed_filters
				#"%sFilter" % classnm[:classnm.find("ResultSet")]
				m_resultsets.append({"name" : classnm, "class" : o, "filter" : e_Filter})
			else:
				colortext.error("Unknown class '%s' found." % classnm)
	return m_filters, m_resultsets, d_filters

def _print_lines(helplines):
	for linepair in helplines:
		colortext.printf(linepair[0], color=linepair[1])

def ShowDatabaseStructure():
	'''Extracts the database structure from the MySQL database and prints it out.'''
	help = []
	
	ddGdb = ddGDatabase()
	help.append(("\n* Database structure *", "white"))
	help.append(("The tables of the database are as follows (orange = primary key, blue = foreign key):\n", "silver"))
	tablenames = sorted([r[0] for r in (ddGdb.execute("SHOW TABLES", cursorClass = StdCursor))])
	for tbl in tablenames:
		help.append((tbl, "green"))
		fieldnames = ddGdb.execute("SHOW COLUMNS FROM %s" % tbl)
		for fld in fieldnames:
			fname = fld["Field"]
			ftype = fld["Type"]
			fdefault = fld["Default"]
			fextra = fld["Extra"]
			fkey = fld["Key"]
			fCanBeNull = fld["Null"]
			padding = " " * (32 - len(fname))
			str = "\t%s%s%s" % (fname, padding, ftype)
			if fkey == "PRI":
				help.append((str, "orange"))
			elif fkey == "MUL":
				help.append((str, "lightblue"))
			else:
				help.append((str, "yellow"))
	
	_print_lines(help)
	
def ShowResultSet():
	'''Explains how to extract sets of data from the database.'''
	help = []
	m_filters, m_resultsets, d_filters = _get_ResultSetFilter_data()
	
	help.append(('''
* ResultSets *''', "white"))

	help.append(('''
ResultSets are used to store data retrieved from the database. They are the result of an SQL query or
stored procedure call. To create a ResultSet, use the following format:''', 'silver'))
	help.append(('''  SomeResultSet(db, <string:SQL>, <tuple:parameters>, <list:AdditionalIDs>)''', "lightblue"))
	help.append(('''The argument db is a datbase connection retrieved from calling common.ddgproject.ddGDatabase().

The SQL string and associated parameters are optional. If they are not supplied and neither is AdditionalIDs
then all records will be returned. If only AdditionalIDs is supplied then only records with those primary keys
in the list are returned. If an SQL string is supplied then that query will be run instead with any associated
parameters. If AdditionalIDs is also supplied then the associated records are included in the result set.

Examples:
1. Return all structures:''', "silver"))
	help.append(('''       sr = StructureResultSet(ddGdb)''', "lightblue"))
	help.append(('''2. Get four specific structures:''', "silver"))
	help.append(('''       sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])''', "lightblue"))
	help.append(('''3. Get all predictions with IDs >= 12395 from the 'testrun' prediction set:''', "silver"))   
	help.append(('''       pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s AND ID>=%s", parameters = ("testrun", 12595))''', "lightblue"))
	help.append(('''
Note in the last example that the arguments to the SQL string are placed in the parameters tuple and %s is
used for both strings and numbers. This is the recommended approach as it avoids any manual mistakes made 
from badly formatting the SQL strings. 

An alternative method of narrowing down results is by using Filters. Each ResultSet has associated Filters that
may be applied to return a filtered set of results. This approach is no more powerful than using SQL (and the 
implementation is typically slower) but it is hopefully easier to use. It also abstracts away the database structure
so will hopefully survive any changes made to the database design. Filters are explained below.''', "silver"))

	help.append(('''
The following ResultSet classes and associated filters are available:\n''', "silver"))
	for rs in sorted(m_resultsets):
		help.append(("\t%s" % rs["name"], "green"))
		print(rs["filter"])
		if rs["filter"]:
			avail_f = [d_filters[rs_f.__name__]["name"] for rs_f in rs["filter"] if d_filters.get(rs_f.__name__)]
			help.append(("\t\t%s" % join(avail_f, ", "), "orange"))
			
	_print_lines(help)
		
def ShowFilter():
	'''Explains how to filter sets of data extracted from the database (ResultSets).'''
	help = []
	m_filters, m_resultsets, d_filters = _get_ResultSetFilter_data()
	
	help.append(('''
* Filters *''', "white"))
	help.append(('''
Database records can be retrieved by SQL queries or stored procedures but also through the use of filters.

---
sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])

sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))
results = sr.getFilteredResults()  

We run an SQL query if the SQL parameter has been specified or if AdditionalIDs is empty.
			# If AdditionalIDs is not empty and the SQL is blank then we do NOT run the SQL query and
			# instead use AdditionalIDs as the record keys.
			
Note that it is less efficient to use union filters than to set the properties using the set methods. However,
speed may not be an issue and they are quicker to write.
 
Filters are used to narrow to
---			
 The following Filters are available:\n''', "silver"))
	
	for filter in sorted(m_filters):
		help.append(("%s" % filter["name"], "green"))
		c = filter["class"]
		
		class_members = c.getMembers()
		
		help.append(("\tAttributes", "orange"))
		for x in sorted(class_members["attributes"]):
			help.append(("\t\t%s" % x[0], "silver"))

		help.append(("\tMethods", "orange"))
		for x in sorted(class_members["methods"]):
			arglist = join(inspect.getargspec(x[1])[0][1:], ", ")
			help.append(("\t\t%s(%s)" % (x[0], arglist), "silver"))
			
		help.append(("\tFunctions", "orange"))
		for x in sorted(class_members["functions"]):
			arglist = join(inspect.getargspec(x[1])[0], ", ")
			help.append(("\t\t%s(%s)" % (x[0], arglist), "silver"))
	
	_print_lines(help)
	
def help():
	'''This help function.'''
	
	help = []
	help.append(("\nThe following help functions are available. To use them, import help and call the functions by name e.g.", "white"))
	help.append(("  ddglib.help.help()", "lightblue"))
	
	s_module = sys.modules[__name__]
	hfns = []
	for m in inspect.getmembers(s_module):
		o = m[1]
		if inspect.isroutine(o) and o.__module__ == "ddglib.help":
			if m[0][0] != "_":
				hfns.append((m[0], m[1].__doc__))
	
	if hfns:
		help.append(("  ddglib.help.%s()" % sorted(hfns)[0][0], "lightblue"))
	help.append(("", "silver"))
	
	for fn in sorted(hfns):
		help.append(("\t%s" % fn[0], "green"))
		help.append(("\t  %s" % fn[1], "silver"))
	
	_print_lines(help)
	

	

