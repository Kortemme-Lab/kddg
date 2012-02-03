#from string import join
from common.ddgproject import FieldNames, StdCursor, ddGDatabase
from common import colortext
from ddglib.filter import *

dbfields = FieldNames()

def show_tables():
	ddGdb = ddGDatabase()
	colortext.printf("\nThe tables of the database are as follows (orange = primary key, blue = foreign key):\n", color="white")
	tablenames = sorted([r[0] for r in (ddGdb.execute("SHOW TABLES", cursorClass = StdCursor))])
	for tbl in tablenames:
		colortext.message(tbl)
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
				colortext.printf(str, color = "orange")
			elif fkey == "MUL":
				colortext.printf(str, color = "lightblue")
			else:
				colortext.printf(str, color = "yellow")
