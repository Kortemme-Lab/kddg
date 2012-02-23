#!/usr/bin/python2.4
# encoding: utf-8
"""
filter.py
Defines generic filter and result set classes for use on databases.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import re
from string import join
from base import kobject

class Filter(kobject):
	
	_ResultSet = None
	
	def __init__(self):
		self.post_SQL_filters = []
		#raise Exception("Calling an abstract class.")
		
	def apply(self, db):
		self.conditions = []
		self.parameters = []
		self.paramsoffset = 0
		
		self._apply()
		conditions = self.conditions
		parameters = self.parameters
		assert(len(conditions) == len(parameters) - self.paramsoffset)
		if conditions:
			result_set = self._ResultSet(db, SQL = "WHERE %s" % join(conditions, " AND "), parameters = tuple(parameters))
		else:
			result_set = self._ResultSet(db)
		
		for postfilter in self.post_SQL_filters:
			# post_SQL_filters functions take a database connection and self._ResultSet object and return a self._ResultSet object 
			result_set &= postfilter(db, result_set)
			
		return result_set

	def __or__(self, disjunct):
		if disjunct.isOfClass(UnionFilter):
			return UnionFilter([self] + disjunct.filters)
		elif Filter in disjunct.getBaseClasses():
			return UnionFilter([self, disjunct])
		else:
			raise Exception("Trying to combine a filter with an object of class '%s'." % disjunct.getClassName())

class UnionFilter(kobject):
	
	def __init__(self, filters = []):
		self.filters = filters
	
	def add(self, filter):
		if Filter in filter.getBaseClasses():
			self.filters.append(filter)
		else:
			raise Exception("Trying to combine a filter with an object of class '%s'." % disjunct.getClassName())
	
	def getFilters(self):
		return self.filters
	
	def __or__(self, disjunct):
		if disjunct.isOfClass(UnionFilter):
			return UnionFilter(self.filters + disjunct.filters)
		elif Filter in disjunct.getBaseClasses():
			return UnionFilter(self.filters + [disjunct])
		else:
			raise Exception("Trying to combine a filter with an object of class '%s'." % disjunct.getClassName())
		
	
class ResultSet(kobject):
	dbname = None
	primary_key = "ID"
	stored_procedures = [] # e.g. ["GetScores"]
	allowed_filters = []
	allowed_restrict_sets = []
	
	stored_procedure_regex = re.compile("CALL (\w+)")
	
	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = [], retrieveAllByDefault = True):
		'''e.g. ResultSet("CALL GetScores", parameters = 16734)
				ResultSet("WHERE ...")'''
		
		self._log = []
		self.db = db
		if parameters:
			if type(parameters) != tuple:
				parameters = (parameters,)
			if not SQL:
				raise Exception("Parameters %s were specified for an empty SQL query." % parameters)
		
		results = []
		if SQL or (retrieveAllByDefault and not AdditionalIDs):
			# We run an SQL query if the SQL parameter has been specified or if AdditionalIDs is empty.
			# If AdditionalIDs is not empty and the SQL is blank then we do NOT run the SQL query and
			# instead use AdditionalIDs as the record keys.  
			sp = ResultSet.stored_procedure_regex.match(SQL)
			if sp:
				if not sp.groups(1) in self.__class__.stored_procedures:
					raise Exception("Cannot call stored procedure %s to create a %s" % (sp.groups(1)[0], self.__class__.__name__))
				else:
					results = db.callproc(SQL, parameters)
					self.log("ResultSet object initialized with stored procedure call %s%s." % (SQL, parameters))
			else:
				if not self.__class__.dbname:
					raise Exception("There is not database table associated with the class %s." % ( self.__class__.__name__))
				SQL = "SELECT %s FROM %s %s" % (self.__class__.primary_key, self.__class__.dbname, SQL)
				results = self.db.execute(SQL, parameters)
				self.log("ResultSet object initialized with SQL query '%s' %% %s." % (SQL, parameters or ""))
		
		AdditionalIDs = set(AdditionalIDs)
		if AdditionalIDs:
			AdditionalIDs = set(AdditionalIDs)
			allIDs = self.db.execute("SELECT %s FROM %s" % (self.__class__.primary_key, self.__class__.dbname), parameters)
			pIDs = set([r[self.__class__.primary_key] for r in allIDs]).intersection(AdditionalIDs)
			if not (len(pIDs) == len(AdditionalIDs)):
				raise Exception("Records associated with the following IDs could not be found: %s." % join(AdditionalIDs.difference(pIDs), ","))
		
		if SQL:
			if not results:
				print("No results were returned from '%s' %% %s." % (SQL, parameters or ""))
			else:
				if not self.__class__.primary_key in results[0]:
					raise Exception("The resulting set from '%s'(%s) must including the primary key field %s of %s." % (SQL, parameters, self.__class__.primary_key, self.__class__.dbname))
				
		self.IDs = set([r[self.__class__.primary_key] for r in results]).union(AdditionalIDs)
		self.initialresults = None
		self.filterchain = []
		self.log("Initial record count: %d" % len(self.IDs))
	
	def log(self, msg):
		self._log.append(msg)
				
	def __repr__(self):
		return join(self._log, "\n")
	
	def getInitialResults(self):
		if not self.initialresults:
			SQL = "SELECT * FROM %s" % self.__class__.dbname
			results = self.db.execute(SQL)
			self.initialresults = [r for r in results if r[self.__class__.primary_key] in self.IDs]
		return self.initialresults

	def addFilter(self, filter, tag = None):
		if not (Filter in filter.getBaseClassesForObject() or filter.isOfClass(UnionFilter)):
			raise Exception("Trying to add a filter which does not derive from the Filter class.")
		
		allowed = False 
		for f in self.__class__.allowed_filters:
			if filter.isOfClass(f):
				allowed = True
				break
		if not allowed:
			if filter.isOfClass(UnionFilter):
				for f in filter.getFilters():
					if not (f.__class__ in self.__class__.allowed_filters):
						raise Exception("Trying to add a %s filter which is not supported by this class (%s)." % (f.getClassName(), self.getClassName()))
			else:
				raise Exception("Trying to add a %s filter which is not supported by this class (%s)." % (filter.getClassName(), self.getClassName()))
		
		self.filterchain.append((filter, tag))
	
	def applyUnionFilter(self, pks, ufilter, tag):
		pkset = set()
		for filter in ufilter.getFilters():
			self.log("  Applying %s:" % (tag or filter.getClassName()))
			pkset = pkset.union(self.applyFilter(pks, filter, tag))
		return pkset
	
	def applyFilter(self, pks, filter, tag):
		self.log("Applying %s:" % (tag or filter.getClassName()))
		raise Exception("This function needs to be implemented for this class (%s)." % self.getClassName())
		
	def getFilteredIDs(self):
		'''Applies the filters, returns a new dict'''
		pks = self.IDs
		#self.log("Primary keys before filtering:")
		#self.log(str(pks))
		for taggedFilter in self.filterchain:
			filter = taggedFilter[0] 
			tag = taggedFilter[1]
			if filter.isOfClass(UnionFilter):
				self.log("Applying %s:"% (tag or "Union Filter"))
				pks = self.applyUnionFilter(pks, filter, tag)
				#self.log("Primary keys after filtering:")
				#self.log(str(pks))
			elif Filter in filter.getBaseClasses():
				self.log("Applying %s:" % (tag or filter.getClassName()))
				pks = self.applyFilter(pks, filter, tag)
				#self.log("Primary keys after filtering:")
				#self.log(str(pks))
			else:
				raise Exception("BLARG!")
		self.log("Filtered record count: %d" % len(pks))
		return pks
	
	def getFilteredResults(self):
		'''Applies the filters, returns a new dict'''
		pks = self.getFilteredIDs()
		SQL = "SELECT * FROM %s" % self.__class__.dbname
		results = self.db.execute(SQL)
		return [r for r in results if r[self.__class__.primary_key] in pks]
	
	def filterBySet(self, resSet):
		if resSet.__class__== self.__class__:
			return self.__class__(self.db, AdditionalIDs = self.IDs.intersection(resSet.IDs),  retrieveAllByDefault = False)
		elif resSet.__class__ in self.__class__.allowed_restrict_sets:
			pks = self._filterBySet(resSet)
			return self.__class__(self.db, AdditionalIDs = pks, retrieveAllByDefault = False)
		else:
			raise Exception("Trying to restrict %s by a %s set which is impossible." % (self.__class__, resSet.__class__))
		
	def __or__(self, resSet):
		if self.hasTheSameClassAs(resSet):
			pks = self.getFilteredIDs().union(resSet.getFilteredIDs())
			disjunction = self.__class__(self.db, AdditionalIDs = pks, retrieveAllByDefault = False)
			disjunction._log = self._log + resSet._log + ["UNION OF RESULT SETS", "Filtered record count: %d" % len(disjunction.IDs)] 
			return disjunction
		else:
			raise Exception("Trying to combine two ResultSets of different types (%s and %s)." % (self.getClassName(), resSet.getClassName()))

	def __and__(self, resSet):
		if self.hasTheSameClassAs(resSet):
			pks = self.getFilteredIDs().intersection(resSet.getFilteredIDs())
			conjunction = self.__class__(self.db, AdditionalIDs = pks, retrieveAllByDefault = False)
			conjunction._log = self._log + resSet._log + ["INTERSECTION OF RESULT SETS", "Filtered record count: %d" % len(conjunction.IDs)] 
			return conjunction
		else:
			raise Exception("Trying to combine two ResultSets of different types (%s and %s)." % (self.getClassName(), resSet.getClassName()))

	def __sub__(self, resSet):
		if self.hasTheSameClassAs(resSet):
			pks = self.getFilteredIDs().difference(resSet.getFilteredIDs())
			conjunction = self.__class__(self.db, AdditionalIDs = pks, retrieveAllByDefault = False)
			conjunction._log = self._log + resSet._log + ["DIFFERENCE OF RESULT SETS", "Filtered record count: %d" % len(conjunction.IDs)] 
			return conjunction
		else:
			raise Exception("Trying to combine two ResultSets of different types (%s and %s)." % (self.getClassName(), resSet.getClassName()))

	def __div__(self, resSet):
		return self.__sub__(resSet)

	def __xor__(self, resSet):
		if self.hasTheSameClassAs(resSet):
			my_pks = self.getFilteredIDs()
			other_pks = resSet.getFilteredIDs()
			allpks = my_pks.union(other_pks)
			commonpks = my_pks.intersection(other_pks)
			pks = allpks.difference(commonpks)
			xconjunction = self.__class__(self.db, AdditionalIDs = pks, retrieveAllByDefault = False)
			xconjunction._log = self._log + resSet._log + ["XOR OF RESULT SETS", "Filtered record count: %d" % len(xconjunction.IDs)] 
			return xconjunction
		else:
			raise Exception("Trying to combine two ResultSets of different types (%s and %s)." % (self.getClassName(), resSet.getClassName()))


