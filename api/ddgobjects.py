import inspect
from klab import colortext
from string import join

class AbstractClass(object):
	'''Helper class to try to emulate abstract classes in Python. For this to work, the abstract class A must call
	     super(A, self).__init__(A)
	   in its __init__ method and each concrete subclass B must call
	     super(B, self).__init__() 
	   in its __init__ method.'''
	
	def __init__(self, abstractClassReference, abstract_metamethods = []):
		'''Checks whether all non-private members (using the convention that private members start with a '_') are implemented.
		   abstract_metamethods is an optional argument to specify non-private/metamethods which are to be considered abstract.'''
		errors = []
		for m in inspect.getmembers(abstractClassReference, inspect.ismethod):
			if m[0][0] != "_" or m[0] in abstract_metamethods:
				if getattr(self, m[0]).im_func == getattr(abstractClassReference, m[0]).im_func:
					errors.append("The abstract method %s of class %s has not been implemented in class %s." % (m[0], abstractClassReference.__name__, self.__class__.__name__))
		if errors:
			raise NotImplementedError("\n%s" % join(errors, "\n"))	
		

class DatasetParser(AbstractClass):
	'''I am just defining this class as an example of the proposed interface. No real OO design is happening yet.'''
	def __init__(self): super(DatasetParser, self).__init__(DatasetParser)
	def compareToDatabase(self, ProThermID, ExperimentScoreID, record): pass
	def parse(self, args = None): pass
	def addDataSet(self): pass

class AbstractDBObject(AbstractClass):
	dict = {}
	
	def __init__(self, ddGdb):
		self.ddGdb = ddGdb
		self.quiet = False
		self.databaseID = None
		super(AbstractDBObject, self).__init__(AbstractDBObject, ['__repr__'])
		
	def __getitem__(self, key):
		return self.dict[key]
	
	def find(self):
		'''If the records already exists in the database then two variables are returned - the record ID in the database and some printable object.
		   Otherwise, false is returned.'''
		pass

	def test(self):
		DBID, details = self.find()
		if DBID != None:
			#if not self.quiet:
			#	colortext.warning("\nWarning: Experiment already exists in the database with Experiment ID=%s." % DBID)
			#	colortext.warning("*** This record ***%s\n" % self)
			#	colortext.warning("*** Database record ***\n%s\n" % details)
			return DBID

	def commit(self):
		DBID, details = self.find()
		if DBID != None:
			#if not self.quiet:
			#	colortext.error("\nExperiment already exists in the database with Experiment ID=%s." % DBID)
			#	colortext.warning("*** This record ***%s\n" % self)
			#	colortext.warning("*** Database record ***\n%s\n" % details)
			return DBID
		
	def remove(self):
		'''Removes the record and associated records from the database using self.databaseID.'''
		pass


class DBObject(AbstractDBObject):
	def __init__(self, ddGdb):
		 super(DBObject, self).__init__(ddGdb) 
		 
	def getDatabaseID(self):
		if not self.databaseID:
			raise Exception("Cannot get the database ID of an uncommitted record.")
		else:
			return self.databaseID

