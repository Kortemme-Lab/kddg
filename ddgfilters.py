import pickle
from string import join
import common.ddgproject
from ddglib.filter import *

dbfields = common.ddgproject.FieldNames()
StdCursor = common.ddgproject.StdCursor

class StructureResultSet(ResultSet):
	dbname = dbfields.Structure
	primary_key = dbfields.PDB_ID
	stored_procedures = []

	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = []):
		super(StructureResultSet, self).__init__(db, SQL, parameters, AdditionalIDs)
	
	def applyFilter(self, pks, filter, tag):
		if filter.isOfClass(StructureFilter):
			structures = filter.apply(self.db)
			return pks.intersection(structures.IDs)

class StructureFilter(Filter):
	'''Note: minresolution and maxresolution are soft bounds i.e. if they are set to 1 and 2.5 respectively then we find structures with resolution r where 1 <= r <= 2.5.
			 If you want to specify 1 <= r < 2.5 instead, passing 2.49999 as the max resolution should usually work.'''
	
	_ResultSet = StructureResultSet
	
	XRay = "X-RAY DIFFRACTION"
	NMR = "SOLUTION NMR"
	_AllowedTechniques = [XRay, NMR]
	
	# Constructors
	
	def __init__(self):
		super(StructureFilter, self).__init__()
		self.nullresolution = None
		self.minresolution = None
		self.maxresolution = None
		self.bfactors_min = None
		self.bfactors_max = None
		self.techniques = None
		self.post_SQL_filters = []
		
	@staticmethod
	def TotalBFactors(min, max):
		sf = StructureFilter()
		sf.setTotalBFactors(min, max)
		return sf

	@staticmethod
	def WithNullResolution(allowNull = True):
		sf = StructureFilter()
		sf.nullresolution = allowNull
		return sf

	@staticmethod
	def Resolution(min, max):
		sf = StructureFilter()
		sf.setResolution(min, max)
		return sf
	
	@staticmethod
	def Techniques(techniques):
		sf = StructureFilter()
		sf.setTechniques(techniques)
		return sf
	
	# Property setters
	
	def setAllowNullResolution(self, allowNull):
		# This function is added to keep in line with the reset of the API
		sf.nullresolution = allowNull
	
	def setTechniques(self, techniques):
		if type(techniques) != list:
			techniques = [techniques]
		for t in techniques:
			if not t in self._AllowedTechniques:
				raise Exception("Technique %s is not recognized." % t) 
		self.techniques = techniques
	
	def setResolution(self, min, max):
		self.minresolution = min
		self.maxresolution = max
	
	def setTotalBFactors(self, min, max):
		self.bfactors_min = min
		self.bfactors_max = max

	# _apply and filter functions
		
	def _apply(self):
		self.post_SQL_filters = []
		if self.nullresolution != None:
			if self.nullresolution:
				self.conditions.append("Resolution IS NULL")
			else:
				self.conditions.append("Resolution IS NOT NULL")
			self.paramsoffset -=  1
		if self.minresolution != None:
			self.conditions.append("Resolution >= %s")
			self.parameters.append(self.minresolution)
		if self.maxresolution != None:
			self.conditions.append("Resolution <= %s")
			self.parameters.append(self.maxresolution)
		if self.techniques:
			tstrs = []
			for t in self.techniques:
				tstrs.append('Techniques LIKE %s')
				self.parameters.append("%%%s%%" % t)
			self.conditions.append(join(tstrs, " OR "))
			self.paramsoffset += (len(self.techniques) - 1)
		
		if self.bfactors_min or self.bfactors_max:
			self.post_SQL_filters.append(self._checkTotalBFactorRange) 
			
	def _checkTotalBFactorRange(self, db, result_set):
		m_min = self.bfactors_min
		m_max = self.bfactors_max
		pkname = self._ResultSet.primary_key
		
		new_IDs = []
		for record in result_set.getInitialResults():
			bfactors =  pickle.loads(record["BFactors"])
			average = bfactors["Total"][0]
			passed = True
			if m_min and average < m_min: 
				passed = False
			if m_max and average > m_max: 
				passed = False
			if passed:
				new_IDs.append(record[pkname])
		return self._ResultSet(db, AdditionalIDs = new_IDs)
		
class ExperimentResultSet(ResultSet):
	dbname = dbfields.Experiment
	primary_key = dbfields.ID
	stored_procedures = ["GetScores"]

class PredictionResultSet(ResultSet):
	dbname = dbfields.Prediction
	primary_key = dbfields.ID
	stored_procedures = []
	
	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = []):
		super(PredictionResultSet, self).__init__(db, SQL, parameters, AdditionalIDs)
		
		self.structure_map = {} 
		for id in self.IDs:
			results = db.execute("SELECT Structure.PDB_ID FROM Prediction INNER JOIN Experiment ON ExperimentID=Experiment.ID INNER JOIN Structure on Experiment.Structure=Structure.PDB_ID WHERE Prediction.ID=%s", parameters=(id,), cursorClass=StdCursor)
			pdbID = results[0][0]
			self.structure_map[pdbID] = self.structure_map.get(pdbID) or []
			self.structure_map[pdbID].append(id)
		
	
	def applyFilter(self, pks, filter, tag):
		if filter.isOfClass(StructureFilter):
			structures = filter.apply(self.db)
			foundIDs = []
			for s in structures.IDs:
				if self.structure_map.get(s):
					foundIDs.extend(self.structure_map[s])
			return pks.intersection(foundIDs)

	
StructureResultSet.allowed_filters = [StructureFilter]	
ExperimentResultSet.allowed_filters = [StructureFilter]	
PredictionResultSet.allowed_filters = [StructureFilter]	