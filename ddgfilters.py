#!/usr/bin/python2.4
# encoding: utf-8
"""
ddgfilters.py
Defines filters and result sets for the ddG database.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""
import pickle
from string import join
import ddgdbapi
from ddglib.filter import *

from ddgdbapi import ddGDatabase
dbfields = ddGDatabase().FieldNames
StdCursor = ddgdbapi.StdCursor

class StructureResultSet(ResultSet):
	dbname = dbfields.Structure
	primary_key = dbfields.Structure.PDB_ID
	stored_procedures = []

	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = [], retrieveAllByDefault = True):
		super(StructureResultSet, self).__init__(db, SQL, parameters, AdditionalIDs, retrieveAllByDefault)
	
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
		self.UniProtACs = []
		self.UniProtIDs = []
		
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
	
	@staticmethod
	def WithUniProtIDs(UniProtIDs, UniProtACs):
		sf = StructureFilter()
		sf.setUniProtIDs(UniProtIDs, UniProtACs)
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

	def setUniProtIDs(self, UniProtACs, UniProtIDs = []):
		self.UniProtACs = dict.fromkeys(UniProtACs, True)
		self.UniProtIDs = dict.fromkeys(UniProtIDs, True)
		
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
		if self.UniProtACs or self.UniProtIDs:
			self.post_SQL_filters.append(self._filterByUniProt) 
	
	def _filterByUniProt(self, db, result_set):
		oldIDs = dict.fromkeys(result_set.IDs, True)
		allowedIDs = self.UniProtIDs
		allowedACs = self.UniProtACs
		
		UniProtKB_AC = dbfields.UniProtKB_AC
		UniProtKB_ID = dbfields.UniProtKB_ID
		
		new_IDs = {}
		mapping = db.callproc("GetPDBUniProtIDMapping")
		for m in mapping:
			pdbID = m["PDB_ID"]
			if (allowedACs.get(m[UniProtKB_AC]) or allowedIDs.get(m[UniProtKB_ID])) and oldIDs.get(pdbID):
				new_IDs[pdbID] = True
		
		return self._ResultSet(db, AdditionalIDs = new_IDs.keys(), retrieveAllByDefault = False)
			
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
		return self._ResultSet(db, AdditionalIDs = new_IDs, retrieveAllByDefault = False)
		
		
class ExperimentResultSet(ResultSet):
	dbname = dbfields.Experiment
	primary_key = dbfields.Experiment.ID
	stored_procedures = ["GetScores"]

	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = [], retrieveAllByDefault = True):
		super(ExperimentResultSet, self).__init__(db, SQL, parameters, AdditionalIDs, retrieveAllByDefault)
		
		self.structure_map = {} 
		for id in self.IDs:
			results = db.execute("SELECT Structure.PDB_ID FROM Experiment INNER JOIN Structure on Experiment.Structure=Structure.PDB_ID WHERE Experiment.ID=%s", parameters=(id,), cursorClass=StdCursor)
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
		elif filter.isOfClass(ExperimentFilter):
			experiments = filter.apply(self.db)
			return pks.intersection(experiments.IDs)		

	def _filterBySet(self, resSet):
		if resSet.isOfClass(StructureResultSet):
			foundIDs = []
			for s in resSet.IDs:
				if self.structure_map.get(s):
					foundIDs.extend(self.structure_map[s])
			return foundIDs
	
	def getStructures(self):
		idstr = join(map(str, list(self.getFilteredIDs())), ",")
		results = self.db.execute("SELECT Structure FROM Experiment WHERE ID IN (%s)" % idstr)
		structureIDs = self.db.execute("SELECT DISTINCT Structure FROM Experiment WHERE ID IN (%s)" % idstr, cursorClass=StdCursor)
		structureIDs = [s[0] for s in structureIDs]
		sr = StructureResultSet.fromIDs(self.db, structureIDs)
		return sr, results 


class ExperimentFilter(Filter):
	'''Note: minresolution and maxresolution are soft bounds i.e. if they are set to 1 and 2.5 respectively then we find structures with resolution r where 1 <= r <= 2.5.
			 If you want to specify 1 <= r < 2.5 instead, passing 2.49999 as the max resolution should usually work.'''
	
	_ResultSet = ExperimentResultSet
	
	ProTherm = "ProTherm-2008-09-08-23581"
	Potapov = "Potapov-2009"
	SenLiu = "SenLiu-ComplexExperimentalDataset"
	DummySource = "DummySource"
	LizKellogg = "LizKellogg:10.1002/prot.22921"
	large = "large"
	small = "small"
	
	
	_AllowedSources = [ProTherm, Potapov, SenLiu]
	
	# Constructors
	
	def __init__(self):
		super(ExperimentFilter, self).__init__()
		self.wildtype_size = None
		self.mutant_size = None
		self.wildtype = None
		self.mutant = None
		self.Source = None
		self.minmutations = None
		self.maxmutations = None
		self.minchains = None
		self.maxchains = None
		self.minstddev = None
		self.maxstddev = None
		
	@staticmethod
	def OnSource(Source):
		sf = ExperimentFilter()
		if not Source in ExperimentFilter._AllowedSources:
			raise Exception("Source %s is not recognized." % Source) 
		sf.setSource(Source)
		return sf
	
	@staticmethod
	def MutationsBetweenAminoAcidSizes(wildtype_size = None, mutant_size = None):
		sf = ExperimentFilter()
		sf.setAminoAcidSizes(wildtype_size, mutant_size)
		return sf
	
	@staticmethod
	def MutationsBetweenAminoAcids(wildtype = None, mutant = None):
		sf = ExperimentFilter()
		sf.setAminoAcids(wildtype, mutant)
		return sf
	
	@staticmethod
	def NumberOfMutations(minmutations = None, maxmutations = None):
		sf = ExperimentFilter()
		sf.setNumberOfMutations(minmutations, maxmutations)
		return sf
	
	@staticmethod
	def NumberOfChains(minchains = None, maxchains = None):
		sf = ExperimentFilter()
		sf.setNumberOfChains(minchains, maxchains)
		return sf
	
	@staticmethod
	def StandardDeviation(minstddev = None, maxstddev = None):
		sf = ExperimentFilter()
		sf.setStandardDeviation(minstddev, maxstddev)
		return sf
	
	# Property setters
	
	def setSource(self, Source):
		self.Source = Source
		
	def setAminoAcidSizes(self, wildtype_size = None, mutant_size = None):
		self.wildtype_size = wildtype_size
		self.mutant_size = mutant_size
	
	def setAminoAcids(self, wildtype = None, mutant = None):
		valid_wildtype = None
		valid_mutant = None
		for aa in ddgdbapi.aas:
			if wildtype == aa[0] or wildtype == aa[1]:
				valid_wildtype = aa[0] 
			if mutant == aa[0] or mutant == aa[1]:
				valid_mutant = aa[0]
		
		if wildtype and (not valid_wildtype):
			raise Exception("The wildtype amino acid '%s' is invalid." % wildtype)
		if mutant and (not valid_mutant):
			raise Exception("The mutant amino acid '%s' is invalid." % mutant)
		
		self.wildtype = valid_wildtype
		self.mutant = valid_mutant
		return

	def setNumberOfMutations(self, minmutations = None, maxmutations = None):
		self.minmutations = minmutations 
		self.maxmutations = maxmutations
	
	def setNumberOfChains(self, minchains = None, maxchains = None):
		self.minchains = minchains 
		self.maxchains = maxchains
	
	def setStandardDeviation(self, minstddev = None, maxstddev = None):
		self.minstddev = minstddev 
		self.maxstddev = maxstddev
	
	# _apply and filter functions
		
	def _apply(self):
		self.post_SQL_filters = []
		if self.Source != None:
			self.conditions.append("Source = %s")
			self.parameters.append(self.Source)		
		if self.wildtype_size or self.mutant_size:
			self.post_SQL_filters.append(self._filterByMutationSize) 
		if self.wildtype or self.mutant:
			self.post_SQL_filters.append(self._filterByMutation) 
		if self.minmutations or self.maxmutations:
			self.post_SQL_filters.append(self._filterByNumberOfMutations) 
		if self.minchains or self.maxchains:
			self.post_SQL_filters.append(self._filterByNumberOfChains) 
		if self.minstddev or self.maxstddev:
			self.post_SQL_filters.append(self._filterByStandardDeviation) 

	def _filterByMutation(self, db, result_set):
		if self.wildtype or self.mutant:
			conditions = []
			if self.wildtype:
				conditions.append('wildtype.Code="%s"' % self.wildtype)
			if self.mutant:
				conditions.append('mutant.Code="%s"' % self.mutant)
			conditions = join(conditions, " AND ")
			results = db.execute('''SELECT ExperimentID, Structure, Source, Chain, ResidueID, WildtypeAA, MutantAA, wildtype.Size, mutant.Size FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE %s;''' % conditions)
			return self._ResultSet(db, AdditionalIDs = [r["ExperimentID"] for r in results], retrieveAllByDefault = False)			
		else:
			return result_set

	def _filterByMutationSize(self, db, result_set):
		if self.wildtype_size or self.mutant_size:
			conditions = []
			if self.wildtype_size:
				conditions.append('wildtype.Size="%s"' % self.wildtype_size)
			if self.mutant_size:
				conditions.append('mutant.Size="%s"' % self.mutant_size)
			conditions = join(conditions, " AND ")
			results = db.execute('''SELECT ExperimentID, Structure, Source, Chain, ResidueID, WildtypeAA, MutantAA, wildtype.Size, mutant.Size FROM Experiment INNER JOIN ExperimentMutation ON ID=ExperimentID INNER JOIN AminoAcid as wildtype ON WildtypeAA = wildtype.Code INNER JOIN AminoAcid as mutant ON MutantAA = mutant.Code WHERE %s;''' % conditions)
			return self._ResultSet(db, AdditionalIDs = [r["ExperimentID"] for r in results], retrieveAllByDefault = False)			
		else:
			return result_set
		#todo: Check that large, none = large, small + large, large

	def _filterByNumberOfMutations(self, db, result_set):
		if self.minmutations or self.maxmutations:
			results = db.execute('''SELECT ExperimentID, COUNT(*) FROM ExperimentMutation GROUP BY ExperimentID;''', cursorClass=StdCursor)
			if self.minmutations and self.maxmutations:
				IDs = [r[0] for r in results if r[1] >= self.minmutations and r[1] <= self.maxmutations] 
			elif self.minmutations and self.maxmutations:
				IDs = [r[0] for r in results if r[1] >= self.minmutations] 
			else:
				IDs = [r[0] for r in results if r[1] <= self.maxmutations] 
			return self._ResultSet(db, AdditionalIDs = IDs, retrieveAllByDefault = False)
		else:
			return result_set

	def _filterByNumberOfChains(self, db, result_set):
		if self.minchains or self.maxchains:
			results = db.execute('''SELECT ExperimentID, COUNT(*) FROM ExperimentChain GROUP BY ExperimentID;''', cursorClass=StdCursor)
			if self.minchains and self.maxchains:
				IDs = [r[0] for r in results if r[1] >= self.minchains and r[1] <= self.maxchains] 
			elif self.minchains and self.maxchains:
				IDs = [r[0] for r in results if r[1] >= self.minchains] 
			else:
				IDs = [r[0] for r in results if r[1] <= self.maxchains] 
			return self._ResultSet(db, AdditionalIDs = IDs, retrieveAllByDefault = False)
		else:
			return result_set

	def _filterByStandardDeviation(self, db, result_set):
		if self.minstddev and self.maxstddev:
			IDs = [ID for ID in result_set.IDs if self.minstddev <= db.getStandardDeviation(ID) <= self.maxstddev]
			return self._ResultSet(db, AdditionalIDs = IDs, retrieveAllByDefault = False)
		elif self.minstddev:
			IDs = [ID for ID in result_set.IDs if self.minstddev <= db.getStandardDeviation(ID)]
			return self._ResultSet(db, AdditionalIDs = IDs, retrieveAllByDefault = False)
		elif self.maxstddev:
			IDs = [ID for ID in result_set.IDs if db.getStandardDeviation(ID) <= self.maxstddev]
			return self._ResultSet(db, AdditionalIDs = IDs, retrieveAllByDefault = False)
		else:
			return result_set
			
class PredictionResultSet(ResultSet):
	dbname = dbfields.Prediction
	primary_key = dbfields.Prediction.ID
	stored_procedures = []
	
	def __init__(self, db, SQL = "", parameters = None, AdditionalIDs = [], retrieveAllByDefault = True):
		super(PredictionResultSet, self).__init__(db, SQL, parameters, AdditionalIDs, retrieveAllByDefault)
		
		self.structure_map = {}
		for id in self.IDs:
			results = db.execute("SELECT Structure.PDB_ID FROM Prediction INNER JOIN Experiment ON ExperimentID=Experiment.ID INNER JOIN Structure on Experiment.Structure=Structure.PDB_ID WHERE Prediction.ID=%s", parameters=(id,), cursorClass=StdCursor)
			pdbID = results[0][0]
			self.structure_map[pdbID] = self.structure_map.get(pdbID) or []
			self.structure_map[pdbID].append(id)
		
		self.experiment_map = {}
		for id in self.IDs:
			results = db.execute("SELECT Experiment.ID FROM Prediction INNER JOIN Experiment ON ExperimentID=Experiment.ID WHERE Prediction.ID=%s", parameters=(id,), cursorClass=StdCursor)
			experimentID = results[0][0]
			self.experiment_map[experimentID] = self.experiment_map.get(experimentID) or []
			self.experiment_map[experimentID].append(id)
		
	
	def applyFilter(self, pks, filter, tag):
		if filter.isOfClass(StructureFilter):
			structures = filter.apply(self.db)
			foundIDs = []
			for s in structures.IDs:
				if self.structure_map.get(s):
					foundIDs.extend(self.structure_map[s])
			return pks.intersection(foundIDs)
		elif filter.isOfClass(ExperimentFilter):
			experiments = filter.apply(self.db)
			foundIDs = []
			for s in experiments.IDs:
				if self.experiment_map.get(s):
					foundIDs.extend(self.experiment_map[s])
			return pks.intersection(foundIDs)
		elif filter.isOfClass(PredictionFilter):
			predictions = filter.apply(self.db)
			return pks.intersection(predictions.IDs)		

	
	def _filterBySet(self, resSet):
		if resSet.isOfClass(StructureResultSet):
			foundIDs = []
			for s in resSet.IDs:
				if self.structure_map.get(s):
					foundIDs.extend(self.structure_map[s])
			return foundIDs
		elif resSet.isOfClass(ExperimentResultSet):
			foundIDs = []
			for s in resSet.IDs:
				if self.structure_map.get(s):
					foundIDs.extend(self.structure_map[s])
			return foundIDs
	
	def getStructures(self):
		idstr = join(map(str, list(self.getFilteredIDs())), ",")
		results = self.db.execute("SELECT Prediction.ID, Experiment.Structure FROM Prediction INNER JOIN Experiment ON ExperimentID = Experiment.ID WHERE Prediction.ID IN (%s)" % idstr)
		structureIDs = self.db.execute("SELECT DISTINCT Experiment.Structure FROM Prediction INNER JOIN Experiment ON ExperimentID = Experiment.ID WHERE Prediction.ID IN (%s)" % idstr, cursorClass=StdCursor)
		structureIDs = [s[0] for s in structureIDs]
		sr = StructureResultSet.fromIDs(self.db, structureIDs)
		return sr, results 

	def getExperiments(self):
		idstr = join(map(str, list(self.getFilteredIDs())), ",")
		results = self.db.execute("SELECT Prediction.ID, ExperimentID FROM Prediction WHERE Prediction.ID IN (%s)" % idstr)
		experimentIDs = self.db.execute("SELECT DISTINCT ExperimentID FROM Prediction WHERE Prediction.ID IN (%s)" % idstr, cursorClass=StdCursor)
		experimentIDs = [s[0] for s in experimentIDs]
		er = ExperimentResultSet.fromIDs(self.db, experimentIDs)
		return er, results 
			

class PredictionFilter(Filter):
	'''Note: minresolution and maxresolution are soft bounds i.e. if they are set to 1 and 2.5 respectively then we find structures with resolution r where 1 <= r <= 2.5.
			 If you want to specify 1 <= r < 2.5 instead, passing 2.49999 as the max resolution should usually work.'''
	
	_ResultSet = PredictionResultSet
	
	#XRay = "X-RAY DIFFRACTION"
	#NMR = "SOLUTION NMR"
	#_AllowedTechniques = [XRay, NMR]
	
	# Constructors
	
	def __init__(self):
		super(Prediction, self).__init__()

		
StructureResultSet.allowed_filters = [StructureFilter]	
ExperimentResultSet.allowed_filters = [StructureFilter, ExperimentFilter]	
PredictionResultSet.allowed_filters = [StructureFilter, ExperimentFilter, PredictionFilter]
ExperimentResultSet.allowed_restrict_sets = [StructureResultSet] 	
PredictionResultSet.allowed_restrict_sets = [StructureResultSet, ExperimentResultSet] 	
