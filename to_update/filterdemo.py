import sys
import time
import profile

#sys.path.insert(0, "..")
#sys.path.insert(0, "ddglib")
from klab import colortext
from klab.deprecated import rosettadb
from klab.debug.profile import ProfileTimer
from ddglib import db_api, ddgdbapi
from ddglib import help as ddg_help
from ddglib.ddgfilters import *

from klab import pdb
import klab.deprecated.rosettahelper

def simpleRunExample(self):
	# Step 1: Open a database connection
	ddGdb = ddgdbapi.ddGDatabase()

	# Step 2: Select database records
	sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])
	
	# Step 3: Add filters
	sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))
	
	# Step 4: Retrieve full database records. 
	# results will be a list each of whose elements is a dict representing a database record.
	results = sr.getFilteredResults()
	
	# Step 5: Optionally print summary
	print("\nSummary: %s\n" % sr)


def help():
	ddg_help.ShowDatabaseStructure()
	ddg_help.ShowResultSet()
	ddg_help.ShowFilter()

def dump_zip(jobnumber):
	ddG_connection = db_api.ddG()
	ddG_connection.dumpData("testzip-%d.zip" % jobnumber, jobnumber)

class JobRunner:
	# Class to contain old code used to kick off jobs

	#e.g.
	# JobRunner.addLinsJobs("lin-3K0NB", "Kellogg:10.1002/prot.22921:protocol16:32231")
	# JobRunner.runLizsSet("lizsettest1", "Kellogg:10.1002/prot.22921:protocol16:32231")

	@staticmethod
	def addLinsJobs(PredictionSet, ProtocolID):
		raise colortext.Exception("Do you really want to run this?")
		colortext.printf("\nAdding Lin's mutations to %s prediction set." % PredictionSet, "lightgreen")
		KeepHETATMLines = False
		FilterTester.openDB()

		# Filter by the DummySource set of experiments
		er1 = ExperimentResultSet(ddGdb)
		ef1 = ExperimentFilter()
		ef1.setSource(ExperimentFilter.DummySource)
		er1.addFilter(ef1)

		# Filter by the particular PDB
		sr = StructureResultSet(ddGdb, 'WHERE PDB_ID="3K0NB_lin"')
		er1 = ExperimentResultSet.fromIDs(ddGdb, er1.getFilteredIDs()).filterBySet(sr)
		FilterTester.test(er1)

		experimentIDs = sorted(list(er1.getFilteredIDs()))
		colortext.message("\nThe number of unique experiments is %d.\n" % len(experimentIDs))
		ddG_connection = db_api.ddG()
		count = 0
		for experimentID in experimentIDs:
			ddG_connection.addPrediction(experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
			count += 1
			if count >= 10:
				colortext.write(".")
				colortext.flush()
				count = 0
		print("")

	@staticmethod
	def runLizsSet(PredictionSet, ProtocolID):
		raise colortext.Exception("Do you really want to run this?")
		colortext.printf("\nAdding Liz's data set to %s prediction set." % PredictionSet, "lightgreen")
		KeepHETATMLines = False
		FilterTester.openDB()

		# Filter by the DummySource set of experiments
		er1 = ExperimentResultSet(ddGdb)
		ef1 = ExperimentFilter()
		ef1.setSource(ExperimentFilter.LizKellogg)
		er1.addFilter(ef1)
		FilterTester.test(er1)

		experimentIDs = sorted(list(er1.getFilteredIDs()))
		colortext.message("\nThe number of unique experiments is %d.\n" % len(experimentIDs))
		ddG_connection = db_api.ddG()
		count = 0
		for experimentID in experimentIDs:
			ddG_connection.addPrediction(experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
			count += 1
			if count >= 10:
				colortext.write(".")
				colortext.flush()
				count = 0
		print("")

	@staticmethod
	def addAllMutationsForAGivenPDB1():
		'''Used to create dummy Experiment records for Lin's DDG run. This should probably be an API function.'''
		ddG_connection = db_api.ddG()
		opdb = common.pdb.PDB("3K0NA_lin.pdb")
		count = 1
		for chainresidueid, wt in sorted(opdb.ProperResidueIDToAAMap().iteritems()):
			chain = chainresidueid[0]
			residueid = chainresidueid[1:].strip()
			allotherAAs = sorted([aa for aa in klab.deprecated.rosettahelper.ROSETTAWEB_SK_AAinv.keys() if aa != wt])
			for otherAA in allotherAAs:
				ms = db_api.MutationSet()
				ms.addMutation(chain, residueid, wt, otherAA)
				print("3K0NA_lin", ms, ms.getChains(), count, 0)
				ddG_connection.createDummyExperiment("3K0NA_lin", ms, ms.getChains(), count, 0, ExperimentSetName = "DummySource")
				count += 1

	@staticmethod
	def addAllMutationsForAGivenPDB2():
		'''Used to create dummy Experiment records for Lin's DDG run. This should probably be an API function.'''
		ddG_connection = db_api.ddG()
		opdb = common.pdb.PDB("3K0On_lin.pdb")
		count = 3098
		for chainresidueid, wt in sorted(opdb.ProperResidueIDToAAMap().iteritems()):
			chain = chainresidueid[0]
			residueid = chainresidueid[1:].strip()
			allotherAAs = sorted([aa for aa in klab.deprecated.rosettahelper.ROSETTAWEB_SK_AAinv.keys() if aa != wt])
			for otherAA in allotherAAs:
				ms = db_api.MutationSet()
				ms.addMutation(chain, residueid, wt, otherAA)
				print("3K0On_lin", ms, ms.getChains(), count, 0)
				ddG_connection.createDummyExperiment("3K0On_lin", ms, ms.getChains(), count, 0, ExperimentSetName = "DummySource")
				count += 1

	@staticmethod
	def addAllMutationsForAGivenPDB3():
		'''Used to create dummy Experiment records for Lin's DDG run. This should probably be an API function.'''
		FilterTester.openDB()
		ddG_connection = db_api.ddG()
		opdb = common.pdb.PDB("pdbs/3K0NB_lin.pdb")

		results = ddGdb.execute('''SELECT SourceID FROM ExperimentScore INNER JOIN Experiment ON ExperimentScore.ExperimentID = Experiment.ID WHERE Source="DummySource"''', cursorClass=ddgdbapi.StdCursor)
		assert(results)
		highestID = max([int(r[0]) for r in results])
		count = highestID + 1

		for chainresidueid, wt in sorted(opdb.ProperResidueIDToAAMap().iteritems()):
			chain = chainresidueid[0]
			residueid = chainresidueid[1:].strip()
			allotherAAs = sorted([aa for aa in klab.deprecated.rosettahelper.ROSETTAWEB_SK_AAinv.keys() if aa != wt])
			for otherAA in allotherAAs:
				ms = db_api.MutationSet()
				ms.addMutation(chain, residueid, wt, otherAA)
				print("3K0NB_lin", ms, ms.getChains(), count, 0, chain, wt, residueid, otherAA)
				ddG_connection.createDummyExperiment("3K0NB_lin", ms, ms.getChains(), count, 0, ExperimentSetName = "DummySource")
				count += 1

class Analyzer:
	# Class to contain old code used to kick off analysis

	@staticmethod
	def testAnalysis():
		ddG_connection = db_api.ddG()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='kellogg16-A' AND Status='done' LIMIT 2000")
		ddG_connection.analyze(pr)

	@staticmethod
	def testAnalysis2():
		ddG_connection = db_api.ddG()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='lizsettest1' AND Status='done' LIMIT 2000")
		ddG_connection.analyze(pr)

		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='lizsettest1' AND Status='done' LIMIT 2000")
		pr.addFilter(ExperimentFilter.MutationsBetweenAminoAcidSizes(ExperimentFilter.large, ExperimentFilter.small))
		FilterTester.test(pr)
		ddG_connection.analyze(pr)

		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet='lizsettest1' AND Status='done' LIMIT 2000")
		pr.addFilter(ExperimentFilter.MutationsBetweenAminoAcidSizes(ExperimentFilter.small, ExperimentFilter.large))
		FilterTester.test(pr)
		ddG_connection.analyze(pr)


class FilterTester:

	@staticmethod
	def profile(command_):
		t1 = time.time()
		profile.run(command_, sort = 'cumulative')
		print("** Total time taken in %s: %0.2f **" % (command_, time.time() - t1))

	@staticmethod
	def test(resultset, expected_size = None):
		print("Applying filters")
		results = resultset.getFilteredResults(just_get_primary_keys = True)
		print(len(results), expected_size)
		assert(len(results) == expected_size)
		print("After application")
		print("\nSummary: %s\n" % resultset)

	@staticmethod
	def openDB():
		if not globals().get("ddGdb"):
			globals()["ddGdb"] = ddgdbapi.ddGDatabase()
			total_number_of_experiments = ddGdb.execute('SELECT COUNT(ID) AS C FROM Experiment')[0]['C']
			globals()["total_number_of_experiments"] = total_number_of_experiments

	# UnionFilter examples

	@staticmethod
	def unionFilterExample1():
		t1 = time.time()
		print("** All structures with null OR non-null resolution **")
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithNullResolution(False) | StructureFilter.WithNullResolution(True))
		FilterTester.test(sr, 848)

	@staticmethod
	def unionFilterExample2():
		print("** All structures with null AND non-null resolution**") 
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithNullResolution(False))
		sr.addFilter(StructureFilter.WithNullResolution(True))
		FilterTester.test(sr, 0)

	# StructureResultSet examples

	@staticmethod
	def allStructures():
		'''Select all Structure records.'''
		print("** All structures **")
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb)
		FilterTester.test(sr, 848)

	@staticmethod
	def getStructuresWithNullResolutionSQL():
		print("** All structures with null resolution **") 
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb, SQL = "WHERE Resolution IS NULL")
		FilterTester.test(sr, 95)

	@staticmethod
	def getStructuresWithNullResolutionFilter():
		print("** All structures with null resolution **") 
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithNullResolution(True))
		FilterTester.test(sr, 95)

	@staticmethod
	def pickSpecific():
		'''Select four specific Structure records and apply a filter.''' 
		print("** 4 specific structures **") 
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb, AdditionalIDs = ['2BQC', '1LAW', '1LHH', '1LHI'])
		sr.addFilter(StructureFilter.TotalBFactors(0,16) | StructureFilter.WithNullResolution(True))
		FilterTester.test(sr, 2)

	@staticmethod
	def getStructuresInResolutionRange():
		print("** All structures with null resolution **") 
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.Resolution(1, 2))
		FilterTester.test(sr, 512)

	@staticmethod
	def getStructuresWithUniProtIDs():
		print("** All structures with null resolution **") 
		FilterTester.openDB()
		sr = StructureResultSet(ddGdb)
		sr.addFilter(StructureFilter.WithUniProtIDs(["P0A7Y4"], ["RNH_ECOLI", "RNP30_RANPI"]))
		FilterTester.test(sr, 15)

	@staticmethod
	def getStructuresFilteredByStructures():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		FilterTester.openDB()
		
		sr1 = StructureResultSet(ddGdb, SQL = "WHERE PDB_ID LIKE %s", parameters = "1A%")
		FilterTester.test(sr1, 53)
		
		sr2 = StructureResultSet(ddGdb, SQL = "WHERE PDB_ID LIKE %s", parameters = "1AY%")
		FilterTester.test(sr2, 2)
		
		sr = sr1.filterBySet(sr2)
		FilterTester.test(sr, 2)


	# ExperimentResultSet examples

	@staticmethod
	def getExperimentsWithSQL():
		'''Select all Structure records.'''
		print("** All structures **") 
		FilterTester.openDB()
		er = ExperimentResultSet(ddGdb, SQL = "WHERE Structure LIKE %s", parameters = "1A%")
		FilterTester.test(er, 287)

		er.addFilter(StructureFilter.Resolution(1, 1.7))
		
		FilterTester.test(er, 1)

	@staticmethod
	def getExperimentsFilteredByStructures():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		FilterTester.openDB()
		
		sr = StructureResultSet(ddGdb, SQL = "WHERE PDB_ID LIKE %s", parameters = "1AY%")
		FilterTester.test(sr, 2)
		
		er = ExperimentResultSet(ddGdb, SQL = "WHERE Structure LIKE %s", parameters = "1A%")
		FilterTester.test(er, 287)
		
		er = er.filterBySet(sr)
		FilterTester.test(er, 30)
		
		er = ExperimentResultSet(ddGdb, SQL = "WHERE Structure LIKE %s", parameters = "1AY%")
		FilterTester.test(er, 30)
		
		#print(er.structure_map.keys())
		
		er.addFilter(StructureFilter.Resolution(1, 1.80))
		FilterTester.test(er, 19)

		er.addFilter(StructureFilter.Resolution(1, 1.70))
		FilterTester.test(er, 0)

	@staticmethod
	def getExperimentsFilteredBySource():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		FilterTester.openDB()
		
		er = ExperimentResultSet(ddGdb)
		FilterTester.test(er, 14151)
		
		er.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
		
		FilterTester.test(er)

	
	@staticmethod
	def getExperimentsFilteredByMutationSize():
		'''Select all Structure records.'''
		print("** Experiments filtered by mutation size **") 
		FilterTester.openDB()

		total_count = 13208
		filters_matrix = [
			(ExperimentFilter.large, ExperimentFilter.large, 3098),
			(ExperimentFilter.large, ExperimentFilter.small, 3401),
			(ExperimentFilter.small, ExperimentFilter.large, 3262),
			(ExperimentFilter.small, ExperimentFilter.small, 3447),
		]
		assert(total_count == sum([f[2] for f in filters_matrix]))

		for f in filters_matrix:
			er = ExperimentResultSet(ddGdb)
			er.addFilter(ExperimentFilter.NumberOfMutations(1, 1))
			FilterTester.test(er, total_count)

			er.addFilter(ExperimentFilter.MutationsBetweenAminoAcidSizes(f[0], f[1]))
			FilterTester.test(er, f[2])

	@staticmethod
	def getExperimentsFilteredByMutationSize_faster():
		'''Another example of speedups using compound filters rather than specific filters.
			Time for getExperimentsFilteredByMutationSize        on my machine: @1.93s
			Time for getExperimentsFilteredByMutationSize_faster on my machine: @1.65s
			For the getExperimentsFilteredByMutationSize run, I disabled the first FilterTester.test call (otherwise it takes @2.8s).
		'''
		print("** Experiments filtered by mutation size **")
		FilterTester.openDB()

		total_count = 13208
		filters_matrix = [
			(ExperimentFilter.large, ExperimentFilter.large, 3098),
			(ExperimentFilter.large, ExperimentFilter.small, 3401),
			(ExperimentFilter.small, ExperimentFilter.large, 3262),
			(ExperimentFilter.small, ExperimentFilter.small, 3447),
		]
		assert(total_count == sum([f[2] for f in filters_matrix]))

		for f in filters_matrix:
			er = ExperimentResultSet(ddGdb)

			ef = ExperimentFilter()
			ef.setNumberOfMutations(1, 1)
			ef.setAminoAcidSizes(f[0], f[1])
			er.addFilter(ef)
			FilterTester.test(er, f[2])

	@staticmethod
	def getExperimentsFilteredByAminoAcids1():
		'''Select all Structure records.'''
		print("** Experiments filtered by residue (from ALA) **") 
		FilterTester.openDB()
		
		er = ExperimentResultSet(ddGdb)
		FilterTester.test(er, total_number_of_experiments)
		
		er.addFilter(ExperimentFilter.MutationsBetweenAminoAcids('ALA', 'G'))
		
		FilterTester.test(er, 144)

	@staticmethod
	def getExperimentsFilteredByAminoAcids2():
		'''Select all Structure records.'''
		print("** Experiments filtered by residue (from ALA) **") 
		FilterTester.openDB()
		
		er = ExperimentResultSet(ddGdb)
		FilterTester.test(er, total_number_of_experiments)
		
		er.addFilter(ExperimentFilter.MutationsBetweenAminoAcids('A', 'GLY'))
		
		FilterTester.test(er, 144)

	@staticmethod
	def getExperimentsFilteredByResolution():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **")
		FilterTester.openDB()

		er = ExperimentResultSet(ddGdb)
		FilterTester.test(er, total_number_of_experiments)

		er.addFilter(StructureFilter.Resolution(1, 2))
		FilterTester.test(er, 973)

	@staticmethod
	def getExperimentsFilteredBySourceAndResolution():
		'''Select all Structure records.'''
		print("** Experiments filtered by structures **") 
		FilterTester.openDB()
		
		er = ExperimentResultSet(ddGdb)
		FilterTester.test(er, total_number_of_experiments)
		
		er.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
		
		FilterTester.test(er)
		
		er.addFilter(StructureFilter.Resolution(1, 2))
		FilterTester.test(er)
		
		
	# PredictionResultSet examples

	@staticmethod
	def getAllPredictions():
		'''Select all Prediction records.'''
		print("** All predictions **")
		FilterTester.openDB()
		pr = PredictionResultSet(ddGdb)
		FilterTester.test(pr, 14373)

	@staticmethod
	def getPredictionsWithSQL():
		'''Select all Structure records.'''
		print("** Specific prediction **")
		FilterTester.openDB()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s AND ID=14061", parameters = "lin-3K0NA")
		FilterTester.test(pr, 1)

	@staticmethod
	def getPredictionsUsingMultipleFilters():
		'''This demonstrates the use of multiple filters.'''
		print("** Multiple filter example **") 
		FilterTester.openDB()
		t1 = time.time()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s", parameters = "AllExperimentsProtocol16")
		print(time.time() - t1)
		t1 = time.time()
		pr.addFilter(StructureFilter.Techniques(StructureFilter.XRay))
		print(time.time() - t1)
		t1 = time.time()
		pr.addFilter(StructureFilter.Resolution(1, 1.5) | StructureFilter.Resolution(3.9, 4))
		print(time.time() - t1)
		t1 = time.time()
		pr.addFilter(StructureFilter.TotalBFactors(0, 10))
		print(time.time() - t1)
		t1 = time.time()
		FilterTester.test(pr, 30)
		print(time.time() - t1)
		t1 = time.time()

	@staticmethod
	def getPredictionsUsingMultipleFilters_Speed():
		'''This demonstrates how slow separate filters are. Separate filters query the entire table whereas single
			filters with multiple criteria drill down further and further, working on subsets of the table.'''
		print("** Multiple filter example **") 
		FilterTester.openDB()

		t1 = time.time()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s", parameters = "AllExperimentsProtocol16")
		pr.addFilter(StructureFilter.Techniques(StructureFilter.XRay))
		pr.addFilter(StructureFilter.Resolution(1, 1.5) | StructureFilter.Resolution(3.9, 4))
		pr.addFilter(StructureFilter.TotalBFactors(0, 10))
		FilterTester.test(pr, 30)
		t2 = time.time()

		print("Time taken: %0.2fs" % (t2 - t1))

		t1 = time.time()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE PredictionSet=%s", parameters = "AllExperimentsProtocol16")
		sf = StructureFilter()
		sf.setTechniques(StructureFilter.XRay)
		sf.setResolution(1, 1.5)
		sf.setTotalBFactors(0, 10)
		pr.addFilter(sf | StructureFilter.Resolution(3.9, 4))
		FilterTester.test(pr, 30)
		t2 = time.time()
		print("Time taken: %0.2fs" % (t2 - t1))

	@staticmethod
	def showResultSetOperations():
		'''Demonstrates how to union, intersect, subtract, and XOR ResultSets.'''

		range_size1 = 14
		range_size2 = 187
		range_size3 = 506
		range_size4 = 8

		print("\n** ResultSet SR1 **\n")
		FilterTester.openDB()
		sr1 = StructureResultSet(ddGdb)
		sr1.addFilter(StructureFilter.Resolution(1, 1.3))
		FilterTester.test(sr1, 14)

		print("\n** ResultSet SR2 **\n")
		sr2 = StructureResultSet(ddGdb)
		sr2.addFilter(StructureFilter.Resolution(2, 2.3))
		FilterTester.test(sr2, 187)
		
		print("\n** ResultSet SR3 **\n")
		sr3 = StructureResultSet(ddGdb)
		sr3.addFilter(StructureFilter.Resolution(1.2, 2))
		FilterTester.test(sr3, 506)
		
		print("\n** ResultSet union - SR1 | SR2 **\n")
		srUnion = sr1 | sr2
		assert(len(srUnion) == range_size1 + range_size2)
		print(join(srUnion._log, "\n"))
		
		print("\n** ResultSet union - SR1 - SR3 **\n")
		srUnion = sr1 - sr3
		print(join(srUnion._log, "\n"))

		print("\n** ResultSet intersection - SR1 & SR3 **\n")
		srIntersection = sr1 & sr3
		print(join(srIntersection._log, "\n"))

		print("\n** ResultSet difference and union sanity check **\n")
		assert(len(srUnion) + len(srIntersection) == range_size1)

		print("\n** ResultSet intersection sanity check **\n")
		sr4 = StructureResultSet(ddGdb)
		sr4.addFilter(StructureFilter.Resolution(1.2, 1.3))
		FilterTester.test(sr4, range_size4)
		
		print("\n** ResultSet difference - SR1 - SR3 **\n")
		
		srDifference = sr1 / sr3
		print(join(srDifference._log, "\n"))

		print("\n** ResultSet exclusive or - SR1 ^ SR3 **\n")
		srXOR = sr1 ^ sr3
		print(join(srXOR._log, "\n"))

		print("\n** ResultSet exclusive or sanity check **\n")
		assert(len(srXOR) == (range_size1 - len(srIntersection)) + (range_size3 - len(srIntersection)))


	@staticmethod
	def showAllEligibleProTherm(PredictionSet, ProtocolID, KeepHETATMLines):
		#inserter = JobInserter()
		colortext.printf("\nAdding ProTherm mutations to %s prediction set." % PredictionSet, "lightgreen")
		#ddGdb = ddgdbapi.ddGDatabase()
		
		MAX_RESOLUTION = 2.1
		MAX_NUMRES_PROTHERM = 350
		MAX_STANDARD_DEVIATION = 1.0

		FilterTester.openDB()
		
		if False:
			t1 = time.time()
			er1 = ExperimentResultSet(ddGdb)
			er1.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
			er1.addFilter(ExperimentFilter.NumberOfMutations(1, 1))
			er1.addFilter(ExperimentFilter.NumberOfChains(1, 1))
			er1.addFilter(ExperimentFilter.StandardDeviation(None, MAX_STANDARD_DEVIATION))
			er1.addFilter(StructureFilter.Resolution(None, MAX_RESOLUTION))
			er1.addFilter(StructureFilter.Techniques(StructureFilter.XRay))
			FilterTester.test(er1)
			t2 = time.time()
			print(t2 - t1)
		
		# This method usually takes around 65% of the time as the method above 
		t1 = time.time()
		ef1 = ExperimentFilter()
		ef1.setSource(ExperimentFilter.ProTherm)
		er1 = ExperimentResultSet(ddGdb)
		er1.addFilter(ExperimentFilter.OnSource(ExperimentFilter.ProTherm))
		FilterTester.test(er1)
		ef1.setNumberOfMutations(1, 1)
		ef1.setNumberOfChains(1, 1)
		ef1.setStandardDeviation(None, MAX_STANDARD_DEVIATION)
		sf1 = StructureFilter()
		sf1.setResolution(None, MAX_RESOLUTION)
		sf1.setTechniques(StructureFilter.XRay)
		er1 = ExperimentResultSet(ddGdb)
		er1.addFilter(ef1)
		er1.addFilter(sf1)
		FilterTester.test(er1)
		t2 = time.time()
		print(t2 - t1)
		
		experimentIDs = sorted(list(er1.getFilteredIDs()))
		colortext.message("\nThe number of unique ProTherm experiments with:\n\t- one mutation;\n\t- structures solved by X-ray diffraction and with <= %d residues;\n\t- a maximum standard deviation in experimental results of <= %0.2f;\n\t- and a resolution of <= %0.2f Angstroms.\nis %d.\n" % (MAX_NUMRES_PROTHERM, MAX_STANDARD_DEVIATION, MAX_RESOLUTION, len(experimentIDs)))
		ddG_connection = db_api.ddG()
		count = 0
		sys.exit(0)
		print("")
		for experimentID in experimentIDs:
			ddG_connection.addPrediction(experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True)
			count += 1
			if count >= 10:
				colortext.write(".")
				colortext.flush()
				count = 0
		print("")


	@staticmethod
	def testPublications():
		ddG_connection = db_api.ddG()

		FilterTester.openDB()
		pr = PredictionResultSet(ddGdb, SQL = "WHERE ID >= 28331 and ID <= 28431")
		print(1)
		er = ExperimentResultSet(ddGdb, SQL = "WHERE ID >= 110906 and ID <= 111006")
		print(2)
		print(len(pr))
		print(len(er))

		ddG_connection.getPublications(pr)
		ddG_connection.getPublications(er)
	

				


if False:
	import analysis
	analyzer = analysis.Analyzer("AllExperimentsProtocol16")
	analyzer.AddPublishedDDGsToAnalysisTables()
	analyzer.plot(analysis.Analyzer.correlation_coefficient, "Kellogg.rr", table_names = ["Kellogg"])
	#			"kellogg.txt", "Kellogg")
	for table_name, a_table in sorted(analyzer.analysis_tables.iteritems()):
		print(a_table)
		print(table_name)
	
#print(analysis.AnalysisPoint.headers)
#print(analysis_tables)
#print(analysis.AnalysisPoint.headers)
#print(analysis_tables["Kellogg"])

#ddG_connection.createPredictionsFromUserDataSet("AllValidPGPK", "AllExperimentsProtocol16", "Kellogg:10.1002/prot.22921:protocol16:32231", False, StoreOutput = False, Description = {}, InputFiles = {}, testonly = False)
	
#ddG_connection = db_api.ddG()
#ddG_connection.addPDBtoDatabase(pdbID = "1FKJ")

ddG_connection = db_api.ddG()

if __name__ == '__main__':
	#help()
    # Tested functions
	tests = [
		FilterTester.unionFilterExample1,
		FilterTester.unionFilterExample2,
		FilterTester.allStructures,
		FilterTester.getStructuresWithNullResolutionSQL,
		FilterTester.getStructuresWithNullResolutionFilter,
		FilterTester.pickSpecific,
		FilterTester.getStructuresInResolutionRange,
		FilterTester.getStructuresWithUniProtIDs,
		FilterTester.getStructuresFilteredByStructures,
		FilterTester.getExperimentsWithSQL,
		FilterTester.getExperimentsFilteredByStructures,
		FilterTester.getExperimentsFilteredByMutationSize,
		FilterTester.getExperimentsFilteredByMutationSize_faster,
		FilterTester.getExperimentsFilteredByAminoAcids1,
		FilterTester.getExperimentsFilteredByAminoAcids2,
		FilterTester.getExperimentsFilteredByResolution,
		FilterTester.getAllPredictions,
		FilterTester.getPredictionsWithSQL,
		FilterTester.getPredictionsUsingMultipleFilters,
		FilterTester.getPredictionsUsingMultipleFilters_Speed,
		FilterTester.showResultSetOperations,
	]
	tests = [FilterTester.getExperimentsFilteredByMutationSize]


	# BROKEN FUNCTIONS
	#FilterTester.getExperimentsFilteredBySource() # Needs update w.r.t. new database schema
	#FilterTester.getExperimentsFilteredBySourceAndResolution() # Needs update w.r.t. new database schema
	#FilterTester.showAllEligibleProTherm("test", "test", False) # Needs update w.r.t. new database schema
	#FilterTester.testPublications # Needs update w.r.t. new database schema

	do_profiling = False

	import gc
	gc.disable()
	if do_profiling:
		FilterTester.profile('FilterTester.getExperimentsFilteredByMutationSize_faster()')
	else:
		for t in tests:
			t1 = time.time()
			t()
			print("** Total time taken in getExperimentsFilteredByMutationSize_faster: %0.2f **" % (time.time() - t1))
	gc.enable()



