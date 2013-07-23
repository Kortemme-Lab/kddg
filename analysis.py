import sys
import os
import re
import pickle
import subprocess
from string import join
import time
import inspect
from tempfile import mkstemp

import tools.colortext as colortext
import ddgdbapi
import tools.deprecated.rosettahelper as rosettahelper
from tools.deprecated.rosettahelper import kJtokcal

import tools.latex as latex

script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class RInterface(object):
	
	@staticmethod
	def _runRScript(RScript):
		rscriptname = rosettahelper.writeTempFile(".", RScript)
		#p = subprocess.Popen(["/opt/R-2.15.1/bin/R","CMD", "BATCH", rscriptname])
		p = subprocess.Popen(["R","CMD", "BATCH", rscriptname])
		while True:
			time.sleep(0.3)
			errcode = p.poll()
			if errcode != None:
				break
		rout = "%s.Rout" % rscriptname
		os.remove(rscriptname)
		#colortext.warning(rosettahelper.readFile(rout)) 
		if errcode != 0:
			if os.path.exists(rout):
				colortext.warning(rosettahelper.readFile(rout)) 
				os.remove(rout)
			raise colortext.Exception("The R script failed with error code %d." % errcode)
		os.remove(rout)
	
	# R methods
	
	@staticmethod
	def _correlation_coefficient_unfixed(inputfname, outfname, filetype):
		'''File suffix: pearsons_r_unfixedDDGs
		   Description: Pearson's r using dataset DDGs
		   '''
		return RInterface.correlation_coefficient(inputfname, outfname, filetype, experiment_field = "DatasetPublishedDDG")
	
	@staticmethod
	def _correlation_coefficient(inputfname, output_filename, filetype, experiment_field = "ExperimentalDDG"):
		'''File suffix: pearsons_r
		   Description: Pearson's r
		   '''
		# todo: The graph size is not being set correctly i.e. the max function call is wrong in the script
		title = "" # "Pearson's r"
		RScript = rosettahelper.readFile(os.path.join(script_path, "R", "pearsons.R")) % vars()
		RInterface._runRScript(RScript)
		
	@staticmethod
	def correlation_coefficient_gplot(inputfname, output_filename, filetype, experiment_field = "ExperimentalDDG"):
		'''File suffix: pearsons_r_gplot
		   Description: Pearson's r
		   Filename: ggplot_pearsons.R
		   Priority: 1
		   '''
		title = "" #Pearson's r"
		RScript = rosettahelper.readFile(os.path.join(script_path, "R", "ggplot_pearsons.R")) % vars()
		RInterface._runRScript(RScript)

	@staticmethod
	def correlation_coefficient_gplot_unfixed(inputfname, outfname, filetype):
		'''File suffix: pearsons_r_unfixedDDGs
		   Description: Pearson's r using dataset DDGs
		   Priority: 2
		   Filename: ggplot_pearsons.R
		   '''
		return RInterface.correlation_coefficient_gplot(inputfname, outfname, filetype, experiment_field = "DatasetPublishedDDG")

	@staticmethod
	def correlation_coefficient_gplot_crop(inputfname, output_filename, filetype, experiment_field = "ExperimentalDDG"):
		'''File suffix: pearsons_r_gplot_cropped
		   Description: Pearson's r filtering with quantile(res, probs = c(0.05,0.95))
		   Filename: ggplot_pearsons_cropped.R
		   Priority: 3
		   '''
		title = "" #Pearson's r filtering with quantile(res, probs = c(0.05,0.95))"
		RScript = rosettahelper.readFile(os.path.join(script_path, "R", "ggplot_pearsons_cropped.R")) % vars()
		RInterface._runRScript(RScript)
	

# Create array of R functions
#RFunctions = {
#	"pearsons_r"				: ("Pearson's r", RInterface.correlation_coefficient),
#	"pearsons_r_unfixedDDGs"	: ("Pearson's r using dataset DDGs", RInterface.correlation_coefficient_unfixed),
#}
re_filesuffix = re.compile(".*File suffix:(.*?)\n")
re_description = re.compile(".*Description:(.*?)\n", re.DOTALL)
re_priority = re.compile(".*Priority:(.*?)\n", re.DOTALL)
re_filename = re.compile(".*Filename:(.*?)\n", re.DOTALL)
RFunctions = {} 
rfunctions = sorted([m[1] for m in inspect.getmembers(RInterface) if inspect.isfunction(m[1]) and m[0][0] != "_"])
for rfunction in rfunctions:
	docstr = (rfunction.__doc__)

	filesuffix = re_filesuffix.match(docstr)
	assert(filesuffix)
	filesuffix = filesuffix.group(1).strip()
	
	description = re_description.match(docstr)
	assert(description)
	description = description.group(1).strip()
	
	priority = re_priority.match(docstr)
	assert(priority)
	priority = priority.group(1).strip()

	filename = re_filename.match(docstr)
	assert(filename)
	filename = filename.group(1).strip()

	assert(not RFunctions.get(filesuffix))
	RFunctions[filesuffix] = (description, rfunction, int(priority), filename)
	

class UserDataSetExperimentalScores(object):
	def __init__(self, ddGdb, UserDataSetID, AnalysisSet):
		scores = {}
		exp_scores = ddGdb.execute('''
SELECT
UserAnalysisSet.*, ExperimentAssayDDG.Value AS ExperimentalDDG
FROM UserAnalysisSet 
INNER JOIN ExperimentAssayDDG ON UserAnalysisSet.ExperimentAssayID=ExperimentAssayDDG.ExperimentAssayID AND	UserAnalysisSet.Type=ExperimentAssayDDG.Type
WHERE UserDataSetID=%s AND Subset=%s ORDER BY Section, RecordNumber''', parameters=(UserDataSetID, AnalysisSet))

		for exp_score in exp_scores:
			Section = exp_score["Section"]
			RecordNumber = exp_score["RecordNumber"]
			
			# Store the experimental DDG values
			scores[Section] = scores.get(Section, {})
			if scores[Section].get(RecordNumber):
				scores[Section][RecordNumber]["ExperimentalDDG"] = scores[Section][RecordNumber].get("ExperimentalDDG", [])
				scores[Section][RecordNumber]["ExperimentalDDG"].append(exp_score["ExperimentalDDG"])
			else:
				scores[Section][RecordNumber] = {
					"ExperimentID" : exp_score["ExperimentID"],
					"PDB_ID" : exp_score["PDB_ID"],
					"ExperimentalDDG" : [exp_score["ExperimentalDDG"]]
				}
			
		for section, sectiondata in scores.iteritems():
			for recordnumber, recorddata in sectiondata.iteritems():
				numrecords = float(len(recorddata["ExperimentalDDG"]))
				recorddata["ExperimentalDDG"] = float(sum(recorddata["ExperimentalDDG"]))/numrecords 
		
		self.UserDataSetID = UserDataSetID
		self.AnalysisSet = AnalysisSet
		self.scores = scores
		
	def iteritems(self):
		return self.scores.iteritems()

class PredictionScores(object):

	def __init__(self, ddGdb, PredictionSet, ddG_score_type = 'kellogg.total'):
		import pickle
		
		# Get the UserDataSet ID and the list of AnalysisSets associated with the Prediction set
		UserDataSetID = ddGdb.execute("SELECT DISTINCT UserDataSetID FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=UserDataSetExperimentID WHERE PredictionSet=%s AND Status='done'", parameters=(PredictionSet,))
		if len(UserDataSetID) == 0:
			raise colortext.Exception("There are no UserDataSets associated with this Prediction set. Unable to proceed automatically.")
		if len(UserDataSetID) > 1:
			raise colortext.Exception("There are %d UserDataSets associated with this Prediction set. Unable to proceed automatically." % len(UserDataSetID))
		UserDataSetID = UserDataSetID[0]["UserDataSetID"]
		
		UserDataSetName = ddGdb.execute("SELECT TextID FROM UserDataSet WHERE ID=%s", parameters=(UserDataSetID,))[0]["TextID"]
		
		NumberOfPredictions = ddGdb.execute("SELECT COUNT(ID) AS NumberOfPredictions FROM Prediction WHERE PredictionSet=%s AND Status='done'", parameters=(PredictionSet,))[0]["NumberOfPredictions"]
	
		AnalysisSets = ddGdb.execute("SELECT DISTINCT Subset FROM UserAnalysisSet WHERE UserDataSetID=%s", parameters=(UserDataSetID,))
		if len(AnalysisSets) == 0:
			raise colortext.Exception("There are no analysis sets associated with UserDataSet %s. Unable to proceed automatically." % UserDataSetName)
		AnalysisSets = sorted([r["Subset"] for r in AnalysisSets])
		
		# Get list of Kellogg record IDs for which we have predictions and experimental DDG values
		dbpredictions = ddGdb.execute("SELECT Prediction.ID as PredictionID, Prediction.ExperimentID AS PExperimentID, UserDataSetExperiment.ExperimentID AS UDSEExperimentID, UserDataSetExperiment.PDBFileID, UserDataSetExperiment.ID AS UserDataSetExperimentID, ddG FROM Prediction INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=UserDataSetExperimentID WHERE PredictionSet=%s AND Status='done'", parameters=(PredictionSet,))
		Predictions = {}
		for p in dbpredictions:
			PDB_ID = p["PDBFileID"]
			Predictions[p["PExperimentID"]] = Predictions.get(p["PExperimentID"], {})
			if Predictions[p["PExperimentID"]].get(PDB_ID):
				raise colortext.Exception("There are multiple predictions for Experiment ID %d with PDB ID %s in the PredictionSet. This case is currently not handled (perhaps averaged values would be acceptable)." % (p["PExperimentID"], PDB_ID))
			#Predictions[p["PExperimentID"]][PDB_ID] = Predictions[p["PExperimentID"]].get(PDB_ID, {})
			
			assert(p["PExperimentID"] == p["UDSEExperimentID"])
			ddG = pickle.loads(p["ddG"])

			if (ddG['version'] == '0.1'):
				pass
			elif (ddG['version'] != '0.23'):
				continue

			# Traverse the score hierarchy
			PredictedDDG = ddG["data"]
			#try:
			for score_path in ddG_score_type.split("."):
				PredictedDDG = PredictedDDG[score_path]
			PredictedDDG = PredictedDDG['ddG']
			#except:
			#	continue
					
			Predictions[p["PExperimentID"]][PDB_ID] = {
				#"ExperimentID" : p["PExperimentID"],
				#"PDB_ID" : p["PDB_ID"],
				"UserDataSetExperimentID" : p["UserDataSetExperimentID"],
				"PredictedDDG" : PredictedDDG
			}

		self.Predictions = Predictions
		self.PredictionSet = PredictionSet
		self.UserDataSetID = UserDataSetID
		self.UserDataSetName = UserDataSetName
		self.AnalysisSets = AnalysisSets
		self.NumberOfPredictions = NumberOfPredictions
	
	def iteritems(self):
		return self.Predictions.iteritems()

class PublishedDatasetScores(object):
	def __init__(self, ddGdb, ShortID):
		'''Creates a dict with the published dataset DDG in kcal/mol (Rosetta convention), the published PDB ID, and the fixed PDB ID (which may be the same as the published).''' 
		scores = {}
		results = ddGdb.execute("SELECT DDGConvention, Section, RecordNumber, PDBFileID, PublishedPDBFileID, PublishedValue as PublishedValueInKcal, AggregateType FROM DataSetDDG INNER JOIN DataSet ON DataSet.ID=DataSetID WHERE ShortID=%s", parameters=(ShortID, ))
		for r in results:
			scores[r["Section"]] = scores.get(r["Section"], {})
			assert(not(scores[r["Section"]].get(r["RecordNumber"])))
			assert(r["DDGConvention"] == "ProTherm" or r["DDGConvention"] == "Rosetta")
			DDG = r["PublishedValueInKcal"] 
			if r["DDGConvention"] == "ProTherm":
				DDG =- DDG
			
			scores[r["Section"]][r["RecordNumber"]] = {
				"PDB_ID" : r["PDBFileID"],
				"PublishedPDB_ID" : r["PublishedPDBFileID"],
				"PublishedDatasetDDG" : DDG,
				"AggregateType" : r["AggregateType"],
			}
		self.scores = scores

class AnalysisPoint(object):
	
	headers = ["ExperimentalDDG", "PredictedDDG", "ExperimentID", "PDB_ID", "section", "recordnumber", "DatasetPublishedDDG"]
	set_of_headers = None
	
	def __init__(self, ExperimentalDDG, PredictedDDG, ExperimentID = None, PDB_ID = None, DatasetPublishedDDG = None, section = None, recordnumber = None):
		if not AnalysisPoint.set_of_headers:
			AnalysisPoint.set_of_headers = set(AnalysisPoint.headers)
			
		self.ExperimentalDDG = ExperimentalDDG
		assert(isinstance(ExperimentalDDG, float))
		
		self.PredictedDDG = PredictedDDG
		assert(isinstance(PredictedDDG, float))
		
		self.DatasetPublishedDDG = DatasetPublishedDDG
		if DatasetPublishedDDG:
			assert(isinstance(DatasetPublishedDDG, float))
		
		self.ExperimentID = ExperimentID
		self.PDB_ID = PDB_ID
		
		self.section = section
		if recordnumber:
			assert(str(recordnumber).isdigit())
			recordnumber = int(recordnumber)
		self.recordnumber = recordnumber
	
	def __setattr__(self, key, val):
		if key not in AnalysisPoint.set_of_headers:
			raise Exception("Cannot set key '%s' to value '%s': key is not allowed." % (key, val))
		else:
			super(AnalysisPoint, self).__setattr__(key, val)
		
	def __repr__(self, delimiter = ","):
		s = []
		for h in self.headers:
			if self.__dict__[h] == 0:
				s.append("0")
			else:
				s.append(str(self.__dict__[h] or ""))
		return join(s, delimiter)

class AnalysisTable(object):

	def __init__(self):
		self.points = []
		
	def add(self, point):
		assert(isinstance(point, AnalysisPoint))
		self.points.append(point)
		
	def __repr__(self, delimiter = ","):
		s = [join(AnalysisPoint.headers, delimiter)]
		for p in self.points:
			s.append(str(p))
		return join(s, "\n")

class AnalysisObject(object):
	
	def __init__(self, dataset, description, filetype, contents):
		self.dataset = dataset
		self.description = description
		self.filetype = filetype
		self.contents = contents
	
	def __repr__(self):
		return "Dataset: %(dataset)s\nAnalysis type: %(description)s\nFiletype: %(filetype)s" % self.__dict__


class Analyzer(object):

	def __init__(self, PredictionSet, quiet = False, ddG_score_type = 'kellogg.total'):
		self.PredictionSet = PredictionSet
		self.quiet = quiet
		self.ddG_score_type = ddG_score_type
		
		self.analysis_tables = None
		self.ddGdb = ddgdbapi.ddGDatabase()
		self.description = []
		self.CreateAnalysisTables()
	
	def CreateAnalysisTables(self):
		ddGdb = self.ddGdb
		PredictionSet = self.PredictionSet
		predictions = PredictionScores(ddGdb, PredictionSet, self.ddG_score_type)
		predicted_scores = predictions.Predictions
		
		s = "Analyzing %d predictions in PredictionSet '%s' for UserDataSet '%s'. " % (predictions.NumberOfPredictions, predictions.PredictionSet.replace("_", "\_"), predictions.UserDataSetName)
		s += "Running analysis over the following analysis sets: '%s'." % (join(predictions.AnalysisSets, "', '"))
		self.description.append(("black", s))
		if not self.quiet:
			colortext.message("Analyzing %d predictions in PredictionSet '%s' for UserDataSet '%s'." % (predictions.NumberOfPredictions, predictions.PredictionSet, predictions.UserDataSetName))
			colortext.message("Running analysis over the following analysis sets: '%s'." % (join(predictions.AnalysisSets, "', '")))
		
		analysis_tables = {}
		# Analyze data for
		for AnalysisSet in predictions.AnalysisSets:
			analysis_table = AnalysisTable()
			
			experiments = UserDataSetExperimentalScores(ddGdb, predictions.UserDataSetID, AnalysisSet)
			
			count = 0
			numMissing = 0
			for section, sectiondata in sorted(experiments.iteritems()):
				for recordnumber, record_data in sorted(sectiondata.iteritems()):
					count += 1
					PDB_ID = record_data["PDB_ID"]
					ExperimentID = record_data["ExperimentID"]
					ExperimentalDDG = record_data["ExperimentalDDG"]
					if predicted_scores.get(ExperimentID) and predicted_scores[ExperimentID].get(PDB_ID):
						PredictedDDG = predicted_scores[ExperimentID][PDB_ID]["PredictedDDG"]
						analysis_table.add(AnalysisPoint(ExperimentalDDG, PredictedDDG, ExperimentID = ExperimentID, PDB_ID = PDB_ID, section = section, recordnumber = recordnumber))
					else:
						numMissing += 1
			if numMissing > 0 and not self.quiet:
				self.description.append(("Bittersweet", "Missing %d predictions out of %d records for analysis set %s." % (numMissing, count, AnalysisSet)))
				colortext.warning("Missing %d predictions out of %d records for analysis set %s." % (numMissing, count, AnalysisSet))
			analysis_tables[AnalysisSet] = analysis_table
		
		self.analysis_tables = analysis_tables 
	
	def AddPublishedDDGsToAnalysisTables(self):
		ddGdb = self.ddGdb
		analysis_tables = self.analysis_tables
		for AnalysisSet, analysis_table in analysis_tables.iteritems():
			published_dataset_scores = PublishedDatasetScores(ddGdb, AnalysisSet).scores
			
			for analysis_point in analysis_table.points:
				if analysis_point.section and analysis_point.recordnumber:
					section = analysis_point.section
					recordnumber = analysis_point.recordnumber
					if published_dataset_scores.get(section) and published_dataset_scores[section].get(recordnumber):
						published_dataset_score = published_dataset_scores[section][recordnumber]["PublishedDatasetDDG"]
						analysis_point.DatasetPublishedDDG = published_dataset_score
					else:
						if not self.quiet:
							colortext.warning("No published dataset score found for %s-%s-%s." % (AnalysisSet, Section, RecordNumber))
	
	@staticmethod
	def _mean(points, numpoints):
		'''Points is expected to be a list (or iterable), numpoints should be an integer.'''
		return sum(points) / numpoints

	def CreateCSVFile(self, analysis_table_name, path = "."):
		analysis_table = self.analysis_tables[analysis_table_name]
		return rosettahelper.writeTempFile(path, str(analysis_table))
	
	def PlotAll(self, createFiles = True, filetype = "pdf"):
		import operator
		filenames = []
		filecontents = []
		for table_name, analysis_table in self.analysis_tables.iteritems():
			for analysisType, data in sorted(RFunctions.iteritems(), key=lambda x: x[1][2]):
				description = data[0]
				RFunction = data[1]
				if createFiles:
					outfname = "%s_%s.%s" % (table_name.lower().replace(" ", "_"), analysisType, filetype)
					self.plot(table_name, RFunction, outfname, filetype = filetype)
					filenames.append(outfname)
				else:
					filecontents.append(self.plot(table_name, RFunction, filetype = filetype))
		if createFiles:
			return filenames
		else:
			return filecontents
	
	def plot(self, table_name, RFunction, output_filename = None, filetype = "pdf"):
		'''Results is expect to be a list of dicts each of which has the keys ExperimentID and ddG.'''
		if (not self.analysis_tables) or (not table_name):
			raise Exception("There are no analysis tables to plot.")
		if not table_name in self.analysis_tables.keys():
			raise Exception("The analysis table '%s' does not exist." % table_name)
		
		gplot = None
		analysis_table = self.analysis_tables[table_name]
		print(table_name)
		print(RFunction)
		if len(analysis_table.points) == 1:
			raise Exception("The analysis table %s set only has one data point. At least two points are required." % table_name)
		else:
			inputfname = self.CreateCSVFile(table_name)
			print(inputfname)
			try:
				if not self.quiet:
					colortext.printf("Running %s." % RFunction)
					if output_filename:
						colortext.printf("Saving graph as %s with filename %s." % (filetype, output_filename))
				
				output_fname = output_filename
				if not output_fname:
					output_fname = rosettahelper.writeTempFile(".", "")
										
				RFunction(inputfname, output_fname, filetype)
				if not output_filename:
					contents = rosettahelper.readBinaryFile(output_fname)
					os.remove(output_fname)
					description = None
					for file_suffix, details in RFunctions.iteritems():
						if details[1] == RFunction:
							description = details[0]
					assert(description)
					gplot = AnalysisObject(table_name, description, filetype, contents)
				else:
					gplot = output_filename
			except Exception, e:
				import traceback
				colortext.error(traceback.format_exc())
				os.remove(inputfname)
				raise Exception(e)
			os.remove(inputfname)
		return gplot

class Reporter(object):
	
	def __init__(self, analyzer):
		self.graphfiles = []
		self.analyzer = analyzer
	
	def clean(self):
		for graphfile in self.graphfiles:
			if os.path.exists(graphfile):
				os.remove(graphfile)
		self.graphfiles = []
		
	def CreateReport(self, outfname = 'test.pdf', description = None, filetype = "pdf"):
		
		latexdoc = latex.LaTeXDocument("$\\Delta\\Delta$G analysis", pdf_title = "DDG analysis")
		
		analysis_objects = self.analyzer.PlotAll(filetype = filetype, createFiles = False)
		
		analysis_by_dataset = {}
		for analysis_object in analysis_objects:
			dataset = analysis_object.dataset
			analysis_by_dataset[dataset] = analysis_by_dataset.get(dataset, [])
			analysis_by_dataset[dataset].append(analysis_object)
			
		count = 0
		try:
			if description:
				latexdoc.addSection("Description")
				for tagged_line in description:
					color = tagged_line[0]
					line = tagged_line[1]
					latexdoc.addLaTeXCode("\\textcolor{%s}{%s}\n" % (color, line))
			
			if analysis_by_dataset:
				latexdoc.addSection("Analysis")
				
			for dataset, analysis_objects in analysis_by_dataset.iteritems():
				latexdoc.addSubsection(dataset)
				
				# The tab solution here is messy. Redo this when we decide what LaTeX floats we want to use
				captiontabs = []
				imagetabs = []
				for analysis_object in analysis_objects:
					count +=1
					filename = "graph%d" + (".%s" % filetype)
					rosettahelper.writeFile(filename % count, analysis_object.contents)
					captiontabs.append(analysis_object.description)
					imagetabs.append(latexdoc.getImageText(filename % count, width_cm = 6))
					self.graphfiles.append(filename)
				
				assert(len(captiontabs) == len(imagetabs))
				captiontabs = [captiontabs[i:i+2] for i in range(0, len(captiontabs), 2)]
				imagetabs = [imagetabs[i:i+2] for i in range(0, len(imagetabs), 2)]
				assert(len(captiontabs) == len(imagetabs))
				
				tabs = []
				for i in range(len(captiontabs)):
					tabs.append(captiontabs[i])
					tabs.append(imagetabs[i])
				
				print(tabs)
				latexdoc.addTabular(tabs, alignment = 'c')
			
			latexdoc.clearPage()
			latexdoc.addSection("R files")
			for analysisType, data in sorted(RFunctions.iteritems(), key=lambda x: x[1][2]):
				description = data[0]
				RFunction = data[1]
				filename = data[3]
				latexdoc.addSubsection(description)
				latexdoc.addTypewriterText(rosettahelper.readFile(os.path.join(script_path, "R", filename)), language="R")
				latexdoc.clearPage()
				
			latexdoc.compile_pdf(outfname)
		except Exception, e:
			self.clean()
			raise
		self.clean()
			

# R notes for Shane
# par is used to specify graphical parameters. mar sets the margin sizes - c(bottom, left, top right)
# read.table is useful for reading delimited files. e.g. read.table("c:/mydata.csv", header=TRUE, sep=",", row.names="id")
# plot creates a plot
#	pch=19 means use a diamond character (web-search for 'R plot symbols' pch)
#	a$V1, a$V2 choose the first and second columns of data table a
#	xlab and ylab set the labels
#	cex	indicates the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50% larger, 0.5 is 50% smaller, etc.
#	cex.lab	defines the magnification of x and y labels relative to cex
#	abline(k,l) adds a straight line where k and l are the intercept and slope values respectively

