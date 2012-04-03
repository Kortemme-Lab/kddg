#!/usr/bin/python2.4
# encoding: utf-8
"""
dbapi.py
High-level functions for interacting with the ddG database.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import os
from string import join
import common.ddgproject as ddgproject
from pdb import PDB, ResidueID2String, checkPDBAgainstMutations, aa1
#from Bio.PDB import *
from rosettahelper import write_file
import common.colortext as colortext  
import traceback
import pickle
import md5
import random
import score 
import analysis
from ddgfilters import PredictionResultSet, ExperimentResultSet, StructureResultSet 

dbfields = ddgproject.FieldNames()

class MutationSet(object):
	def __init__(self):
		self.mutations = []
	
	def addMutation(self, chainID, residueID, wildtypeAA, mutantAA):
		self.mutations.append((chainID, residueID, wildtypeAA, mutantAA))
	
	def getChains(self):
		return sorted(list(set([m[0] for m in self.mutations])))

class ddG(object):
	'''This class is responsible for inserting prediction jobs to the database.''' 
	
	def __init__(self):
		self.ddGDB = ddgproject.ddGDatabase()
		self.ddGDataDB = ddgproject.ddGPredictionDataDatabase()
	
	def __del__(self):
		self.ddGDB.close()
		self.ddGDataDB.close()
			
	def _createResfile(self, pdb, mutations):
		'''The mutations here are in the original PDB numbering. pdb is assumed to use Rosetta numbering.
			We use the pdb mapping from PDB numbering to Rosetta numbering to generate the resfile.
		'''
		resfile = []
		for mutation in mutations:
			chain = mutation[0]
			resid = mutation[1]
			wt = mutation[2]
			mt = mutation[3]
			
			# Check that the expected wildtype exists in the PDB 
			readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
			assert(wt == aa1[readwt])
			resid = resid.strip()
			resfile.append("%(resid)s %(chain)s PIKAA %(mt)s" % vars())
		if resfile:
			resfile = ["NATAA", "start"] + resfile
			return join(resfile, "\n")
		else:
			raise Exception("An error occurred creating a resfile for the ddG job.")
	
	def getData(self, predictionID):
		results = self.ddGDataDB.execute("SELECT * FROM PredictionData WHERE ID=%s", parameters = (predictionID,))
		if results:
			assert(len(results) == 1)
			return results[0]["Data"]
	
	def getPublications(self, result_set):
		if result_set:
			structures = None
			experiments = None
			
			if result_set.isOfClass(ExperimentResultSet):
				experiments = result_set 
			elif ExperimentResultSet in result_set.__class__.allowed_restrict_sets:
				experiments, experiment_map = result_set.getExperiments()			
			
			if result_set.isOfClass(StructureResultSet):
				structures = result_set 
			elif StructureResultSet in result_set.__class__.allowed_restrict_sets:
				structures, structure_map = result_set.getStructures()
			
			if structures:
				colortext.printf("\nRelated publications for structures:", "lightgreen")
				for id in sorted(structures.IDs):
					pubs = self.ddGDB.callproc("GetPublications", parameters=(id,))
					print(id)
					for pub in pubs:
						print("\t%s: %s" % (pub["Type"], pub["PublicationID"]))
					 
			if experiments:
				colortext.printf("\nRelated publications for experiments:", "lightgreen")
				for id in sorted(experiments.IDs):
					pubs = self.ddGDB.callproc("GetExperimentPublications", parameters=(id,))
					print(id)
					for pub in pubs:
						print("\t%s: %s" % (pub["Type"], pub["SourceLocation.ID"]))
						
				experimentsets = [e[0] for e in self.ddGDB.execute("SELECT DISTINCT Source FROM Experiment WHERE ID IN (%s)" % join(map(str, list(experiments.IDs)), ","), cursorClass = ddgproject.StdCursor)]
				
				if experimentsets:
					colortext.printf("\nRelated publications for experiment-set sources:", "lightgreen")
					for id in sorted(experimentsets):
						print(id)
						pubs = self.ddGDB.execute("SELECT ID, Type FROM SourceLocation WHERE SourceID=%s", parameters = (id,))
						for pub in pubs:
							print("\t%s: %s" % (pub["Type"], pub["ID"]))
		else:
			raise Exception("Empty result set.")
				

		
	def dumpData(self, outfile, predictionID):
		write_file(outfile, self.getData(predictionID))
	
	def analyze(self, prediction_result_set, outpath = os.getcwd()):
		PredictionIDs = sorted(list(prediction_result_set.getFilteredIDs()))
		colortext.printf("Analyzing %d records:" % len(PredictionIDs), "lightgreen")
		#results = self.ddGDB.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % join(map(str, PredictionIDs), ","))
		
		#for r in results:
		#	r["ddG"] = pickle.loads(r["ddG"])
		#	predicted_score = r["ddG"]["data"]["ddG"]
		#	experimental_scores = [expscore["ddG"] for expscore in self.ddGDB.callproc("GetScores", parameters = r["ExperimentID"])]
		#	mean_experimental_score = float(sum(experimental_scores)) / float(len(experimental_scores))
	
		results = self.ddGDB.execute("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % join(map(str, PredictionIDs), ","))
		
		analysis.plot(analysis._R_mean_unsigned_error, analysis._createMAEFile, results, "my_plot1.pdf", average_fn = analysis._mean)
		analysis.plot(analysis._R_correlation_coefficient, analysis._createAveragedInputFile, results, "my_plot2.pdf", average_fn = analysis._mean)
		colortext.printf("Done", "lightgreen")
		
			
		
		#score.ddgTestScore
	
	def addPDBtoDatabase(self, filepath = None, pdbID = None, protein = None, source = None, UniProtAC = None, UniProtID = None):
		if filepath:
			if not os.path.exists(filepath):
				raise Exception("The file %s does not exist." % filepath)
			filename = os.path.split(filepath)[-1]
			rootname, extension = os.path.splitext(filename)
			if not extension.lower() == ".pdb":
				raise Exception("Aborting: The file does not have a .pdb extension.")
		elif pdbID:
			rootname = pdbID
		try:
			Structure = ddgproject.PDBStructure(rootname, protein = protein, source = source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID)
			Structure.getPDBContents()
			
			sql = ("SELECT PDB_ID FROM Structure WHERE %s=" % dbfields.PDB_ID) + "%s"
			results = self.ddGDB.execute(sql, parameters = (rootname,))
			if results:
				ddgproject.getUniProtMapping(pdbID, storeInDatabase = True)
				raise Exception("There is already a structure in the database with the ID %s." % rootname)
			Structure.commit(self.ddGDB)
		except Exception, e:
			colortext.error(str(e))
			colortext.error(traceback.format_exc())
			raise Exception("An exception occurred committing %s to the database." % filepath)

	def createDummyExperiment(self, pdbID, mutationset, chains, sourceID, ddG, ExperimentSetName = "DummySource"):
		Experiment = ddgproject.ExperimentSet(pdbID, ExperimentSetName)
		for m in mutationset.mutations:
			Experiment.addMutation(m[0], m[1], m[2], m[3])
		for c in chains:
			Experiment.addChain(c)
		Experiment.addExperimentalScore(sourceID, ddG, pdbID)
		Experiment.commit(self.ddGDB)
			
			
	def addPrediction(self, experimentID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = False, Description = {}, InputFiles = {}):
		'''This function inserts a prediction into the database.
			The parameters define:
				the experiment we are running the prediction for;
				the name of the set of predictions for later grouping;
				the short description of the Command to be used for prediction;
				whether HETATM lines are to be kept or not.
			We strip the PDB based on the chains used for the experiment and KeepHETATMLines.
			We then add the prediction record, including the stripped PDB and the inverse mapping
			from Rosetta residue numbering to PDB residue numbering.''' 
			
		parameters = (experimentID,)
		
		try:
			sql = ("SELECT %(PDB_ID)s, %(Content)s FROM %(Experiment)s INNER JOIN %(Structure)s WHERE %(Experiment)s.%(Structure)s=%(PDB_ID)s AND %(Experiment)s.%(ID)s=" % dbfields) + "%s"
			results = self.ddGDB.execute(sql, parameters = parameters)
			if len(results) != 1:
				raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
			
			# Get the related PDB ID and file
			result = results[0]
			pdbID = result[dbfields.PDB_ID]
			contents = result[dbfields.Content]
			
			pdb = PDB(contents.split("\n"))
			
			# Check that the mutated positions exist and that the wild-type matches the PDB
			mutations = [result for result in self.ddGDB.callproc("GetMutations", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			checkPDBAgainstMutations(pdbID, pdb, mutations)
			
			# Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
			chains = [result[0] for result in self.ddGDB.callproc("GetChains", parameters = parameters, cursorClass = ddgproject.StdCursor)]
			pdb.stripForDDG(chains, KeepHETATMLines)
			
			# - Post stripping checks -
			# Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
			# then check again that the mutated positions exist and that the wild-type matches the PDB
			remappedMutations = pdb.remapMutations(mutations, pdbID)
			remappedMutations = [[m[0], ResidueID2String(m[1]), m[2], m[3]] for m in remappedMutations]
			
			resfile = self._createResfile(pdb, remappedMutations)
			
			# Check to make sure that we haven't stripped all the ATOM lines
			if not pdb.GetAllATOMLines():
				raise colortext.Exception("No ATOM lines remain in the stripped PDB file of %s." % pdbID)				
			
			# Check to make sure that CSE and MSE are not present in the PDB
			badresidues = pdb.CheckForPresenceOf(["CSE", "MSE"])
			if badresidues:
				raise colortext.Exception("Found residues [%s] in the stripped PDB file of %s. These should be changed to run this job under Rosetta." % (join(badresidues, ", "), pdbID))
			
			# Turn the lines array back into a valid PDB file				
			strippedPDB = join(pdb.lines, "\n")
		except Exception, e:
			colortext.error("\nError: '%s'.\n" % (str(e)))
			colortext.error(traceback.format_exc())
			raise colortext.Exception("An exception occurred retrieving the experimental data for Experiment ID #%s." % experimentID)
		
		InputFiles["RESFILE"] = resfile
		
		ExtraParameters = {}
		InputFiles = pickle.dumps(InputFiles)
		Description = pickle.dumps(Description)
		ExtraParameters = pickle.dumps(ExtraParameters)
		
		params = {
			dbfields.ExperimentID : experimentID,
			dbfields.PredictionSet : PredictionSet,
			dbfields.ProtocolID : ProtocolID,
			dbfields.KeptHETATMLines : KeepHETATMLines,
			dbfields.StrippedPDB : strippedPDB,
			dbfields.ResidueMapping : pickle.dumps(pdb.get_ddGInverseResmap()),
			dbfields.InputFiles : InputFiles,
			dbfields.Description : Description,
			dbfields.Status : dbfields.queued,
			dbfields.ExtraParameters : ExtraParameters,
			dbfields.StoreOutput : StoreOutput,
		}
		
		self.ddGDB.insertDict('Prediction', params)
		
		# Add cryptID string
		predictionID = self.ddGDB.getLastRowID()
		entryDate = self.ddGDB.execute("SELECT EntryDate FROM Prediction WHERE ID=%s", parameters = (predictionID,))[0]["EntryDate"]	
		rdmstring = join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16), '')
		cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
		cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
		entryDate = self.ddGDB.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))	
		