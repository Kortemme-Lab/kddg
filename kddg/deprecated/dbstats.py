# This file exists as a placeholder for old code used to generate statistics. Yes, it's bad practice but it should
# save time later. The data will be pulled from the database rather than the raw files and used to create graphs.

# THESE FUNCTIONS SHOULD BE MERGED WITH analysis.py AND monomer_analysis.py TO FORM AN ANALYSIS LAYER OF THE API

def getEmptyMutationMatrix(self):
	smallAminoAcids = set(['A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'])
	largeAminoAcids = set(["E", "F", "H", "I", "K", "L", "M", "Q", "R", "W", "Y"])
	allAminoAcids = smallAminoAcids.union(largeAminoAcids) 
	
	mutationMatrix = {}
	for aa1 in allAminoAcids:
		mutationMatrix[aa1] = {}
		for aa2 in allAminoAcids:
			mutationMatrix[aa1][aa2] = 0
	return mutationMatrix

def printMutationMatrix(self, mutationMatrix, divisor = 1.0):
	smallAminoAcids = set(['A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'])
	largeAminoAcids = set(["E", "F", "H", "I", "K", "L", "M", "Q", "R", "W", "Y"])
	allAminoAcids = smallAminoAcids.union(largeAminoAcids) 
	allAminoAcids = sorted(list(allAminoAcids))
	
	print(",%s" % join(allAminoAcids,","))
	for wt in allAminoAcids:
		print("%s,%s" % (wt, join(map(str, [float(mutationMatrix[wt][mt])/float(divisor) for mt in allAminoAcids]),",")))
	
def getMutationTypesFromOtherDataSets(self):
	PTids = set()
	
	smallAminoAcids = set(['A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'])
	largeAminoAcids = set(["E", "F", "H", "I", "K", "L", "M", "Q", "R", "W", "Y"])
	allAminoAcids = smallAminoAcids.union(largeAminoAcids) 
	
	print("Guerois")
	mutationMatrix = self.getEmptyMutationMatrix()
	mutationsToType = dict.fromkeys(list(allAminoAcids), 0)
	mutationsFromType = dict.fromkeys(list(allAminoAcids), 0)
	singleMutationsByType = {"s" : {"s" : 0, "l" : 0}, "l" : {"s" : 0, "l" : 0}}
	lines = rosettahelper.readFileLines("../rawdata/guerois/guerois-annotated.csv")
	assert(lines[0].split("\t")[7] == 'Mutations')
	numrecords = 0
	for line in lines:
		if line[0] != "#" and line.strip():
			mutation = line.split("\t")[7].strip()
			if mutation.find(",") == -1:
				assert(mutation[0] in allAminoAcids)
				assert(mutation[-1] in allAminoAcids)
				assert(mutation[1:-1].isdigit() or mutation[1:-2].isdigit())
				fromsize = None
				tosize = None
				if mutation[0] in smallAminoAcids:
					fromsize = "s"
				else:
					fromsize = "l"
				if mutation[-1] in smallAminoAcids:
					tosize = "s"
				else:
					tosize = "l"
				mutationsFromType[mutation[0]] += 1
				mutationsToType[mutation[-1]] += 1
				mutationMatrix[mutation[0]][mutation[-1]] += 1
				singleMutationsByType[fromsize][tosize] += 1
				numrecords += 1
	#print(singleMutationsByType)
	#print("mutationsFromType:")
	#for k, v in sorted(mutationsFromType.iteritems()):
	#	print("%s, %d" % (k, v))
	#print(sum(mutationsFromType.values()))
	#print("mutationsToType:")
	#for k, v in sorted(mutationsToType.iteritems()):
	#	print("%s, %d" % (k, v))
	#print(sum(mutationsToType.values()))
	self.printMutationMatrix(mutationMatrix)
	print("")
	self.printMutationMatrix(mutationMatrix, divisor = float(numrecords)/100.0)
	

	print("Liz")
	mutationMatrix = self.getEmptyMutationMatrix()
	mutationsToType = dict.fromkeys(list(allAminoAcids), 0)
	mutationsFromType = dict.fromkeys(list(allAminoAcids), 0)
	singleMutationsByType = {"s" : {"s" : 0, "l" : 0}, "l" : {"s" : 0, "l" : 0}}
	lines = rosettahelper.readFileLines("../rawdata/liz_kellogg/ProteinsPaper-annotated.csv")
	assert(lines[0].split("\t")[7] == 'WildTypeAA')
	assert(lines[0].split("\t")[9] == 'MutantAA')
	numrecords = 0
	for line in lines:
		if line[0] != "#" and line.strip():
			wt = line.split("\t")[7].strip()
			mt = line.split("\t")[9].strip()
			assert(wt in allAminoAcids)
			assert(mt in allAminoAcids)
			fromsize = None
			tosize = None
			if wt in smallAminoAcids:
				fromsize = "s"
			else:
				fromsize = "l"
			if mt in smallAminoAcids:
				tosize = "s"
			else:
				tosize = "l"
			mutationsFromType[wt] += 1
			mutationsToType[mt] += 1
			mutationMatrix[wt][mt] += 1
			singleMutationsByType[fromsize][tosize] += 1
			numrecords += 1
	#print(singleMutationsByType)
	#print("mutationsFromType:")
	#for k, v in sorted(mutationsFromType.iteritems()):
	#	print("%s, %d" % (k, v))
	#print(sum(mutationsFromType.values()))
	#print("mutationsToType:")
	#for k, v in sorted(mutationsToType.iteritems()):
	#	print("%s, %d" % (k, v))
	#print(sum(mutationsToType.values()))
	self.printMutationMatrix(mutationMatrix)
	print("")
	self.printMutationMatrix(mutationMatrix, divisor = float(numrecords)/100.0)
	
	print("Potapov")
	mutationMatrix = self.getEmptyMutationMatrix()
	mutationsToType = dict.fromkeys(list(allAminoAcids), 0)
	mutationsFromType = dict.fromkeys(list(allAminoAcids), 0)
	singleMutationsByType = {"s" : {"s" : 0, "l" : 0}, "l" : {"s" : 0, "l" : 0}}
	lines = rosettahelper.readFileLines("../rawdata/potapov/mutants-annotated.csv")
	assert(lines[0].split("\t")[6] == 'WildTypeAA')
	assert(lines[0].split("\t")[8] == 'MutantAA')
	numrecords = 0
	for line in lines:
		if line[0] != "#" and line.strip():
			wt = line.split("\t")[6].strip()
			mt = line.split("\t")[8].strip()
			if mt != "LA":
				assert(wt in allAminoAcids)
				assert(mt in allAminoAcids)
				fromsize = None
				tosize = None
				if wt in smallAminoAcids:
					fromsize = "s"
				else:
					fromsize = "l"
				if mt in smallAminoAcids:
					tosize = "s"
				else:
					tosize = "l"
				mutationsFromType[wt] += 1
				mutationsToType[mt] += 1
				mutationMatrix[wt][mt] += 1
				singleMutationsByType[fromsize][tosize] += 1
				numrecords += 1
	#print(singleMutationsByType)
	#print("mutationsFromType:")
	#for k, v in sorted(mutationsFromType.iteritems()):
	#	print("%s, %d" % (k, v))
	#print(sum(mutationsFromType.values()))
	#print("mutationsToType:")
	#for k, v in sorted(mutationsToType.iteritems()):
	#	print("%s, %d" % (k, v))
	#print(sum(mutationsToType.values()))
	self.printMutationMatrix(mutationMatrix)
	print("")
	self.printMutationMatrix(mutationMatrix, divisor = float(numrecords)/100.0)
	
	

def getDataForRosettaCon(self):
	# Get data to determine whether or not to store parsed data
	ddGdb = self.ddGdb
	FieldNames = ddGdb.FieldNames
	
	publicationSources = {}
	for r in ddGdb.execute('SELECT * FROM %s' % FieldNames.Source._name):
		publicationID = r[FieldNames.Source.ID]
		publicationSources[publicationID] = r
		publicationSources[publicationID]["DDGValueLocations"] = [location["Location"] for location in ddGdb.execute('SELECT * FROM SourceDDGValueLocation WHERE SourceID=%s', parameters=(publicationID,))]

	ExistingDBIDs = {}
	for r in ddGdb.execute('SELECT RecordNumber, DataSetDDGSource.ExperimentAssayID FROM DataSetDDG INNER JOIN DataSetDDGSource ON DataSetDDGID=DataSetDDG.ID INNER JOIN ExperimentAssayDDG ON DataSetDDGSource.ExperimentAssayID=ExperimentAssayDDG.ExperimentAssayID WHERE DataSetID="ProTherm_v25616_2011/12/21"', cursorClass = dbi.StdCursor):
		ExistingDBIDs[int(r[0])] = int(r[1])
	
	PublicationsToCheck = [r["SourceID"] for r in ddGdb.execute('SELECT SourceID FROM ProThermUnits WHERE DDGConvention IS NULL')] 
	PublicationsToCheckProThermRecords = dict.fromkeys(PublicationsToCheck)
	for p in PublicationsToCheckProThermRecords.keys():
		PublicationsToCheckProThermRecords[p] = []
	PublicationsToCheck = set(PublicationsToCheck)

	setOfSuccessfullyParsedIDs = set() 
	AllPDBIDs = {}
	ID = None
	mutation = {}
	
	experiments = {}
	count = 0
	recordcount = 0
	chains = {}
	
	mutationsByType = {}
	
	smallAminoAcids = set(['A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'])
	largeAminoAcids = set(["E", "F", "H", "I", "K", "L", "M", "Q", "R", "W", "Y"])
	allAminoAcids = smallAminoAcids.union(largeAminoAcids) 
	singleMutationsByType = {"s" : {"s" : 0, "l" : 0}, "l" : {"s" : 0, "l" : 0}}
	
	mutationsToType = dict.fromkeys(list(allAminoAcids), 0)
	mutationsFromType = dict.fromkeys(list(allAminoAcids), 0)
	
	# Variables for the experimental conditions and publications data
	PMIDlist = {}
	
	ptReader = ProThermReader(os.path.join("..", "rawdata", "ProTherm", "ProTherm25616.dat"))
	#ptReader.ExistingScores = ExistingScores
	if not ptReader.test():
		return

	ptReader.ExistingDBIDs = ExistingDBIDs
	
	check_just_these_cases = []
	check_these_cases = check_just_these_cases or ptReader.list_of_available_keys
	
	colortext.message("Parsing ProTherm")
	requiredProThermIDs = self.readProThermIDsFromOtherDataSets()
	
	newPublicationsToCheck = set()
	
	mutationMatrix = self.getEmptyMutationMatrix()
	nummatrixrecords = 0
	
	colortext.printf("|" + ("*" * (int(len(ptReader.list_of_available_keys)/1000) - 2)) + "|")
	for ID in check_these_cases:
		
		#Progress meter
		if ID % 1000 == 0:
			colortext.write(".", "green")
			colortext.flush()
		
		# Skip bad cases in ProTherm
		if ID in ptReader.skipTheseCases:
			continue
		
		record = ptReader.readRecord(ID)
		store = True
		
		if record["ddG"] == None and record["ddG_H2O"] == None:
			# "No ddG value for record %d. Skipping." % ID
			continue
		
		# *** Experiment records ***
		
		# Get PDB ID
		if record["PDB_wild"]:
			pdbID = record["PDB_wild"].upper()
		
		# Parse chain
		chainID = None
		if record["MUTATED_CHAIN"] and len(record["MUTATED_CHAIN"]) == 1:
			chainID = record["MUTATED_CHAIN"]
			assert(len(chainID) == 1)
			chains[chainID] = True
		else:
			colortext.error("Error processing chain: ID %d, %s" %  (ID, record["MUTATED_CHAIN"]))
			store = False
		
		# Parse mutant
		mutantlist = []
		if record["PDB_mutant"]:
			mutantlist = record["PDB_mutant"].split(",")
			for mutantID in mutantlist:
				mutantID = mutantID.strip()
				if not len(mutantID) == 4:
					raise Exception('Error parsing mutant list "%s" in record %d: ' % (mutantlist, ID))
				AllPDBIDs[mutantID] = True
						
		# Parse mutations
		mutations = None
		try:
			mutations = ptReader.getMutations(ID, record)
			if not mutations:
				store = False
		except Exception, e:
			colortext.error(str(e))
			colortext.error(traceback.format_exc())
			colortext.warning("An exception occurred parsing the mutation '%s' in record %d." % (record["MUTATION"], ID))
			continue
			raise Exception("An exception occurred parsing the mutation '%s' in record %d." % (record["MUTATION"], ID))
		
		# We have enough information to create the Experiment records
		if mutations and store:
			nummut = len(mutations)
			mutationsByType[nummut] = mutationsByType.get(nummut, {"ddG" : 0, "ddG_H2O" : 0})
			fromsize = None
			tosize = None
			if nummut == 1:
				if mutations[0].WildTypeAA in smallAminoAcids:
					fromsize = "s"
				else:
					fromsize = "l"
				if mutations[0].MutantAA in smallAminoAcids:
					tosize = "s"
				else:
					tosize = "l"
			if record["ddG"] != None:
				mutationsByType[nummut]["ddG"] += 1
				if nummut == 1:
					singleMutationsByType[fromsize][tosize] += 1
					mutationsFromType[mutations[0].WildTypeAA] += 1
					mutationsToType[mutations[0].MutantAA] += 1
					mutationMatrix[mutations[0].WildTypeAA][mutations[0].MutantAA] += 1
					nummatrixrecords += 1
				recordcount += 1
			if record["ddG_H2O"] != None:
				if nummut == 1:
					singleMutationsByType[fromsize][tosize] += 1
					mutationsFromType[mutations[0].WildTypeAA] += 1
					mutationsToType[mutations[0].MutantAA] += 1
					mutationMatrix[mutations[0].WildTypeAA][mutations[0].MutantAA] += 1
					nummatrixrecords += 1
				mutationsByType[nummut]["ddG_H2O"] += 1
				recordcount += 1
		# *** ExperimentAssay records ***
		# Parse references
		record["dbReferencePK"] = None
		try:
			referenceID = ptReader.getReference(ID, record)
			record["dbReferencePK"] = "PMID:%s" % referenceID
			if not referenceID:
				store = False
		except Exception, e:
			colortext.error(str(e))
			colortext.error(traceback.format_exc())
			raise Exception("An exception occurred parsing the reference '%s' in record %d." % (record["REFERENCE"], ID))
		
		# Parse ddG and ddG_H2O
		ddG = None
		ddG_H2O = None
		dbExperimentDDGAssay = None
		dbExperimentDDGH2OAssay = None
		assert(record.get("dbReferencePK"))
		if record["ddG"]:
			try:
				ddG = ptReader.getDDGInKcal(ID, record, useRosettaConvention = True)
			except Exception, e:
				colortext.error(str(e))
				colortext.error(traceback.format_exc())
				#ptReader.printRecord(ID)
				#raise Exception("An exception occurred parsing the ddG '%s' in record %d." % (record["ddG"], ID))
		if record["ddG_H2O"]:
			try:
				ddG_H2O = ptReader.getDDGH2OInKcal(ID, record, useRosettaConvention = True)
			except Exception, e:
				colortext.error(str(e))
				colortext.error(traceback.format_exc())
				#ptReader.printRecord(ID)
				#raise Exception("An exception occurred parsing the ddG_H2O '%s' in record %d." % (record["ddG_H2O"], ID))
			#if ExistingScores.get(ID) and (abs(ddG - ExistingScores[ID]) > 0.000001):
			#	colortext.error("ProTherm record %s exists as ExperimentScore entry %s but the values disagree (%s and %s respectively)." % (ID, ExistingDBIDs[ID], ddG, ExistingScores[ID]))
		
		if ddG == None and ddG_H2O == None:
			continue
	
	for k, v in sorted(mutationsByType.iteritems()):
		print("#mutations: %d" % k)
		print("#ddG: %d" % v["ddG"])
		print("#ddG_H2O: %d" % v["ddG_H2O"])
		print("")
	
	print("ProTherm:")
	print(recordcount)
	print(singleMutationsByType)
	print("mutationsFromType:")
	for k, v in sorted(mutationsFromType.iteritems()):
		print("%s, %d" % (k, v))
	print(sum(mutationsFromType.values()))
	print("mutationsToType:")
	for k, v in sorted(mutationsToType.iteritems()):
		print("%s, %d" % (k, v))
	print(sum(mutationsToType.values()))
	
	print("printMutationMatrix")
	self.printMutationMatrix(mutationMatrix)
	self.printMutationMatrix(mutationMatrix, divisor = float(nummatrixrecords)/100.0)
	
class ExperimentSet(DBObject):
	
	def __init__(self, ddGdb, pdbid, source, interface = None):
		#todo: delete
		raise Exception("Out of date.")
		self.ddGdb = ddGdb
		FieldNames_ = ddGdb.FlatFieldNames
		self.dict = {
			FieldNames_.Structure	: pdbid,
			FieldNames_.Source		: source,
			FieldNames_.ScoreVariance: None,
			"Interface"				: interface,
			"Mutants"				: {},
			"Mutations"				: [],
			"ExperimentChains"		: [],
			"ExperimentScores"		: [],
			"StdDeviation"			: None,
			"WithinStdDeviation"	: None
		}
	
	def addMutant(self, mutant):
		self.dict["Mutants"][mutant] = True

	def addMutation(self, chainID, residueID, wildtypeAA, mutantAA, ID = None, SecondaryStructureLocation=None):
		errors = []
		residueID = ("%s" % residueID).strip()
		if not chainID in AllowedChainLetters:
			errors.append("The chain '%s' is invalid." % chainID)
		if not wildtypeAA in AllowedAminoAcids:
			errors.append("The wildtype amino acid '%s' is invalid." % wildtypeAA)
		if not mutantAA in AllowedAminoAcids:
			errors.append("The mutant amino acid '%s' is invalid." % mutantAA)
		if not residueID.isdigit():
			if not residueID[0:-1].isdigit():
				errors.append("The residue '%s' is invalid." % residueID)
			elif residueID[-1] not in AllowedInsertionCodes:
				errors.append("The insertion code '%s' of residue '%s' is invalid." % (residue[-1], residueID))
		if errors:
			ID = ID or ""
			if ID:
				ID = ", ID %s" % ID
			errors = join(['\t%s\n' % e for e in errors], "")
			raise Exception("An exception occurred processing a mutation in the dataset %s%s.\n%s" % (self.dict[FieldNames_.Source], ID, errors))
		self.dict["Mutations"].append({
			FieldNames_.Chain 		: chainID,
			FieldNames_.ResidueID	: residueID,
			FieldNames_.WildTypeAA	: wildtypeAA,
			FieldNames_.MutantAA	: mutantAA
			})
	
	def addChain(self, chainID, ID = ""):
		if not chainID in AllowedChainLetters:
			raise Exception("An exception occurred processing a chain in the dataset %s%s.\n\tThe chain '%s' is invalid." % (self.dict[FieldNames_.Source], ID, errors, chainID))
		self.dict["ExperimentChains"].append(chainID)
	
	def getChains(self):
		return self.dict["ExperimentChains"]
	
	def setMutantIfUnset(self, mutant):
		if not self.dict[FieldNames_.Mutant]:
			self.dict[FieldNames_.Mutant] = mutant
	
	def addExperimentalScore(self, sourceID, ddG, pdbID, numMeasurements = 1):
		if pdbID != self.dict[FieldNames_.Structure]:
			raise colortext.Exception("Adding experimental score related to PDB structure %s to an experiment whose PDB structure should be %s." % (pdbID, self.dict[FieldNames_.Structure]))
		self.dict["ExperimentScores"].append({
			FieldNames_.SourceID				: sourceID,
			FieldNames_.ddG						: ddG,
			FieldNames_.NumberOfMeasurements	: numMeasurements
			})
	
	def mergeScores(self, maxStdDeviation = 1.0):
		d = self.dict
		
		n = len(d["ExperimentScores"])
		if n > 1:
			n = float(n)
			sum = 0
			for experimentalResult in d["ExperimentScores"]:
				if experimentalResult[FieldNames_.NumberOfMeasurements] != 1:
					raise Exception("Cannot merge scores when individual scores are from more than one measurement. Need to add logic to do proper weighting.")
				sum += experimentalResult[FieldNames_.ddG]
			mean = sum / n
			squaredsum = 0
			for experimentalResult in d["ExperimentScores"]:
				diff = (experimentalResult[FieldNames_.ddG] - mean)
				squaredsum += diff * diff
			variance = squaredsum / n
			d[FieldNames_.ScoreVariance] = variance
			stddev = sqrt(variance)
			d["StdDeviation"] = stddev 
			d["WithinStdDeviation"] = stddev <= maxStdDeviation
		else:
			d[FieldNames_.ScoreVariance] = 0
			d["WithinStdDeviation"] = True	
	
	def isEligible(self):
		d = self.dict
		if d["WithinStdDeviation"] == None:
			raise Exception("Standard deviation not yet computed.")
		else:
			return d["WithinStdDeviation"]
	
	def __repr__(self):
		raise Exception('''This is unlikely to work as I have not tested it in a while. In particular, ddG is not a string anymore.''')
		d = self.dict
		str = []
		str.append("%s: %s" % (FieldNames_.Structure, d[FieldNames_.Structure]))
		str.append("%ss: %s" % (FieldNames_.Mutant, join(d["Mutants"].keys(), ', ')))
		str.append("%s: %s" % (FieldNames_.Source, d[FieldNames_.Source]))
		str.append("Chains: %s" % (join([chain for chain in d["ExperimentChains"]], ", ")))
		str.append("Mutations:")
		for mutation in d["Mutations"]:
			str.append("\t%s%s: %s -> %s" % (mutation[FieldNames_.Chain], mutation[FieldNames_.ResidueID], mutation[FieldNames_.WildTypeAA], mutation[FieldNames_.MutantAA]))
		str.append("Experimental Scores:")
		for score in d["ExperimentScores"]:
			n = score[FieldNames_.NumberOfMeasurements]
			if n > 1:
				str.append("\t%s\t%0.2f (%d measurements)" % (score[FieldNames_.SourceID], score[FieldNames_.ddG], score[FieldNames_.NumberOfMeasurements]))
			else:
				str.append("\t%s\t%0.2f" % (score[FieldNames_.SourceID], score[FieldNames_.ddG]))
		return join(str, "\n")
			
	def commit(self, testonly = False, pdbPath = None, mutationAllowedToBeStoredDespiteMissingCoordinates = False):
		'''Commits the set of experiments associated with the mutation to the database. Returns the unique ID of the associated Experiment record.'''
		d = self.dict
		failed = False
		
		for score in d["ExperimentScores"]:
			scoresPresent = True
			results = self.ddGdb.locked_execute("SELECT Source, SourceID FROM Experiment INNER JOIN ExperimentScore ON Experiment.ID=ExperimentID WHERE Source=%s AND SourceID=%s", parameters = (d[FieldNames_.Source], score[FieldNames_.SourceID]))
			if results:
				return
		
		if len(d["ExperimentScores"]) == 0:
			raise Exception("This experiment has no associated scores.")
		
		if not d[FieldNames_.ScoreVariance]:
			self.mergeScores()
		
		if d["Mutants"]:
			for mutant in d["Mutants"].keys():
				MutantStructure = PDBStructure(self.ddGdb, mutant)
				results = self.ddGdb.execute("SELECT PDB_ID FROM Structure WHERE PDB_ID = %s", parameters = (MutantStructure.dict[self.ddGdb.FlatFieldNames.PDB_ID])) 
				if not results:
					MutantStructure.commit()
		
		# Sanity check that the chain information is correct (ProTherm has issues)
		pdbID = d[FieldNames_.Structure] 
		associatedRecords = sorted([score[FieldNames_.SourceID] for score in d["ExperimentScores"]])
		associatedRecordsStr = "%s (records: %s)" % (d[FieldNames_.Source], join(map(str, sorted([score[FieldNames_.SourceID] for score in d["ExperimentScores"]])),", "))
		chainsInPDB = PDBChains.get(d[FieldNames_.Structure])
		if not chainsInPDB:
			raise Exception("The chains for %s were not read in properly." % associatedRecordsStr)
		for c in self.dict["ExperimentChains"]:
			if not c in chainsInPDB:
				if len(chainsInPDB) == 1 and len(self.dict["ExperimentChains"]) == 1:
					colortext.warning("%s: Chain '%s' of %s does not exist in the PDB %s. Chain %s exists. Use that chain instead." % (pdbID, c, associatedRecordsStr, pdbID, chainsInPDB[0]))
					self.ddGdb.addChainWarning(pdbID, associatedRecords, c)
					failed = True
				else:
					self.ddGdb.addChainError(pdbID, c)
					raise colortext.Exception("Error committing experiment:\n%s: Chain '%s' of %s does not exist in the PDB %s. Chains %s exist." % (pdbID, c, associatedRecordsStr, pdbID, join(chainsInPDB, ", ")))
				
		# Sanity check that the wildtypes of all mutations are correct
		if pdbPath:
			WildTypeStructure = PDBStructure(self.ddGdb, pdbID, filepath = os.path.join(pdbPath, "%s.pdb" % pdbID))
		else:
			WildTypeStructure = PDBStructure(self.ddGdb, pdbID)
		contents = WildTypeStructure.getPDBContents()
		pdb = PDB(contents.split("\n"))
		
		badResidues = ["CSE", "MSE"]
		foundRes = pdb.CheckForPresenceOf(badResidues)
		if foundRes:
			colortext.warning("The PDB %s contains residues which could affect computation (%s)." % (pdbID, join(foundRes, ", ")))
			failed = True
			for res in foundRes:
				colortext.warning("The PDB %s contains %s. Check." % (pdbID, res))
		for mutation in d["Mutations"]:
			foundMatch = False
			for resid, wtaa in sorted(pdb.ProperResidueIDToAAMap().iteritems()):
				c = resid[0]
				resnum = resid[1:].strip()
				if mutation[FieldNames_.Chain] == c and mutation[FieldNames_.ResidueID] == resnum and mutation[FieldNames_.WildTypeAA] == wtaa:
					foundMatch = True
			if not foundMatch and not(mutationAllowedToBeStoredDespiteMissingCoordinates):
				#raise colortext.Exception("%s: Could not find a match for mutation %s %s:%s -> %s in %s." % (pdbID, mutation[FieldNames_.Chain], mutation[FieldNames_.ResidueID], mutation[FieldNames_.WildTypeAA], mutation[FieldNames_.MutantAA], associatedRecordsStr ))
				colortext.error("%s: Could not find a match for mutation %s %s:%s -> %s in %s." % (pdbID, mutation[FieldNames_.Chain], mutation[FieldNames_.ResidueID], mutation[FieldNames_.WildTypeAA], mutation[FieldNames_.MutantAA], associatedRecordsStr ))
				failed = True
				

				#raise Exception(colortext.make_error("%s: Could not find a match for mutation %s %s:%s -> %s in %s." % (pdbID, mutation[FieldNames_.Chain], mutation[FieldNames_.ResidueID], mutation[FieldNames_.WildTypeAA], mutation[FieldNames_.MutantAA], associatedRecordsStr )))
				
		# To disable adding new experiments:	return here
		if failed:
			return False
		
		SQL = 'INSERT INTO Experiment (Structure, Source) VALUES (%s, %s);'
		vals = (d[FieldNames_.Structure], d[FieldNames_.Source]) 
		#print(SQL % vals)
		if not testonly:
			self.ddGdb.locked_execute(SQL, parameters = vals)
			self.databaseID = self.ddGdb.getLastRowID()
			ExperimentID = self.databaseID
			#print(ExperimentID)
		else:
			ExperimentID = None 
			
		for chain in d["ExperimentChains"]:
			SQL = 'INSERT INTO ExperimentChain (ExperimentID, Chain) VALUES (%s, %s);'
			vals = (ExperimentID, chain) 
			#print(SQL % vals)
			if not testonly:
				self.ddGdb.locked_execute(SQL, parameters = vals)
		
		interface = d["Interface"]
		if interface:
			SQL = 'INSERT INTO ExperimentInterface (ExperimentID, Interface) VALUES (%s, %s);'
			vals = (ExperimentID, interface) 
			#print(SQL % vals)
			if not testonly:
				self.ddGdb.locked_execute(SQL, parameters = vals)
		
		for mutant in d["Mutants"].keys():
			SQL = 'INSERT INTO ExperimentMutant (ExperimentID, Mutant) VALUES (%s, %s);'
			vals = (ExperimentID, mutant) 
			#print(SQL % vals)
			if not testonly:
				self.ddGdb.locked_execute(SQL, parameters = vals)
		
		for mutation in d["Mutations"]:
			SQL = 'INSERT INTO ExperimentMutation (ExperimentID, Chain, ResidueID, WildTypeAA, MutantAA) VALUES (%s, %s, %s, %s, %s);'
			vals = (ExperimentID, mutation[FieldNames_.Chain], mutation[FieldNames_.ResidueID], mutation[FieldNames_.WildTypeAA], mutation[FieldNames_.MutantAA]) 
			#print(SQL % vals)
			if not testonly:
				self.ddGdb.locked_execute(SQL, parameters = vals)
			
		for score in d["ExperimentScores"]:
			SQL = 'INSERT INTO ExperimentScore (ExperimentID, SourceID, ddG, NumberOfMeasurements) VALUES (%s, %s, %s, %s);'
			vals = (ExperimentID, score[FieldNames_.SourceID], score[FieldNames_.ddG], score[FieldNames_.NumberOfMeasurements]) 
			#print(SQL % vals)
			if not testonly:
				self.ddGdb.locked_execute(SQL, parameters = vals)
		
		if not testonly:
			return self.databaseID
		else:
			return None
		