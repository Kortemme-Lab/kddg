# This file exists as a placeholder for old code used to generate statistics. Yes, it's bad practice but it should
# save time later. The data will be pulled from the database rather than the raw files and used to create graphs.

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
	for r in ddGdb.execute('SELECT RecordNumber, DataSetDDGSource.ExperimentAssayID FROM DataSetDDG INNER JOIN DataSetDDGSource ON DataSetDDGID=DataSetDDG.ID INNER JOIN ExperimentAssayDDG ON DataSetDDGSource.ExperimentAssayID=ExperimentAssayDDG.ExperimentAssayID WHERE DataSetID="ProTherm_v25616_2011/12/21"', cursorClass = ddgdbapi.StdCursor):
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
		