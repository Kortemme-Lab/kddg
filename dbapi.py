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
import ddgdbapi
from tools.pdb import PDB, ResidueID2String, checkPDBAgainstMutations, aa1
#from Bio.PDB import *
from tools.fs.io import write_file
from tools import colortext
import traceback
import pickle
import md5
import random
import score 
#import analysis
from ddgfilters import PredictionResultSet, ExperimentResultSet, StructureResultSet

#todo: dbfields = ddgdbapi.FieldNames()

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
        self.ddGDB = ddgdbapi.ddGDatabase()
        self.ddGDataDB = ddgdbapi.ddGPredictionDataDatabase()

    def __del__(self):
        pass
        #self.ddGDB.close()
        #self.ddGDataDB.close()

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

    def _createMutfile(self, pdb, mutations):
        '''The mutations here are in the original PDB numbering. pdb is assumed to use Rosetta numbering.
            We use the pdb mapping from PDB numbering to Rosetta numbering to generate the mutfile.
        '''
        mutfile = []
        for mutation in mutations:
            chain = mutation[0]
            resid = mutation[1]
            wt = mutation[2]
            mt = mutation[3]

            # Check that the expected wildtype exists in the PDB
            readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
            assert(wt == aa1[readwt])
            resid = resid.strip()
            mutfile.append("%(wt)s %(resid)s %(mt)s" % vars())
        if mutfile:
            mutfile = ["total %d" % len(mutations), "%d" % len(mutations)] + mutfile
            return join(mutfile, "\n")
        else:
            raise Exception("An error occurred creating a mutfile for the ddG job.")

    def getData(self, predictionID):
        results = self.ddGDataDB.execute_select("SELECT * FROM PredictionData WHERE ID=%s", parameters = (predictionID,))
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

                experimentsets = [e[0] for e in self.ddGDB.execute_select("SELECT DISTINCT Source FROM Experiment WHERE ID IN (%s)" % join(map(str, list(experiments.IDs)), ","), cursorClass = ddgdbapi.StdCursor)]

                if experimentsets:
                    colortext.printf("\nRelated publications for experiment-set sources:", "lightgreen")
                    for id in sorted(experimentsets):
                        print(id)
                        pubs = self.ddGDB.execute_select("SELECT ID, Type FROM SourceLocation WHERE SourceID=%s", parameters = (id,))
                        for pub in pubs:
                            print("\t%s: %s" % (pub["Type"], pub["ID"]))
        else:
            raise Exception("Empty result set.")



    def dumpData(self, outfile, predictionID):
        write_file(outfile, self.getData(predictionID))

    def analyze(self, prediction_result_set, outpath = os.getcwd()):
        PredictionIDs = sorted(list(prediction_result_set.getFilteredIDs()))
        colortext.printf("Analyzing %d records:" % len(PredictionIDs), "lightgreen")
        #results = self.ddGDB.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % join(map(str, PredictionIDs), ","))

        #for r in results:
        #	r["ddG"] = pickle.loads(r["ddG"])
        #	predicted_score = r["ddG"]["data"]["ddG"]
        #	experimental_scores = [expscore["ddG"] for expscore in self.ddGDB.callproc("GetScores", parameters = r["ExperimentID"])]
        #	mean_experimental_score = float(sum(experimental_scores)) / float(len(experimental_scores))

        results = self.ddGDB.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % join(map(str, PredictionIDs), ","))

        analysis.plot(analysis._R_mean_unsigned_error, analysis._createMAEFile, results, "my_plot1.pdf", average_fn = analysis._mean)
        analysis.plot(analysis._R_correlation_coefficient, analysis._createAveragedInputFile, results, "my_plot2.pdf", average_fn = analysis._mean)
        colortext.printf("Done", "lightgreen")



        #score.ddgTestScore

    def addPDBtoDatabase(self, filepath = None, pdbID = None, protein = None, source = None, UniProtAC = None, UniProtID = None, testonly = False):
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
            Structure = ddgdbapi.PDBStructure(rootname, protein = protein, source = source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID, testonly = testonly)
            #Structure.getPDBContents(self.ddGDB)
            sql = ("SELECT PDB_ID FROM Structure WHERE %s=" % dbfields.PDB_ID) + "%s"
            results = self.ddGDB.execute_select(sql, parameters = (rootname,))
            if results:
                #ddgdbapi.getUniProtMapping(pdbID, storeInDatabase = True)
                raise Exception("There is already a structure in the database with the ID %s." % rootname)
            Structure.commit(self.ddGDB, testonly = testonly)
        except Exception, e:
            colortext.error(str(e))
            colortext.error(traceback.format_exc())
            raise Exception("An exception occurred committing %s to the database." % filepath)

    def createDummyExperiment(self, pdbID, mutationset, chains, sourceID, ddG, ExperimentSetName = "DummySource"):
        #todo
        raise Exception("Out of date function.")
        Experiment = ddgdbapi.ExperimentSet(pdbID, ExperimentSetName)
        for m in mutationset.mutations:
            Experiment.addMutation(m[0], m[1], m[2], m[3])
        for c in chains:
            Experiment.addChain(c)
        Experiment.addExperimentalScore(sourceID, ddG, pdbID)
        Experiment.commit(self.ddGDB)

    def charge_PredictionSet_by_number_of_residues(self, PredictionSet):
        '''This function assigns a cost for a prediction equal to the number of residues in the chains.'''
        from tools.bio.rcsb import parseFASTAs

        ddGdb = self.ddGDB
        predictions = ddGdb.execute_select("SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet=%s", parameters=(PredictionSet,))

        PDB_chain_lengths ={}
        for prediction in predictions:
            chain_records = ddGdb.execute_select('SELECT PDBFileID, Chain FROM Experiment INNER JOIN ExperimentChain ON ExperimentID=Experiment.ID WHERE ExperimentID=%s', parameters=(prediction['ExperimentID']))
            num_residues = 0
            for chain_record in chain_records:
                key = (chain_record['PDBFileID'], chain_record['Chain'])

                if PDB_chain_lengths.get(key) == None:
                    fasta = ddGdb.execute_select("SELECT FASTA FROM PDBFile WHERE ID=%s", parameters = (chain_record['PDBFileID'],))
                    assert(len(fasta) == 1)
                    fasta = fasta[0]['FASTA']
                    f = parseFASTAs(fasta)
                    PDB_chain_lengths[key] = len(f[chain_record['PDBFileID']][chain_record['Chain']])
                chain_length = PDB_chain_lengths[key]
                num_residues += chain_length

            print("UPDATE Prediction SET Cost=%0.2f WHERE ID=%d" % (num_residues, prediction['ID']))

            predictions = ddGdb.execute("UPDATE Prediction SET Cost=%s WHERE ID=%s", parameters=(num_residues, prediction['ID'],))

    def create_PredictionSet(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40):
        Status = 'halted'
        if not halted:
            Status = 'active'
        d = {
            'ID'        : PredictionSetID,
            'Status'    : Status,
            'Priority'  : Priority,
            'BatchSize' : BatchSize,
        }
        self.ddGDB.insertDictIfNew("PredictionSet", d, ['ID'])

    def createPredictionsFromUserDataSet(self, userdatasetTextID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = False, Description = {}, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False):

        assert(self.ddGDB.execute_select("SELECT ID FROM PredictionSet WHERE ID=%s", parameters=(PredictionSet,)))

        #results = self.ddGDB.execute_select("SELECT * FROM UserDataSet WHERE TextID=%s", parameters=(userdatasetTextID,))
        results = self.ddGDB.execute_select("SELECT UserDataSetExperiment.* FROM UserDataSetExperiment INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%s", parameters=(userdatasetTextID,))
        if not results:
            return False

        if not(quiet):
            colortext.message("Creating predictions for UserDataSet %s using protocol %s" % (userdatasetTextID, ProtocolID))
            colortext.message("%d records found in the UserDataSet" % len(results))

        count = 0
        showprogress = not(quiet) and len(results) > 300
        if showprogress:
            print("|" + ("*" * (int(len(results)/100)-2)) + "|")
        for r in results:
            self.addPrediction(r["ExperimentID"], r["ID"], PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = r["PDBFileID"], StoreOutput = StoreOutput, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = testonly)
            count += 1
            if showprogress:
                if count > 100:
                    colortext.write(".", "cyan", flush = True)
                    count = 0
            if shortrun and count > 4:
                break
        print("")
        return(True)

    def addPrediction(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False):
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
        assert(ReverseMutation == False) # todo: allow this later
        try:
            predictionPDB_ID = None

            sql = "SELECT PDBFileID, Content FROM Experiment INNER JOIN PDBFile WHERE Experiment.PDBFileID=PDBFile.ID AND Experiment.ID=%s"
            results = self.ddGDB.execute_select(sql, parameters = parameters)
            if len(results) != 1:
                raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
            experimentPDB_ID = results[0]["PDBFileID"]

            if PDB_ID:
                #sql = "SELECT ID, Content FROM PDBFile WHERE ID=%s"
                results = self.ddGDB.execute_select("SELECT ID, Content FROM PDBFile WHERE ID=%s", parameters=(PDB_ID))
                if len(results) != 1:
                    raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
                predictionPDB_ID = results[0]["ID"]
            else:
                predictionPDB_ID = experimentPDB_ID

            # Get the related PDB ID and file
            assert(len(results) == 1)
            result = results[0]
            pdbID = result["ID"]
            contents = result["Content"]

            pdb = PDB(contents.split("\n"))

            # Check that the mutated positions exist and that the wild-type matches the PDB
            mutations = self.ddGDB.call_select_proc("GetMutations", parameters = parameters)

            # todo: Hack. This should be removed when PDB homologs are dealt with properly.
            for mutation in mutations:
                if experimentPDB_ID == "1AJ3" and predictionPDB_ID == "1U5P":
                    assert(int(mutation['ResidueID']) < 1000)
                    mutation['ResidueID'] = str(int(mutation['ResidueID']) + 1762)

            checkPDBAgainstMutations(pdbID, pdb, mutations)

            # Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
            chains = [result['Chain'] for result in self.ddGDB.call_select_proc("GetChains", parameters = parameters)]
            pdb.stripForDDG(chains, KeepHETATMLines, numberOfModels = 1)

            # - Post stripping checks -
            # Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
            # then check again that the mutated positions exist and that the wild-type matches the PDB
            remappedMutations = pdb.remapMutations(mutations, pdbID)
            remappedMutations = [[m[0], ResidueID2String(m[1]), m[2], m[3]] for m in remappedMutations]

            #resfile = self._createResfile(pdb, remappedMutations)
            mutfile = self._createMutfile(pdb, remappedMutations)

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
            colortext.error("Error in %s, %s: " % (experimentID, UserDataSetExperimentID))
            colortext.warning(str(e))
            return
            colortext.error("\nError: '%s'.\n" % (str(e)))
            colortext.error(traceback.format_exc())
            raise colortext.Exception("An exception occurred retrieving the experimental data for Experiment ID #%s." % experimentID)

        #InputFiles["RESFILE"] = resfile
        InputFiles["MUTFILE"] = mutfile

        ExtraParameters = {}
        InputFiles = pickle.dumps(InputFiles)
        Description = pickle.dumps(Description)
        ExtraParameters = pickle.dumps(ExtraParameters)

        PredictionFieldNames = self.ddGDB.FieldNames.Prediction
        params = {
            PredictionFieldNames.ExperimentID		: experimentID,
            PredictionFieldNames.UserDataSetExperimentID : UserDataSetExperimentID,
            PredictionFieldNames.PredictionSet		: PredictionSet,
            PredictionFieldNames.ProtocolID			: ProtocolID,
            PredictionFieldNames.KeptHETATMLines	: KeepHETATMLines,
            PredictionFieldNames.StrippedPDB		: strippedPDB,
            PredictionFieldNames.ResidueMapping		: pickle.dumps(pdb.get_ddGInverseResmap()),
            PredictionFieldNames.InputFiles			: InputFiles,
            PredictionFieldNames.Description		: Description,
            PredictionFieldNames.Status 			: "queued",
            PredictionFieldNames.ExtraParameters	: ExtraParameters,
            PredictionFieldNames.StoreOutput		: StoreOutput,
        }
        if not testonly:
            self.ddGDB.insertDict('Prediction', params)

            # Add cryptID string
            predictionID = self.ddGDB.getLastRowID()
            entryDate = self.ddGDB.execute_select("SELECT EntryDate FROM Prediction WHERE ID=%s", parameters = (predictionID,))[0]["EntryDate"]
            rdmstring = join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16), '')
            cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
            cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
            entryDate = self.ddGDB.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))
