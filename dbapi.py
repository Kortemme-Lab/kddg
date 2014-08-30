#!/usr/bin/python2.4
# encoding: utf-8
"""
dbapi.py
High-level functions for interacting with the ddG database.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import os
import string
import shutil
import glob
import traceback
import pickle
import md5
import random
import datetime
import zipfile

try:
    import matplotlib
    # A non-interactive backend to generate PNGs. matplotlib.use('PS') is used for PS files. If used, this command must be run before importing matplotlib.pyplot.
    matplotlib.use("AGG")
    import matplotlib.pyplot as plt
    import textwrap
except ImportError:
    plt=None


import score
import ddgdbapi
from tools.bio.pdb import PDB
from tools.bio.basics import residue_type_3to1_map as aa1
from tools.bio.basics import Mutation
from tools.bio.alignment import ScaffoldModelChainMapper
#from Bio.PDB import *
from tools.fs.io import write_file, read_file
from tools.process import Popen
from tools import colortext
#import analysis
#from ddgfilters import PredictionResultSet, ExperimentResultSet, StructureResultSet

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

    def __init__(self, passwd = None, username = 'kortemmelab'):
        self.ddGDB = ddgdbapi.ddGDatabase(passwd = passwd, username = username)
        self.ddGDataDB = ddgdbapi.ddGPredictionDataDatabase(passwd = passwd, username = username)

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
            return '\n'.join(resfile)
        else:
            raise Exception("An error occurred creating a resfile for the ddG job.")

    def _createMutfile(self, pdb, mutations):
        '''The mutations here are in the original PDB numbering. pdb is assumed to use Rosetta numbering.
            We use the pdb mapping from PDB numbering to Rosetta numbering to generate the mutfile.
        '''
        mutfile = []

        for mutation in mutations:
            chain = mutation.Chain
            resid = PDB.ResidueID2String(mutation.ResidueID)
            wt = mutation.WildTypeAA
            mt = mutation.MutantAA

            # Check that the expected wildtype exists in the PDB
            readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
            assert(wt == aa1[readwt])
            resid = resid.strip()
            mutfile.append("%(wt)s %(resid)s %(mt)s" % vars())
        if mutfile:
            mutfile = ["total %d" % len(mutations), "%d" % len(mutations)] + mutfile
            return '\n'.join(mutfile)
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

                experimentsets = [e[0] for e in self.ddGDB.execute_select("SELECT DISTINCT Source FROM Experiment WHERE ID IN (%s)" % ','.join(map(str, list(experiments.IDs))), cursorClass = ddgdbapi.StdCursor)]

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

        results = self.ddGDB.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % ','.join(map(str, PredictionIDs)))

        analysis.plot(analysis._R_mean_unsigned_error, analysis._createMAEFile, results, "my_plot1.pdf", average_fn = analysis._mean)
        analysis.plot(analysis._R_correlation_coefficient, analysis._createAveragedInputFile, results, "my_plot2.pdf", average_fn = analysis._mean)
        colortext.printf("Done", "lightgreen")



        #score.ddgTestScore

    def addPDBtoDatabase(self, filepath = None, pdbID = None, protein = None, source = None, UniProtAC = None, UniProtID = None, testonly = False):
        #todo: use either this or add_pdb_file but not both
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
            Structure = ddgdbapi.PDBStructure(self.ddGDB, rootname, protein = protein, source = source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID, testonly = testonly)
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

    def add_pdb_file(self, filepath, pdb_id):
        #todo: use either this or addPDBtoDatabase but not both

        existing_pdb = self.ddGDB.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
        if not existing_pdb:
            pdb_contents = read_file(filepath)
            p = PDB(pdb_contents)

            fasta = []
            for c, sequence in p.atom_sequences.iteritems():
                fasta.append('>%s:%s|PDBID|CHAIN|SEQUENCE' % (pdb_id.replace(':', '_'), c))
                fasta.append(str(sequence))
            fasta = '\n'.join(fasta)

            d = {
                'ID' : pdb_id,
                'FileSource' : 'Biosensor project',
                'Content' : read_file(filepath),
                'FASTA' : fasta,
                'Resolution' : None,
                'Techniques' : 'Rosetta model',
                'BFactors' : '',
                'Publication' : None
            }
            self.ddGDB.insertDictIfNew('PDBFile', d, ['ID'])

    def createDummyExperiment(self, pdbID, mutationset, chains, sourceID, ddG, ExperimentSetName = "DummySource"):
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        raise Exception("Out of date function.")
        Experiment = ddgdbapi.ExperimentSet(pdbID, ExperimentSetName)
        for m in mutationset.mutations:
            Experiment.addMutation(m[0], m[1], m[2], m[3])
        for c in chains:
            Experiment.addChain(c)
        Experiment.addExperimentalScore(sourceID, ddG, pdbID)
        Experiment.commit(self.ddGDB)

    def createDummyExperiment_ankyrin_repeat(self, pdbID, mutations, chain):
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        experiment = ddgdbapi.ExperimentDefinition(self.ddGDB, pdbID, interface = None)
        experiment.addChain(chain)
        for m in mutations:
            experiment.addMutation(m)
        experiment.commit(False)

    def add_mutant(self, pdb_ID, mutant_mutations):
        '''Use this function to add one set of mutations ON THE SAME CHAIN (i.e. corresponding to one mutant) to the database.
           todo: generalize this to allow different chains
        '''
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        chains = set([m.Chain for m in mutant_mutations])
        assert(len(chains) == 1)
        colortext.warning("Adding mutation: %s." % ', '.join(map(str, mutant_mutations)))
        self.createDummyExperiment_ankyrin_repeat(pdb_ID, mutant_mutations, chains.pop())


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
        if halted:
            Status = 'halted'
        else:
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

    def add_predictions_by_pdb_id(self, pdb_ID, PredictionSet, ProtocolID, status = 'active', priority = 5, KeepHETATMLines = False, strip_other_chains = True):
        ''' This function adds predictions for all Experiments corresponding to pdb_ID to the specified prediction set.
            This is useful for custom runs e.g. when we are using the DDG scheduler for design rather than for benchmarking.
        '''
        colortext.printf("\nAdding any mutations for this structure which have not been queued/run in the %s prediction set." % PredictionSet, "lightgreen")

        d = {
            'ID' : PredictionSet,
            'Status' : status,
            'Priority' : 9,
            'BatchSize' : 40,
            'EntryDate' : datetime.datetime.now(),
        }
        self.ddGDB.insertDictIfNew('PredictionSet', d, ['ID'])

        # Update the priority and activity if necessary
        self.ddGDB.execute('UPDATE PredictionSet SET Status=%s AND Priority=%s WHERE ID=%s', parameters = (status, priority, PredictionSet))

        # Determine the set of experiments to add
        ExperimentIDs = set([r['ID'] for r in self.ddGDB.execute_select('SELECT ID FROM Experiment WHERE PDBFileID=%s', parameters=(pdb_ID,))])
        ExperimentIDsInPredictionSet = set([r['ExperimentID'] for r in self.ddGDB.execute_select('SELECT ExperimentID FROM Prediction WHERE PredictionSet=%s', parameters=(PredictionSet,))])
        experiment_IDs_to_add = sorted(ExperimentIDs.difference(ExperimentIDsInPredictionSet))

        if experiment_IDs_to_add:
            colortext.printf("\nAdding %d jobs to the prediction set." % len(experiment_IDs_to_add), "lightgreen")
            count = 0
            for experiment_ID in experiment_IDs_to_add:
                colortext.write('.', "lightgreen")
                self.addPrediction(experiment_ID, None, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True, strip_other_chains = strip_other_chains)
                count +=1
                #if count >= 1:
                #    sys.exit(0)
        else:
            colortext.printf("\nAll jobs are already in the queue or have been run.", "lightgreen")
        print('')

    def addPrediction(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False, strip_other_chains = True):
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
            pdbID = results[0]["PDBFileID"]

            if PDB_ID:
                #sql = "SELECT ID, Content FROM PDBFile WHERE ID=%s"
                results = self.ddGDB.execute_select("SELECT ID, Content FROM PDBFile WHERE ID=%s", parameters=(PDB_ID))
                if len(results) != 1:
                    raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
                predictionPDB_ID = results[0]["ID"]
                pdbID = results[0]["ID"]
            else:
                predictionPDB_ID = experimentPDB_ID

            # Get the related PDB ID and file
            assert(len(results) == 1)
            result = results[0]
            contents = result["Content"]

            pdb = PDB(contents.split("\n"))

            # Check that the mutated positions exist and that the wild-type matches the PDB
            mutations = self.ddGDB.call_select_proc("GetMutations", parameters = parameters)

            # todo: Hack. This should be removed when PDB homologs are dealt with properly.
            mutation_objects = []
            for mutation in mutations:
                if experimentPDB_ID == "1AJ3" and predictionPDB_ID == "1U5P":
                    assert(int(mutation['ResidueID']) < 1000)
                    mutation['ResidueID'] = str(int(mutation['ResidueID']) + 1762)
                mutation_objects.append(Mutation(mutation['WildTypeAA'], mutation['ResidueID'], mutation['MutantAA'], mutation['Chain']))

            #todo: a
            #checkPDBAgainstMutations(pdbID, pdb, mutations)
            pdb.validate_mutations(mutation_objects)

            #for mutation in mutations:
            #    if experimentPDB_ID == "ub_OTU":
            #        mutation['ResidueID'] = str(int(mutation['ResidueID']) + 172)

            # Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
            chains = [result['Chain'] for result in self.ddGDB.call_select_proc("GetChains", parameters = parameters)]
            if strip_other_chains:
                pdb.stripForDDG(chains, KeepHETATMLines, numberOfModels = 1)
            else:
                pdb.stripForDDG(True, KeepHETATMLines, numberOfModels = 1)

            #print('\n'.join(pdb.lines))
            # - Post stripping checks -
            # Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
            # then check again that the mutated positions exist and that the wild-type matches the PDB
            colortext.warning('mutations %s' % (str(mutations)))

            remappedMutations = pdb.remapMutations(mutations, pdbID)

            #resfile = self._createResfile(pdb, remappedMutations)
            mutfile = self._createMutfile(pdb, remappedMutations)

            # Check to make sure that we haven't stripped all the ATOM lines
            if not pdb.GetAllATOMLines():
                raise colortext.Exception("No ATOM lines remain in the stripped PDB file of %s." % pdbID)

            # Check to make sure that CSE and MSE are not present in the PDB
            badresidues = pdb.CheckForPresenceOf(["CSE", "MSE"])
            if badresidues:
                raise colortext.Exception("Found residues [%s] in the stripped PDB file of %s. These should be changed to run this job under Rosetta." % (', '.join(badresidues), pdbID))

            # Turn the lines array back into a valid PDB file
            strippedPDB = '\n'.join(pdb.lines)
        except Exception, e:
            colortext.error("Error in %s, %s: .\n%s" % (experimentID, UserDataSetExperimentID, traceback.format_exc()))
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
            rdmstring = ''.join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16))
            cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
            cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
            entryDate = self.ddGDB.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))

    def get_flattened_prediction_results(self, PredictionSet):
        #todo: add this as a stored procedure
        return self.ddGDB.execute_select('''
SELECT Prediction.ID AS PredictionID, Prediction.ExperimentID, Experiment.PDBFileID, ExperimentMutations.FlattenedMutations, Prediction.Scores, TIMEDIFF(Prediction.EndDate, Prediction.StartDate) AS TimeTaken FROM Prediction INNER JOIN
(
  SELECT ExperimentID, GROUP_CONCAT(Mutation SEPARATOR ', ') AS FlattenedMutations FROM
  (
    SELECT ExperimentID, CONCAT(Chain, ' ', WildTypeAA, ResidueID, MutantAA) As Mutation FROM ExperimentMutation
  ) AS FlattenedMutation
  GROUP BY ExperimentID
) AS ExperimentMutations
ON Prediction.ExperimentID=ExperimentMutations.ExperimentID
INNER JOIN Experiment ON Prediction.ExperimentID=Experiment.ID
WHERE Prediction.PredictionSet=%s AND Prediction.Scores IS NOT NULL
ORDER BY Prediction.ExperimentID''', parameters=(PredictionSet,))


    #### todo: these following functions should be refactored and renamed. In particular, the graphing functions should
    ####       be moved into the tools repository

    def create_abacus_graph_for_a_single_structure(self, PredictionSet, scoring_method, scoring_type, graph_title = None, PredictionIDs = None, graph_filename = None, cached_results = None, num_datapoints = 0):
        '''This function is meant to
            The num_datapoints is mainly for debugging - tuning the resolution/DPI to fit the number of datapoints.'''
        import json

        results = cached_results
        if not results:
            results = self.get_flattened_prediction_results(PredictionSet)

        pdb_ids = set()
        for r in results:
            pdb_ids.add(r['PDBFileID'])

        if len(pdb_ids) != 1:
            raise Exception('This function is only meant to be called when the PredictionSet or the set of results contains records for a single structure. The set of results contains %d structures.' % len(pdb_ids))

        sortable_results = {}
        for r in results:
            if (not PredictionIDs) or (r['PredictionID'] in PredictionIDs):
                sortable_results[(json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], r['ExperimentID'])] = r
        count = 0

        set_of_mutations = set()
        for k, r in sorted(sortable_results.iteritems()):
            #if r['FlattenedMutations'].find('A E141L') != -1 and r['FlattenedMutations'].find('A S142A') != -1 and r['FlattenedMutations'].find('A L78Y') != -1:
            #    print('%f, %s' % (k[0], r['FlattenedMutations']))
            #if r['FlattenedMutations'].find('A W103M') != -1 and r['FlattenedMutations'].find('A F70Y') != -1:
            #    if r['FlattenedMutations'].find('A E141L') == -1 and r['FlattenedMutations'].find('A S142A') == -1 and r['FlattenedMutations'].find('A L78Y') == -1:
            #        print('%f, %s' % (k[0], r['FlattenedMutations']))

            if r['FlattenedMutations'].find('A W103M') != -1 and r['FlattenedMutations'].find('A F70Y') != -1:
                if r['FlattenedMutations'].find('A E141L') == -1 and r['FlattenedMutations'].find('A S142A') == -1 and r['FlattenedMutations'].find('A L78Y') == -1:
                    #print('%f, %s' % (k[0], r['FlattenedMutations']))
                    count += 1
            #A E141L, A S142A

            mutations = [m for m in map(string.strip, r['FlattenedMutations'].split(',')) if m]
            for m in mutations:
                set_of_mutations.add((int(m.split()[1][1:-1]), m))
            #if r['FlattenedMutations'].find('A L78Y') == -1:
            #    print('%f, %s' % (k[0], r['FlattenedMutations']))
            #    #count += 1

        pruned_data = []
        for k, r in sorted(sortable_results.iteritems()):
            line = []
            print(json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], r['FlattenedMutations'])
            for m in sorted(set_of_mutations):
                if r['FlattenedMutations'].find(m[1]) != -1:
                    line.append(1)
                else:
                    line.append(0)
            pruned_data.append((json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], line))

        labels = [m[1].split()[1] for m in sorted(set_of_mutations)]
        graph_title = graph_title or r'$\Delta\Delta$G predictions for %s (%s.%s)' % (PredictionSet, scoring_method.replace(',0A', '.0$\AA$').replace('_', ' '), scoring_type)

        pruned_data = pruned_data[0:num_datapoints or len(pruned_data)]
        colortext.message('Creating graph with %d datapoints...' % len(pruned_data))

        number_of_non_zero_datapoints = 0
        for p in pruned_data:
            if 1 in p[1]:
                number_of_non_zero_datapoints += 1
                if number_of_non_zero_datapoints > 1:
                    break
        if number_of_non_zero_datapoints < 2:
            raise Exception('The dataset must contain at least two non-zero points.')

        if graph_filename:
            return self.write_abacus_graph(graph_filename, graph_title, labels, pruned_data, scoring_method, scoring_type)
        else:
            return self.create_abacus_graph(graph_title, labels, pruned_data, scoring_method, scoring_type)

    def write_abacus_graph(self, graph_filename, graph_title, labels, data, scoring_method, scoring_type):
        byte_stream = self.create_abacus_graph(graph_title, labels, data, scoring_method, scoring_type)
        write_file(graph_filename, byte_stream.getvalue(), 'wb')

    def create_abacus_graph(self, graph_title, labels, data, scoring_method, scoring_type):
        '''Even though this is technically a scatterplot, I call this an abacus graph because it is basically a bunch of beads on lines.'''

        from io import BytesIO
        if plt:
            assert(data)

            image_dpi = 300.0
            horizontal_margin = 400.0 # an estimate of the horizontal space not used by the graph
            horizontal_spacing = 100.0 # an estimate of the horizontal space between points on the same line
            vertical_margin = 100.0 # an estimate of the vertical space not used by the graph
            vertical_spacing = 50.0 # the rough amount of pixels between abacus lines
            point_size = 50 # the size of datapoints in points^2.
            y_offset = 1.0

            points_per_line = set([len(line[1]) for line in data])
            assert(len(points_per_line) == 1)
            points_per_line = points_per_line.pop()
            assert(len(labels) == points_per_line)

            number_of_lines = float(len(data))
            number_of_labels = float(len(labels))
            height_in_inches = max(600/image_dpi, (vertical_margin + (vertical_spacing * number_of_lines)) / image_dpi) # Use a minimum of 600 pixels in height. This avoids graphs with a small number of lines (<=10) not to become squashed.
            width_in_inches = max(700/image_dpi, (horizontal_margin + (horizontal_spacing * points_per_line)) / image_dpi) # Use a minimum of 600 pixels in width. This avoids graphs with a small number of labels (e.g. 1) not to become squashed.

            graph_color_scheme = matplotlib.cm.jet

            #y_offset = (1.75 * data_length) / 128
            #image_dpi = (400 * data_length) / 128
            #image_dpi = 400
            #point_sizes = {1 : 100, 64: 75, 128: 50, 192: 25, 256: 10}
            #index = round(data_length / 64.0) * 64
            #point_size = point_sizes.get(index, 10)


            fig = plt.figure(figsize=(width_in_inches, height_in_inches)) # figsize is specified in inches - w, h
            fig.set_dpi(image_dpi)

            # Create three identically-sized lists. Each triple is of an x-coordinate, a y-coordinate, and the DDG value
            # and corresponds to a 1 in the matrix i.e. we should draw a point/abacus bead at these coordinates.
            x_values = []
            y_values = []
            x_coordinate_skip = 3 # x-axis distance between two points
            y_coordinate_skip = 7 # y-axis distance between two points
            ddg_values = []
            y = 0
            for line in data:
                x = 0
                y += y_coordinate_skip
                w = line[0]
                #plt.text(30, y, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=9)
                for point in line[1]:
                    x += x_coordinate_skip
                    if point == 1:
                        x_values.append(x)
                        y_values.append(y)
                        ddg_values.append(line[0])

            # Draw the scatter plot
            plt.scatter(x_values, y_values, c=ddg_values, s=point_size, cmap=graph_color_scheme, edgecolors='none', zorder=99)

            # Define the limits of the cartesian coordinates. Add extra space on the right for the DDG values.
            extra_space = 1.3
            plt.axis((0, (points_per_line + 1 + extra_space) * x_coordinate_skip, -15, (y_coordinate_skip * number_of_lines) + 15))

            plt.tick_params(
                axis='both',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='off',      # ticks along the bottom edge are off
                left='off',      # ticks along the left edge are off
                labelleft='off', # labels along the bottom edge are off
                top='off',         # ticks along the top edge are off
                labelbottom='off') # labels along the bottom edge are off

            # Add the mutation labels at the botom of the diagram
            x = 1.9
            for i in range(len(labels)):
                l = labels[i]
                plt.text(x, -5 + ((i % 2) * -5), l, fontdict=None, withdash=True, fontsize=6)
                x += x_coordinate_skip

            added_zero_line = False
            last_y_value = 0
            y = 0
            for line in data:
                x = 0
                y += 7
                plt.plot([1, 25], [y, y], color='#999999', linestyle='-', linewidth=0.1)

                # Add a DDG value on every third line
                if y % 21 == 7:
                    plt.text(((points_per_line + 0.6) * x_coordinate_skip) , y-y_offset, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)
                if not added_zero_line:
                    if line[0] > 0:
                        plt.plot([1, 25], [0.5 + ((y + last_y_value) / 2), 0.5 + ((y + last_y_value) / 2)], color='k', linestyle='-', linewidth=1)
                        added_zero_line = True
                    else:
                        last_y_value = y

            plt.text(((points_per_line + 0.6) * x_coordinate_skip), y-y_offset, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)

            # Set the colorbar font size and then add a colorbar
            #cbar.ax.tick_params(labelsize=6)
            #plt.colorbar(use_gridspec=True)

            #ax = fig.add_subplot(111)

            # Add a title. Note: doing this after the colorbar code below messes up the alignment.
            # Adjust the wrap length to the width of the graph
            wrap_length = 40 + max(0, (points_per_line - 3) * 14)
            graph_title = "\n".join(textwrap.wrap(graph_title, wrap_length))

            plt.title(graph_title, fontdict={'fontsize' : 6})

            from mpl_toolkits.axes_grid1 import make_axes_locatable
            ax = fig.add_subplot(111)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05) # pad specifies the padding between the right of the graph and the left of the colorbar, size seems to specify the width of the colobar relative to the graph

            CS3 = plt.contourf([[0,0],[0,0]], ddg_values, cmap=graph_color_scheme)
            #plt.colorbar(CS3)
            #cbar = fig.colorbar(CS3, format='%.2f')
            cbar = fig.colorbar(CS3, format='%.2f', cax=cax)
            cbar.set_label('$\Delta\Delta$G',size=6)
            cbar.ax.tick_params(labelsize=5)

             # Use the tight_layout command to tighten up the spaces. The pad, w_pad, and h_pad parameters are specified in fraction of fontsize.
            plt.tight_layout(pad=0.5)

            #quadmesh = ax.pcolormesh(theta,phi,data)
            #cb = fig.colorbar(quadmesh,ax=ax, shrink=.5, pad=.2, aspect=10)
            #cax = ax.imshow(ddg_values, interpolation='nearest', cmap=matplotlib.cm.coolwarm)
            #cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
            #surf = ax.contourf(X,Y,Z, 8, cmap=cm.jet)
            #cbar = fig.colorbar(surf, use_gridspec=True, shrink=0.5, aspect=20, fraction=.12,pad=.02)
            #cbar.set_label('Activation',size=18)

            byte_stream = BytesIO()
            plt.savefig(byte_stream, dpi=image_dpi, format="png")

            plt.close(fig)

            return byte_stream
        else:
            return None

    def extract_data(output_dir, PredictionID):
        assert(os.path.exists(output_dir))
        archive = self.getData(PredictionID)
        write_file(os.path.join(output_dir, '%d.zip' % PredictionID), archive, 'wb')
        p = Popen(output_dir, ['unzip', '%d.zip' % PredictionID])
        os.remove(os.path.join(output_dir, '%d.zip' % PredictionID))
        if p.errorcode != 0:
            raise colortext.Exception(p.stderr)
        else:
            colortext.warning(p.stdout)

    def test_results(output_dir, PredictionSet):
        PredictionIDs = []
        results = get_flattened_prediction_results(PredictionSet)
        mutation_lists = {}
        for r in results:
            PredictionIDs.append(r['PredictionID'])
            mutation_lists[r['PredictionID']] = r['FlattenedMutations']
        RandomPredictionIDs = [PredictionIDs[random.randint(0, len(PredictionIDs) - 1)] for k in range(10)]
        RandomPredictionIDs = [54090L, 53875L, 54085L, 54079L, 54008L, 53853L, 53952L, 54056L, 53935L, 53893L]

        # Retrieve and unzip results
        if not(os.path.exists(output_dir)):
            os.mkdir(output_dir)
        for PredictionID in PredictionIDs:#RandomPredictionIDs:
            if not(os.path.exists(os.path.join(output_dir, str(PredictionID)))):
                colortext.message('Retrieving archive for Prediction %d.' % PredictionID)
                self.extract_data(output_dir, PredictionID)

        # Get the sequences of the wildtype and mutant structures
        count = 0
        for PredictionID in PredictionIDs:#RandomPredictionIDs:
            wildtype_sequences = set()
            mutation_sequences = set()
            working_dir = os.path.join(os.path.join(output_dir, str(PredictionID)))
            for f in glob.glob(os.path.join(working_dir, '*.pdb')):
                if os.path.split(f)[1].startswith('mut_'):
                    p = PDB.from_filepath(f)
                    assert(len(p.atom_sequences) == 1)
                    sequence = str(p.atom_sequences.values()[0])
                    mutation_sequences.add(sequence)
                elif os.path.split(f)[1].startswith('repacked_wt_'):
                    p = PDB.from_filepath(f)
                    assert(len(p.atom_sequences) == 1)
                    sequence = str(p.atom_sequences.values()[0])
                    wildtype_sequences.add(sequence)

            assert(len(wildtype_sequences) == 1)
            assert(len(mutation_sequences) == 1)
            wildtype_sequence = wildtype_sequences.pop()
            mutation_sequence = mutation_sequences.pop()

            colortext.message('Prediction %d. Mutations: %s' % (PredictionID, mutation_lists[PredictionID]))
            assert(len(wildtype_sequence) == len(mutation_sequence))
            s = ''
            t = ''
            for x in range(len(wildtype_sequence)):
                if wildtype_sequence[x] != mutation_sequence[x]:
                    s += colortext.make(wildtype_sequence[x], color="green")
                    t += colortext.make(mutation_sequence[x], color="yellow")
                else:
                    s += wildtype_sequence[x]
                    t += mutation_sequence[x]
            print(s)
            print(t)

    def create_pymol_session(download_dir, PredictionID, task_number, keep_files = True):
        '''Create a PyMOL session for a pair of structures.'''

        # Retrieve and unzip results
        if not(os.path.exists(download_dir)):
            os.mkdir(download_dir)
        working_dir = os.path.join(os.path.join(download_dir, str(PredictionID)))
        if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
            self.extract_data(download_dir, PredictionID)
        if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
            raise Exception('Could not extract the models for task #%d of Prediction #%d.' % (task_number, PredictionID))

        # Retrieve the two structures corresponding to the task_number
        files = sorted(glob.glob(os.path.join(working_dir, '*_round_%d.pdb' % task_number)), reverse = True)
        assert(os.path.split(files[0])[1].startswith('repacked_wt_'))
        assert(os.path.split(files[1])[1].startswith('mut_'))

        # Creator the alignment object and write the PSE file
        chain_mapper = ScaffoldModelChainMapper.from_file_paths(files[0], files[1])

        # Remove the downloaded files
        if not keep_files:
            shutil.rmtree(download_dir)
        return chain_mapper.generate_pymol_session()

    def create_pymol_session_in_memory(self, PredictionID, task_number):
        '''Create a PyMOL session for a pair of structures.
        '''
        # todo: This should replace create_pymol_session since creating a PyMOL session should be a separate task from extracting an archive. write_pymol_session should be changed to use this function

        from io import BytesIO

        # Retrieve and unzip results
        archive = self.getData(PredictionID)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            # Get the name of the files from the zip
            wildtype_filename = os.path.join(str(PredictionID), 'repacked_wt_round_%d.pdb' % task_number)
            mutant_filename = None
            for filepath in sorted(zipped_content.namelist()):
                filename = os.path.split(filepath)[1]
                if filename.startswith('mut_') and filename.endswith('_round_%d.pdb' % task_number):
                    mutant_filename = os.path.join(str(PredictionID), filename)
                    break

            PyMOL_session = None
            file_list = zipped_content.namelist()

            # If both files exist in the zip, extract their contents in memory and create a PyMOL session pair (PSE, script)
            if (mutant_filename in file_list) and (wildtype_filename in file_list):
                wildtype_pdb = zipped_content.open(wildtype_filename, 'r').read()
                mutant_pdb = zipped_content.open(mutant_filename, 'U').read()
                chain_mapper = ScaffoldModelChainMapper.from_file_contents(wildtype_pdb, mutant_pdb)
                PyMOL_session = chain_mapper.generate_pymol_session(pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol')

            zipped_content.close()
            return PyMOL_session

        except Exception, e:
            zipped_content.close()
            raise Exception(str(e))

    def write_pymol_session(download_dir, PredictionID, task_number, keep_files = True):
        PSE_file = create_pymol_session(download_dir, PredictionID, task_number, keep_files = keep_files)
        write_file(output_filepath, PSE_file[0], 'wb')

