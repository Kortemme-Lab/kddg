#!/usr/bin/python2.4
# encoding: utf-8
"""
monomer_api.py
High-level functions for interacting with the protein stability sections of the ddG database.

Classes:
MonomericStabilityDDGInterface - an class used to interface with the database. Call get_interface to get a user API based on this class.
AnalysisBreakdown - an class used to run analyses on the data

Note: I moved this code from db_api.py during a large refactor and have not tested it yet.
      A lot of functionality is currently broken but all the pieces are there. See Trac ticket #1375.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

from io import BytesIO
import os
import zipfile
import traceback
import copy
import pprint

import numpy
from sqlalchemy import and_

from klab import colortext
from klab.bio.alignment import ScaffoldModelChainMapper

from klab.benchmarking.analysis.ddg_monomeric_stability_analysis import DBBenchmarkRun as MonomericStabilityBenchmarkRun
from klab.benchmarking.analysis.ddg_binding_affinity_analysis import DBBenchmarkRun as BindingAffinityBenchmarkRun

import db_schema as dbmodel
from api_layers import *
from db_api import ddG, PartialDataException


def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a MonomericStabilityDDGInterface interface object. It hides the
    private methods from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(MonomericStabilityDDGInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)


class MonomericStabilityDDGInterface(ddG):


    def __init__(self, passwd = None, username = 'kortemmelab', hostname = None, rosetta_scripts_path = None, rosetta_database_path = None):
        super(MonomericStabilityDDGInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)
        self.prediction_data_path = self.DDG_db.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionDataPath"')[0]['Value']


    #########################################################################################
    ## Broken API layer
    ##
    ## This section contains useful functions which need to be updated to work with the new
    ## schema or code
    #########################################################################################


    #== Deprecated functions =================================================================

    @deprecated
    def get_prediction_experiment_chains(self, predictionset): raise Exception('This function has been deprecated. Use get_pdb_chains_used_for_prediction_set instead.')


    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================


    @informational_pdb
    def get_pdb_chains_for_prediction(self, prediction_id):
        raise Exception('This needs to be implemented.')



    #todo: remove this and replace with get_pdb_mutations_for_mutagenesis
    @informational_pdb
    def get_pdb_mutations_for_experiment(self, experiment_id):
        '''Returns the PDB mutations for a mutagenesis experiment as well as the PDB residue information.'''
        pdb_mutations = []

#            SELECT ExperimentMutation.*, Experiment.PDBFileID, PDBResidue.ResidueType,
#            PDBResidue.BFactorMean, PDBResidue.BFactorDeviation,
#            PDBResidue.ComplexExposure, PDBResidue.ComplexDSSP, PDBResidue.MonomericExposure, PDBResidue.MonomericDSSP
#            FROM
#            Experiment INNER JOIN ExperimentMutation ON Experiment.ID = ExperimentMutation.ExperimentID
#            INNER JOIN
#            PDBResidue ON Experiment.PDBFileID = PDBResidue.PDBFileID AND ExperimentMutation.Chain = PDBResidue.Chain AND ExperimentMutation.ResidueID = PDBResidue.ResidueID AND ExperimentMutation.WildTypeAA = PDBResidue.ResidueAA
#            WHERE Experiment.ID=%s ORDER BY Chain, ResidueID''', parameters=(experiment_id,)):


        for pdb_mutation in self.DDG_db.execute_select('''
            SELECT ExperimentMutation.*, Experiment.PDBFileID
            FROM
            Experiment INNER JOIN ExperimentMutation ON Experiment.ID = ExperimentMutation.ExperimentID
            WHERE Experiment.ID=%s ORDER BY Chain, ResidueID''', parameters=(experiment_id,)):

                pdb_mutation['ResidueType'] = None
                pdb_mutation['BFactorMean'] = None
                pdb_mutation['BFactorDeviation'] = None
                pdb_mutation['ComplexExposure'] = None
                pdb_mutation['ComplexDSSP'] = None
                pdb_mutation['MonomericExposure'] = None
                pdb_mutation['MonomericDSSP'] = None
                pdb_mutations.append(pdb_mutation)
        return pdb_mutations


    #todo: remove this and replace with get_user_dataset_experiment_details
    def get_experiment_details(self, experiment_id):
        e = self.DDG_db.execute_select('SELECT * FROM Experiment WHERE ID=%s', parameters=(experiment_id,))
        if len(e) != 1:
            raise colortext.Exception('Experiment %d does not exist.' % (experiment_id, ))
        e = e[0]

        pdb_mutations = self.get_pdb_mutations_for_experiment(experiment_id)
        assert(len(pdb_mutations) > 0)
        pdb_id = set([m['PDBFileID'] for m in pdb_mutations])
        assert(len(pdb_id) == 1)
        pdb_id = pdb_id.pop()
        return dict(
            Mutagenesis = dict(
                ExperimentID = e['ID'],
            ),
            Structure = dict(
                PDBFileID = pdb_id,
            ),
            PDBMutations = pdb_mutations,
        )


    @informational_job
    def get_job_details(self, prediction_id, include_files = True, truncate_content = None):
        prediction_record = self.DDG_db.execute_select('SELECT * FROM Prediction WHERE ID=%s', parameters=(prediction_id,))
        if not prediction_record:
            raise Exception('No details could be found for prediction #%d in the database.' % prediction_id)
        prediction_record = prediction_record[0]
        prediction_record['Files'] = {}
        if include_files:
            prediction_record['Files'] = self.get_job_files(prediction_id, truncate_content = truncate_content)

        # mutfile_content = self.create_mutfile(prediction_id)

        # Read the UserPPDataSetExperiment details
        user_dataset_experiment_id = prediction_record['UserDataSetExperimentID']
        if user_dataset_experiment_id:
            ude_details = self.get_user_dataset_experiment_details(user_dataset_experiment_id)
            assert(ude_details['Mutagenesis']['PPMutagenesisID'] == prediction_record['PPMutagenesisID'])
            for k, v in ude_details.iteritems():
                assert(k not in prediction_record)
                prediction_record[k] = v
        else:
            # todo: Remove this later
            e_details = self.get_experiment_details(prediction_record['ExperimentID'])
            for k, v in e_details.iteritems():
                assert(k not in prediction_record)
                prediction_record[k] = v
        return prediction_record


    @informational_job
    def get_predictions_experimental_details(self, prediction_id, userdatset_experiment_ids_to_subset_ddgs = None, include_files = False, reference_ids = set(), include_experimental_data = True):

        details = self.get_job_details(prediction_id, include_files = include_files)

        # Add the DDG values for the related analysis sets
        if include_experimental_data:
            userdatset_experiment_ids_to_subset_ddgs = userdatset_experiment_ids_to_subset_ddgs or self.get_experimental_ddgs_by_analysis_set(user_dataset_experiment_id, reference_ids = reference_ids)
            assert('DDG' not in details)
            details['DDG'] = userdatset_experiment_ids_to_subset_ddgs[user_dataset_experiment_id]
        else:
            details['DDG'] = None

        return details


    ###########################################################################################
    ## Prediction creation/management layer
    ##
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    @job_creator
    def add_prediction_set(self, prediction_set_id, halted = True, priority = 5, batch_size = 40, allow_existing_prediction_set = False):
        return super(MonomericStabilityDDGInterface, self).add_prediction_set(prediction_set_id, halted = halted, priority = priority, batch_size = batch_size, allow_existing_prediction_set = allow_existing_prediction_set, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False)


    @job_creator
    def add_job(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, keep_hetatm_lines, PDB_ID = None, ReverseMutation = False, InputFiles = {}, testonly = False, strip_other_chains = True):
        '''This function inserts a prediction into the database.
            The parameters define:
                the experiment we are running the prediction for;
                the name of the set of predictions for later grouping;
                the short description of the Command to be used for prediction;
                whether HETATM lines are to be kept or not.
            We strip the PDB based on the chains used for the experiment and keep_hetatm_lines.
            We then add the prediction record, including the stripped PDB and the inverse mapping
            from Rosetta residue numbering to PDB residue numbering.'''
        raise Exception('This function needs to be rewritten.')
        raise Exception('Make sure to call charge by residue function, _charge_prediction_set_by_residue_count')

        parameters = (experimentID,)
        assert(ReverseMutation == False) # todo: allow this later
        try:
            predictionPDB_ID = None

            sql = "SELECT PDBFileID, Content FROM Experiment INNER JOIN PDBFile WHERE Experiment.PDBFileID=PDBFile.ID AND Experiment.ID=%s"
            results = self.DDG_db.execute_select(sql, parameters = parameters)
            if len(results) != 1:
                raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
            experimentPDB_ID = results[0]["PDBFileID"]
            pdbID = results[0]["PDBFileID"]

            if PDB_ID:
                #sql = "SELECT ID, Content FROM PDBFile WHERE ID=%s"
                results = self.DDG_db.execute_select("SELECT ID, Content FROM PDBFile WHERE ID=%s", parameters=(PDB_ID))
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
            mutations = self.DDG_db.call_select_proc("GetMutations", parameters = parameters)

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
            chains = [result['Chain'] for result in self.DDG_db.call_select_proc("GetChains", parameters = parameters)]
            if strip_other_chains:
                pdb.stripForDDG(chains, keep_hetatm_lines, numberOfModels = 1)
            else:
                pdb.stripForDDG(True, keep_hetatm_lines, numberOfModels = 1)

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
        ExtraParameters = pickle.dumps(ExtraParameters)

        PredictionFieldNames = self.DDG_db.FieldNames.Prediction
        params = {
            PredictionFieldNames.ExperimentID		: experimentID,
            PredictionFieldNames.UserDataSetExperimentID : UserDataSetExperimentID,
            PredictionFieldNames.PredictionSet		: PredictionSet,
            PredictionFieldNames.ProtocolID			: ProtocolID,
            PredictionFieldNames.KeptHETATMLines	: keep_hetatm_lines,
            PredictionFieldNames.StrippedPDB		: strippedPDB,
            PredictionFieldNames.ResidueMapping		: pickle.dumps(pdb.get_ddGInverseResmap()),
            PredictionFieldNames.InputFiles			: InputFiles,
            PredictionFieldNames.Status 			: "queued",
            PredictionFieldNames.ExtraParameters	: ExtraParameters,
        }
        if not testonly:
            self.DDG_db.insertDict('Prediction', params)

            # Add cryptID string
            predictionID = self.DDG_db.getLastRowID()
            entryDate = self.DDG_db.execute_select("SELECT EntryDate FROM Prediction WHERE ID=%s", parameters = (predictionID,))[0]["EntryDate"]
            rdmstring = ''.join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16))
            cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
            cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
            entryDate = self.DDG_db.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))
            return predictionID


    @job_creator
    def add_jobs_by_pdb_id(self, pdb_ID, PredictionSet, ProtocolID, status = 'active', priority = 5, keep_hetatm_lines = False, strip_other_chains = True):
        raise Exception('This function needs to be rewritten.')
        colortext.printf("\nAdding any mutations for this structure which have not been queued/run in the %s prediction set." % PredictionSet, "lightgreen")

        d = {
            'ID' : PredictionSet,
            'Status' : status,
            'Priority' : priority,
            'BatchSize' : 40,
            'EntryDate' : datetime.datetime.now(),
        }
        self.DDG_db.insertDictIfNew('PredictionSet', d, ['ID'])

        # Update the priority and activity if necessary
        self.DDG_db.execute('UPDATE PredictionSet SET Status=%s, Priority=%s WHERE ID=%s', parameters = (status, priority, PredictionSet))

        # Determine the set of experiments to add
        ExperimentIDs = set([r['ID'] for r in self.DDG_db.execute_select('SELECT ID FROM Experiment WHERE PDBFileID=%s', parameters=(pdb_ID,))])
        ExperimentIDsInPredictionSet = set([r['ExperimentID'] for r in self.DDG_db.execute_select('SELECT ExperimentID FROM Prediction WHERE PredictionSet=%s', parameters=(PredictionSet,))])
        experiment_IDs_to_add = sorted(ExperimentIDs.difference(ExperimentIDsInPredictionSet))

        if experiment_IDs_to_add:
            colortext.printf("\nAdding %d jobs to the prediction set." % len(experiment_IDs_to_add), "lightgreen")
            count = 0
            for experiment_ID in experiment_IDs_to_add:
                colortext.write('.', "lightgreen")
                self.addPrediction(experiment_ID, None, PredictionSet, ProtocolID, keep_hetatm_lines, strip_other_chains = strip_other_chains)
                count +=1
        else:
            colortext.printf("\nAll jobs are already in the queue or have been run.", "lightgreen")
        print('')


    @job_creator
    def add_prediction_run(self, user_dataset_name, prediction_set_id, protocol_id, keep_hetatm_lines, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False):
        raise Exception('This function needs to be rewritten.')

        assert(self.DDG_db.execute_select("SELECT ID FROM PredictionSet WHERE ID=%s", parameters=(prediction_set_id,)))

        #results = self.DDG_db.execute_select("SELECT * FROM UserDataSet WHERE TextID=%s", parameters=(user_dataset_name,))
        results = self.DDG_db.execute_select("SELECT UserDataSetExperiment.* FROM UserDataSetExperiment INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%s", parameters=(user_dataset_name,))
        if not results:
            return False

        if not(quiet):
            colortext.message("Creating predictions for UserDataSet %s using protocol %s" % (user_dataset_name, protocol_id))
            colortext.message("%d records found in the UserDataSet" % len(results))

        count = 0
        showprogress = not(quiet) and len(results) > 300
        if showprogress:
            print("|" + ("*" * (int(len(results)/100)-2)) + "|")
        for r in results:

            existing_results = self.DDG_db.execute_select("SELECT * FROM Prediction WHERE PredictionSet=%s AND UserDataSetExperimentID=%s", parameters=(prediction_set_id, r["ID"]))
            if len(existing_results) > 0:
                #colortext.warning('There already exist records for this UserDataSetExperimentID. You probably do not want to proceed. Skipping this entry.')
                continue

            PredictionID = self.addPrediction(r["ExperimentID"], r["ID"], prediction_set_id, protocol_id, keep_hetatm_lines, PDB_ID = r["PDBFileID"], ReverseMutation = False, InputFiles = {}, testonly = testonly)
            count += 1
            if showprogress:
                if count > 100:
                    colortext.write(".", "cyan", flush = True)
                    count = 0
            if shortrun and count > 4:
                break
        print("")
        return(True)


    @job_creator
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set):
        raise Exception('not implemented yet')
        #assert(existing_prediction_set exists and has records)
        #assert(new_prediction_set is empty)
        #for each prediction record, add the record and all associated predictionfile records,


    #== Input file generation API ===========================================================
    #
    # This part of the API is responsible for creating input files for predictions


    @job_input
    def create_resfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    @job_input
    def create_mutfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    #== Job execution/completion API ===========================================================
    #
    # This part of the API is responsible for starting jobs and setting them as failed or
    # completed


    @job_execution
    def get_job(self, prediction_set):
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_execution
    def start_job(self, prediction_id, prediction_set):
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_execution
    def get_max_number_of_cluster_jobs(self, prediction_set_id, priority):
        return self.DDG_db.execute_select('SELECT Value FROM _DBCONSTANTS WHERE VariableName="MaxStabilityClusterJobs"')['Value']


    @job_completion
    def fail_job(self, prediction_id, prediction_set, maxvmem, ddgtime):
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_completion
    def parse_prediction_scores(self, stdout):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore records.'''
        self._parse_ddg_monomer_scores_per_structure(stdout)


    @job_completion
    def store_scores(self, prediction_set, prediction_id, scores):
        '''Stores a list of dicts suitable for database storage e.g. PredictionStructureScore records.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @job_completion
    def complete_job(self, prediction_id, prediction_set, scores, maxvmem, ddgtime):
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @staticmethod
    def _parse_ddg_monomer_scores_per_structure(stdout):
        '''Returns a dict mapping the DDG scores from a ddg_monomer run to a list of structure numbers.'''

        # Parse the stdout into two mappings (one for wildtype structures, one for mutant structures) mapping
        # structure IDs to a dict containing the score components
        wildtype_scores = {}
        mutant_scores = {}
        s1 = 'score before mutation: residue'
        s1_len = len(s1)
        s2 = 'score after mutation: residue'
        s2_len = len(s2)
        for line in stdout.split('\n'):
            idx = line.find(s1)
            if idx != -1:
                idx += s1_len
                mtchs = re.match('.*?(\d+) %s' % s1, line)
                structure_id = int(mtchs.group(1))
                assert(structure_id not in wildtype_scores)
                tokens = line[idx:].split()
                d = {'total' : float(tokens[0])}
                for x in range(1, len(tokens), 2):
                    component_name = tokens[x].replace(':', '')
                    assert(rosetta_weights.get(component_name))
                    component_value = float(tokens[x + 1])
                    d[component_name] = component_value
                wildtype_scores[structure_id] = d
            else:
                idx = line.find(s2)
                if idx != -1:
                    idx += s2_len
                    mtchs = re.match('.*?(\d+) %s' % s2, line)
                    structure_id = int(mtchs.group(1))
                    assert(structure_id not in mutant_scores)
                    tokens = line[idx:].split()
                    d = {'total' : float(tokens[1])}
                    for x in range(2, len(tokens), 2):
                        component_name = tokens[x].replace(':', '')
                        assert(rosetta_weights.get(component_name))
                        component_value = float(tokens[x + 1])
                        d[component_name] = component_value
                    mutant_scores[structure_id] = d

        # Sanity checks
        num_structures = max(wildtype_scores.keys())
        expected_keys = set(range(1, num_structures + 1))
        assert(expected_keys == set(wildtype_scores.keys()))
        assert(expected_keys == set(mutant_scores.keys()))

        # Create a list of lists - MutantScoreOrder - of structure IDs e.g. [[5,1,34], [23], [12,3], ...] which is ordered
        # by increasing energy so that each sublist contains structure IDs of equal energy and if structures have the same
        # energy then their IDs are in the same sublist
        d = {}
        for structure_id, scores in sorted(mutant_scores.iteritems()):
            d[scores['total']] = d.get(scores['total'], [])
            d[scores['total']].append(structure_id)
        MutantScoreOrder = []
        for score, structure_ids in sorted(d.iteritems()):
            MutantScoreOrder.append(structure_ids)

        # Sanity check - make sure that MutantScoreOrder is really ordered such that each set of structure IDs contains
        # structures of the same energy and of a lower energy than the following set of structure IDs in the list
        for x in range(len(MutantScoreOrder) - 1):
            s1 = set([mutant_scores[n]['total'] for n in MutantScoreOrder[x]])
            assert(len(s1) == 1)
            if x + 1 < len(MutantScoreOrder):
                s2 = set([mutant_scores[n]['total'] for n in MutantScoreOrder[x + 1]])
                assert(len(s2) == 1)
                assert(s1.pop() < s2.pop())

        return dict(
            WildType = wildtype_scores,
            Mutant = mutant_scores,
            MutantScoreOrder = MutantScoreOrder,
        )


    ###########################################################################################
    ## Prediction results layer
    ##
    ## This part of the API for returning data about completed predictions.
    ###########################################################################################


    @job_results
    def get_ddg_scores_per_structure(self, prediction_id):
        # At present, we only use ddg_monomer
        raise Exception('Reimplement using the database records.')


    ###########################################################################################
    ## Analysis layer
    ##
    ## This part of the API is responsible for running analysis on completed predictions
    ###########################################################################################


    @analysis_api
    def determine_best_pair(self, prediction_id, score_method_id):
        # Iterates over the (wildtype, mutant) pairs in the PredictionStructureScore table and returns the structure ID
        # for the pair with the lowest energy mutant
        # Note: There are multiple ways to select the best pair. For example, if multiple mutants have the same minimal total
        # score, we could have multiple wildtype structures to choose from. In this case, we choose a pair where the wildtype
        # structure has the minimal total score.

        lowest_wt_score = self.DDG_db.execute_select('SELECT StructureID, total FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="WildType" ORDER BY total LIMIT 1', parameters=(prediction_id, score_method_id))
        lowest_mutant_score = self.DDG_db.execute_select('SELECT StructureID, total FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" ORDER BY total LIMIT 1', parameters=(prediction_id, score_method_id))
        if lowest_wt_score and lowest_mutant_score:
            return lowest_wt_score[0]['StructureID'], lowest_mutant_score[0]['StructureID']
        return None, None

        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        mutant_complexes = []
        wildtype_complexes = []
        for structure_id, scores in scores.iteritems():
            if scores.get('MutantComplex'):
                mutant_complexes.append((scores['MutantComplex']['total'], structure_id))
            if scores.get('WildTypeComplex'):
                wildtype_complexes.append((scores['WildTypeComplex']['total'], structure_id))
        wildtype_complexes = sorted(wildtype_complexes)
        mutant_complexes = sorted(mutant_complexes)
        if wildtype_complexes and mutant_complexes:
            return wildtype_complexes[0][1], mutant_complexes[0][1]
        return None, None



        lowest_mutant_score = self.DDG_db.execute_select('SELECT total FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" ORDER BY total LIMIT 1', parameters=(prediction_id, score_method_id))
        if lowest_mutant_score:
            lowest_mutant_score = lowest_mutant_score[0]['total']
            mutant_structure_ids = [r['StructureID'] for r in self.DDG_db.execute_select('SELECT StructureID FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" AND total=%s', parameters=(prediction_id, score_method_id, lowest_mutant_score))]
            if len(mutant_structure_ids) > 1:
                return self.DDG_db.execute_select(('SELECT StructureID FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="WildType" AND StructureID IN (' + ','.join(map(str, mutant_structure_ids)) + ') ORDER BY total LIMIT 1'), parameters=(prediction_id, score_method_id ))[0]['StructureID']
            else:
                return mutant_structure_ids[0]
        return None


    @analysis_api
    def get_prediction_data(self, prediction_id, score_method_id, main_ddg_analysis_type, top_x = 3, expectn = None, extract_data_for_case_if_missing = True, root_directory = None, dataframe_type = "Stability"):
        assert(dataframe_type == "Binding affinity") # todo: stability case needs to be written
        try:
            top_x_ddg = self.get_top_x_ddg_affinity(prediction_id, score_method_id, top_x = top_x, expectn = expectn)
        except Exception, e:
            colortext.pcyan(str(e))
            colortext.warning(traceback.format_exc())
            if extract_data_for_case_if_missing:
                self.extract_data_for_case(prediction_id, root_directory = root_directory, force = True, score_method_id = score_method_id)
            try:
                top_x_ddg = self.get_top_x_ddg_affinity(prediction_id, score_method_id, top_x = top_x, expectn = expectn)
            except PartialDataException, e:
                raise
            except Exception, e:
                raise
        top_x_ddg_stability = self.get_top_x_ddg_stability(prediction_id, score_method_id, top_x = top_x, expectn = expectn)
        return {
            main_ddg_analysis_type : top_x_ddg,
            'DDGStability_Top%d' % top_x : top_x_ddg_stability,
        }


    @analysis_api
    def get_top_x_ddg_affinity(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''This function was taken from the PPI API.'''
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
            return None
        try:
            #colortext.warning(prediction_id)
            #pprint.pprint(scores)
            wt_total_scores = [(scores[struct_num]['WildTypeComplex']['total'], struct_num) for struct_num in scores]
            wt_total_scores.sort()
            top_x_wt_struct_nums = [t[1] for t in wt_total_scores[:top_x]]
            #print(wt_total_scores)

            mut_total_scores = [(scores[struct_num]['MutantComplex']['total'], struct_num) for struct_num in scores]
            mut_total_scores.sort()
            top_x_mut_struct_nums = [t[1] for t in mut_total_scores[:top_x]]
            #print(mut_total_scores)

            top_x_score = numpy.average([
                (scores[mut_struct_num]['MutantComplex']['total'] - scores[mut_struct_num]['MutantLPartner']['total'] - scores[mut_struct_num]['MutantRPartner']['total']) -
                (scores[wt_struct_num]['WildTypeComplex']['total'] - scores[wt_struct_num]['WildTypeLPartner']['total'] - scores[wt_struct_num]['WildTypeRPartner']['total'])
                for wt_struct_num, mut_struct_num in zip(top_x_wt_struct_nums, top_x_mut_struct_nums)
            ])
            return top_x_score
        except Exception, e:
            print(e)
            colortext.warning(traceback.format_exc())
            raise PartialDataException('The case is missing some data.')



    @analysis_api
    def get_top_x_ddg_stability(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''Returns the TopX value for the prediction only considering the complex scores. This computation may work as a
           measure of a stability DDG value.'''
        scores = self.get_prediction_scores(prediction_id, expectn = expectn).get(score_method_id)
        if scores == None:
            return None

        wt_total_scores = [(scores[struct_num]['WildTypeComplex']['total'], struct_num) for struct_num in scores]
        wt_total_scores.sort()
        top_x_wt_struct_nums = [t[1] for t in wt_total_scores[:top_x]]

        mut_total_scores = [(scores[struct_num]['MutantComplex']['total'], struct_num) for struct_num in scores]
        mut_total_scores.sort()
        top_x_mut_struct_nums = [t[1] for t in mut_total_scores[:top_x]]

        return numpy.average([scores[mut_struct_num]['MutantComplex']['total'] - scores[wt_struct_num]['WildTypeComplex']['total']
                           for wt_struct_num, mut_struct_num in zip(top_x_wt_struct_nums, top_x_mut_struct_nums)])


    @analysis_api
    def get_top_x_ddg(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''Returns the TopX value for the prediction. Typically, this is the mean value of the top X predictions for a
           case computed using the associated Score records in the database.'''

        # Make sure that we have as many cases as we expect
        if expectn != None:
            scores = self.get_prediction_scores(prediction_id)[score_method_id]
            num_cases = 0
            for k in scores.keys():
                if type(k) == type(1L):
                    num_cases += 1
            if num_cases != expectn:
                raise Exception('Expected scores for {0} runs; found {1}.'.format(expectn, num_cases))

        # Call the stored procedure which takes the top_x-lowest wildtype scores and gets their average then subtracts this from the average of the top_x-lowest mutant scores
        try:
            r = self.DDG_db.call_select_proc('MonomericStabilityTopX', parameters=(55808, score_method_id, 3), quiet=False)
            assert(len(r) == 1)
            return r[0]['TopX']
        except Exception, e:
            raise Exception('An error occurred determining the Top{0} score for prediction #{1} using score method {2}: "{3}"\n{4}'.format(top_x, prediction_id, score_method_id, str(e), traceback.print_exc()))


    @analysis_api
    def get_analysis_dataframe(self, prediction_set_id,
            experimental_data_exists = True,
            create_binding_affinity_dataframe = False, # Hack added for PUBS class
            prediction_set_series_name = None, prediction_set_description = None, prediction_set_credit = None,
            prediction_set_color = None, prediction_set_alpha = None,
            use_existing_benchmark_data = True,
            include_derived_mutations = False,
            use_single_reported_value = False,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            score_method_id = None,
            expectn = None,
            allow_failures = False,
            extract_data_for_case_if_missing = True,
            ):

        #todo: rename function since we return BenchmarkRun objects

        assert(score_method_id)
        assert(experimental_data_exists == False) # todo: I am implementing the case needed for the PUBS class - the general case still needs to be implemented. Change BindingAffinityBenchmarkRun below accordingly (BindingAffinityBenchmarkRun for affinity, MonomericStabilityBenchmarkRun for stability)

        # We allow different dataframe types as sometimes there will be no associated experimental data
        dataframe_type = "Stability"
        if (not experimental_data_exists) and create_binding_affinity_dataframe:
            # A hacky case used for the PUBS year 1 results which were monomeric stability predictions rescored to be used as binding affinity predictions
            dataframe_type = "Binding affinity"
        assert(create_binding_affinity_dataframe == False or (not experimental_data_exists)) # Part of the PUBS hack. We cannot create this dataframe since we do not have associated experimental data

        parameters = copy.copy(locals())
        del parameters['self']
        del parameters['create_binding_affinity_dataframe']
        return super(MonomericStabilityDDGInterface, self)._get_analysis_dataframe(BindingAffinityBenchmarkRun, **parameters)


    ################################################################################################
    ## Application layer
    ## These functions combine the database and prediction data with useful klab
    ################################################################################################


    #== PyMOL API ===========================================================


    @app_pymol
    def create_pymol_session_in_memory(self, prediction_id, wt_task_number, mutant_task_number, pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol'):

        # Retrieve and unzip results
        archive = self.get_job_data(prediction_id)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            # Get the name of the files from the zip
            wildtype_filename = os.path.join(str(prediction_id), 'repacked_wt_round_%d.pdb' % wt_task_number)
            mutant_filename = None
            for filepath in sorted(zipped_content.namelist()):
                filename = os.path.split(filepath)[1]
                if filename.startswith('mut_') and filename.endswith('_round_%d.pdb' % mutant_task_number):
                    mutant_filename = os.path.join(str(prediction_id), filename)
                    break

            PyMOL_session = None
            file_list = zipped_content.namelist()

            # If both files exist in the zip, extract their contents in memory and create a PyMOL session pair (PSE, script)
            if (mutant_filename in file_list) and (wildtype_filename in file_list):
                wildtype_pdb = zipped_content.open(wildtype_filename, 'r').read()
                mutant_pdb = zipped_content.open(mutant_filename, 'U').read()

                # todo: this should be structure_1_name = 'Wildtype', structure_2_name = 'Mutant' but the underlying PyMOL script needs to be parameterized
                chain_mapper = ScaffoldModelChainMapper.from_file_contents(wildtype_pdb, mutant_pdb, structure_1_name = 'Scaffold', structure_2_name = 'Model')
                PyMOL_session = chain_mapper.generate_pymol_session(pymol_executable = pymol_executable)

            zipped_content.close()
            return PyMOL_session

        except Exception, e:
            zipped_content.close()
            raise Exception(str(e))






    ################################################################################################
    ## Subclass-specific API layer
    ## These are functions written specifically for this class which are not necessarily available
    ## in sibling classes
    ################################################################################################


    @analysis_api
    def get_predictionset_data(self, predictionset, userdataset_textid, cached_pdb_details = None, only_single = False):
        '''
            A helper function for analysis / generating graphs.
            Arguments:
                predictionset - the name of a PredictionSet
                cached_pdb_details - a cached copy of the pdb_details returned by this function. Generating this dict is the slowest step so caching is recommended.
                only_single - restrict to single mutations
            Returns:
                amino_acids - details about amino acid types
                pdb_details - details about PDB file techniques, resolution, chain lengths, and whether it is a transmembrane protein
                predictions - a mapping: Prediction IDs -> ExperimentID, UserDataSetExperimentID, Experiment and Prediction (this is the one used for the prediction) PDB IDs, scores. If single mutation then also mutation details, DSSP, and exposure.
                analysis_datasets - a mapping: analysis subset (e.g. "Guerois") -> Prediction IDs -> (prediction) PDB_ID, ExperimentID, ExperimentDDG (mean of experimental values)
        '''
        UserDataSetID = self.DDG_db.execute_select("SELECT ID FROM UserDataSet WHERE TextID=%s", parameters=(userdataset_textid,))
        assert(UserDataSetID)
        UserDataSetID = UserDataSetID[0]['ID']

        amino_acids = self.get_amino_acids_for_analysis()
        prediction_chains = self._get_pdb_chains_used_for_prediction_set(predictionset)

        # Get the list of mutation predictions
        if only_single:
            prediction_records = self.DDG_db.execute_select('''
                SELECT a.ID AS PredictionID FROM
                (
                SELECT Prediction.ID, Prediction.ExperimentID, UserDataSetExperimentID, COUNT(Prediction.ID) AS NumMutations
                FROM Prediction
                INNER JOIN ExperimentMutation ON ExperimentMutation.ExperimentID=Prediction.ExperimentID
                WHERE PredictionSet = %s
                GROUP BY Prediction.ID
                ) AS a
                WHERE a.NumMutations=1''', parameters=(predictionset,))
        else:
            prediction_records = self.DDG_db.execute_select('''
                SELECT a.ID AS PredictionID FROM
                (
                SELECT Prediction.ID, Prediction.ExperimentID, UserDataSetExperimentID, COUNT(Prediction.ID) AS NumMutations
                FROM Prediction
                INNER JOIN ExperimentMutation ON ExperimentMutation.ExperimentID=Prediction.ExperimentID
                WHERE PredictionSet = %s
                GROUP BY Prediction.ID
                ) AS a''', parameters=(predictionset,))
        allowed_prediction_ids = set([m['PredictionID'] for m in prediction_records])

        kellogg_score_method_id = self.DDG_db.execute_select('''SELECT ID FROM ScoreMethod WHERE MethodName='Global' AND MethodType='Protocol 16' ''')
        assert(len(kellogg_score_method_id) == 1)
        kellogg_score_method_id = kellogg_score_method_id[0]['ID']

        noah_8A_positional_score_method_id = self.DDG_db.execute_select('''SELECT ID FROM ScoreMethod WHERE MethodName='Local' AND MethodType='Position' ''')
        assert(len(noah_8A_positional_score_method_id) == 1)
        noah_8A_positional_score_method_id = noah_8A_positional_score_method_id[0]['ID']

        # Hack - add on the mutations from the datasets which were represented as single mutants (in the original datasets) but which are double mutants
        # See ExperimentAssays ( 917, 918, 919, 920, 922, 7314, 932, 933, 936, 937, 938, 2076, 7304, 7305, 7307, 7308, 7309, 7310, 7312, 7315, 7316, 7317, 7320 )
        # or Experiments (111145, 110303, 110284, 110287, 110299, 110300, 110285, 110286, 110289, 114180, 114175, 114177, 114171, 110304, 110305, 114179, 114168, 114170, 114172, 114173, 114178, 114167)
        # or PubMed IDs 7479708, 9079363, and 9878405)
        badly_entered_predictions = self.DDG_db.execute_select('''
                    SELECT Prediction.ID AS PredictionID FROM Prediction
                    INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID
                    WHERE PredictionSet=%s
                    AND ExperimentID IN (111145, 110303, 110284, 110287, 110299, 110300, 110285, 110286, 110289, 114180, 114175, 114177, 114171, 110304, 110305, 114179, 114168, 114170, 114172, 114173, 114178, 114167)''', parameters=(predictionset,))
        badly_entered_predictions = set([r['PredictionID'] for r in badly_entered_predictions])
        allowed_prediction_ids = allowed_prediction_ids.union(badly_entered_predictions)


        # Read in the PredictionStructureScore records
        kellogg_structure_score_query = self.DDG_db.execute_select('''
            SELECT PredictionID, ScoreType, StructureID, total FROM PredictionStructureScore
            WHERE PredictionID >=%s
            AND PredictionID <=%s
            AND (ScoreType = 'Mutant' OR ScoreType = 'WildType')
            AND ScoreMethodID=%s
            ''', parameters=(min(allowed_prediction_ids), max(allowed_prediction_ids), kellogg_score_method_id))
        kellogg_structure_scores = {}
        for kss in kellogg_structure_score_query:
            PredictionID = kss['PredictionID']
            kellogg_structure_scores[PredictionID] = kellogg_structure_scores.get(PredictionID, {'WildType' : {}, 'Mutant' : {}})
            kellogg_structure_scores[PredictionID][kss['ScoreType']][kss['StructureID']] = kss['total']

        # Get the Prediction records for the mutation predictions and the list of PDB IDs
        num_predictions = 0
        predictions = {}
        experiment_to_prediction_map = {}
        prediction_ids = set()
        pdb_ids = set()
        failures = 0
        prediction_results = self.DDG_db.execute_select('''
            SELECT Prediction.ID AS PredictionID, Prediction.ExperimentID, UserDataSetExperimentID, Experiment.PDBFileID AS ePDB, UserDataSetExperiment.PDBFileID AS pPDB, Scores
            FROM Prediction
            INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID
            INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=Prediction.UserDataSetExperimentID
            WHERE PredictionSet=%s''', parameters=(predictionset,))
        for p in prediction_results:
            id = p['PredictionID']
            experiment_id = p['ExperimentID']
            if id not in allowed_prediction_ids:
                continue
            num_predictions += 1
            experiment_to_prediction_map[(experiment_id, p['pPDB'])] = id
            assert(id not in predictions)
            prediction_ids.add(id)
            pdb_ids.add(p['ePDB'])
            pdb_ids.add(p['pPDB'])
            if p['Scores']:
                scores = json.loads(p['Scores'])

                # Retrieve the scores from the PredictionStructureScore records
                # todo: this is only done for the Kellogg scores at present
                kellogg_output_score = scores['data']['kellogg']['total']['ddG']
                individual_scores = kellogg_structure_scores.get(id)
                assert(individual_scores)
                assert(len(individual_scores['Mutant']) == 50) # todo: parameterize
                assert(len(individual_scores['WildType']) == 50) # todo: parameterize
                sorted_mutant_scores = sorted(individual_scores['Mutant'].values())
                sorted_wildtype_scores = sorted(individual_scores['WildType'].values())

                # Compute three scores - the best pair, the average of the best 3 pairs, and the average of the best 5 pairs.
                kellogg_top1 = sorted_mutant_scores[0] - sorted_wildtype_scores[0]
                kellogg_top3 = (sum(sorted_mutant_scores[:3]) - sum(sorted_wildtype_scores[:3]))/3.0
                kellogg_top5 = (sum(sorted_mutant_scores[:5]) - sum(sorted_wildtype_scores[:5]))/5.0

                print(kellogg_output_score, kellogg_top1, kellogg_top3, kellogg_top5)
                assert(abs(kellogg_output_score - kellogg_top1) < 0.01)
                p['Kellogg_top1'] = kellogg_top1
                p['Kellogg_top3'] = kellogg_top3
                p['Kellogg_top5'] = kellogg_top5

                if scores['data'].get('noah_8,0A'):
                    p['Noah'] = scores['data']['noah_8,0A']['positional']['ddG']
                else:
                    p['Noah'] = None
            else:
                p['Kellogg_top1'] = None
                p['Kellogg_top3'] = None
                p['Kellogg_top5'] = None
                p['Noah'] = None
                failures += 1
            del p['PredictionID']
            del p['Scores']
            predictions[id] = p
        assert(len(experiment_to_prediction_map) == num_predictions)


        # Get the PDB chain for each prediction
        missing_count = 0
        for pc in prediction_chains:
            if pc['ID'] not in allowed_prediction_ids:
                continue
            prediction = predictions.get(pc['ID'])
            if prediction:
                assert(None == prediction.get('Chain'))
                prediction['Chain'] = pc['Chain']
            else:
                raise Exception('Missing chain data')


        # Get the mutation details for each single mutation prediction
        mutation_details_1 = self.DDG_db.execute_select('''
SELECT a.ID AS PredictionID, UserDataSetExperiment.PDBFileID as pPDB, ExperimentMutation.Chain, ExperimentMutation.ResidueID, ExperimentMutation.WildTypeAA, ExperimentMutation.MutantAA,
PDBResidue.MonomericExposure, PDBResidue.MonomericDSSP
FROM
(
SELECT Prediction.ID, Prediction.ExperimentID, UserDataSetExperimentID, COUNT(Prediction.ID) AS NumMutations
FROM Prediction
INNER JOIN ExperimentMutation ON ExperimentMutation.ExperimentID=Prediction.ExperimentID
WHERE PredictionSet = %s
GROUP BY Prediction.ID
) AS a
INNER JOIN ExperimentMutation ON a.ExperimentID=ExperimentMutation.ExperimentID
INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=a.UserDataSetExperimentID
INNER JOIN PDBResidue
  ON (PDBResidue.PDBFileID=UserDataSetExperiment.PDBFileID
  AND PDBResidue.Chain=ExperimentMutation.Chain
  AND TRIM(PDBResidue.ResidueID)=TRIM(ExperimentMutation.ResidueID))
WHERE a.NumMutations=1''', parameters=(predictionset,))
        # Hack for 1U5P. Note: TRIM removes warnings e.g. "Warning: Truncated incorrect INTEGER value: '1722 '".
        mutation_details_2 = self.DDG_db.execute_select('''
SELECT a.ID AS PredictionID, UserDataSetExperiment.PDBFileID as pPDB, ExperimentMutation.Chain, ExperimentMutation.ResidueID, ExperimentMutation.WildTypeAA, ExperimentMutation.MutantAA,
PDBResidue.MonomericExposure, PDBResidue.MonomericDSSP
FROM
(
SELECT Prediction.ID, Prediction.ExperimentID, UserDataSetExperimentID, COUNT(Prediction.ID) AS NumMutations
FROM Prediction
INNER JOIN ExperimentMutation ON ExperimentMutation.ExperimentID=Prediction.ExperimentID
WHERE PredictionSet = %s
GROUP BY Prediction.ID
) AS a
INNER JOIN ExperimentMutation ON a.ExperimentID=ExperimentMutation.ExperimentID
INNER JOIN UserDataSetExperiment ON UserDataSetExperiment.ID=a.UserDataSetExperimentID
INNER JOIN PDBResidue
  ON (PDBResidue.PDBFileID=UserDataSetExperiment.PDBFileID
  AND PDBResidue.Chain=ExperimentMutation.Chain
  AND CAST(TRIM(PDBResidue.ResidueID) AS UNSIGNED) - 1762 = CAST(TRIM(ExperimentMutation.ResidueID) AS UNSIGNED))
WHERE a.NumMutations=1 AND UserDataSetExperiment.PDBFileID="1U5P" ''', parameters=(predictionset,), quiet=True)
        mutation_details = mutation_details_1 + mutation_details_2
        all_prediction_ids = set([m['PredictionID'] for m in prediction_records])
        found_prediction_ids = set([m['PredictionID'] for m in mutation_details])
        #assert(len(found_prediction_ids) == len(all_prediction_ids))
        for m in mutation_details:
            prediction = predictions[m['PredictionID']]
            prediction['DSSP'] = dssp_elision.get(m['MonomericDSSP'])
            prediction['Exposure'] = m['MonomericExposure']
            prediction['WTAA'] = m['WildTypeAA']
            prediction['MutantAA'] = m['MutantAA']
            prediction['ResidueID'] = m['ResidueID']

        # Add missing fields for the set of badly_entered_predictions and multiple mutations
        for prediction_id, d in sorted(predictions.iteritems()):
            if 'DSSP' not in d.keys():
                # todo: if there is a G or P in any mutation, mark this record ["GP"] = True
                # use this below to separate the GP mutations, rather than checking the wtaa and mutaa there
                prediction['DSSP'] = None
                prediction['Exposure'] = None
                prediction['WTAA'] = None
                prediction['MutantAA'] = None
                prediction['ResidueID'] = None

        # We can derive the following data:
        #   TM, Resolution, XRay per PDB
        #   for prediction_id, d in predictions.iteritems():
        #     assoc_pdb = pdb_details[d['pPDB']]
        #     d['TM'] = assoc_pdb['TM'] == 1
        #     d['XRay'] = assoc_pdb['XRay']
        #     d['Resolution'] = assoc_pdb['Resolution']
        # Derive GP (if wt or mutant is glycine or proline)
        # Derive WTPolarity, WTAromaticity, MutantPolarity, MutantAromaticity
        # Derive SL, LS, SS, LL
        # Derive ChainLength: prediction['ChainLength'] = pdbs[pc['PDBFileID']]['chains'][pc['Chain']]

        # AnalysisSets are defined on UserDataSets. The main 'datasets' are defined in the Subset field of the AnalysisSet records associated
        # with the UserDataSet i.e. AnalysisSets for a UserDataSet := "SELECT DISTINCT Subset FROM UserAnalysisSet WHERE UserDataSetID=x".
        from analysis import UserDataSetExperimentalScores
        analysis_data = {}
        analysis_subsets = [r['Subset'] for r in self.DDG_db.execute_select("SELECT DISTINCT Subset FROM UserAnalysisSet WHERE UserDataSetID=%s", parameters=(UserDataSetID,))]
        for analysis_subset in analysis_subsets:
            analysis_data[analysis_subset] = {}
            adata = analysis_data[analysis_subset]
            adata['Missing'] = []

            UDS_scores = UserDataSetExperimentalScores(self.DDG_db, 1, analysis_subset)
            count = 0
            for section, sectiondata in sorted(UDS_scores.iteritems()):
                for recordnumber, record_data in sorted(sectiondata.iteritems()):
                    PDB_ID = record_data["PDB_ID"]
                    ExperimentID = record_data["ExperimentID"]
                    ExperimentalDDG = record_data["ExperimentalDDG"]
                    prediction_id = experiment_to_prediction_map.get((record_data['ExperimentID'], record_data['PDB_ID']))
                    if prediction_id != None:
                        adata[prediction_id] = record_data # PDB_ID, ExperimentID, ExperimentalDDG
                    else:
                        adata['Missing'].append((record_data['ExperimentID'], record_data['PDB_ID']))
                    count += 1

        return AnalysisBreakdown(
            amino_acids,
            self.get_pdb_details_for_analysis(pdb_ids, cached_pdb_details = cached_pdb_details),
            predictions,
            analysis_data, # a mapping: analysis subset (e.g. "Guerois") -> Prediction IDs -> (prediction) PDB_ID, ExperimentID, ExperimentDDG (mean of experimental values)
        )


    ### Dataset stats functions


    @analysis_api
    def get_analysis_set_overlap_by_Experiment(self, restrict_to_subsets = set(), UserDataSetID = 1):
        ''' Returns the overlap between analysis sets of a UserDataSet where overlap is determined by the set of ExperimentIDs.
            Caveat: This assumes that the Experiments do not overlap. While this is mostly true at present, there are probably
                    still some duplicates.
            Returns a symmetric matrix (as a pandas dataframe) with the pairwise overlaps.
            Usage: self.get_dataset_overlap_by_Experiment(restrict_to_subsets = ['CuratedProTherm', 'Guerois', 'Kellogg', 'Potapov']).'''

        import numpy
        import pandas as pd

        # Read the list of experiments from the database
        analysis_set_experiments = {}
        results = self.DDG_db.execute_select('SELECT * FROM UserAnalysisSet WHERE UserDataSetID=%s', parameters=(UserDataSetID,))
        for r in results:
            subset = r['Subset']
            if (len(restrict_to_subsets) == 0) or (subset in restrict_to_subsets):
                analysis_set_experiments[subset] = analysis_set_experiments.get(subset, set())
                analysis_set_experiments[subset].add(r['ExperimentID'])

        analysis_sets = sorted(analysis_set_experiments.keys())
        m = []
        for x in analysis_sets:
            mx = []
            for y in analysis_sets:
                mx.append(len(analysis_set_experiments[x].intersection(analysis_set_experiments[y])))
            m.append(mx)
        df = pd.DataFrame(m, index = analysis_sets, columns = analysis_sets)
        return df


    @analysis_api
    def get_analysis_set_disjoint_by_Experiment(self, primary_subset, other_subsets = set(), UserDataSetID = 1):
        ''' Returns the overlap between analysis sets of a UserDataSet where overlap is determined by the set of ExperimentIDs.
            Caveat: This assumes that the Experiments do not overlap. While this is mostly true at present, there are probably
                    still some duplicates.
            Returns a symmetric matrix (as a pandas dataframe) with the pairwise overlaps.
            Usage: self.get_dataset_overlap_by_Experiment(other_subsets = ['CuratedProTherm', 'Guerois', 'Kellogg', 'Potapov']).'''

        import numpy
        import pandas as pd

        # Read the list of experiments from the database
        analysis_set_experiments = {}
        primary_analysis_set_experiments = set()
        results = self.DDG_db.execute_select('SELECT * FROM UserAnalysisSet WHERE UserDataSetID=%s', parameters=(UserDataSetID,))
        for r in results:
            subset = r['Subset']
            if subset != primary_subset:
                if (len(other_subsets) == 0) or (subset in other_subsets):
                    analysis_set_experiments[subset] = analysis_set_experiments.get(subset, set())
                    analysis_set_experiments[subset].add(r['ExperimentID'])
            else:
                primary_analysis_set_experiments.add(r['ExperimentID'])

        # Create the subsets
        analysis_set_common_experiments = {}
        analysis_set_disjoint_experiments = {}
        analysis_sets = sorted(analysis_set_experiments.keys())
        m = []
        for x in analysis_sets:
            analysis_set_common_experiments[x] = analysis_set_experiments[x].intersection(primary_analysis_set_experiments)
            analysis_set_disjoint_experiments[x] = analysis_set_experiments[x].difference(primary_analysis_set_experiments)
            assert(len(analysis_set_common_experiments[x]) + len(analysis_set_disjoint_experiments[x]) == len(analysis_set_experiments[x]))

        assert(primary_subset not in analysis_set_experiments.keys())
        all_common_analysis_set_experiments = {}
        for x in analysis_set_experiments.keys():
            all_common_analysis_set_experiments[x] = analysis_set_common_experiments[x]
        all_common_analysis_set_experiments[primary_subset] = primary_analysis_set_experiments
        all_common_analysis_sets = sorted(all_common_analysis_set_experiments.keys())
        m = []
        for x in all_common_analysis_sets:
            mx = []
            for y in all_common_analysis_sets:
                mx.append(len(all_common_analysis_set_experiments[x].intersection(all_common_analysis_set_experiments[y])))
            m.append(mx)
        df = pd.DataFrame(m, index = all_common_analysis_sets, columns = all_common_analysis_sets)

        other_sets = analysis_set_experiments.keys()
        if len(other_sets) == 3:
            print(other_sets)
            print('Common to all three other sets in the intersection: %d' % len(analysis_set_common_experiments[other_sets[0]].intersection(analysis_set_common_experiments[other_sets[1]]).intersection(analysis_set_common_experiments[other_sets[2]])))
            print('Common to all three other sets in the disjoint set: %d' % len(analysis_set_disjoint_experiments[other_sets[0]].intersection(analysis_set_disjoint_experiments[other_sets[1]]).intersection(analysis_set_disjoint_experiments[other_sets[2]])))
            print('ere')

        print('\nCommon to primary set\n')
        print(df)

        all_disjoint_analysis_sets = sorted(analysis_set_disjoint_experiments.keys())
        m = []
        for x in all_disjoint_analysis_sets:
            mx = []
            for y in all_disjoint_analysis_sets:
                mx.append(len(analysis_set_disjoint_experiments[x].intersection(analysis_set_disjoint_experiments[y])))
            m.append(mx)
        df = pd.DataFrame(m, index = all_disjoint_analysis_sets, columns = all_disjoint_analysis_sets)
        print('\nDisjoint from primary set\n')
        print(df)

        return df


    @analysis_api
    def get_analysis_set_overlap_by_Experiment_as_radii(self, max_radius, restrict_to_subsets = set(), UserDataSetID = 1):
        '''Todo: look at where this was called and figure out what I was doing. I think this was used to decide on the relative sizes of overlaps between datasets in a Venn diagram.'''
        import numpy
        df = self.get_analysis_set_overlap_by_Experiment(restrict_to_subsets = restrict_to_subsets, UserDataSetID = UserDataSetID)

        # Determine the relative sizes for the radii (a = pi.r^2 so pi cancels out).
        radii_ratios = df.apply(lambda x: numpy.sqrt(x))

        # Get the max value and determine the scalar
        max_value = max(numpy.hstack(radii_ratios.values))
        scalar = float(max_radius) / max_value

        # Returned the scaled matrix
        return radii_ratios.apply(lambda x: x * scalar)








    ################################################################################################
    ## Private API layer
    ## These are helper functions used internally by the class but which are not intended for export
    ################################################################################################


    ###########################################################################################
    ## Subclass layer
    ##
    ## These functions need to be implemented by subclasses
    ###########################################################################################

    # Concrete functions


    def _get_sqa_prediction_table(self): return dbmodel.Prediction
    def _get_sqa_prediction_structure_scores_table(self): return dbmodel.PredictionStructureScore

    def _get_prediction_table(self): return 'Prediction'
    def _get_prediction_structure_scores_table(self): return 'PredictionStructureScore'
    def _get_prediction_type(self): return 'ProteinStability'
    def _get_prediction_dataset_type(self): return 'Protein stability'
    def _get_prediction_type_description(self): return 'monomeric stability'
    def _get_user_dataset_experiment_table(self): return 'UserDataSetExperiment'
    def _get_user_dataset_experiment_tag_table(self): raise Exception('To be added.')
    def _get_allowed_score_types(self): return set(['DDG', 'WildType', 'Mutant'])


    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================

    @informational_pdb
    def _get_pdb_chains_used_for_prediction_set(self, prediction_set):
        return self.DDG_db.execute_select('''
            SELECT Prediction.ID, Experiment.PDBFileID, Chain
            FROM Prediction
            INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID
            INNER JOIN ExperimentChain ON ExperimentChain.ExperimentID=Prediction.ExperimentID
            WHERE PredictionSet=%s''', parameters=(prediction_set,))


    ###########################################################################################
    ## Prediction layer
    ##
    ## This part of the API is responsible for inserting prediction jobs in the database via
    ## the trickle-down proteomics paradigm.
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    def _charge_prediction_set_by_residue_count(self, PredictionSet):
        '''This function assigns a cost for a prediction equal to the number of residues in the chains.'''
        raise Exception('This function needs to be rewritten.')
        from klab.bio.rcsb import parseFASTAs

        DDG_db = self.DDG_db
        predictions = DDG_db.execute_select("SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet=%s", parameters=(PredictionSet,))

        PDB_chain_lengths ={}
        for prediction in predictions:
            chain_records = DDG_db.execute_select('SELECT PDBFileID, Chain FROM Experiment INNER JOIN ExperimentChain ON ExperimentID=Experiment.ID WHERE ExperimentID=%s', parameters=(prediction['ExperimentID']))
            num_residues = 0
            for chain_record in chain_records:
                key = (chain_record['PDBFileID'], chain_record['Chain'])

                if PDB_chain_lengths.get(key) == None:
                    fasta = DDG_db.execute_select("SELECT FASTA FROM PDBFile WHERE ID=%s", parameters = (chain_record['PDBFileID'],))
                    assert(len(fasta) == 1)
                    fasta = fasta[0]['FASTA']
                    f = parseFASTAs(fasta)
                    PDB_chain_lengths[key] = len(f[chain_record['PDBFileID']][chain_record['Chain']])
                chain_length = PDB_chain_lengths[key]
                num_residues += chain_length

            print("UPDATE Prediction SET Cost=%0.2f WHERE ID=%d" % (num_residues, prediction['ID']))

            predictions = DDG_db.execute("UPDATE Prediction SET Cost=%s WHERE ID=%s", parameters=(num_residues, prediction['ID'],))

