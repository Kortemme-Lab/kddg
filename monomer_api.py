#!/usr/bin/python2.4
# encoding: utf-8
"""
monomer_api.py
High-level functions for interacting with the protein stability sections of the ddG database.

Classes:
MonomericStabilityDDGInterface - an class used to interface with the database
AnalysisBreakdown - an class used to run analyses on the data

Note: I moved this code from db_api.py during a large refactor and have not tested it yet.
      A lot of functionality is currently broken but all the pieces are there. See Trac ticket #1375.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

from io import BytesIO
import os
import zipfile

from api_layers import *
from db_api import ddG
from tools import colortext
from tools.bio.alignment import ScaffoldModelChainMapper


def get_interface(passwd, username = 'kortemmelab'):
    '''This is the function that should be used to get a MonomericStabilityDDGInterface interface object. It hides the
    private methods from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(MonomericStabilityDDGInterface, passwd = passwd, username = username)


class MonomericStabilityDDGInterface(ddG):


    def __init__(self, passwd = None, username = 'kortemmelab'):
        super(MonomericStabilityDDGInterface, self).__init__(passwd = passwd, username = username)
        self.prediction_data_path = self.DDG_db.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionDataPath"')[0]['Value']



    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================


    @informational_pdb
    def get_pdb_chains_for_prediction(self, prediction_id):
        '''Returns the PDB file ID and a list of chains for the prediction.'''
        raise Exception('This needs to be implemented.')


    @informational_pdb
    def get_pdb_chains_used_for_prediction_set(self, prediction_set):
        return self.DDG_db.execute_select('''
            SELECT Prediction.ID, Experiment.PDBFileID, Chain
            FROM Prediction
            INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID
            INNER JOIN ExperimentChain ON ExperimentChain.ExperimentID=Prediction.ExperimentID
            WHERE PredictionSet=%s''', parameters=(prediction_set,))


    ###########################################################################################
    ## Prediction creation/management layer
    ##
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    @job_creator
    def add_prediction_set(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40, allow_existing_prediction_set = False):
        '''Adds a new PredictionSet (a construct used to group Predictions) to the database.'''
        return super(MonomericStabilityDDGInterface, self).add_prediction_set(PredictionSetID, halted = halted, Priority = Priority, BatchSize = BatchSize, allow_existing_prediction_set = allow_existing_prediction_set, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False)


    @job_creator
    def add_job(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False, strip_other_chains = True):
        '''This function inserts a prediction into the database.
            The parameters define:
                the experiment we are running the prediction for;
                the name of the set of predictions for later grouping;
                the short description of the Command to be used for prediction;
                whether HETATM lines are to be kept or not.
            We strip the PDB based on the chains used for the experiment and KeepHETATMLines.
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

        PredictionFieldNames = self.DDG_db.FieldNames.Prediction
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
    def add_jobs_by_pdb_id(self, pdb_ID, PredictionSet, ProtocolID, status = 'active', priority = 5, KeepHETATMLines = False, strip_other_chains = True):
        ''' This function adds predictions for all Experiments corresponding to pdb_ID to the specified prediction set.
            This is useful for custom runs e.g. when we are using the DDG scheduler for design rather than for benchmarking.
        '''
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
                self.addPrediction(experiment_ID, None, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True, strip_other_chains = strip_other_chains)
                count +=1
        else:
            colortext.printf("\nAll jobs are already in the queue or have been run.", "lightgreen")
        print('')


    @job_creator
    def add_prediction_run(self, userdatasetTextID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = False, Description = {}, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False):
        raise Exception('This function needs to be rewritten.')

        assert(self.DDG_db.execute_select("SELECT ID FROM PredictionSet WHERE ID=%s", parameters=(PredictionSet,)))

        #results = self.DDG_db.execute_select("SELECT * FROM UserDataSet WHERE TextID=%s", parameters=(userdatasetTextID,))
        results = self.DDG_db.execute_select("SELECT UserDataSetExperiment.* FROM UserDataSetExperiment INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%s", parameters=(userdatasetTextID,))
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

            existing_results = self.DDG_db.execute_select("SELECT * FROM Prediction WHERE PredictionSet=%s AND UserDataSetExperimentID=%s", parameters=(PredictionSet, r["ID"]))
            if len(existing_results) > 0:
                #colortext.warning('There already exist records for this UserDataSetExperimentID. You probably do not want to proceed. Skipping this entry.')
                continue

            PredictionID = self.addPrediction(r["ExperimentID"], r["ID"], PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = r["PDBFileID"], StoreOutput = StoreOutput, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = testonly)
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
        '''Abstract function.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')
        #returns None if no queued jobs exist or if the PredictionSet is halted otherwise return details necessary to run the job


    @job_execution
    def start_job(self, prediction_id, prediction_set):
        '''Abstract function.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')
         # sets the job status to 'active'


    @job_completion
    def fail_job(self, prediction_id, prediction_set, maxvmem, ddgtime):
        '''Abstract function.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')
        # sets a job to 'failed'. prediction_set is redundant but acts as a sanity check


    @job_completion
    def parse_prediction_scores(self, stdout):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        self._parse_ddg_monomer_scores_per_structure(stdout)


    @job_completion
    def store_scores(self, scores, prediction_set, prediction_id):
        '''Stores a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @job_completion
    def complete_job(self, prediction_id, prediction_set, scores, maxvmem, ddgtime):
        '''Abstract function.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')
        # sets a job to 'completed' and stores scores using self.store_scores. prediction_set is redundant but acts as a sanity check


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
    def determine_best_pair(self, prediction_id, score_method_id = 1):
        # Iterates over the (wildtype, mutant) pairs in the PredictionStructureScore table and returns the structure ID
        # for the pair with the lowest energy mutant
        # Note: There are multiple ways to select the best pair. For example, if multiple mutants have the same minimal total
        # score, we could have multiple wildtype structures to choose from. In this case, we choose a pair where the wildtype
        # structure has the minimal total score.
        lowest_mutant_score = self.DDG_db.execute_select('SELECT total FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" ORDER BY total LIMIT 1', parameters=(prediction_id, score_method_id))
        if lowest_mutant_score:
            lowest_mutant_score = lowest_mutant_score[0]['total']
            mutant_structure_ids = [r['StructureID'] for r in self.DDG_db.execute_select('SELECT StructureID FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" AND total=%s', parameters=(prediction_id, score_method_id, lowest_mutant_score))]
            if len(mutant_structure_ids) > 1:
                return self.DDG_db.execute_select(('SELECT StructureID FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="WildType" AND StructureID IN (' + ','.join(map(str, mutant_structure_ids)) + ') ORDER BY total LIMIT 1'), parameters=(prediction_id, score_method_id ))[0]['StructureID']
            else:
                return mutant_structure_ids[0]
        return None


    ################################################################################################
    ## Application layer
    ## These functions combine the database and prediction data with useful tools
    ################################################################################################


    #== PyMOL API ===========================================================


    @app_pymol
    def create_pymol_session_in_memory(self, prediction_id, task_number, pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol'):
        '''Create a PyMOL session for a pair of structures.'''

        # Retrieve and unzip results
        archive = self.get_job_data(prediction_id)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            # Get the name of the files from the zip
            wildtype_filename = os.path.join(str(prediction_id), 'repacked_wt_round_%d.pdb' % task_number)
            mutant_filename = None
            for filepath in sorted(zipped_content.namelist()):
                filename = os.path.split(filepath)[1]
                if filename.startswith('mut_') and filename.endswith('_round_%d.pdb' % task_number):
                    mutant_filename = os.path.join(str(prediction_id), filename)
                    break

            PyMOL_session = None
            file_list = zipped_content.namelist()

            # If both files exist in the zip, extract their contents in memory and create a PyMOL session pair (PSE, script)
            if (mutant_filename in file_list) and (wildtype_filename in file_list):
                wildtype_pdb = zipped_content.open(wildtype_filename, 'r').read()
                mutant_pdb = zipped_content.open(mutant_filename, 'U').read()
                chain_mapper = ScaffoldModelChainMapper.from_file_contents(wildtype_pdb, mutant_pdb)
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
        prediction_chains = self.get_pdb_chains_used_for_prediction_set(predictionset)

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


    def _get_prediction_table(self): return 'Prediction'
    def _get_prediction_type(self): return 'ProteinStability'
    def _get_prediction_type_description(self): return 'monomeric stability'


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
        from tools.bio.rcsb import parseFASTAs

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


























class AnalysisBreakdown(object):

    def __init__(self, amino_acids, pdb_details, predictions, analysis_datasets):
        self.amino_acids = amino_acids
        self.pdb_details = pdb_details

        # split the predictions over mutations to/from glycine and proline and other predictions
        GP = set(['G', 'P'])
        single_mutation_GP_predictions = {}
        single_mutation_no_GP_predictions = {}
        multiple_mutation_predictions = {}
        for p, details in predictions.iteritems():
            if not details.get('WTAA'):
                multiple_mutation_predictions[p] = details
            elif (details['WTAA'] in GP or details['MutantAA'] in GP):
                single_mutation_GP_predictions[p] = details
            else:
                single_mutation_no_GP_predictions[p] = details

        self.predictions = predictions
        self.single_mutation_GP_predictions = single_mutation_GP_predictions
        self.single_mutation_no_GP_predictions = single_mutation_no_GP_predictions
        self.multiple_mutation_predictions = multiple_mutation_predictions
        print('%d total predictions: %d Single No GP, %d Single GP, %d multiple' % (len(predictions), len(single_mutation_no_GP_predictions), len(single_mutation_GP_predictions), len(multiple_mutation_predictions)))
        self.analysis_datasets = analysis_datasets


    def analyze_subset_all(self, analysis_subset, scoring_method, prediction_details_map = {}):
        colortext.message('ANALYZING SUBSET analyze_subset_all of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.predictions, prediction_details_map)

    def analyze_subset_single_no_GP(self, analysis_subset, scoring_method):
        colortext.message('ANALYZING SUBSET single_mutation_no_GP_predictions of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.single_mutation_no_GP_predictions)

    def analyze_subset_single_GP(self, analysis_subset, scoring_method):
        colortext.message('ANALYZING SUBSET single_mutation_GP_predictions of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.single_mutation_GP_predictions)

    def analyze_subset_multiple(self, analysis_subset, scoring_method):
        colortext.message('ANALYZING SUBSET multiple_mutation_predictions of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.multiple_mutation_predictions)

    def _analyze_subset_sub(self, analysis_subset, scoring_method, predictions, prediction_details_map = {}):
        'Analyzes a subset using the main datapoints.'

        analysis_dataset = self.analysis_datasets[analysis_subset]
        xvalues = []
        yvalues = []
        print('ID,Experimental,Predicted')

        # for benchmarks paper
        ddgapio = MonomericStabilityDDGInterface().DDG_db

        for prediction_id, details in sorted(predictions.iteritems()):
            if prediction_id in analysis_dataset:
                ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                predicted_score = details[scoring_method]
                if predicted_score != None:

                    # for the Benchmarks paper...
                    if False:
                        mutations = ddgapio.execute_select('''
SELECT ExperimentMutation.*
FROM Prediction
INNER JOIN ExperimentMutation ON ExperimentMutation.ExperimentID=Prediction.ExperimentID
WHERE Prediction.ID = %s''', parameters=(prediction_id,))
                        if len(mutations) == 1:
                            if abs(ExperimentalDDG - predicted_score) < 0.5:
                                if abs(ExperimentalDDG) > 2:
                                    if (mutations[0]['WildTypeAA'] == 'D'):
                                        print(mutations[0]['WildTypeAA'], mutations[0]['MutantAA'], ExperimentalDDG, predicted_score, prediction_id)


                    xvalues.append(ExperimentalDDG)
                    yvalues.append(predicted_score)
                    if prediction_details_map:
                        assert(prediction_details_map.get(int(prediction_id)))
                        #print('%f,%f,%s,%s' % (ExperimentalDDG, predicted_score, prediction_id, ','.join(map(str, prediction_details_map[int(prediction_id)]))))
                    else:
                        pass
                        #print('%f,%f,%s,%s' % (ExperimentalDDG, predicted_score, prediction_id))

        print('*' * 30)
        print(min(xvalues), max(xvalues), min(yvalues), max(yvalues))
        colortext.message('Analyzing %d values for dataset %s using scoring method %s.' % (len(xvalues), analysis_subset, scoring_method))
        stats = None
        if (len(xvalues) >= 8):
            stats = get_xy_dataset_correlations(xvalues, yvalues)
            #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
            colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
            colortext.warning('''
                              <span class="DDG_correlation_score ">%0.2f</span> /
                              <span class="DDG_stability_classification_score ">%0.2f</span> /
                              <span class="DDG_MAE_score">%0.2f</span>''' % (stats['pearsonr'][0], stats['fraction_correct'], stats['MAE']))
            pprint.pprint(stats)
        else:
            colortext.warning('Not enough data.')
        return stats


    def analyze_subset_by_specific_resolutions(self, analysis_subset, scoring_method, bins = [1.5, 2.0, 2.5]):
        ''' Analyzes a subset using specific PDB resolution bins.
            The bins argument defines the resolution bins.
            The first bin is <x where x is the smallest value in bins.
            The last bin is >x where x is the largest value in bins.
        '''

        bins = sorted(set(bins))
        assert(bins[0] > 0)
        limits = [0] + bins + [max(bins[-1], 100)]
        xvalues = {None : []}
        yvalues = {None : []}
        for x in range(0, len(limits) - 1):
            xvalues[(limits[x], limits[x + 1])] = []
            yvalues[(limits[x], limits[x + 1])] = []

        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    resolution = self.pdb_details[details['pPDB']]['Resolution']
                    if resolution:
                        for x in range(0, len(limits) - 1):
                            if limits[x] <= resolution < limits[x + 1]:
                                xvalues[(limits[x], limits[x + 1])].append(ExperimentalDDG)
                                yvalues[(limits[x], limits[x + 1])].append(predicted_score)
                                break
                    else:
                        xvalues[None].append(ExperimentalDDG)
                        yvalues[None].append(predicted_score)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for bin %s.' % (len(xvalues[k]), str(k)))
                results[k] = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                #pprint.pprint(results[k])
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results

    def analyze_subset_by_binned_resolutions(self, analysis_subset, scoring_method, num_bins = 9):
        ''' Analyzes a subset using PDB resolution bins. This function attempts to break up the result set into num_bins
            bins of somewhat equal size.
            Additionally, there is a special bin for PDB files with null resolution.'''

        assert(num_bins > 1)

        count = 0
        xyvalues = {None: []}
        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    resolution = self.pdb_details[details['pPDB']]['Resolution']
                    xyvalues[resolution] = xyvalues.get(resolution, [])
                    xyvalues[resolution].append((ExperimentalDDG, predicted_score))
                    count += 1

        # determine the ideal number per bin, ignoring the bin with null resolution
        ideal_per_bin = (count - len(xyvalues.get(None))) / num_bins

        xvalues = {None : [p[0] for p in xyvalues[None]]}
        yvalues = {None : [p[1] for p in xyvalues[None]]}
        current_bin = []
        bin_start = 0
        for bin_end in sorted(xyvalues.keys()):
            if bin_end != None:
                current_bin += xyvalues[bin_end]
                if (len(current_bin) > ideal_per_bin) or (bin_end == sorted(xyvalues.keys())[-1]):
                    xvalues[(bin_start, bin_end)] = [p[0] for p in current_bin]
                    yvalues[(bin_start, bin_end)] = [p[1] for p in current_bin]
                    bin_start = bin_end
                    current_bin = []
        assert(sum([len(xvalues[k]) for k in xvalues]) == count)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for bin %s.' % (len(xvalues[k]), str(k)))
                pprint.pprint(get_xy_dataset_correlations(xvalues[k], yvalues[k]))
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results


    def analyze_subset_by_binned_chain_length(self, analysis_subset, scoring_method, num_bins = 9):
        ''' Analyzes a subset using PDB resolution bins. This function attempts to break up the result set into num_bins
            bins of somewhat equal size.
            Additionally, there is a special bin for PDB files with null resolution.'''

        assert(num_bins > 1)

        pdb_details = self.pdb_details
        count = 0
        xyvalues = {}
        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    chain = details['Chain']
                    ppdb = details['pPDB']
                    chain_length = pdb_details[ppdb]['chains'][chain]
                    xyvalues[chain_length] = xyvalues.get(chain_length, [])
                    xyvalues[chain_length].append((ExperimentalDDG, predicted_score))
                    count += 1

        # determine the ideal number per bin, ignoring the bin with null resolution
        ideal_per_bin = (count) / num_bins

        xvalues = {}#None : [p[0] for p in xyvalues[None]]}
        yvalues = {}#sNone : [p[1] for p in xyvalues[None]]}
        current_bin = []
        bin_start = max(xyvalues.keys())
        print(sorted(xyvalues.keys()))
        print(sorted(xyvalues.keys(), reverse = True))
        for bin_end in sorted(xyvalues.keys(), reverse = True):
            current_bin += xyvalues[bin_end]
            if (len(current_bin) > ideal_per_bin) or (bin_end == sorted(xyvalues.keys(), reverse = True)[-1]):
                xvalues[(bin_start, bin_end)] = [p[0] for p in current_bin]
                yvalues[(bin_start, bin_end)] = [p[1] for p in current_bin]
                bin_start = bin_end
                current_bin = []
        print(sum([len(xvalues[k]) for k in xvalues]), count)
        assert(sum([len(xvalues[k]) for k in xvalues]) == count)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for bin %s.' % (len(xvalues[k]), str(k)))
                pprint.pprint(get_xy_dataset_correlations(xvalues[k], yvalues[k]))
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results

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

        #print(self.pdb_details[details['pPDB']])
        #        print(self.pdb_details[details['pPDB']]['chains'][details['Chain']])
        #        #{u'XRay': True, u'chains': {u'A': 680}, u'Technique': u'X-RAY DIFFRACTION', u'Resolution': 2.4, u'TM': 1}


    def analyze_subset_by_exposure(self, analysis_subset, scoring_method, cut_off = 0.25):
        ''' Analyzes a subset using exposure of the wildtype residue as calculated by DSSP.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.
        '''

        assert(0 <= cut_off <= 1.0)
        xvalues = {None : [], 'B' : [], 'E' : []}
        yvalues = {None : [], 'B' : [], 'E' : []}

        main_subsets = [self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        all_results = dict.fromkeys(main_subset_names)
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {'B' : [], 'E' : []}
                yvalues = {'B' : [], 'E' : []}

                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):
                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            exposure = details.get('Exposure')
                            if exposure != None:
                                if exposure <= cut_off:
                                    xvalues['B'].append(ExperimentalDDG)
                                    yvalues['B'].append(predicted_score)
                                else:
                                    xvalues['E'].append(ExperimentalDDG)
                                    yvalues['E'].append(predicted_score)
                            else:
                                print(details)
                                xvalues[None].append(ExperimentalDDG)
                                yvalues[None].append(predicted_score)

                results = {}
                colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for exposure type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
                        #pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results
        return all_results


    def analyze_subset_by_wildtype_charge(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        xvalues = {None : [], 'C' : [], 'P' : [], 'H' : []}
        yvalues = {None : [], 'C' : [], 'P' : [], 'H' : []}

        amino_acids = self.amino_acids
        CAA = [aa for aa in amino_acids if amino_acids[aa]['Polarity'] == 'C']
        PAA = [aa for aa in amino_acids if amino_acids[aa]['Polarity'] == 'P']
        HAA = [aa for aa in amino_acids if amino_acids[aa]['Polarity'] == 'H']

        main_subsets = [self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        all_results = dict.fromkeys(main_subset_names)
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {'Change' : [], 'Polar/Charged' : [], 'Hydrophobic/Non-polar' : []}
                yvalues = {'Change' : [], 'Polar/Charged' : [], 'Hydrophobic/Non-polar' : []}

                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):

                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            wtaa = details.get('WTAA')
                            if wtaa:
                                mutaa = details.get('MutantAA')
                                if ((wtaa in CAA or wtaa in PAA) and (mutaa in HAA)) or ((mutaa in CAA or mutaa in PAA) and (wtaa in HAA)):
                                    # change in charge
                                    xvalues['Change'].append(ExperimentalDDG)
                                    yvalues['Change'].append(predicted_score)
                                elif (wtaa in CAA or wtaa in PAA) and (mutaa in CAA or mutaa in PAA):
                                    xvalues['Polar/Charged'].append(ExperimentalDDG)
                                    yvalues['Polar/Charged'].append(predicted_score)
                                elif (wtaa in HAA) and (mutaa in HAA):
                                    xvalues['Hydrophobic/Non-polar'].append(ExperimentalDDG)
                                    yvalues['Hydrophobic/Non-polar'].append(predicted_score)
                                else:
                                     raise Exception('Should not reach here.')

                results = {}
                colortext.message('Analyzing %s=%d values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), sum([len(xvalues[k]) for k in sorted(xvalues.keys())]), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for polarity type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
                        #pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results
        return all_results


    def analyze_subset_by_aromaticity(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        xvalues = {None : [], '-' : [], 'L' : [], 'R' : []}
        yvalues = {None : [], '-' : [], 'L' : [], 'R' : []}

        amino_acids = self.amino_acids
        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    wtaa = details.get('WTAA')
                    if wtaa:
                        aromaticity = amino_acids[wtaa]['Aromaticity']
                        xvalues[aromaticity].append(ExperimentalDDG)
                        yvalues[aromaticity].append(predicted_score)
                    else:
                        xvalues[None].append(ExperimentalDDG)
                        yvalues[None].append(predicted_score)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for aromaticity type %s.' % (len(xvalues[k]), str(k)))
                results[k] = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                #pprint.pprint(results[k])
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results


    def analyze_subset_by_mutation_size(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        #xvalues = {None : [], 'XX' : [], 'SL' : [], 'LS' : []}
        #yvalues = {None : [], 'XX' : [], 'SL' : [], 'LS' : []}

        amino_acids = self.amino_acids
        SAA = [aa for aa in amino_acids if amino_acids[aa]['Size'] == 'small']
        LAA = [aa for aa in amino_acids if amino_acids[aa]['Size'] == 'large']

        amino_acid_volumes = {}
        for aa, details in amino_acids.iteritems():
            amino_acid_volumes[aa] = details['van_der_Waals_volume']
        assert(len(amino_acid_volumes) == 20)

        main_subsets = [self.predictions, self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['all_mutations', 'single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        main_subsets = [self.predictions]
        main_subset_names = ['all_mutations']
        all_results = dict.fromkeys(main_subset_names)
        failed_cases = 0
        non_cases = set()
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {'XX' : [], 'SL' : [], 'LS' : [], 'Failed' : []}
                yvalues = {'XX' : [], 'SL' : [], 'LS' : [], 'Failed' : []}

                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):
                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            wtaa = details.get('WTAA')
                            if wtaa:
                                mutaa = details.get('MutantAA')

                                if details.get('MutationIsReversed') != None:
                                    # todo: this is not currently an issue but it will be once we include reverse mutations
                                    if details['MutationIsReversed']:
                                        # Note: For reverse mutations, we need to switch the order since we only store the forward mutation
                                        wtaa, mutaa = mutaa, wtaa

                                if wtaa == mutaa:
                                    colortext.warning('Error in analysis: Record mutating %s to %s in Prediction #%s.' % (wtaa, mutaa, prediction_id))
                                    error = True
                                elif amino_acid_volumes[wtaa] < amino_acid_volumes[mutaa]:
                                    xvalues['SL'].append(ExperimentalDDG)
                                    yvalues['SL'].append(predicted_score)
                                elif amino_acid_volumes[wtaa] > amino_acid_volumes[mutaa]:
                                    xvalues['LS'].append(ExperimentalDDG)
                                    yvalues['LS'].append(predicted_score)
                                else:
                                    assert(amino_acid_volumes[wtaa] == amino_acid_volumes[mutaa])
                                    xvalues['XX'].append(ExperimentalDDG)
                                    yvalues['XX'].append(predicted_score)

                                if False:
                                    if wtaa in SAA and mutaa in LAA:
                                        xvalues['SL'].append(ExperimentalDDG)
                                        yvalues['SL'].append(predicted_score)
                                    elif wtaa in LAA and mutaa in SAA:
                                        xvalues['LS'].append(ExperimentalDDG)
                                        yvalues['LS'].append(predicted_score)
                                    else:
                                        non_cases.add('%s->%s' % (wtaa, mutaa))
                                        xvalues['XX'].append(ExperimentalDDG)
                                        yvalues['XX'].append(predicted_score)
                            else:
                                failed_cases += 1
                                #raise Exception('Should not reach here.')
                                xvalues['Failed'].append(ExperimentalDDG)
                                yvalues['Failed'].append(predicted_score)

                results = {}
                colortext.message('Analyzing %s=%d values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), sum([len(xvalues[k]) for k in sorted(xvalues.keys())]), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for mutation size type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE']))

                        colortext.warning('''
                                          <span class="DDG_correlation_score ">%0.2f</span> /
                                          <span class="DDG_stability_classification_score ">%0.2f</span> /
                                          <span class="DDG_MAE_score">%0.2f</span>''' % (stats['pearsonr'][0], stats['fraction_correct'], stats['MAE']))

                        pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results

        print(sorted(non_cases))
        return all_results


    def analyze_subset_by_secondary_structure(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        xvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}
        yvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}

        main_subsets = [self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        all_results = dict.fromkeys(main_subset_names)
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}
                yvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}
                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):
                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            dssp = details.get('DSSP')
                            xvalues[dssp].append(ExperimentalDDG)
                            yvalues[dssp].append(predicted_score)

                results = {}
                colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for secondary structure type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
                        #pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results
        return all_results


#class MonomericStabilityDDGUserInterface(MonomericStabilityDDGInterface): pass

#import inspect
#for m in inspect.getmembers(MonomericStabilityDDGInterface, predicate=inspect.ismethod):
#    print(m)
#    MonomericStabilityDDGUserInterface.__dict__[m[0]] = m[1]
