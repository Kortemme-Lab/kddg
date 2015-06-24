#!/usr/bin/python2.4
# encoding: utf-8
"""
ppi_api.py
High-level functions for interacting with the protein-protein interaction sections of the ddG database.

Classes:
BindingAffinityDDGInterface - an class used to interface with the database. Call get_interface to get a user API based on this class.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

import pprint
from io import BytesIO
import os
import zipfile

from api_layers import *
from db_api import ddG
from tools import colortext
from tools.bio.alignment import ScaffoldModelChainMapper
from tools.fs.fsio import read_file


def get_interface(passwd, username = 'kortemmelab'):
    '''This is the function that should be used to get a BindingAffinityDDGInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(BindingAffinityDDGInterface, passwd = passwd, username = username)


class BindingAffinityDDGInterface(ddG):
    '''This is the internal API class that should be NOT used to interface with the database.'''


    def __init__(self, passwd = None, username = 'kortemmelab'):
        super(BindingAffinityDDGInterface, self).__init__(passwd = passwd, username = username)
        self.prediction_data_path = self.DDG_db.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionPPIDataPath"')[0]['Value']


    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================


    @informational_pdb
    def get_pdb_chains_for_prediction(self, prediction_id):
        # look up the complex associated with the dataset record for the list of chains
        raise Exception('This needs to be implemented.')


    ##### Public API: Rosetta-related functions


    @job_input
    def create_resfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    @job_input
    def create_mutfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    ###########################################################################################
    ## Prediction layer
    ##
    ## This part of the API is responsible for inserting prediction jobs in the database via
    ## the trickle-down proteomics paradigm.
    ###########################################################################################


    #== Job creation/management API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via the
    # trickle-down proteomics paradigm.


    #   PredictionSet interface

    @job_creator
    def add_prediction_set(self, prediction_set_id, halted = True, priority = 5, batch_size = 40, allow_existing_prediction_set = False):
        return super(BindingAffinityDDGInterface, self).add_prediction_set(prediction_set_id, halted = halted, priority = priority, batch_size = batch_size, allow_existing_prediction_set = allow_existing_prediction_set, contains_protein_stability_predictions = False, contains_binding_affinity_predictions = True)


    @job_creator
    def add_job_command_lines(self, *args, **kwargs):
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_creator
    def add_prediction_run(self, prediction_set_id, user_dataset_name, protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, input_files = {}, quiet = False, test_only = False, only_single_mutations = False, short_run = False):
        '''Adds all jobs corresponding to a user dataset e.g. add_prediction_run("my first run", "AllBindingAffinityData", tagged_subset = "ZEMu").
           If keep_hetatm_lines is False then all HETATM records for the PDB prediction chains will be removed. Otherwise, they are kept.
           input_files is a global parameter for the run which is generally empty. Any files added here will be associated to all predictions in the run.

           Returns False if no predictions were added to the run else return True if all predictions (and there were some) were added to the run.'''

        # Check preconditions
        assert(not(input_files)) # todo: do something with input_files when we use that here - call self._add_file_content, associate the filenames with the FileContent IDs, and pass that dict to add_job which will create PredictionPPIFile records
        assert(only_single_mutations == False) # todo: support this later? it may make more sense to just define new UserDataSets
        self._add_prediction_run_preconditions(prediction_set_id, user_dataset_name, tagged_subset)

        # Get the list of user dataset experiment records
        user_dataset_experiments = self.get_user_dataset_experiment_ids(user_dataset_name, tagged_subset = tagged_subset)
        assert(set([u['IsComplex'] for u in user_dataset_experiments]) == set([1,]))
        if not user_dataset_experiments:
            return False

        # Count the number of individual PDB files
        pdb_file_ids = set([u['PDBFileID'] for u in user_dataset_experiments])
        if not quiet:
            tagged_subset_str = ''
            if tagged_subset:
                tagged_subset_str = 'subset "%s" of ' % tagged_subset
            colortext.message('Adding %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))

        import sys; sys.exit(0)

        #SetNumber
        #PPMutagenesisID
        #PPComplexID
        #PDBFileID

        count = 0
        showprogress = not(quiet) and len(results) > 300
        if showprogress:
            print("|" + ("*" * (int(len(results)/100)-2)) + "|")
        for r in results:

            existing_results = self.DDG_db.execute_select("SELECT * FROM Prediction WHERE PredictionSet=%s AND UserDataSetExperimentID=%s", parameters=(prediction_set_id, r["ID"]))
            if len(existing_results) > 0:
                #colortext.warning('There already exist records for this UserDataSetExperimentID. You probably do not want to proceed. Skipping this entry.')
                continue

            PredictionID = self.add_job(r["ExperimentID"], r["ID"], prediction_set_id, protocol_id, keep_hetatm_lines, PDB_ID = r["PDBFileID"], ReverseMutation = False, input_files = input_files, test_only = test_only)
            count += 1
            if showprogress:
                if count > 100:
                    colortext.write(".", "cyan", flush = True)
                    count = 0
            if short_run and count > 4:
                break
        print("")
        return(True)


    @job_creator
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set):
        raise Exception('not implemented yet')
        #assert(existing_prediction_set exists and has records)
        #assert(new_prediction_set is empty)
        #for each prediction record, add the record and all associated predictionfile records,


    @job_creator
    def add_job(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, keep_hetatm_lines, PDB_ID = None, ReverseMutation = False, input_files = {}, test_only = False, strip_other_chains = True):
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

        # todo: do something with input_files when we use that here - see add_prediction_run
        assert(not(input_files))

        ExtraParameters = {}
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
            PredictionFieldNames.Status 			: "queued",
            PredictionFieldNames.ExtraParameters	: ExtraParameters,
        }
        if not test_only:
            self.DDG_db.insertDict('Prediction', params)

            # Add cryptID string
            predictionID = self.DDG_db.getLastRowID()
            entryDate = self.DDG_db.execute_select("SELECT EntryDate FROM Prediction WHERE ID=%s", parameters = (predictionID,))[0]["EntryDate"]
            rdmstring = ''.join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16))
            cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
            cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
            entryDate = self.DDG_db.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))
            return predictionID


    @job_execution
    def get_max_number_of_cluster_jobs(self, prediction_set_id, priority):
        return self.DDG_db.execute_select('SELECT Value FROM _DBCONSTANTS WHERE VariableName="MaxStabilityClusterJobs"')['Value']


    @job_completion
    def parse_prediction_scores(self, stdout):
        '''Returns a list of dicts suitable for database storage e.g. PredictionPPIStructureScore records.'''
        raise Exception('not implemented yet')


    @job_completion
    def store_scores(self, scores, prediction_set, prediction_id):
        '''Stores a list of dicts suitable for database storage e.g. PredictionPPIStructureScore records.'''
        raise Exception('not implemented yet')



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


    def _get_prediction_table(self): return 'PredictionPPI'
    def _get_prediction_type(self): return 'BindingAffinity'
    def _get_prediction_dataset_type(self): return 'Binding affinity'
    def _get_prediction_type_description(self): return 'binding affinity'
    def _get_user_dataset_experiment_table(self): return 'UserPPDataSetExperiment'
    def _get_user_dataset_experiment_tag_table(self): return 'UserPPDataSetExperimentTag'

    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================

    @informational_pdb
    def _get_pdb_chains_used_for_prediction_set(self, prediction_set):
        raise Exception('not implemented yet')
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



