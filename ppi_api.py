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
from tools.bio.pdb import PDB
from tools.bio.alignment import ScaffoldModelChainMapper
from tools.bio.basics import ChainMutation
from tools.fs.fsio import read_file
from tools.rosetta.input_files import Mutfile

def get_interface(passwd, username = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a BindingAffinityDDGInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(BindingAffinityDDGInterface, passwd = passwd, username = username, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)


class BindingAffinityDDGInterface(ddG):
    '''This is the internal API class that should be NOT used to interface with the database.'''


    def __init__(self, passwd = None, username = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None):
        super(BindingAffinityDDGInterface, self).__init__(passwd = passwd, username = username, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)
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
    def add_prediction_run(self, prediction_set_id, user_dataset_name, extra_rosetta_command_flags = None, protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, input_files = {}, quiet = False, test_only = False, only_single_mutations = False, short_run = False, test_run_first = True):
        '''Adds all jobs corresponding to a user dataset e.g. add_prediction_run("my first run", "AllBindingAffinityData", tagged_subset = "ZEMu").
           If keep_hetatm_lines is False then all HETATM records for the PDB prediction chains will be removed. Otherwise, they are kept.
           input_files is a global parameter for the run which is generally empty. Any files added here will be associated to all predictions in the run.

           The extra_rosetta_command_flags parameter e.g. "-ignore_zero_occupancy false" is used to determine the mapping
           from PDB to Rosetta numbering. These flags should correspond to those used in the protocol otherwise errors could occur.

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
        tagged_subset_str = ''
        if not quiet:
            if tagged_subset:
                tagged_subset_str = 'subset "%s" of ' % tagged_subset

        # Test all predictions before creating records
        if test_only or test_run_first:
            if not quiet:
                colortext.message('Testing %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
            # Progress counter setup
            count, records_per_dot = 0, 50
            showprogress = not(quiet) and len(user_dataset_experiments) > 300
            if showprogress: print("|" + ("*" * (int(len(user_dataset_experiments)/records_per_dot)-2)) + "|")
            for ude in user_dataset_experiments:
                # If the mutagenesis already exists in the prediction set, do not test it again
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
                if len(existing_results) > 0: continue

                # Test the prediction setup
                print('%d/%d' % (count, len(user_dataset_experiments)))
                prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = True)
                # Progress counter
                count += 1
                if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
                if short_run and count > 4: break
            if not quiet: print('')
            return True

        if test_only:
            return

        return

        # Progress counter setup
        if not quiet:
            colortext.message('Adding %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
        count, records_per_dot = 0, 50
        showprogress = not(quiet) and len(user_dataset_experiments) > 300
        if showprogress: print("|" + ("*" * (int(len(user_dataset_experiments)/records_per_dot)-2)) + "|")

        # Add the individual predictions
        for ude in user_dataset_experiments:

            # If the mutagenesis already exists in the prediction set, do not add it again
            existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
            if len(existing_results) > 0:
                continue

            # Add the prediction
            prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only)

            # Progress counter
            count += 1
            if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
            if short_run and count > 4: break

        if not quiet: print('')
        return True


    @job_creator
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set):
        raise Exception('not implemented yet')
        #assert(existing_prediction_set exists and has records)
        #assert(new_prediction_set is empty)
        #for each prediction record, add the record and all associated predictionfile records,


    @job_creator
    def add_job_by_user_dataset_record(self, prediction_set_id, user_dataset_name, user_dataset_experiment_id, protocol_id, extra_rosetta_command_flags = None, keep_hetatm_lines = False, input_files = {}, test_only = False):
        '''Add a prediction job based on a user dataset record. This is typically called during add_prediction_run rather than directly by the user.
           user_dataset_name is implied by user_dataset_experiment_id but we include it for sanity checking errors in data-entry.

           The extra_rosetta_command_flags variable is used to add additional flags e.g. "-ignore_zero_occupancy false". These should be added if they are used in the protocol.'''

        try:
            user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
        except:
            raise colortext.Exception('The user dataset "%s" does not exist for this API.' % user_dataset_name)

        ude = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(user_dataset_experiment_id, user_dataset_id))
        if not len(ude) == 1:
            raise colortext.Exception('User dataset experiment %d does not exist for/correspond to this user dataset.' % user_dataset_experiment_id)
        ude = ude[0]
        #colortext.message(pprint.pformat(ude))
        return self._add_job(prediction_set_id, protocol_id, ude['PPMutagenesisID'], ude['PPComplexID'], ude['PDBFileID'], ude['SetNumber'], extra_rosetta_command_flags = extra_rosetta_command_flags, user_dataset_experiment_id = user_dataset_experiment_id, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only)


    @job_creator
    def add_job(self, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, keep_hetatm_lines = False, input_files = {}, test_only = False):
        '''This function inserts a prediction into the database.
            The parameters define:
                - the prediction set id used to group this prediction with other predictions for analysis;
                - the protocol to be used to run the prediction;
                - the set of mutations and PDB complex associated with the mutagenesis experiment;
                - whether HETATM lines are to be kept or not.
                - additional Rosetta flags e.g. "-ignore_zero_occupancy false" used to determine the mapping from PDB to Rosetta numbering. These flags should correspond to those used in the protocol otherwise errors could occur.
            We strip the PDB based on the chains defined by the complex and keep_hetatm_lines and store the PDB in the database.
            Next, the mapping from Rosetta numbering to PDB numbering is determined and stored in the database.
            Then, the appropriate input files e.g. resfiles or mutfiles are generated and stored in the database.
            Finally, we add the prediction record and associate it with the generated files.'''
        return self._add_job(prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only)


    @informational_pdb
    def get_chain_sets_for_mutatagenesis(self, mutagenesis_id, complex_id = None):
        '''Gets a list of possibilities for the associated complex and calls get_chains_for_mutatagenesis on each.
             e.g. returns {('1KI1', 0) : {'L' : ['A','B'], 'R' : ['C']}, ('12AB', 2) : {'L' : ['L','H'], 'R' : ['A']}, ...}
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.'''

        pp_mutagenesis = self.DDG_db.execute_select("SELECT * FROM PPMutagenesis WHERE ID=%s", parameters = (mutagenesis_id,))
        # Sanity checks
        assert(len(pp_mutagenesis) == 1)
        if complex_id:
            assert(pp_mutagenesis[0]['PPComplexID'] == complex_id)
        else:
            complex_id = pp_mutagenesis[0]['PPComplexID']

        d = {}
        for pdb_set in self.DDG_db.execute_select("SELECT * FROM PPIPDBSet WHERE PPComplexID=%s AND IsComplex=1", parameters = (complex_id,)):
            pdb_set_number = pdb_set['SetNumber']
            pdb_file_ids = self.DDG_db.execute_select("SELECT DISTINCT PDBFileID FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND SetNumber=%s", parameters = (complex_id, pdb_set_number))
            assert(len(pdb_file_ids) == 1)
            pdb_file_id = pdb_file_ids[0]['PDBFileID']
            d[(pdb_file_id, pdb_set_number)] = self.get_chains_for_mutatagenesis(mutagenesis_id, pdb_file_id, pdb_set_number)
        return d


    @informational_pdb
    def get_chains_for_mutatagenesis(self, mutagenesis_id, pdb_file_id, pdb_set_number, complex_id = None):
        '''Returns a dictionary mapping 'L' to the list of left chains and 'R' to the list of right chains.
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.
        '''

        pp_mutagenesis = self.DDG_db.execute_select("SELECT * FROM PPMutagenesis WHERE ID=%s", parameters = (mutagenesis_id,))

        # Sanity checks
        assert(len(pp_mutagenesis) == 1)
        if complex_id:
            assert(pp_mutagenesis[0]['PPComplexID'] == complex_id)
            pdb_set = self.DDG_db.execute_select("SELECT * FROM PPIPDBSet WHERE PPComplexID=%s AND SetNumber=%s", parameters = (complex_id, pdb_set_number))
            assert(len(pdb_set) == 1 and pdb_set[0]['IsComplex'] == 1) # complex structure check
        else:
            complex_id = pp_mutagenesis[0]['PPComplexID']

        complex_chains = dict(L = [], R = [])
        crecords = self.DDG_db.execute_select("SELECT * FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND SetNumber=%s ORDER BY ChainIndex", parameters = (complex_id, pdb_set_number))
        for c in crecords:
            assert(c['PDBFileID'] == pdb_file_id) # complex structure check
            complex_chains[c['Side']].append(c['Chain'])
        assert(complex_chains['L'] and complex_chains['R'])
        assert(len(set(complex_chains['L']).intersection(set(complex_chains['R']))) == 0) # in one unbound case, the same chain appears twice on one side (2CLR_DE|1CD8_AA, may be an error since this was published as 1CD8_AB but 1CD8 has no chain B) but it seems reasonable to assume that a chain should only appear on one side
        return complex_chains


    def _add_job(self, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, user_dataset_experiment_id = None, keep_hetatm_lines = False, input_files = {}, test_only = False):
        '''This is the internal function which adds a prediction job to the database. We distinguish it from add_job as
           prediction jobs added using that function should have no associated user dataset experiment ID.

           The extra_rosetta_command_flags variable is used to add additional flags e.g. "-ignore_zero_occupancy false". These should be added if they are used in the protocol.
           '''

        # todo: do something with input_files when we use that here - see add_prediction_run
        assert(not(input_files))

        if not(self.rosetta_scripts_path and self.rosetta_database_path):
            raise Exception('The rosetta_scripts_path and rosetta_database_path API variables need to be set when adding predictions as we use the Features Reporter to create the mapping between PDB residues and Rosetta residues.')
        else:
            if not(os.path.exists(self.rosetta_scripts_path)):
                raise Exception('The path "%s" to the RosettaScripts executable does not exist.' % self.rosetta_scripts_path)
            if not(os.path.exists(self.rosetta_database_path)):
                raise Exception('The path "%s" to the Rosetta database does not exist.' % self.rosetta_database_path)

        # Information for debugging
        pp_complex = self.DDG_db.execute_select("SELECT * FROM PPComplex WHERE ID=%s", parameters = (pp_complex_id,))

        # Determine the list of PDB chains that will be kept
        pdb_chains = self.get_chains_for_mutatagenesis(pp_mutagenesis_id, pdb_file_id, pp_complex_pdb_set_number, complex_id = pp_complex_id)
        pdb_chains_to_keep = set(pdb_chains['L'] + pdb_chains['R'])

        pdb_residues_to_rosetta_cache = {}
        cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, self.rosetta_database_path, extra_command_flags)

        # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
        pdb_file = self.DDG_db.execute_select("SELECT * FROM PDBFile WHERE ID=%s", parameters = (pdb_file_id,))
        p = PDB(pdb_file[0]['Content'])
        p.strip_to_chains(list(pdb_chains_to_keep))
        if not keep_hetatm_lines:
            p.strip_HETATMs()
        stripped_p = PDB('\n'.join(p.lines))

        # Check for CSE and MSE
        try:
            if 'CSE' in p.residue_types:
                raise Exception('This case contains a CSE residue which may (or may not) cause an issue.')
            elif 'MSE' in p.residue_types:
                raise Exception('This case contains an MSE residue which may (or may not) cause an issue.')
                # It looks like MSE (and CSE?) may now be handled - https://www.rosettacommons.org/content/pdb-files-rosetta-format
        except Exception, e:
            colortext.error('%s: %s, chains %s' % (str(e), str(stripped_p.pdb_id), str(pdb_chains_to_keep)))

        # Assert that there are no empty sequences
        assert(sorted(stripped_p.atom_sequences.keys()) == sorted(pdb_chains_to_keep))
        for chain_id, sequence in stripped_p.atom_sequences.iteritems():
            assert(len(sequence) > 0)

        # Get the PDB mutations and check that they make sense in the context of the stripped PDB file
        # Note: the schema assumes that at most one set of mutations can be specified per PDB file per complex per mutagenesis. We may want to relax that in future by adding the SetNumber to the PPMutagenesisPDBMutation table
        complex_mutations = self.DDG_db.execute_select('SELECT * FROM PPMutagenesisMutation WHERE PPMutagenesisID=%s', parameters=(pp_mutagenesis_id,))
        pdb_complex_mutations = self.DDG_db.execute_select('SELECT * FROM PPMutagenesisPDBMutation WHERE PPMutagenesisID=%s AND PPComplexID=%s AND PDBFileID=%s', parameters=(pp_mutagenesis_id, pp_complex_id, pdb_file_id))
        assert(len(complex_mutations) == len(pdb_complex_mutations))
        mutations = [ChainMutation(m['WildTypeAA'], m['ResidueID'], m['MutantAA'], Chain = m['Chain']) for m in pdb_complex_mutations]
        try:
            stripped_p.validate_mutations(mutations)
        except Exception, e:
            colortext.error('%s: %s' % (str(e), str(mutations)))
            #colortext.warning('PPMutagenesisID=%d, ComplexID=%d, PDBFileID=%s, SetNumber=%d, UserDatasetExperimentID=%d' % (pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, user_dataset_experiment_id))
            #colortext.warning('SKEMPI record: %s' % self.DDG_db.execute_select('SELECT * FROM PPMutagenesis WHERE ID=%s', parameters=(pp_mutagenesis_id,))[0]['SKEMPI_KEY'])
            #colortext.warning('PDB chains to keep: %s' % str(pdb_chains_to_keep))
            #colortext.warning('PPIPDBPartnerChain records: %s' % pprint.pformat(self.DDG_db.execute_select('SELECT PPIPDBPartnerChain.* FROM PPIPDBPartnerChain INNER JOIN PPIPDBSet ON PPIPDBSet.PPComplexID=PPIPDBPartnerChain.PPComplexID AND PPIPDBSet.SetNumber=PPIPDBPartnerChain.SetNumber WHERE PPIPDBPartnerChain.PPComplexID=%s AND IsComplex=1 ORDER BY PPIPDBPartnerChain.SetNumber, PPIPDBPartnerChain.ChainIndex', parameters=(pp_complex_id,))))

        # Determine the mapping from the stripped PDB to Rosetta numbering
        # Note: we assume that this stripped PDB will be the input to the Rosetta protocol and that
        stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, self.rosetta_database_path, extra_command_flags = extra_rosetta_command_flags)
        rosetta_mutations = stripped_p.map_pdb_residues_to_rosetta_residues(mutations)

        # Make mutfile
        mf = Mutfile.from_mutagenesis(rosetta_mutations)
        colortext.warning(mf)

        # Make JSON mapping
        atom_to_rosetta_residue_map = stripped_p.get_atom_sequence_to_rosetta_json_map()
        rosetta_to_atom_residue_map = stripped_p.get_rosetta_sequence_to_atom_json_map()

        #pprint.pprint(stripped_p.rosetta_to_atom_sequence_maps)
        #pprint.pprint(stripped_p.get_atom_sequence_to_rosetta_map())
        return


        # get_atom_sequence_to_rosetta_json_map

        print('')
        import sys
        sys.exit(0)
        #stripped_p.map_pdb_residues_to_rosetta_residues(mutations)

        # Assert that there are no empty sequences in the Rosetta-processed PDB file
        return
        assert(sorted(stripped_p.atom_sequences.keys()) == sorted(pdb_chains_to_keep))
        for chain_id, sequence in stripped_p.atom_sequences.iteritems():
            assert(len(sequence) > 0)


        return
        sys.exit(0)

        # to create a stripped pdb, resfile/mutfile, and residue mapping

        # The next main piece of work is rewritting the PDB object to return the residue mapping using the features reporter.
        # This will involve replacing the stripForDDG function and will result in removing the dependency of PDB on BioPython.


        parameters = (experimentID,)
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


        ExtraParameters = {}
        ExtraParameters = pickle.dumps(ExtraParameters)

        prediction_id = None
        try:
            raise Exception('test')
            # Check to see whether the user dataset experiment/protocol combination already exists for this prediction set in which case we early-out
            qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%ss AND ProtocolID=%ss' % self._get_prediction_table()
            existing_record = self.DDG_db.execute_select(qry, parameters=(prediction_set_id, user_dataset_experiment_id, protocol_id))
            if existing_record:
                assert(len(existing_record) == 1)
                return existing_record[0]['ID']

            prediction_record = dict(
                PredictionSet = prediction_set_id,
                PPMutagenesisID = pp_mutagenesis_id,
                UserPPDataSetExperimentID = user_dataset_experiment_id,
                ProtocolID = protocol_id,
                TemporaryProtocolField = None,
                Status = 'queued',
                Cost = read_number_of_residues(),
                KeptHETATMLines = keep_hetatm_lines,
            )
            if not test_only:
                self.DDG_db.insertDictIfNew(self._get_prediction_table(), prediction_record, ['PredictionSet', 'UserPPDataSetExperimentID', 'ProtocolID'])
                existing_record = self.DDG_db.execute_select(qry, parameters=(prediction_set_id, user_dataset_experiment_id, protocol_id))
                assert(len(existing_record) == 1)
                prediction_id = existing_record[0]['ID']

                # add the stripped pdb, resfile/mutfile, and residue mapping associations

                prediction_file_record = dict(
                    PredictionPPIID = prediction_id,
                    FileContentID = None,
                    Filename = None,
                    Filetype = None,
                    FileRole = None,
                    Stage = 'Input'
                )
                #self.DDG_db.insertDictIfNew('PredictionPPIFile', prediction_file_record, ['PredictionPPIID', 'Filename', 'Stage'])
        except:
            if prediction_id:
                pass # delete all associated file records and then the prediction record
            prediction_id = None

        return prediction_id


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



