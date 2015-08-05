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
import json
import zipfile
import traceback

from api_layers import *
from db_api import ddG
from tools import colortext
from tools.bio.pdb import PDB
from tools.bio.basics import ChainMutation
from tools.fs.fsio import read_file
from tools.rosetta.input_files import Mutfile, Resfile

def get_interface_with_config_file(host_config_name = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None):
    # Uses ~/.my.cnf to get authentication information
    ### Example .my.cnf (host_config_name will equal guybrush2):
    ### [clientguybrush2]
    ### user=kyleb
    ### password=notmyrealpass
    ### host=guybrush.ucsf.edu
    my_cnf_path = os.path.expanduser(os.path.join('~', '.my.cnf'))
    if not os.path.isfile( os.path.expanduser(my_cnf_path) ):
        raise Exception("A .my.cnf file must exist at: " + my_cnf_path)

    # These three variables must be set in a section of .my.cnf named host_config_name
    user = None
    password = None
    host = None
    with open(my_cnf_path, 'r') as f:
        parsing_config_section = False
        for line in f:
            if line.strip() == '[client%s]' % host_config_name:
                parsing_config_section = True
            elif line.strip() == '':
                parsing_config_section = False
            elif parsing_config_section:
                if '=' in line:
                    key, val = line.strip().split('=')
                    if key == 'user':
                        user = val
                    elif key == 'password':
                        password = val
                    elif key == 'host':
                        host = val
                else:
                    parsing_config_section = False

    if not user or not password or not host:
        raise Exception("Couldn't find host(%s), username(%s), or password in section %s in %s" % (host, user, host_config_name, my_cnf_path) )

    return get_interface(password, username=user, hostname=host, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)

def get_interface(passwd, username = 'kortemmelab', hostname='kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a BindingAffinityDDGInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(BindingAffinityDDGInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)


class BindingAffinityDDGInterface(ddG):
    '''This is the internal API class that should be NOT used to interface with the database.'''


    def __init__(self, passwd = None, username = 'kortemmelab', hostname = None, rosetta_scripts_path = None, rosetta_database_path = None):
        super(BindingAffinityDDGInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)
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


    @informational_pdb
    def get_chain_sets_for_mutatagenesis(self, mutagenesis_id, complex_id = None):
        '''Gets a list of possibilities for the associated complex and calls get_chains_for_mutatagenesis on each.
             e.g. returns {('1KI1', 0) : {'L' : ['A','B'], 'R' : ['C']}, ('12AB', 2) : {'L' : ['L','H'], 'R' : ['A']}, ...}
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.

           This is a useful method for listing the possible complexes to use in a prediction or to determine whether one
           may be missing. and we need to update the database.'''

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


    @informational_pdb
    def get_pdb_mutations_for_mutagenesis(self, mutagenesis_id, pdb_file_id, set_number, complex_id = None):
        '''Returns the PDB mutations for a mutagenesis experiment as well as the PDB residue information.'''
        pdb_mutations = []
        for pdb_mutation in self.DDG_db.execute_select('''
            SELECT PPMutagenesisPDBMutation.*, PDBResidue.ResidueType, PDBResidue.BFactorMean,
            PDBResidue.BFactorDeviation, PDBResidue.SecondaryStructurePosition, PDBResidue.AccessibleSurfaceArea, ComplexExposure, ComplexDSSP
            FROM
            PPMutagenesisPDBMutation
            INNER JOIN
            PDBResidue ON PPMutagenesisPDBMutation.PDBFileID = PDBResidue.PDBFileID AND PPMutagenesisPDBMutation.Chain = PDBResidue.Chain AND PPMutagenesisPDBMutation.ResidueID = PDBResidue.ResidueID AND PPMutagenesisPDBMutation.WildTypeAA = PDBResidue.ResidueAA
            WHERE PPMutagenesisID=%s AND PDBResidue.PDBFileID=%s AND SetNumber=%s ORDER BY Chain, ResidueID''', parameters=(mutagenesis_id, pdb_file_id, set_number)):
                if complex_id:
                    assert(pdb_mutation['PPComplexID'] == complex_id)
                pdb_mutations.append(pdb_mutation)
        return pdb_mutations


    @informational_job
    def get_complex_details(self, complex_id):
        results = self.DDG_db.execute_select('SELECT * FROM PPComplex WHERE ID=%s', parameters=(complex_id, ))
        if len(results) == 1:
            return results[0]
        return None


    @informational_job
    def get_user_dataset_experiment_details(self, user_dataset_experiment_id, user_dataset_id = None):
        if user_dataset_id:
            ude = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(user_dataset_experiment_id, user_dataset_id))
            if len(ude) != 1:
                raise colortext.Exception('User dataset experiment %d does not exist for/correspond to the user dataset %s.' % (user_dataset_experiment_id, user_dataset_id))
        else:
            ude = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s', parameters=(user_dataset_experiment_id,))
            if len(ude) != 1:
                raise colortext.Exception('User dataset experiment %d does not exist.' % (user_dataset_experiment_id, ))
        ude = ude[0]
        user_dataset_id = ude['UserDataSetID']
        assert(ude['IsComplex'] == 1)

        pdb_mutations = self.get_pdb_mutations_for_mutagenesis(ude['PPMutagenesisID'], ude['PDBFileID'], ude['SetNumber'], complex_id = ude['PPComplexID'])
        return dict(
            Mutagenesis = dict(
                PPMutagenesisID = ude['PPMutagenesisID'],
            ),
            Complex = self.get_complex_details(ude['PPComplexID']),
            Structure = dict(
                PDBFileID = ude['PDBFileID'],
                SetNumber = ude['SetNumber'],
                Partners = self.get_chains_for_mutatagenesis(ude['PPMutagenesisID'], ude['PDBFileID'], ude['SetNumber'], complex_id = ude['PPComplexID']),
            ),
            PDBMutations = pdb_mutations,
        )


    def _export_dataset(self, dataset_id):
        '''Returns a dict containing the dataset information.'''
        dataset_record = self.DDG_db.execute_select('SELECT * FROM DataSet WHERE ID=%s', parameters=(dataset_id,))
        if not dataset_record:
            raise Exception('Dataset %s does not exist in the database.' % dataset_id)
        dataset_record = dataset_record[0]
        if dataset_record['DatasetType'] != 'Binding affinity' and dataset_record['DatasetType'] != 'Protein stability and binding affinity':
            raise Exception('The dataset %s does not contain any binding affinity data..' % dataset_id)

        prediction_record['Files'] = {}
        if include_files:
            prediction_record['Files'] = self.get_job_files(prediction_id, truncate_content = truncate_content)

        # Read the UserPPDataSetExperiment details
        user_dataset_experiment_id = prediction_record['UserPPDataSetExperimentID']
        ude_details = self.get_user_dataset_experiment_details(user_dataset_experiment_id)
        assert(ude_details['Mutagenesis']['PPMutagenesisID'] == prediction_record['PPMutagenesisID'])
        for k, v in ude_details.iteritems():
            assert(k not in prediction_record)
            prediction_record[k] = v
        return prediction_record


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
    def add_development_protocol_command_lines(self, prediction_set_id, protocol_name, application, template_command_line, rosetta_script_file = None):
        dev_protocol_id = self._get_dev_protocol_id(protocol_name)
        if not dev_protocol_id:
            dev_protocol_id = self._create_dev_protocol(protocol_name, application, template_command_line)

        rosetta_script = None
        if rosetta_script_file:
            with open(rosetta_script_file, 'r') as f:
                rosetta_script = f.read()

        prediction_ids = self.get_prediction_ids(prediction_set_id)

        # All functions within the next with block should use the same database cursor.
        # The commands then function as parts of a transaction which is rolled back if errors occur within the block
        # or else is committed.
        self.DDG_db._get_connection()
        con = self.DDG_db.connection
        with con:
            cur = con.cursor()
            for prediction_id in prediction_ids:
                qry = 'UPDATE PredictionPPI SET DevelopmentProtocolID=%s WHERE ID=%s'
                cur.execute(qry, (dev_protocol_id, prediction_id))

                if rosetta_script:
                # Add the Rosetta script to the database, not using cursor
                    self._add_prediction_file(prediction_id, rosetta_script, os.path.basename(rosetta_script_file), 'RosettaScript', 'RosettaScript', 'Input', rm_trailing_line_whitespace = True)

                
    @job_creator
    def add_job(self, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = None):
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
        return self._add_job(prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)


    @job_creator
    def add_job_by_user_dataset_record(self, prediction_set_id, user_dataset_name, user_dataset_experiment_id, protocol_id, extra_rosetta_command_flags = None, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = None):
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
        return self._add_job(prediction_set_id, protocol_id, ude['PPMutagenesisID'], ude['PPComplexID'], ude['PDBFileID'], ude['SetNumber'], extra_rosetta_command_flags = extra_rosetta_command_flags, user_dataset_experiment_id = user_dataset_experiment_id, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = test_only, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)


    @job_creator
    def add_prediction_run(self, prediction_set_id, user_dataset_name, extra_rosetta_command_flags = None, protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, input_files = {}, quiet = False, test_only = False, only_single_mutations = False, short_run = False, test_run_first = True, show_full_errors = False):
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

        # Create a cache to speed up job insertion
        pdb_residues_to_rosetta_cache = {}

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
                if protocol_id:
                    existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
                else:
                    existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID IS NULL", parameters=(prediction_set_id, ude['ID']))
                if len(existing_results) > 0: continue

                # Test the prediction setup
                prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags, keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = True, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)
                # Progress counter
                count += 1
                if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
                if short_run and count > 4: break
            if not quiet: print('')

        if test_only:
            return

        # Progress counter setup
        failed_jobs = {}
        if not quiet:
            colortext.message('Adding %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
        count, records_per_dot = 0, 50
        showprogress = not(quiet) and len(user_dataset_experiments) > 300
        if showprogress: print("|" + ("*" * (int(len(user_dataset_experiments)/records_per_dot)-2)) + "|")

        # Add the individual predictions
        for ude in user_dataset_experiments:

            # If the mutagenesis already exists in the prediction set, do not add it again
            if protocol_id:
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
            else:
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID IS NULL", parameters=(prediction_set_id, ude['ID']))
            if len(existing_results) == 0:
                # Add the prediction
                try:
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
                    prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = False, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)
                except Exception, e:
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
                    ude_record = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(ude['ID'], user_dataset_id))
                    ude_record = ude_record[0]
                    assert(ude_record['ID'] == ude['ID'])
                    colortext.error('Adding the prediction for UserPPDataSetExperimentID %(ID)d failed (%(PDBFileID)s).' % ude_record)
                    failed_jobs[ude_record['PDBFileID']] = failed_jobs.get(ude_record['PDBFileID'], 0)
                    failed_jobs[ude_record['PDBFileID']] += 1
                    if show_full_errors:
                        print(e)
                        print(traceback.format_exc())

            # Progress counter
            count += 1
            if showprogress and count % records_per_dot == 0: colortext.write(".", "green", flush = True)
            if short_run and count > 4: break

        if failed_jobs:
            colortext.error('Some jobs failed to run:\n%s' % pprint.pformat(failed_jobs))
        if not quiet: print('')
        return True



    def _create_pdb_residues_to_rosetta_cache_mp(self, pdb_residues_to_rosetta_cache, pdb_file_id, pdb_chains_to_keep, extra_rosetta_command_flags, keep_hetatm_lines):
        # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
        raise Exception('Shane should finish this')
        assert(type(pdb_residues_to_rosetta_cache) == None)# use the manager dictproxy)
        pdb_file = self.DDG_db.execute_select("SELECT * FROM PDBFile WHERE ID=%s", parameters = (pdb_file_id,))
        p = PDB(pdb_file[0]['Content'])
        p.strip_to_chains(list(pdb_chains_to_keep))
        if not keep_hetatm_lines:
            p.strip_HETATMs()
        stripped_p = PDB('\n'.join(p.lines))

        stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, self.rosetta_database_path, extra_command_flags = extra_rosetta_command_flags)
        atom_to_rosetta_residue_map = stripped_p.get_atom_sequence_to_rosetta_json_map()
        rosetta_to_atom_residue_map = stripped_p.get_rosetta_sequence_to_atom_json_map()
        cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, self.rosetta_database_path, extra_rosetta_command_flags)
        pdb_residues_to_rosetta_cache[cache_key] = dict(
            stripped_p = stripped_p,
            atom_to_rosetta_residue_map = atom_to_rosetta_residue_map,
            rosetta_to_atom_residue_map = rosetta_to_atom_residue_map)


    @job_creator
    def add_prediction_run_mp(self, prediction_set_id, user_dataset_name, extra_rosetta_command_flags = None, protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, input_files = {}, quiet = False, only_single_mutations = False, short_run = False, show_full_errors = False):
        '''This is a multiprocessing version of add_prediction_run and should be used in favor of that function as it runs faster.

           It takes advantage of parallelism at two points - creating the stripped PDB files and mutfiles for input and
           inserting the jobs (MD5 is run multiple times for each job).
           It was simple/quicker to write this as a 2-step method with a bottleneck in the middle i.e. it waits until all
           stripped PDB files are generated before adding the jobs.
           This could be made even more parallel by removing the bottleneck i.e. the process which strips the PDBs could
           then call _add_job immediately rather than waiting for the other calls to _create_pdb_residues_to_rosetta_cache_mp
           to complete.
           '''

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

        # Create a cache to speed up job insertion
        pdb_residues_to_rosetta_cache = manager dictproxy

        # Create the stripped PDBs and residue maps in parallel using the multiprocessing module
    #todo: write this function on Monday - get_user_dataset_pdb_partner_chains should return a set (<list of {'id' : pdb_file_id, 'L' : <list of chain ids>, , 'R' : <list of chain ids>} dicts>)
        pdb_partner_chains = self.get_user_dataset_pdb_partner_chains(user_dataset_name, tagged_subset = tagged_subset)
        for ppc in pdb_partner_chains:

            apply_async self._create_pdb_residues_to_rosetta_cache_mp(pdb_residues_to_rosetta_cache, ppc['id'], set(ppc['L'] + ppc['R']), extra_rosetta_command_flags, keep_hetatm_lines)
        .join()

        # Progress counter setup
        failed_jobs = {}
        if not quiet:
            colortext.message('Adding %d predictions spanning %d PDB files for %suser dataset "%s" using protocol %s.' % (len(user_dataset_experiments), len(pdb_file_ids), tagged_subset_str, user_dataset_name, str(protocol_id or 'N/A')))
        count, records_per_dot = 0, 50
        showprogress = not(quiet) and len(user_dataset_experiments) > 300
        if showprogress: print("|" + ("*" * (int(len(user_dataset_experiments)/records_per_dot)-2)) + "|")

        # Add the individual predictions
        for ude in user_dataset_experiments:

            # If the mutagenesis already exists in the prediction set, do not add it again
            if protocol_id:
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID=%s", parameters=(prediction_set_id, ude['ID'], protocol_id))
            else:
                existing_results = self.DDG_db.execute_select("SELECT * FROM PredictionPPI WHERE PredictionSet=%s AND UserPPDataSetExperimentID=%s AND ProtocolID IS NULL", parameters=(prediction_set_id, ude['ID']))
            if len(existing_results) == 0:
                # Add the prediction
                try:
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
                    prediction_id = self.add_job_by_user_dataset_record(prediction_set_id, user_dataset_name, ude['ID'], protocol_id, extra_rosetta_command_flags = extra_rosetta_command_flags,  keep_hetatm_lines = keep_hetatm_lines, input_files = input_files, test_only = False, pdb_residues_to_rosetta_cache = pdb_residues_to_rosetta_cache)
                except Exception, e:
                    user_dataset_id = self.get_defined_user_datasets()[user_dataset_name]['ID']
                    ude_record = self.DDG_db.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE ID=%s AND UserDataSetID=%s', parameters=(ude['ID'], user_dataset_id))
                    ude_record = ude_record[0]
                    assert(ude_record['ID'] == ude['ID'])
                    colortext.error('Adding the prediction for UserPPDataSetExperimentID %(ID)d failed (%(PDBFileID)s).' % ude_record)
                    failed_jobs[ude_record['PDBFileID']] = failed_jobs.get(ude_record['PDBFileID'], 0)
                    failed_jobs[ude_record['PDBFileID']] += 1
                    if show_full_errors:
                        print(e)
                        print(traceback.format_exc())

            # Progress counter
            count += 1
            if showprogress and count % records_per_dot == 0: colortext.write(".", "cyan", flush = True)
            if short_run and count > 4: break

        if failed_jobs:
            colortext.error('Some jobs failed to run:\n%s' % pprint.pformat(failed_jobs))
        if not quiet: print('')
        return True


    @job_creator
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set):
        raise Exception('not implemented yet')
        #assert(existing_prediction_set exists and has records)
        #assert(new_prediction_set is empty)
        #for each prediction record, add the record and all associated predictionfile records,


    def _add_job(self, prediction_set_id, protocol_id, pp_mutagenesis_id, pp_complex_id, pdb_file_id, pp_complex_pdb_set_number, extra_rosetta_command_flags = None, user_dataset_experiment_id = None, keep_hetatm_lines = False, input_files = {}, test_only = False, pdb_residues_to_rosetta_cache = {}):
        '''This is the internal function which adds a prediction job to the database. We distinguish it from add_job as
           prediction jobs added using that function should have no associated user dataset experiment ID.

           pdb_residues_to_rosetta_cache can be used to speed up job insertion. When the same PDB/chains combination is used again, this cache uses the old mapping rather than run RosettaScripts again.

           The extra_rosetta_command_flags variable is used to add additional flags e.g. "-ignore_zero_occupancy false".
           These are used to generate a mapping from PDB to Rosetta numbering so they should be set according to how they
           are set in the protocol. In particular, include any flags which have an effect on what residues are present.
           '-ignore_zero_occupancy false' and '-ignore_unrecognized_res' are typically used.
           '''

        # todo: do something with input_files when we use that here - see add_prediction_run
        assert(not(input_files))

        # Preliminaries
        if not self.rosetta_scripts_path or not os.path.exists(self.rosetta_scripts_path):
            raise Exception('The path "%s" to the RosettaScripts executable does not exist.' % self.rosetta_scripts_path)

        cache_maps = False
        if isinstance(pdb_residues_to_rosetta_cache, dict):
            cache_maps = True

        # Information for debugging
        pp_complex = self.DDG_db.execute_select("SELECT * FROM PPComplex WHERE ID=%s", parameters = (pp_complex_id,))

        # Determine the list of PDB chains that will be kept
        pdb_chains = self.get_chains_for_mutatagenesis(pp_mutagenesis_id, pdb_file_id, pp_complex_pdb_set_number, complex_id = pp_complex_id)
        pdb_chains_to_keep = set(pdb_chains['L'] + pdb_chains['R'])
        if self.rosetta_database_path:
            cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, self.rosetta_database_path, extra_rosetta_command_flags)
        else:
            cache_key = (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)), self.rosetta_scripts_path, extra_rosetta_command_flags)

        if cache_maps and pdb_residues_to_rosetta_cache.get(cache_key):
            stripped_p = pdb_residues_to_rosetta_cache[cache_key]['stripped_p']
        else:
            # Retrieve the PDB file content, strip out the unused chains, and create a PDB object
            pdb_file = self.DDG_db.execute_select("SELECT * FROM PDBFile WHERE ID=%s", parameters = (pdb_file_id,))
            p = PDB(pdb_file[0]['Content'])
            p.strip_to_chains(list(pdb_chains_to_keep))
            if not keep_hetatm_lines:
                p.strip_HETATMs()
            stripped_p = PDB('\n'.join(p.lines))

        # Determine PDB chains to move
        pdb_chains_to_move_str = ','.join(sorted(set(pdb_chains['R'])))

        # Check for CSE and MSE
        try:
            if 'CSE' in stripped_p.residue_types:
                raise Exception('This case contains a CSE residue which may (or may not) cause an issue.')
            elif 'MSE' in stripped_p.residue_types:
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


        # Make JSON mappings
        if cache_maps and pdb_residues_to_rosetta_cache.get(cache_key):
            atom_to_rosetta_residue_map = pdb_residues_to_rosetta_cache[cache_key]['atom_to_rosetta_residue_map']
            rosetta_to_atom_residue_map = pdb_residues_to_rosetta_cache[cache_key]['rosetta_to_atom_residue_map']
        else:
            if self.rosetta_database_path:
                stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, self.rosetta_database_path, extra_command_flags = extra_rosetta_command_flags)
            else:
                stripped_p.construct_pdb_to_rosetta_residue_map(self.rosetta_scripts_path, extra_command_flags = extra_rosetta_command_flags)
            atom_to_rosetta_residue_map = stripped_p.get_atom_sequence_to_rosetta_json_map()
            rosetta_to_atom_residue_map = stripped_p.get_rosetta_sequence_to_atom_json_map()

        # Make mutfile and resfile
        rosetta_mutations = stripped_p.map_pdb_residues_to_rosetta_residues(mutations)
        mf = Mutfile.from_mutagenesis(rosetta_mutations)
        mutfile_name = 'mutations.mutfile'
        rf = Resfile.from_mutageneses(mutations)
        resfile_name = 'mutations.resfile'

        if cache_maps and (not pdb_residues_to_rosetta_cache.get(cache_key)):
            pdb_residues_to_rosetta_cache[cache_key] = dict(
                stripped_p = stripped_p,
                atom_to_rosetta_residue_map = atom_to_rosetta_residue_map,
                rosetta_to_atom_residue_map = rosetta_to_atom_residue_map)

        # Assert that there are no empty sequences in the Rosetta-processed PDB file
        total_num_residues = 0
        d = json.loads(rosetta_to_atom_residue_map)
        stripped_p_chains = stripped_p.atom_sequences.keys()
        for chain_id in stripped_p_chains:
            num_chain_residues = len([z for z in d.values() if z[0] == chain_id])
            total_num_residues += num_chain_residues
            assert(num_chain_residues > 0)

        pdb_filename = '%s_%s.pdb' % (pdb_file_id, ''.join(sorted(pdb_chains_to_keep)))
        
        # Create parameter substitution dictionary
        parameter_sub_dict = {
            '%%input_pdb%%' : pdb_filename,
            '%%chainstomove%%' : pdb_chains_to_move_str,
            '%%pathtoresfile%%' : resfile_name,
            '%%pathtomutfile%%' : mutfile_name,
        }
            
        if test_only:
            return

        # All functions within the next with block should use the same database cursor.
        # The commands then function as parts of a transaction which is rolled back if errors occur within the block
        # or else is committed.
        self.DDG_db._get_connection()
        con = self.DDG_db.connection
        with con:
            cur = con.cursor()
            
            if protocol_id:
                qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID=%%s' % self._get_prediction_table()
                cur.execute(qry, (prediction_set_id, user_dataset_experiment_id, protocol_id))
            else:
                qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID IS NULL' % self._get_prediction_table()
                cur.execute(qry, (prediction_set_id, user_dataset_experiment_id))
            existing_record = cur.fetchall()
            assert(len(existing_record) == 0)

            prediction_record = dict(
                PredictionSet = prediction_set_id,
                PPMutagenesisID = pp_mutagenesis_id,
                UserPPDataSetExperimentID = user_dataset_experiment_id,
                ProtocolID = protocol_id,
                JSONParameters = json.dumps(parameter_sub_dict),
                DevelopmentProtocolID = None,
                ExtraParameters = extra_rosetta_command_flags,
                Status = 'queued',
                Cost = total_num_residues,
                KeptHETATMLines = keep_hetatm_lines,
            )
            sql, params = self.DDG_db.create_insert_dict_string(self._get_prediction_table(), prediction_record, ['PredictionSet', 'UserPPDataSetExperimentID', 'ProtocolID'])
            cur.execute(sql, params)
            prediction_id = cur.lastrowid

            # Add the stripped PDB file
            self._add_prediction_file(prediction_id, '\n'.join(stripped_p.lines), pdb_filename, 'PDB', 'StrippedPDB', 'Input', db_cursor = cur, rm_trailing_line_whitespace = True, forced_mime_type = 'chemical/x-pdb')

            # Add the mutfile
            self._add_prediction_file(prediction_id, str(mf), mutfile_name, 'Mutfile', 'Mutfile', 'Input', db_cursor = cur, rm_trailing_line_whitespace = True)
            # Add the resfile
            self._add_prediction_file(prediction_id, str(rf), resfile_name, 'Resfile', 'Resfile', 'Input', db_cursor = cur, rm_trailing_line_whitespace = True)

            # Add the residue mappings
            self._add_prediction_file(prediction_id, rosetta_to_atom_residue_map, 'rosetta2pdb.resmap.json', 'RosettaPDBMapping', 'Rosetta residue->PDB residue map', 'Input', db_cursor = cur, forced_mime_type = "application/json")
            self._add_prediction_file(prediction_id, atom_to_rosetta_residue_map, 'pdb2rosetta.resmap.json', 'RosettaPDBMapping', 'PDB residue->Rosetta residue map', 'Input', db_cursor = cur, forced_mime_type = "application/json")

        if protocol_id:
            qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID=%%s' % self._get_prediction_table()
            existing_record = self.DDG_db.execute_select(qry, parameters=(prediction_set_id, user_dataset_experiment_id, protocol_id))
        else:
            qry = 'SELECT * FROM %s WHERE PredictionSet=%%s AND UserPPDataSetExperimentID=%%s AND ProtocolID IS NULL' % self._get_prediction_table()
            existing_record = self.DDG_db.execute_select(qry, parameters=(prediction_set_id, user_dataset_experiment_id))
        assert(len(existing_record) == 1)
        prediction_id = existing_record[0]['ID']
        return prediction_id


    #== Job execution/completion API ===========================================================
    #
    # This part of the API is responsible for starting jobs and setting them as failed or
    # completed


    @informational_job
    def get_job_details(self, prediction_id, include_files = True, truncate_content = None):
        prediction_record = self.DDG_db.execute_select('SELECT * FROM PredictionPPI WHERE ID=%s', parameters=(prediction_id,))
        if not prediction_record:
            raise Exception('No details could be found for prediction #%d in the database.' % prediction_id)
        prediction_record = prediction_record[0]
        prediction_record['Files'] = {}
        if include_files:
            prediction_record['Files'] = self.get_job_files(prediction_id, truncate_content = truncate_content)

        mutfile_content = self.create_mutfile(prediction_id)


        # Read the UserPPDataSetExperiment details
        user_dataset_experiment_id = prediction_record['UserPPDataSetExperimentID']
        ude_details = self.get_user_dataset_experiment_details(user_dataset_experiment_id)
        assert(ude_details['Mutagenesis']['PPMutagenesisID'] == prediction_record['PPMutagenesisID'])
        for k, v in ude_details.iteritems():
            assert(k not in prediction_record)
            prediction_record[k] = v
        return prediction_record


    @job_execution
    def set_job_temporary_protocol_field(self, prediction_id, prediction_set_id, temporary_protocol_field):
        raise Exception('not implemented yet')


    @job_execution
    def start_job(self, prediction_id, prediction_set_id):
        '''Sets the job status to "active". prediction_set must be passed and is used as a sanity check.'''
        prediction_record = self.DDG_db.execute_select('SELECT * FROM PredictionPPI WHERE ID=%s AND PredictionSet=%s', parameters=(prediction_id, prediction_set_id))
        if prediction_record['Protocol'] == None:
            print('empty Protocol')
            if prediction_record['DevelopmentProtocolID'] == None:
                raise Exception('Neither the Protocol nor the DevelopmentProtocolID is set for this job - it cannot be started without this information.')

        raise Exception('This function needs to be implemented by subclasses of the API.')


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


    @job_completion
    def complete_job(self, prediction_id, prediction_set, scores, maxvmem, ddgtime, files = []):
        '''Sets a job to 'completed' and stores scores. prediction_set must be passed and is used as a sanity check.'''

        raise Exception('This function needs to be implemented by subclasses of the API.')


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
    def _get_prediction_structure_scores_table(self): return 'PredictionPPIStructureScore'
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


    @informational_job
    def get_development_protocol(self, development_protocol_id):
        results = self.DDG_db.execute_select('SELECT * FROM DevelopmentProtocol WHERE ID = %s', parameters=(development_protocol_id) )
        assert( len(results) == 1 )
        return results[0]


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


    def _get_dev_protocol_id(self, name):
        dev_protocol_ids = self.DDG_db.execute_select("SELECT ID FROM DevelopmentProtocol WHERE Name=%s", parameters = (name))
        if len(dev_protocol_ids) == 0:
            return None
        elif len(dev_protocol_ids) == 1:
            return int(dev_protocol_ids[0]['ID'])
        else:
            raise Exception("DevelopmentProtocol table was originally set up so that names are unique; this has obviously changed")

        
    def _create_dev_protocol(self, name, application, template_command_line):
        dev_prot_record = {
            'Name' : name,
            'Application' : application,
            'TemplateCommandLine' : template_command_line,
        }
        sql, params = self.DDG_db.create_insert_dict_string('DevelopmentProtocol', dev_prot_record)
        self.DDG_db.execute(sql, params)

