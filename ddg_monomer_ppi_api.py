import os
import json
import sqlite3
import interface_calc
import sys
import shutil
import multiprocessing
import random
import time
import getpass
import tempfile
import cPickle as pickle
import copy
import zipfile
import datetime

from sqlalchemy import and_, update
from sqlalchemy.orm import load_only, Load

import klab.cluster_template.parse_settings as parse_settings
from klab.Reporter import Reporter
from klab.MultiWorker import MultiWorker
from klab.fs.zip_util import zip_file_with_gzip, unzip_file
from klab.fs.io import sanitize_filename
from klab.cluster_template.write_run_file import process as write_run_file

from ppi_api import BindingAffinityDDGInterface, get_interface
from api_layers import *
from db_api import ddG
import db_schema as dbmodel
from util import fill_empty_score_dict

DeclarativeBase = dbmodel.DeclarativeBase

# Constants for cluster runs
rosetta_scripts_xml_file = os.path.join('ddglib', 'score_partners.xml')
output_db3 = 'output.db3'
setup_run_with_multiprocessing = True
tmpdir_location = '/dbscratch/%s/tmp' % getpass.getuser()

def total_seconds(td):
    '''
    Included in python 2.7 but here for backwards-compatibility for old Python versions
    '''
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

def get_interface(passwd, username = 'kortemmelab', hostname = 'guybrush.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, port = 3306, file_content_buffer_size = None):
    '''This is the function that should be used to get a DDGMonomerInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(
        DDGMonomerInterface,
        passwd = passwd, username = username, hostname = hostname,
        rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path,
        port = port,
        file_content_buffer_size = file_content_buffer_size,
    )

class DDGMonomerInterface(BindingAffinityDDGInterface):
    def __init__(
            self,
            passwd = None, username = None, hostname = None,
            rosetta_scripts_path = None, rosetta_database_path = None,
            port = 3306,
            file_content_buffer_size = None
    ):
        super(DDGMonomerInterface, self).__init__(
            passwd = passwd, username = username, hostname = hostname,
            rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path,
            port = port,
            file_content_buffer_size = file_content_buffer_size
        )
        self.rescore_args = {} # Stores arguments for rescoring step, to be dumped later
        self.all_score_types = ['WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex']
        self.all_score_types_index = {}
        for i, score_type in enumerate(self.all_score_types):
            self.all_score_types_index[score_type] = i
        self.master_scores_list = []
        self.ddg_output_path_cache = {} # Stores paths to ddG job output directories, or unzipped job output directories
        self.unzipped_ddg_output_paths = [] # Stores paths to unzipped ddG job output directories (that need to be cleared at the end of this object's life, or before)




    def create_cluster_run_rescore_dir(self, output_dir, passed_job_name = None):
        first_task_count = 0
        dir_count = 1
        rescore_args_keys = sorted( self.rescore_args.keys() )

        while first_task_count + 1 < len(rescore_args_keys):
            settings = parse_settings.get_dict()
            job_name = '%s-%d' % (passed_job_name, dir_count) or '%s-%s_rescore_ddg_monomer-%d' % (time.strftime("%y%m%d"), getpass.getuser(), dir_count)
            dir_count += 1
            job_output_dir = os.path.join(output_dir, job_name)
            if not os.path.isdir(job_output_dir):
                try:
                    os.makedirs(job_output_dir)
                except OSError:
                    pass

            num_jobs = 0
            max_prediction_id = 0
            break_loop = False
            last_task_count = first_task_count
            for ddg_output_path in rescore_args_keys[first_task_count:]:
                last_task_count += 1
                for prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, extra_flags in self.rescore_args[ddg_output_path]:
                    num_jobs += 1
                    if prediction_id > max_prediction_id:
                        max_prediction_id = prediction_id
                    if num_jobs >= 50000:
                        break_loop = True
                if break_loop:
                    break

            settings['scriptname'] = 'rescore_ddg_monomer'
            settings['tasks_per_process'] = 50
            settings['numjobs'] = num_jobs
            settings['mem_free'] = '1.4G'
            settings['appname'] = 'rosetta_scripts'
            settings['output_dir'] = job_output_dir
            job_dict = {}
            job_data_dir = os.path.join(job_output_dir, 'data')
            if not os.path.isdir(job_data_dir):
                try:
                    os.makedirs(job_data_dir)
                except OSError:
                    pass

            data_protocol_path = os.path.join(job_data_dir, os.path.basename(rosetta_scripts_xml_file))
            if not os.path.isfile(data_protocol_path):
                shutil.copy(rosetta_scripts_xml_file, data_protocol_path)
            rel_protocol_path = os.path.relpath(data_protocol_path, job_output_dir)
            settings['rosetta_args_list'] = ['-inout:dbms:database_name', output_db3]

            r = Reporter('copying/zipping ddG output PDBs', entries='PDBs')
            print 'Total number of tasks: %d' % num_jobs
            r.set_total_count(num_jobs)


            if setup_run_with_multiprocessing:
                worker = MultiWorker(process_cluster_rescore_helper, n_cpu=min(multiprocessing.cpu_count(), 16), reporter=r)

            for ddg_output_path in rescore_args_keys[first_task_count:last_task_count]:
                for prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, extra_flags in self.rescore_args[ddg_output_path]:
                    if setup_run_with_multiprocessing:
                        worker.addJob( (ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path, max_prediction_id, extra_flags ) )
                    else:
                        task_name, arg_dict = process_cluster_rescore_helper( ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path, max_prediction_id, extra_flags )
                        job_dict[task_name] = arg_dict
                        r.increment_report()

            if setup_run_with_multiprocessing:
                worker.finishJobs()
                for task_name, arg_dict in worker.data:
                    job_dict[task_name] = arg_dict
            else:
                r.done()
            write_run_file(settings, job_dict = job_dict)
            first_task_count = last_task_count


    def add_rescore_cluster_run(self, ddg_output_path, chains_to_move, score_method_id, prediction_id):
        structs_with_both_rounds = self.find_structs_with_both_rounds(ddg_output_path)

        score_method_details = self.get_score_method_details()[score_method_id]
        method_name = score_method_details['MethodName']
        author = score_method_details['Authors']

        if method_name.startswith('Rescore-') and author == 'Kyle Barlow':
            score_fxn = method_name[8:].lower()
            if score_fxn.startswith('beta'):
                extra_flags = ['-' + score_fxn]
            else:
                extra_flags = []
        else:
            score_fxn = 'interface'
            extra_flags = []

        if ddg_output_path not in self.rescore_args:
            self.rescore_args[ddg_output_path] = []

        for round_num in structs_with_both_rounds:
            wt_pdb, mutant_pdb = structs_with_both_rounds[round_num]
            self.rescore_args[ddg_output_path].append([
                prediction_id,
                os.path.abspath(wt_pdb),
                chains_to_move,
                score_fxn,
                round_num,
                'wt',
                extra_flags,
            ])
            self.rescore_args[ddg_output_path].append([
                prediction_id,
                os.path.abspath(mutant_pdb),
                chains_to_move,
                score_fxn,
                round_num,
                'mut',
                extra_flags,
            ])

    def update_prediction_id_status(self, prediction_set_id, root_directory, verbose = True):
        # Looks for output file in root directory and reads for job status
        unfinished_prediction_ids = self.get_unfinished_prediction_ids(prediction_set_id)
        if len(unfinished_prediction_ids) == 0:
            return
        if verbose:
            r = Reporter('parsing output files (if they exist) for queued predictions')
            r.set_total_count( len(unfinished_prediction_ids) )

        # Constants
        output_file_substring = 'run-2.py.o'
        date_format_string = '%Y-%m-%d %H:%M:%S'
        starting_time_string = 'Starting time:'
        ending_time_string = 'Ending time:'
        tsession = self.importer.session
        for prediction_id in unfinished_prediction_ids:
            ddg_output_dir = self.find_ddg_output_directory(prediction_id, root_directory)
            if ddg_output_dir:
                output_files = [os.path.join(ddg_output_dir, f) for f in os.listdir(ddg_output_dir) if output_file_substring in f]
                if len(output_files) == 0:
                    if verbose:
                        print '\nDirectory %s missing output file\n' % ddg_output_dir
                        r.decrement_total_count()
                elif len(output_files) > 1:
                    raise Exception('ERROR: found more than one output file (defined as having string %s in name) in directory %s' % (output_file_substring, ddg_output_dir))
                else:
                    # Parse output file
                    starting_time = None
                    ending_time = None
                    return_code = None
                    virtual_memory_usage = None
                    elapsed_time = None
                    nstruct = None
                    status = 'active'
                    ddg_iterations_str = '-ddg::iterations '
                    with open(output_files[0], 'r') as f:
                        for line in f:
                            if line.startswith(starting_time_string):
                                starting_time = line[len(starting_time_string):].strip()
                            elif line.startswith(ending_time_string):
                                ending_time = line[len(ending_time_string):].strip()
                            elif 'return code' in line:
                                return_code = int(line.strip().split()[-1].strip())
                            elif 'virtual memory usage' in line:
                                line = line.strip()
                                assert( line.endswith('G') )
                                line = line[:-1]
                                virtual_memory_usage = float(line.strip().split()[-1].strip())
                            elif ddg_iterations_str in line:
                                nstruct = int( line[line.find(ddg_iterations_str)+len(ddg_iterations_str):line.find("'", line.find(ddg_iterations_str)+len(ddg_iterations_str))] )
                    if starting_time and ending_time:
                        starting_time_dt = datetime.datetime.strptime(starting_time, date_format_string)
                        ending_time_dt = datetime.datetime.strptime(ending_time, date_format_string)
                        if nstruct != None:
                            elapsed_time = total_seconds(ending_time_dt - starting_time_dt) / (60.0 * float(nstruct)) # 60 seconds in a minute, and average time per structure


                    if return_code != None:
                        if return_code == 0:
                            status = 'done'
                        else:
                            status = 'failed'

                    prediction_record = tsession.query(self.PredictionTable).filter(self.PredictionTable.ID == prediction_id).one()
                    prediction_record.StartDate = starting_time
                    prediction_record.EndDate = ending_time
                    prediction_record.Status = status
                    prediction_record.maxvmem = virtual_memory_usage
                    prediction_record.DDGTime = elapsed_time
                    prediction_record.ERRORS = str(return_code)
                    tsession.flush()
                    tsession.commit()


                    if verbose:
                        r.increment_report()
            elif verbose:
                r.decrement_total_count()

        if verbose:
            r.done()

    def find_ddg_output_directory(self, prediction_id, root_directory):
        # Returns the ddG output directory path in root directory if it exists
        # If not, checks to see if zipped output directory exists in root directory, and if so,
        #    returns path to unzipped temporary directory
        if prediction_id in self.ddg_output_path_cache:
            return self.ddg_output_path_cache[prediction_id]

        ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
        if not os.path.isdir( ddg_output_path ):
            ddg_output_zip = os.path.join(root_directory, '%d.zip' % prediction_id)
            if os.path.isfile(ddg_output_zip):
                tmp_dir = unzip_to_tmp_dir(prediction_id, ddg_output_zip)
                self.unzipped_ddg_output_paths.append( tmp_dir )
                ddg_output_path = os.path.join(tmp_dir, '%d-ddg' % prediction_id)

        if not os.path.isdir( ddg_output_path):
            return None
        else:
            self.ddg_output_path_cache[prediction_id] = ddg_output_path
            return ddg_output_path

    def cleanup_tmp_ddg_output_directories(self):
        # Removes any temporary unzipped ddG monomer output directories
        if len(self.unzipped_ddg_output_paths) == 0:
            return

        r = Reporter('Removing temporary directories')
        r.set_total_count( len(self.unzipped_ddg_output_paths) )
        for tmp_dir in self.unzipped_ddg_output_paths:
            shutil.rmtree(tmp_dir)
            r.increment_report()
        self.unzipped_ddg_output_paths = []
        r.done()

    @job_completion
    def extract_data(self, prediction_set_id, root_directory = None, force = False, score_method_id = None, max_prediction_ids_to_process=None, setup_cluster_run = False):
        '''Extracts the data for the prediction set run and stores it into the database.

           For all PredictionIDs associated with the PredictionSet:
             - looks for a subdirectory of root_directory with the same name as the ID e.g. /some/folder/21412
             - call extract_data_for_case

           Note: we do not use a transaction at this level. We could but it may end up being a very large transaction
           depending on the dataset size. It seems to make more sense to me to use transactions at the single prediction
           level i.e. in extract_data_for_case

           root_directory defaults to /kortemmelab/shared/DDG/ppijobs.
           If force is True then existing records should be overridden.
        '''
        root_directory = root_directory or self.prediction_data_path

        ### Find all prediction_ids that need to have updated states
        if root_directory != self.prediction_data_path:
            # Only looks for output files if run hasn't been zipped yet
            self.update_prediction_id_status(prediction_set_id, root_directory)

        ### Find all prediction_ids with partial score data and remove this data

        ### Find all prediction_ids with missing scores and setup rescoring or rescore on the fly
        prediction_ids = self.get_prediction_ids_without_scores(prediction_set_id, score_method_id = score_method_id)

        random.shuffle( prediction_ids )
        print '%d prediction_ids are yet to be scored' % len(prediction_ids)
        if setup_cluster_run:
            job_name = '%s-%s' % (time.strftime("%y%m%d"), sanitize_filename(prediction_set_id))
            r = Reporter('fetching job details for prediction_ids', entries='job details')
            r.set_total_count(len(prediction_ids))
            for prediction_id in prediction_ids:
                ddg_output_path = self.find_ddg_output_directory(prediction_id, root_directory)
                job_details = self.get_job_details(prediction_id)
                substitution_parameters = json.loads(job_details['JSONParameters'])
                chains_to_move = substitution_parameters['%%chainstomove%%']
                if ddg_output_path:
                    self.add_rescore_cluster_run(ddg_output_path, chains_to_move, score_method_id, prediction_id)
                r.increment_report()
            r.done()
            self.create_cluster_run_rescore_dir( os.path.join(tmpdir_location, 'cluster_run'), passed_job_name = job_name )
        else:
            r = Reporter('rescoring prediction directories with multiprocessing', entries='prediction ids')
            if max_prediction_ids_to_process:
                r.set_total_count( max_prediction_ids_to_process )
            else:
                r.set_total_count( len(prediction_ids) )

            missing_data_count = 0
            for prediction_id in prediction_ids:
                ddg_output_path = self.find_ddg_output_directory(prediction_id, root_directory)
                if ddg_output_path:
                    self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)
                    r.increment_report()
                else:
                    r.decrement_total_count()
                    missing_data_count += 1
                if max_prediction_ids_to_process and r.n >= max_prediction_ids_to_process:
                    print 'Breaking early; processed %d prediction ids' % r.n
                    break
            r.done()
            if missing_data_count > 0:
                print 'Missing data folders/zips for %d prediction_ids' % missing_data_count

        self.cleanup_tmp_ddg_output_directories() # Remove temporary directories (if created)

    def find_structs_with_both_rounds(self, ddg_output_path):
        '''Searchs directory ddg_output_path to find ddg_monomer output structures for all rounds with both wt and mut structures'''
        output_pdbs = [x for x in os.listdir(ddg_output_path) if '.pdb' in x]
        mut_output_pdbs = { int(x.split('.')[0].split('_')[3]) : x for x in output_pdbs if x.startswith('mut')}
        wt_output_pdbs = { int(x.split('.')[0].split('_')[3]) : x for x in output_pdbs if '_wt_' in x}
        all_round_nums = set()
        structs_with_both_rounds = {}
        for round_num in mut_output_pdbs.keys():
            all_round_nums.add( round_num )
        for round_num in wt_output_pdbs.keys():
            all_round_nums.add( round_num )
        for round_num in all_round_nums:
            if round_num in wt_output_pdbs and round_num in mut_output_pdbs:
                structs_with_both_rounds[round_num] = (
                    os.path.join(ddg_output_path, wt_output_pdbs[round_num]),
                    os.path.join(ddg_output_path, mut_output_pdbs[round_num]),
                )
        return structs_with_both_rounds

    @job_completion
    def parse_prediction_scores(self, prediction_id, root_directory = None, ddg_output_path = None, chains_to_move = None, score_method_id = None, prediction_structure_scores_table = None, prediction_id_field = None, use_multiprocessing = True, verbose = False):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        root_directory = root_directory or self.prediction_data_path
        if not ddg_output_path:
            ddg_output_path = self.find_ddg_output_directory(prediction_id, root_directory)
        structs_with_both_rounds = self.find_structs_with_both_rounds(ddg_output_path)

        scores = []
        output_db3s = {}
        score_method_details = self.get_score_method_details()[score_method_id]
        method_name = score_method_details['MethodName']
        author = score_method_details['Authors']
        job_details = self.get_job_details(prediction_id)
        if not chains_to_move:
            substitution_parameters = json.loads(job_details['JSONParameters'])
            chains_to_move = substitution_parameters['%%chainstomove%%']

        if method_name.startswith('Rescore-') and author == 'Kyle Barlow':
            score_fxn = method_name[8:].lower()
            if score_fxn.startswith('beta'):
                extra_flags = ['-' + score_fxn]
            else:
                extra_flags = []
        else:
            score_fxn = 'interface'
            extra_flags = []

        file_tuples = [] # List of names, contents
        for file_info in job_details['Files']['Input']:
            filename = file_info['Filename']
            if 'params' in filename:
                params_path = os.path.join('/tmp', filename)
                if not os.path.isfile(params_path):
                    with open(params_path, 'w') as f:
                        f.write( file_info['Content'] )
                extra_flags.append('-extra_res_fa')
                extra_flags.append(params_path)

        if len(structs_with_both_rounds) > 0:
            if use_multiprocessing and verbose:
                print '\nOpening rescoring pool for prediction', prediction_id
            elif verbose:
                print 'Individual rescoring each structure for prediction', prediction_id

            if verbose:
                print '%d tasks to run' % (len(structs_with_both_rounds) * 2)
            if use_multiprocessing:
                p = multiprocessing.Pool()
            def finish_rescore(tup):
                round_num, struct_type, output_db3 = tup
                if round_num != None and struct_type != None and output_db3 != None:
                    if round_num not in output_db3s:
                        output_db3s[round_num] = {}
                    assert( struct_type not in output_db3s[round_num] )
                    output_db3s[round_num][struct_type] = output_db3
            for round_num in structs_with_both_rounds:
                wt_pdb, mutant_pdb = structs_with_both_rounds[round_num]
                kwargs = {
                    'rosetta_database_path' : self.rosetta_database_path,
                    'score_fxn' : score_fxn,
                    'round_num' : round_num,
                    'struct_type' : 'wt',
                    'extra_flags' : extra_flags,
                }
                if use_multiprocessing:
                    p.apply_async( interface_calc.rescore_ddg_monomer_pdb, (
                        os.path.abspath(wt_pdb),
                        self.rosetta_scripts_path,
                        chains_to_move,
                        ), kwargs, callback=finish_rescore )
                else:
                    finish_rescore( interface_calc.rescore_ddg_monomer_pdb(
                        os.path.abspath(wt_pdb),
                        self.rosetta_scripts_path,
                        chains_to_move,
                        **kwargs
                    ) )

                kwargs = {
                    'rosetta_database_path' : self.rosetta_database_path,
                    'score_fxn' : score_fxn,
                    'round_num' : round_num,
                    'struct_type' : 'mut',
                    'extra_flags' : extra_flags,
                }
                if use_multiprocessing:
                    p.apply_async( interface_calc.rescore_ddg_monomer_pdb, (
                        os.path.abspath(mutant_pdb),
                        self.rosetta_scripts_path,
                        chains_to_move,
                    ), kwargs, callback=finish_rescore )
                else:
                    finish_rescore( interface_calc.rescore_ddg_monomer_pdb(
                        os.path.abspath(mutant_pdb),
                        self.rosetta_scripts_path,
                        chains_to_move,
                        **kwargs
                    ) )

            if use_multiprocessing:
                p.close()
                p.join()
            if use_multiprocessing and verbose:
                print 'Rescoring pool finished for prediction', prediction_id
                print
            elif verbose:
                print 'Rescoring finished for prediction', prediction_id
                print
        else:
            print '\nlen(structs_with_both_rounds) == 0 for prediction_id:', prediction_id
            print

        for round_num in output_db3s:
            if 'wt' in output_db3s[round_num]:
                wt_output_db3 = output_db3s[round_num]['wt']
            else:
                wt_output_db3 = None

            if 'mut' in output_db3s[round_num]:
                mut_output_db3 = output_db3s[round_num]['mut']
            else:
                mut_output_db3 = None

            if wt_output_db3 and mut_output_db3:
                scores_list = make_scores_list(self.get_score_dict(prediction_id=prediction_id, score_method_id=score_method_id, structure_id=round_num, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field), wt_output_db3, mut_output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
                scores.extend(scores_list)

            if wt_output_db3:
                shutil.rmtree( os.path.dirname(wt_output_db3) )
            if mut_output_db3:
                shutil.rmtree( os.path.dirname(mut_output_db3) )
        return scores

    def add_scores_from_cluster_rescore(self, output_dirs, prediction_structure_scores_table, prediction_id_field, score_method_id):
        structs_with_some_scores = self.get_structs_with_some_scores()

        if isinstance(output_dirs, basestring):
            output_dirs = [output_dirs]

        available_db3_files = {}
        available_db3_files_set = set()
        r = Reporter('looking for output .db3 files in output directories', entries = 'output dirs')
        r.set_total_count( len(output_dirs) )
        for output_dir in output_dirs:
            job_dict_path = os.path.join(os.path.join(output_dir, 'data'), 'job_dict.pickle')

            with open(job_dict_path, 'r') as f:
                job_dict = pickle.load(f)

            for task_name in job_dict:
                prediction_id, round_num, struct_type = split_task_name_path(task_name)
                if struct_type == 'wt':
                    # We will process both wt and mut types below
                    continue
                mut_task_name = task_name
                wt_task_name = mut_task_name.replace('mut', 'wt')
                wt_task_dir = os.path.join(output_dir, wt_task_name)
                mut_task_dir = os.path.join(output_dir, mut_task_name)
                wt_db3_file = os.path.join(wt_task_dir, 'output.db3.gz')
                mut_db3_file = os.path.join(mut_task_dir, 'output.db3.gz')
                if os.path.isfile(wt_db3_file) and os.path.isfile(mut_db3_file):
                    available_db3_files[(prediction_id, round_num)] = (mut_db3_file, wt_db3_file)
                    available_db3_files_set.add( (prediction_id, round_num) )
                else:
                    # Check other output directories for db3 files
                    if not os.path.isfile(wt_db3_file):
                        wt_db3_file = None
                    if not os.path.isfile(mut_db3_file):
                        mut_db3_file = None
                    for inner_output_dir in sorted(output_dirs, key=lambda k: random.random()):
                        if not wt_db3_file:
                            inner_wt_task_dir = os.path.join(inner_output_dir, wt_task_name)
                            inner_wt_db3_file = os.path.join(wt_task_dir, 'output.db3.gz')
                            if os.path.isfile(inner_wt_db3_file):
                                wt_db3_file = inner_wt_db3_file
                        if not mut_db3_file:
                            inner_mut_task_dir = os.path.join(inner_output_dir, mut_task_name)
                            inner_mut_db3_file = os.path.join(mut_task_dir, 'output.db3.gz')
                            if os.path.isfile(inner_mut_db3_file):
                                mut_db3_file = inner_mut_db3_file
                        if wt_db3_file and mut_db3_file:
                            available_db3_files[(prediction_id, round_num)] = (mut_db3_file, wt_db3_file)
                            available_db3_files_set.add( (prediction_id, round_num) )
                            break
            r.increment_report()
        r.done()

        db3_files_to_process = available_db3_files_set.difference(structs_with_some_scores)
        print 'Found %d scores in db3 files needed to add to database' % len(db3_files_to_process)
        r = Reporter('parsing output db3 files and saving scores in database', entries='scores')
        r.set_total_count( len(db3_files_to_process) )

        p = multiprocessing.Pool(min(multiprocessing.cpu_count(), 16)) # Multiprocessing
        def save_scores_helper(return_tuple):
            args, kwargs = return_tuple
            prediction_set_id, prediction_id, scores_list = args
            self.master_scores_list.extend(scores_list)
            if len(self.master_scores_list) >= 150: # Save scores in batches
                self.store_scores_for_many_predictions(None, self.master_scores_list, safe=False, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
                self.master_scores_list = []
            #### self.store_scores(None, prediction_id, scores_list, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field) # Save each score one at a time
            r.increment_report()

        for prediction_id, round_num in db3_files_to_process:
            empty_score_dict = self.get_score_dict(prediction_id=prediction_id, score_method_id=score_method_id, structure_id=round_num, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
            mut_output_db3, wt_output_db3 = available_db3_files[(prediction_id, round_num)]
            p.apply_async( read_db3_scores_helper, (empty_score_dict, prediction_id, round_num, wt_output_db3, mut_output_db3, score_method_id, prediction_structure_scores_table, prediction_id_field), callback=save_scores_helper) # Multiprocessing
            #### save_scores_helper( read_db3_scores_helper(empty_score_dict, prediction_id, round_num, wt_output_db3, mut_output_db3, score_method_id, prediction_structure_scores_table, prediction_id_field) ) # Non-multiprocessing
        p.close() # Multiprocessing
        p.join() # Multiprocessing
        if len(self.master_scores_list) > 0:
            self.store_scores_for_many_predictions(None, self.master_scores_list, safe=False, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
        r.done()

def read_db3_scores_helper(empty_score_dict, prediction_id, round_num, wt_output_db3, mut_output_db3, score_method_id, prediction_structure_scores_table, prediction_id_field):
    wt_tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
    mut_tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
    new_wt_output_db3_path = os.path.join(wt_tmp_dir, os.path.basename(wt_output_db3))
    shutil.copy(wt_output_db3, new_wt_output_db3_path)
    wt_output_db3 = unzip_file(new_wt_output_db3_path)

    new_mut_output_db3_path = os.path.join(mut_tmp_dir, os.path.basename(mut_output_db3))
    shutil.copy(mut_output_db3, new_mut_output_db3_path)
    mut_output_db3 = unzip_file(new_mut_output_db3_path)

    scores_list = make_scores_list(empty_score_dict, wt_output_db3, mut_output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)

    shutil.rmtree(wt_tmp_dir)
    shutil.rmtree(mut_tmp_dir)
    return ( (None, prediction_id, scores_list), {'prediction_structure_scores_table' : prediction_structure_scores_table, 'prediction_id_field' : prediction_id_field} )

def make_task_name_path(max_prediction_id, prediction_id, round_num, struct_type):
    max_num_digits = len( str(max_prediction_id) )
    format_str = '%%0%dd' % max_num_digits
    prediction_id_str = format_str % prediction_id
    task_name = prediction_id_str[0]
    for digit in prediction_id_str[1:]:
        task_name = os.path.join(task_name, digit)
    end_task_name = '%d/%d/%s' % (prediction_id, round_num, struct_type)
    return os.path.join(task_name, end_task_name)

def split_task_name_path(task_name_path):
    if '_' in task_name_path:
        data = task_name_path.split('_')
    else:
        data = task_name_path.split('/')
    return (long(data[-3]), int(data[-2]), data[-1])

def process_cluster_rescore_helper(ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path, max_prediction_id, extra_flags):
    task_name = make_task_name_path(max_prediction_id, prediction_id, round_num, struct_type)
    prediction_id_data_dir = os.path.join(job_data_dir, task_name)
    pdb_name = os.path.basename(pdb_path)
    if pdb_name.endswith('.gz'):
        pdb_name_zipped = pdb_name
    else:
        pdb_name_zipped = pdb_name + '.gz'
    data_pdb_path = os.path.join(prediction_id_data_dir, pdb_name_zipped)

    if os.path.isdir(prediction_id_data_dir) and not os.path.isfile(data_pdb_path):
        shutil.rmtree(prediction_id_data_dir)

    if not os.path.isdir(prediction_id_data_dir):
        try:
            os.makedirs(prediction_id_data_dir)
        except OSError:
            pass

        new_pdb = os.path.join(prediction_id_data_dir, pdb_name)
        shutil.copy(pdb_path, new_pdb)
        if not new_pdb.endswith('.gz'):
            new_pdb = zip_file_with_gzip(new_pdb)
        assert( os.path.abspath(data_pdb_path) == os.path.abspath(new_pdb) )

    rel_pdb_path = os.path.relpath(data_pdb_path, job_output_dir)
    arg_dict = {
        '-parser:script_vars' : ['chainstomove=%s' % chains_to_move,  'currentscorefxn=%s' % score_fxn],
        '-s' : rel_pdb_path,
        '-parser:protocol' : rel_protocol_path,
        'FLAGLIST' : extra_flags,
    }
    return (task_name, arg_dict)

def make_scores_list(empty_score_dict, wt_output_db3, mut_output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = None, prediction_id_field = None):
    # "left" structure has struct_id=1
    wtl_score = add_scores_from_db3_file(wt_output_db3, 1, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'WildTypeLPartner', score_method_id, prediction_structure_scores_table, prediction_id_field))
    wtr_score = add_scores_from_db3_file(wt_output_db3, 2, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'WildTypeRPartner', score_method_id, prediction_structure_scores_table, prediction_id_field))
    wtc_score = add_scores_from_db3_file(wt_output_db3, 3, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'WildTypeComplex', score_method_id, prediction_structure_scores_table, prediction_id_field))
    ml_score = add_scores_from_db3_file(mut_output_db3, 1, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'MutantLPartner', score_method_id, prediction_structure_scores_table, prediction_id_field))
    mr_score = add_scores_from_db3_file(mut_output_db3, 2, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'MutantRPartner', score_method_id, prediction_structure_scores_table, prediction_id_field))
    mc_score = add_scores_from_db3_file(mut_output_db3, 3, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'MutantComplex', score_method_id, prediction_structure_scores_table, prediction_id_field))
    return [wtl_score, wtr_score, wtc_score, ml_score, mr_score, mc_score]

def add_scores_from_db3_file(db3_file, struct_id, round_num, score_dict):
    if db3_file.endswith('.gz'):
        tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
        new_db3_path = os.path.join(tmp_dir, os.path.basename(db3_file))
        shutil.copy(db3_file, new_db3_path)
        db3_file = unzip_file(new_db3_path)
    else:
        tmp_dir = None

    conn = sqlite3.connect(db3_file)
    c = conn.cursor()
    score_types = set()
    for row in c.execute('SELECT score_type_name FROM score_types WHERE batch_id=1'):
        score_types.add(row[0])
    for empty_score_type in score_dict:
        if empty_score_type == 'total':
            score_type = 'total_score'
        else:
            score_type = empty_score_type
        score_rows = [row for row in c.execute('SELECT structure_scores.score_value from structure_scores INNER JOIN score_types ON score_types.batch_id=structure_scores.batch_id AND score_types.score_type_id=structure_scores.score_type_id WHERE structure_scores.struct_id=%d AND score_type_name="%s"' % (struct_id, score_type))]
        if len(score_rows) == 1:
            score_dict[empty_score_type] = float( score_rows[0][0] )
        elif len(score_rows) > 1:
            raise Exception('Matched too many score rows')
    score_dict['StructureID'] = round_num
    conn.close()
    if tmp_dir:
        shutil.rmtree(tmp_dir)
    return score_dict

def unzip_to_tmp_dir(prediction_id, zip_file):
    tmp_dir = tempfile.mkdtemp(prefix='unzip_to_tmp_')
    unzip_path = os.path.join(tmp_dir, '%d-ddg' % prediction_id)
    os.makedirs(unzip_path)
    with zipfile.ZipFile(zip_file, 'r') as job_zip:
        job_zip.extractall(unzip_path)
    return tmp_dir
