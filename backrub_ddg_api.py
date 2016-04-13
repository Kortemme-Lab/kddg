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
import numpy as np

from sqlalchemy import and_, update
from sqlalchemy.orm import load_only, Load

import klab.cluster_template.parse_settings as parse_settings
from klab.Reporter import Reporter
from klab.MultiWorker import MultiWorker
from klab.fs.zip_util import zip_file_with_gzip, unzip_file
from klab.fs.io import sanitize_filename
from klab.cluster_template.write_run_file import process as write_run_file

from ddg_monomer_ppi_api import DDGMonomerInterface, read_db3_scores_helper, total_seconds
from api_layers import *
from db_api import ddG
import db_schema as dbmodel
from util import fill_empty_score_dict

from zip_to_shared import find_all_files

DeclarativeBase = dbmodel.DeclarativeBase

# Constants for cluster runs
output_db3 = 'output.db3'
date_format_string = '%Y-%m-%d %H:%M:%S'
setup_run_with_multiprocessing = True

def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, port = 3306, file_content_buffer_size = None):
    '''This is the function that should be used to get a BackrubDDGInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(
        BackrubDDGInterface,
        passwd = passwd, username = username, hostname = hostname,
        rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path,
        port = port,
        file_content_buffer_size = file_content_buffer_size,
    )

class BackrubDDGInterface(DDGMonomerInterface):
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
        self.structs_with_all_scores = set()
        self.prediction_id_status_cache = {}

    def add_scores_from_cluster_rescore(self, output_dirs, prediction_structure_scores_table, prediction_id_field, score_method_id):

        if isinstance(output_dirs, basestring):
            output_dirs = [output_dirs]

        available_db3_files = {}
        available_db3_files_set = set()
        r = Reporter('looking for output .db3 files in output directories', entries = 'output dirs')
        r.set_total_count( len(output_dirs) )
        output_log_files = set()
        len_output_dirs = {}
        for output_dir in output_dirs:
            wt_job_dict_path = os.path.join(os.path.join(output_dir, 'data-6'), 'job_dict-6.pickle')
            mut_job_dict_path = os.path.join(os.path.join(output_dir, 'data-7'), 'job_dict-7.pickle')

            with open(wt_job_dict_path, 'r') as f:
                wt_job_dict = pickle.load(f)
            with open(mut_job_dict_path, 'r') as f:
                mut_job_dict = pickle.load(f)

            wt_keys = sorted(wt_job_dict.keys())
            mut_keys = sorted(mut_job_dict.keys())
            task_id = 1
            for wt_key, mut_key in zip(wt_keys, mut_keys):
                assert( wt_key.split('/')[0] == mut_key.split('/')[0] ) # Make sure that sorting has matched predictions
                assert( wt_key.split('/')[1].split('-')[0] == mut_key.split('/')[1].split('-')[0] ) # Make sure that sorting has matched structure numbers
                assert( wt_key.endswith( 'rescore_sep_wt' ) )
                assert( mut_key.endswith( 'rescore_sep_mutant' ) )
                prediction_id = long( wt_key.split('/')[0] )
                round_num = long( wt_key.split('/')[1].split('-')[0] )
                wt_task_dir = os.path.join(output_dir, wt_key)
                mut_task_dir = os.path.join(output_dir, mut_key)
                wt_db3_file = os.path.join(wt_task_dir, 'output.db3.gz')
                mut_db3_file = os.path.join(mut_task_dir, 'output.db3.gz')
                if os.path.isdir( os.path.join(output_dir, str(prediction_id)) ):
                    l = len( find_all_files( os.path.join(output_dir, str(prediction_id)) ) )
                    if l not in len_output_dirs:
                        len_output_dirs[l] = set()
                    len_output_dirs[l].add( (output_dir, prediction_id) )
                if os.path.isfile(wt_db3_file) and os.path.isfile(mut_db3_file) and not self.struct_has_all_scores(prediction_id, round_num):
                    available_db3_files[(prediction_id, round_num)] = (mut_db3_file, wt_db3_file)
                    available_db3_files_set.add( (prediction_id, round_num) )

                # Find output log files
                minimize_dir = os.path.join(output_dir, '%d/%04d-0000-minimize/' % (prediction_id, round_num))
                if os.path.isdir( minimize_dir ):
                    log_files = [f for f in os.listdir(minimize_dir) if '.o' in f and f.endswith('.%d' % task_id)]
                    if len( log_files ) == 0:
                        pass
                        # print 'Missing log file in:', minimize_dir
                    elif len( log_files ) == 1:
                        log_file = log_files[0]
                        prediction_set_id = log_file[:log_file.find('_run.py')] # Remove past _run.py
                        output_log_files.add( (prediction_set_id, os.path.join(minimize_dir, log_file), prediction_id, round_num) )
                    else:
                        raise Exception('Found too many log files in: ' + str(minimize_dir))
                task_id += 1
            r.increment_report()
        r.done()

        db3_files_to_process = available_db3_files_set
        print 'Found %d scores in db3 files needed to add to database' % len(db3_files_to_process)
        output_dirs_to_zip = len_output_dirs[sorted( len_output_dirs.keys() )[-1]] # Take only folders with the most structures for zipping
        self.parse_db3_files_to_process(db3_files_to_process, available_db3_files, output_log_files, score_method_ids, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field, output_dirs_to_zip = output_dirs_to_zip )

    def parse_db3_files_to_process(self, db3_files_to_process, available_db3_files, output_log_files, score_method_id, prediction_structure_scores_table = None, prediction_id_field = None, output_dirs_to_zip = [], use_multiprocessing = True, verbose = True):
        if not prediction_structure_scores_table:
            prediction_structure_scores_table = self.prediction_structure_scores_table
        if not prediction_id_field:
            prediction_id_field = self.prediction_table + 'ID'
        if verbose:
            r = Reporter('parsing output files', entries='output files')
            r.set_total_count( len(output_log_files) )
        for prediction_set_id, output_log_file, prediction_id, round_num in output_log_files:
            self.update_prediction_id_status(prediction_set_id, output_log_file, prediction_id, round_num, verbose = verbose)
            if verbose:
                r.increment_report()
        if verbose:
            r.done()

        self.finalize_update_prediction_id_status(verbose = verbose)

        if verbose:
            r = Reporter('parsing output db3 files and saving scores in database', entries='scores')
            r.set_total_count( len(db3_files_to_process) )
        if use_multiprocessing:
            p = multiprocessing.Pool(min(multiprocessing.cpu_count(), 16)) # Multiprocessing
            self.master_scores_list = []
        def save_scores_helper(return_tuple):
            args, kwargs = return_tuple
            prediction_set_id, prediction_id, scores_list = args
            if use_multiprocessing:
                self.master_scores_list.extend(scores_list)
                if len(self.master_scores_list) >= 150: # Save scores in batches
                        self.store_scores_for_many_predictions(None, self.master_scores_list, safe=False, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
                        self.master_scores_list = []
            else:
                self.store_scores(None, prediction_id, scores_list, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field) # Save each score one at a time
            if verbose:
                r.increment_report()

        for prediction_id, round_num in db3_files_to_process:
            empty_score_dict = self.get_score_dict(prediction_id=prediction_id, score_method_id=score_method_id, structure_id=round_num, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
            mut_output_db3, wt_output_db3 = available_db3_files[(prediction_id, round_num)]
            if use_multiprocessing:
                p.apply_async( read_db3_scores_helper, (empty_score_dict, prediction_id, round_num, wt_output_db3, mut_output_db3, score_method_id, prediction_structure_scores_table, prediction_id_field), callback=save_scores_helper) # Multiprocessing
            else:
                save_scores_helper( read_db3_scores_helper(empty_score_dict, prediction_id, round_num, wt_output_db3, mut_output_db3, score_method_id, prediction_structure_scores_table, prediction_id_field) ) # Non-multiprocessing
        if use_multiprocessing:
            p.close() # Multiprocessing
            p.join() # Multiprocessing
        if use_multiprocessing and len(self.master_scores_list) > 0:
            self.store_scores_for_many_predictions(None, self.master_scores_list, safe=False, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
        if verbose:
            r.done()

        if len(output_dirs_to_zip) > 0:
            if verbose:
                r = Reporter('zipping prediction id output directories', entries='directories')
                r.set_total_count( len(output_dirs_to_zip) )
            for output_dir, prediction_id in output_dirs_to_zip:
                if self.zip_prediction_id(prediction_id, output_dir):
                    if verbose:
                        r.increment_report()
                else:
                    if verbose:
                        r.decrement_total_count()
            if verbose:
                r.done()

    def update_prediction_id_status(self, prediction_set_id, output_log_file, prediction_id, round_num, verbose = True):
        # Looks for output file in root directory and reads for job status
        unfinished_prediction_ids = self.get_unfinished_prediction_ids(prediction_set_id)
        if len(unfinished_prediction_ids) == 0 or prediction_id not in unfinished_prediction_ids:
            return

        # Constants
        starting_time_string = 'Starting time:'
        ending_time_string = 'Ending time:'

        # Parse output file
        starting_time = None
        ending_time = None
        return_code = None
        return_codes = []
        max_virtual_memory_usage = None
        elapsed_time = None
        status = 'active'
        with open(output_log_file, 'r') as f:
            for line in f:
                if starting_time == None and line.startswith(starting_time_string):
                    starting_time = line[len(starting_time_string):].strip()
                elif line.startswith(ending_time_string):
                    ending_time = line[len(ending_time_string):].strip()
                elif 'return code' in line:
                    return_codes.append( int(line.strip().split()[-1].strip()) )
                elif 'virtual memory usage' in line:
                    line = line.strip()
                    if line.endswith('G'):
                        byte_factor = 1.0
                    elif line.endswith('M'):
                        byte_factor = 1000.0
                    else:
                        raise Exception()
                    line = line[:-1]
                    virtual_memory_usage = float(line.strip().split()[-1].strip()) / byte_factor
                    if max_virtual_memory_usage == None or virtual_memory_usage > max_virtual_memory_usage:
                        max_virtual_memory_usage = virtual_memory_usage
        if starting_time and ending_time:
            starting_time_dt = datetime.datetime.strptime(starting_time, date_format_string)
            ending_time_dt = datetime.datetime.strptime(ending_time, date_format_string)
            elapsed_time = total_seconds(ending_time_dt - starting_time_dt) / 60.0 # 60 seconds in a minute
        else:
            starting_time_dt = None
            ending_time_dt = None

        if len(return_codes) > 0:
            status = 'done'
            for return_code in return_codes:
                if return_code != 0:
                    status = 'failed'

        if starting_time_dt != None:
            if prediction_id not in self.prediction_id_status_cache:
                self.prediction_id_status_cache[prediction_id] = {
                    'starting_time_dt' : [],
                    'ending_time_dt' : [],
                    'status' : [],
                    'max_virtual_memory_usage' : [],
                    'elapsed_time' : [],
                    'return_code' : [],
                }

            self.prediction_id_status_cache[prediction_id]['starting_time_dt'].append(
                starting_time_dt
            )
            self.prediction_id_status_cache[prediction_id]['ending_time_dt'].append(
                ending_time_dt
            )
            self.prediction_id_status_cache[prediction_id]['status'].append(
                status
            )
            self.prediction_id_status_cache[prediction_id]['max_virtual_memory_usage'].append(
                max_virtual_memory_usage
            )
            self.prediction_id_status_cache[prediction_id]['elapsed_time'].append(
                elapsed_time
            )
            self.prediction_id_status_cache[prediction_id]['return_code'].append(
                return_code
            )


    def finalize_update_prediction_id_status(self, verbose = True):
        if len(self.prediction_id_status_cache) == 0:
            return

        if verbose:
            r = Reporter('updating prediction id statii in database', entries='prediction IDs')
            r.set_total_count( len(self.prediction_id_status_cache) )

        tsession = self.importer.session
        for prediction_id in self.prediction_id_status_cache:
            prediction_record = tsession.query(self.PredictionTable).filter(self.PredictionTable.ID == prediction_id).one()
            cache = self.prediction_id_status_cache[prediction_id]
            prediction_record.StartDate = datetime.datetime.strftime(min(cache['starting_time_dt']), date_format_string)
            prediction_record.EndDate = datetime.datetime.strftime(max(cache['ending_time_dt']), date_format_string)
            if 'failed' in cache['status']:
                prediction_record.Status = 'failed'
            else:
                prediction_record.Status = 'done'
            prediction_record.maxvmem = max( cache['max_virtual_memory_usage'])
            prediction_record.DDGTime = np.mean( cache['elapsed_time'] )
            prediction_record.ERRORS = str( max(cache['return_code']) )
            tsession.flush()
            tsession.commit()
            if verbose:
                r.increment_report()
        if verbose:
            r.done()


    def struct_has_all_scores(self, prediction_id, round_num, expect_n = 6, verbose = True, prediction_id_field = None, prediction_structure_scores_table = None):
        if (prediction_id, round_num) in self.structs_with_all_scores:
            return True

        prediction_id_field = prediction_id_field or self.prediction_table + 'ID'
        prediction_structure_scores_table = prediction_structure_scores_table or self.prediction_table + 'StructureScore'

        prediction_ids_and_structs_score_count = {}
        for row in self.DDG_db.execute_select("SELECT %s, ScoreType, StructureID FROM %s WHERE ScoreType IN ('WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex') AND %s=%d" % (prediction_id_field, prediction_structure_scores_table, prediction_id_field, prediction_id)):
            assert( prediction_id == long(row[prediction_id_field]) )
            score_type = row['ScoreType']
            structure_id = int(row['StructureID'])
            if (prediction_id, structure_id) not in prediction_ids_and_structs_score_count:
                prediction_ids_and_structs_score_count[(prediction_id, structure_id)] = 0
            prediction_ids_and_structs_score_count[(prediction_id, structure_id)] += 1
        structs_with_all_scores = set()
        for prediction_id, structure_id in prediction_ids_and_structs_score_count:
            if prediction_ids_and_structs_score_count[(prediction_id, structure_id)] == expect_n:
                self.structs_with_all_scores.add( (prediction_id, structure_id) )
            if verbose and prediction_ids_and_structs_score_count[(prediction_id, structure_id)] != expect_n:
                print 'Missing data:', prediction_id, structure_id, prediction_ids_and_structs_score_count[(prediction_id, structure_id)]
        return (prediction_id, round_num) in self.structs_with_all_scores

    def zip_prediction_id(self, prediction_id, root_dir, delete_if_exists = False, remove_output_dir = True):
        zip_path = '/kortemmelab/shared/DDG/ppijobs/%d.zip' % prediction_id
        ddg_dir = os.path.join(root_dir, str(prediction_id))
        all_files = find_all_files(ddg_dir)
        if os.path.isfile(zip_path):
            if delete_if_exists:
                os.remove(zip_path)
            else:
                if remove_output_dir:
                    shutil.rmtree( ddg_dir )
                return False

        filtered_files = []
        for f in all_files:
            if '.pdb.gz' in f:
                if 'rescore_sep_' in f:
                    filtered_files.append( f )
            else:
                filtered_files.append(f)
        with zipfile.ZipFile(zip_path, 'w') as job_zip:
            for f in filtered_files:
                f_path = os.path.join(ddg_dir, f)
                job_zip.write(f_path, arcname=f)
        # Delete prediction ID directory if it contains expected number of files
        if remove_output_dir:
            shutil.rmtree( ddg_dir )
        return True
