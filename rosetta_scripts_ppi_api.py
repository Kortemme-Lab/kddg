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
import re

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

def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, port = 3306):
    '''This is the function that should be used to get a RosettaScriptsInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(RosettaScriptsInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port)

class RosettaScriptsInterface(BindingAffinityDDGInterface):
    def __init__(self, passwd = None, username = None, hostname = None, rosetta_scripts_path = None, rosetta_database_path = None, port = 3306):
        super(RosettaScriptsInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path, port = port)

    @job_completion
    def extract_data(self, prediction_set_id, root_directory = None, force = False, score_method_id = None):
        root_directory = root_directory or self.prediction_data_path

        ### Find all prediction_ids that need to have updated states
        self.update_prediction_id_status(prediction_set_id, root_directory)

        ### Find all prediction_ids with partial score data and remove this data

        ### Find all prediction_ids with missing scores and setup rescoring or rescore on the fly
        prediction_ids = self.get_prediction_ids_without_scores(prediction_set_id, score_method_id = score_method_id)

        self.add_scores_from_db3s(root_directory, score_method_id = score_method_id)

    def add_scores_from_db3s(self, output_dir, expect_n = None, round_num = 1, score_method_id = None):
        prediction_structure_scores_table = self.prediction_table + 'StructureScore'
        prediction_id_field = self.prediction_table + 'ID'

        structs_with_some_scores = self.get_structs_with_some_scores(expect_n = expect_n)

        available_db3_files = {}
        available_db3_files_set = set()
        # Take the last step data directory job pickle
        data_dirs = sorted( [d for d in os.listdir(output_dir) if d.startswith('data')] )
        m = re.match( '(?:data-)(\d+)', data_dirs[-1] )
        if m:
            job_dict_file = 'job_dict-%d.pickle' % int(m.group(1))
        else:
            job_dict_file = 'job_dict.pickle'
        job_dict_path = os.path.join(os.path.join(output_dir, data_dirs[-1]), job_dict_file)

        with open(job_dict_path, 'r') as f:
            job_dict = pickle.load(f)

        for task_name in job_dict:
            task_name = str(task_name)
            prediction_id = long(task_name)
            db3_file = os.path.join(output_dir, os.path.join(task_name, 'output.db3.gz'))
            if os.path.isfile(db3_file):
                available_db3_files[ (prediction_id, round_num) ] = db3_file
                available_db3_files_set.add( (prediction_id, round_num) )
            else:
                print 'Missing db3 file:', db3_file

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
                self.store_scores_for_many_predictions(None, self.master_scores_list, safe=False, prediction_structure_scores_table = prediction_structure_scores_table , prediction_id_field = prediction_id_field)
                self.master_scores_list = []
            r.increment_report()

        for prediction_id, round_num in db3_files_to_process:
            empty_score_dict = self.get_score_dict(prediction_id=prediction_id, score_method_id = score_method_id, structure_id=round_num, prediction_structure_scores_table = prediction_structure_scores_table , prediction_id_field = prediction_id_field)

            output_db3 = available_db3_files[(prediction_id, round_num)]
            p.apply_async( read_db3_scores_helper, (empty_score_dict, prediction_id, round_num, output_db3, score_method_id, prediction_structure_scores_table, prediction_id_field), callback=save_scores_helper) # Multiprocessing
        p.close() # Multiprocessing
        p.join() # Multiprocessing
        if len(self.master_scores_list) > 0:
            self.store_scores_for_many_predictions(None, self.master_scores_list, safe=False, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
        r.done()

def read_db3_scores_helper(empty_score_dict, prediction_id, round_num, output_db3,score_method_id, prediction_structure_scores_table, prediction_id_field):
    tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
    new_output_db3_path = os.path.join(tmp_dir, os.path.basename(output_db3))
    shutil.copy(output_db3, new_output_db3_path)
    output_db3 = unzip_file(new_output_db3_path)

    scores_list = make_scores_list(empty_score_dict, output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)

    shutil.rmtree(tmp_dir)
    return ( (None, prediction_id, scores_list), {'prediction_structure_scores_table' : prediction_structure_scores_table, 'prediction_id_field' : prediction_id_field} )

def make_scores_list(empty_score_dict, db3_file, prediction_id, round_num, score_method_id, prediction_structure_scores_table = None, prediction_id_field = None):
    # Fill empty score dict
    score_dict = fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'DDG', score_method_id, prediction_structure_scores_table, prediction_id_field)

    if db3_file.endswith('.gz'):
        tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
        new_db3_path = os.path.join(tmp_dir, os.path.basename(db3_file))
        shutil.copy(db3_file, new_db3_path)
        db3_file = unzip_file(new_db3_path)
    else:
        tmp_dir = None

    # There should only be one structure in this db3 file output
    assert( len(c.execute('SELECT struct_id FROM structures')) == 1)
    struct_id = 1

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
