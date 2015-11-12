from ddg_monomer_ppi_api import DDGMonomerInterface
from api_layers import *
from db_api import ddG
import os
import json
import sqlite3
import interface_calc
import sys
import shutil
import multiprocessing
import random


def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a DDGMonomerInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(RescoreddGDatabase, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)

class RescoreddGDatabase(DDGMonomerInterface):
    def __init__(self, passwd = None, username = None, hostname = None, rosetta_scripts_path = None, rosetta_database_path = None):
        super(RescoreddGDatabase, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)
        self.rescore_args = []

    def add_rescore_cluster_run(self, ddg_output_path, chains_to_move, score_method_id):
        structs_with_both_rounds = self.find_structs_with_both_rounds(ddg_output_path)

        score_method_details = self.get_score_method_details()[score_method_id]
        method_name = score_method_details['MethodName']
        author = score_method_details['Authors']
        if method_name.startswith('Rescore-') and author == 'Kyle Barlow':
            score_fxn = method_name[8:].lower()
        else:
            score_fxn = 'interface'

        for round_num in structs_with_both_rounds:
            wt_pdb, mutant_pdb = structs_with_both_rounds[round_num]
            self.rescore_args.append([
                os.path.abspath(wt_pdb),
                chains_to_move,
                score_fxn,
                round_num,
                'wt'
            ])
            self.rescore_args.append([
                os.path.abspath(mutant_pdb),
                chains_to_move,
                score_fxn,
                round_num,
                'mut'
            ])
        print self.rescore_args
    @job_completion
    def store_scores(self, prediction_set, prediction_id, scores):
        '''Stores scores for one prediction.
           scores should be a list of dicts suitable for database storage e.g. PredictionStructureScore or
           PredictionPPIStructureScore records.
           This function uses a transaction so if any of the insertions fail then they are all rolled back.
           '''
        self._check_prediction(prediction_id, prediction_set)
        self._check_scores_for_main_fields(scores, prediction_id)
        self._check_score_fields(scores)

        try:
            con = self.DDG_db.connection
            with con:
                db_cursor = con.cursor()
                for score in scores:
                    sql, params, record_exists = self.DDG_db.create_insert_dict_string(self._get_prediction_structure_scores_table(), score, PKfields = [self._get_prediction_id_field(), 'ScoreMethodID', 'ScoreType', 'StructureID'], check_existing = True)
                    if not record_exists:
                        db_cursor.execute(sql, params)
        except Exception, e:
            raise colortext.Exception('Failed to insert scores for Prediction #{0}: "{1}".\n{2}'.format(prediction_id, str(e), traceback.format_exc()))
