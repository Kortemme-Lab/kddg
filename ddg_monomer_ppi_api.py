from ppi_api import BindingAffinityDDGInterface
from api_layers import *
from db_api import ddG
import os
import json
import sqlite3
import interface_calc
import sys

def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a DDGMonomerInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(DDGMonomerInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)

class DDGMonomerInterface(BindingAffinityDDGInterface):
    @job_completion
    def extract_data(self, prediction_set_id, root_directory = None, force = False, score_method_id = None):
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
        print prediction_set_id
        root_directory = root_directory or self.prediction_data_path
        prediction_ids = self.get_prediction_ids(prediction_set_id)
        for prediction_id in prediction_ids:
            ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
            if os.path.isdir( ddg_output_path ):
                print self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)
                sys.exit(0)        


    @job_completion
    def parse_prediction_scores(self, prediction_id, root_directory = None, score_method_id = None):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        root_directory = root_directory or self.prediction_data_path
        ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
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
        scores = []
        if len(structs_with_both_rounds) > 0:
            for round_num in structs_with_both_rounds:
                wt_pdb, mutant_pdb = structs_with_both_rounds[round_num]

                output_db3 = self.rescore_ddg_monomer_pdb(wt_pdb, prediction_id)
                # "left" structure has struct_id=1
                wtl_score = self.add_scores_from_db3_file(output_db3, 1, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeLPartner', score_method_id = score_method_id))
                wtr_score = self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeRPartner', score_method_id = score_method_id)
                wtc_score = self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeComplex', score_method_id = score_method_id)
                ml_score = self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantLPartner', score_method_id = score_method_id)
                mr_score = self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantRPartner', score_method_id = score_method_id)
                mc_score = self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantComplex', score_method_id = score_method_id)
                ddg_score = self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'DDG', score_method_id = score_method_id)
                shutil.rmtree( os.path.dirname(output_db3) )

                scores.extend([wtl_score, wtr_score, wtc_score, ml_score, mr_score, mc_score])
        return scores
    
    def rescore_ddg_monomer_pdb(self, pdb_file, prediction_id):
        job_details = self.get_job_details(prediction_id)
        substitution_parameters = json.loads(job_details['JSONParameters'])
        output_db3 = interface_calc.rescore_ddg_monomer_pdb(
            os.path.abspath(pdb_file),
            self.rosetta_scripts_path,
            substitution_parameters['%%chainstomove%%'],
            rosetta_database_path = self.rosetta_database_path,
        )
        return output_db3

    def add_scores_from_db3_file(self, db3_file, struct_id, score_dict):
        conn = sqlite3.connect(db3_file)
        c = conn.cursor()
        print db3_file
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
        print score_dict
        sys.exit(0)
        return score_dict

