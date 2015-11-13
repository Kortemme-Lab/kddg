from ppi_api import BindingAffinityDDGInterface
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
from ddglib.ppi_api import get_interface
import klab.cluster_template.parse_settings as parse_settings
from klab.Reporter import Reporter
from klab.MultiWorker import MultiWorker
from klab.fs.zip_util import zip_file_with_gzip
from klab.cluster_template.write_run_file import process as write_run_file
import time
import getpass

# Constants for cluster runs
rosetta_scripts_xml_file = os.path.join('ddglib', 'score_partners.xml')
output_db3 = 'output.db3'
setup_run_with_multiprocessing = True

def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a DDGMonomerInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(DDGMonomerInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)

class DDGMonomerInterface(BindingAffinityDDGInterface):
    def __init__(self, passwd = None, username = None, hostname = None, rosetta_scripts_path = None, rosetta_database_path = None):
        super(DDGMonomerInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)
        self.rescore_args = {} # Stores arguments for rescoring step, to be dumped later

    def get_prediction_ids_with_scores(self, prediction_set_id, score_method_id = None):
        '''Returns a set of all prediction_ids that already have an associated score in prediction_set_id
        '''
        if score_method_id:
            q = "SELECT DISTINCT PredictionPPIID FROM PredictionPPIStructureScore INNER JOIN PredictionPPI ON PredictionPPI.ID=PredictionPPIStructureScore.PredictionPPIID WHERE PredictionPPI.PredictionSet='%s' AND PredictionPPIStructureScore.ScoreMethodID=%d" % (prediction_set_id, score_method_id)
        else:
            q = "SELECT DISTINCT PredictionPPIID FROM PredictionPPIStructureScore INNER JOIN PredictionPPI ON PredictionPPI.ID=PredictionPPIStructureScore.PredictionPPIID WHERE PredictionPPI.PredictionSet='%s'" % (prediction_set_id)
        return_set = set()
        for r in self.DDG_db.execute_select(q):
            return_set.add( r['PredictionPPIID'] )
        return return_set

    def create_cluster_run_rescore_dir(self, output_dir):
        settings = parse_settings.get_dict()
        job_name = '%s-%s_rescore_ddg_monomer' % (time.strftime("%y%m%d"), getpass.getuser())
        job_output_dir = os.path.join(output_dir, job_name)
        if not os.path.isdir(job_output_dir):
            os.makedirs(job_output_dir)

        num_jobs = 0
        for ddg_output_path in self.rescore_args:
            for arg_tuple in self.rescore_args[ddg_output_path]:
                num_jobs += 1
            
        settings['scriptname'] = 'rescore_ddg_monomer'
        settings['tasks_per_process'] = 50
        settings['numjobs'] = num_jobs
        settings['mem_free'] = '1.4G'
        settings['appname'] = 'rosetta_scripts'
        settings['output_dir'] = job_output_dir
        job_dict = {}
        job_data_dir = os.path.join(job_output_dir, 'data')
        if not os.path.isdir(job_data_dir):
            os.makedirs(job_data_dir)
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

        for ddg_output_path in self.rescore_args:
            for prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type in self.rescore_args[ddg_output_path]:
                if setup_run_with_multiprocessing:
                    worker.addJob( (ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path ) )
                else:
                    task_name, arg_dict = process_cluster_rescore_helper( ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path )
                    job_dict[task_name] = arg_dict
                    r.increment_report()

        if setup_run_with_multiprocessing:
            worker.finishJobs()
            for task_name, arg_dict in worker.data:
                job_dict[task_name] = arg_dict
        else:
            r.done()
        write_run_file(settings, job_dict = job_dict)

    def add_rescore_cluster_run(self, ddg_output_path, chains_to_move, score_method_id, prediction_id):
        structs_with_both_rounds = self.find_structs_with_both_rounds(ddg_output_path)

        score_method_details = self.get_score_method_details()[score_method_id]
        method_name = score_method_details['MethodName']
        author = score_method_details['Authors']
        if method_name.startswith('Rescore-') and author == 'Kyle Barlow':
            score_fxn = method_name[8:].lower()
        else:
            score_fxn = 'interface'

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
                'wt'
            ])
            self.rescore_args[ddg_output_path].append([
                prediction_id,
                os.path.abspath(mutant_pdb),
                chains_to_move,
                score_fxn,
                round_num,
                'mut'
            ])

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
        root_directory = root_directory or self.prediction_data_path
        all_prediction_ids = [x for x in self.get_prediction_ids(prediction_set_id)]
        all_prediction_ids_set = set()
        for prediction_id in all_prediction_ids:
            all_prediction_ids_set.add( prediction_id )
        scored_prediction_ids_set = self.get_prediction_ids_with_scores(prediction_set_id, score_method_id = score_method_id)
        prediction_ids = [x for x in all_prediction_ids_set.difference(scored_prediction_ids_set)]

        random.shuffle( prediction_ids )
        print '%d prediction_ids to process' % len(prediction_ids)
        for prediction_id in prediction_ids:
            ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
            if os.path.isdir( ddg_output_path ):
                self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)

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
    def parse_prediction_scores(self, prediction_id, root_directory = None, ddg_output_path = None, chains_to_move = None, score_method_id = None):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        root_directory = root_directory or self.prediction_data_path
        if not ddg_output_path:
            ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
        structs_with_both_rounds = self.find_structs_with_both_rounds(ddg_output_path)

        scores = []
        output_db3s = {}
        score_method_details = self.get_score_method_details()[score_method_id]
        method_name = score_method_details['MethodName']
        author = score_method_details['Authors']
        if not chains_to_move:
            job_details = self.get_job_details(prediction_id)
            substitution_parameters = json.loads(job_details['JSONParameters'])
            chains_to_move = substitution_parameters['%%chainstomove%%']
        if method_name.startswith('Rescore-') and author == 'Kyle Barlow':
            score_fxn = method_name[8:].lower()
        else:
            score_fxn = 'interface'

        if len(structs_with_both_rounds) > 0:
            print 'Opening rescoring pool for prediction', prediction_id
            print '%d tasks to run' % (len(structs_with_both_rounds) * 2)
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
                }
                p.apply_async( interface_calc.rescore_ddg_monomer_pdb, (
                    os.path.abspath(wt_pdb),
                    self.rosetta_scripts_path,
                    chains_to_move,
                    ), kwargs, callback=finish_rescore )

                kwargs = {
                    'rosetta_database_path' : self.rosetta_database_path,
                    'score_fxn' : score_fxn,
                    'round_num' : round_num,
                    'struct_type' : 'mut',
                }
                p.apply_async( interface_calc.rescore_ddg_monomer_pdb, (
                    os.path.abspath(mutant_pdb),
                    self.rosetta_scripts_path,
                    chains_to_move,
                    ), kwargs, callback=finish_rescore )
            p.close()
            p.join()
            print 'Rescoring pool finished for prediction', prediction_id
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
                # "left" structure has struct_id=1
                wtl_score = self.add_scores_from_db3_file(wt_output_db3, 1, round_num, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeLPartner', score_method_id = score_method_id))
                wtr_score = self.add_scores_from_db3_file(wt_output_db3, 2, round_num, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeRPartner', score_method_id = score_method_id))
                wtc_score = self.add_scores_from_db3_file(wt_output_db3, 3, round_num, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'WildTypeComplex', score_method_id = score_method_id))
                ml_score = self.add_scores_from_db3_file(mut_output_db3, 1, round_num, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantLPartner', score_method_id = score_method_id))
                mr_score = self.add_scores_from_db3_file(mut_output_db3, 2, round_num, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantRPartner', score_method_id = score_method_id))
                mc_score = self.add_scores_from_db3_file(mut_output_db3, 3, round_num, self.get_score_dict(prediction_id = prediction_id, structure_id = round_num, score_type = 'MutantComplex', score_method_id = score_method_id))
                scores.extend([wtl_score, wtr_score, wtc_score, ml_score, mr_score, mc_score])

            if wt_output_db3:
                shutil.rmtree( os.path.dirname(wt_output_db3) )
            if mut_output_db3:
                shutil.rmtree( os.path.dirname(mut_output_db3) )
        return scores


    def rescore_ddg_monomer_pdb(self, pdb_file, prediction_id, score_method_id):

        return output_db3

    def add_scores_from_db3_file(self, db3_file, struct_id, round_num, score_dict):
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
        return score_dict

def process_cluster_rescore_helper(ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path):
    task_name = '%d_%d_%s' % (prediction_id, round_num, struct_type)
    prediction_id_data_dir = os.path.join(job_data_dir, task_name)
    pdb_name = os.path.basename(pdb_path)
    if pdb_name.endswith('.gz'):
        pdb_name_zipped = pdb_name
    else:
        pdb_name_zipped = pdb_name + '.gz'
    data_pdb_path = os.path.join(prediction_id_data_dir, pdb_name_zipped)
    if not os.path.isdir(prediction_id_data_dir):
        os.mkdir(prediction_id_data_dir)
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
    }
    return (task_name, arg_dict)
