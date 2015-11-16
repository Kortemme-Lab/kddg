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
from klab.fs.zip_util import zip_file_with_gzip, unzip_file
from klab.fs.io import sanitize_filename
from klab.cluster_template.write_run_file import process as write_run_file
import time
import getpass
import tempfile
import cPickle as pickle
import copy

# Constants for cluster runs
rosetta_scripts_xml_file = os.path.join('ddglib', 'score_partners.xml')
output_db3 = 'output.db3'
setup_run_with_multiprocessing = True
tmpdir_location = '/dbscratch/%s/tmp' % getpass.getuser()

def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
    '''This is the function that should be used to get a DDGMonomerInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(DDGMonomerInterface, passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)

class DDGMonomerInterface(BindingAffinityDDGInterface):
    def __init__(self, passwd = None, username = None, hostname = None, rosetta_scripts_path = None, rosetta_database_path = None):
        super(DDGMonomerInterface, self).__init__(passwd = passwd, username = username, hostname = hostname, rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path)
        self.rescore_args = {} # Stores arguments for rescoring step, to be dumped later
        self.all_score_types = ['WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex']
        self.all_score_types_index = {}
        for i, score_type in enumerate(self.all_score_types):
            self.all_score_types_index[score_type] = i
        self.master_scores_list = []

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

    def get_prediction_ids_without_scores(self, prediction_set_id, score_method_id = None):
        all_prediction_ids = [x for x in self.get_prediction_ids(prediction_set_id)]
        all_prediction_ids_set = set()
        for prediction_id in all_prediction_ids:
            all_prediction_ids_set.add( prediction_id )
        scored_prediction_ids_set = self.get_prediction_ids_with_scores(prediction_set_id, score_method_id = score_method_id)
        return [x for x in all_prediction_ids_set.difference(scored_prediction_ids_set)]
        
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
                for prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type in self.rescore_args[ddg_output_path]:
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
                for prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type in self.rescore_args[ddg_output_path]:
                    if setup_run_with_multiprocessing:
                        worker.addJob( (ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path, max_prediction_id ) )
                    else:
                        task_name, arg_dict = process_cluster_rescore_helper( ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path, max_prediction_id )
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
        prediction_ids = self.get_prediction_ids_without_scores(prediction_set_id, score_method_id = score_method_id)

        random.shuffle( prediction_ids )
        print '%d prediction_ids to process' % len(prediction_ids)
        if setup_cluster_run:
            job_name = sanitize_filename(prediction_set_id)
            r = Reporter('fetching job details for prediction_ids', entries='job details')
            r.set_total_count(len(prediction_ids))
            for prediction_id in prediction_ids:
                ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
                job_details = self.get_job_details(prediction_id)
                substitution_parameters = json.loads(job_details['JSONParameters'])
                chains_to_move = substitution_parameters['%%chainstomove%%']
                self.add_rescore_cluster_run(ddg_output_path, chains_to_move, score_method_id, prediction_id)
                r.increment_report()
            r.done()
            self.create_cluster_run_rescore_dir( os.path.join(tmpdir_location, 'cluster_run'), passed_job_name = job_name )
        else:
            processed_count = 0
            for prediction_id in prediction_ids:
                ddg_output_path = os.path.join(root_directory, '%d-ddg' % prediction_id)
                if os.path.isdir( ddg_output_path ):
                    self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)
                    processed_count += 1
                if max_prediction_ids_to_process and processed_count >= max_prediction_ids_to_process:
                    print 'Breaking early; processed %d prediction ids' % processed_count
                    break

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
    def parse_prediction_scores(self, prediction_id, root_directory = None, ddg_output_path = None, chains_to_move = None, score_method_id = None, prediction_structure_scores_table = None, prediction_id_field = None):
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
                scores_list = make_scores_list(self.get_score_dict(prediction_id=prediction_id, score_method_id=score_method_id, structure_id=round_num, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field), wt_output_db3, mut_output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
                scores.extend(scores_list)

            if wt_output_db3:
                shutil.rmtree( os.path.dirname(wt_output_db3) )
            if mut_output_db3:
                shutil.rmtree( os.path.dirname(mut_output_db3) )
        return scores

    def add_scores_from_cluster_rescore(self, output_dir, prediction_structure_scores_table, prediction_id_field, score_method_id):
        job_dict_path = os.path.join(os.path.join(output_dir, 'data'), 'job_dict.pickle')

        with open(job_dict_path, 'r') as f:
            job_dict = pickle.load(f)

        settings = parse_settings.get_dict()
        rosetta_scripts_path = settings['local_rosetta_installation_path'] + '/source/bin/' + 'rosetta_scripts' + settings['local_rosetta_binary_type']

        DDGdb = self.DDG_db

        prediction_ids_and_structs_score_count = {}
        for row in DDGdb.execute_select("SELECT %s, ScoreType, StructureID FROM %s WHERE ScoreType IN ('WildTypeLPartner', 'WildTypeRPartner', 'WildTypeComplex', 'MutantLPartner', 'MutantRPartner', 'MutantComplex') AND ScoreMethodID=%d" % (prediction_id_field, prediction_structure_scores_table, score_method_id)):
            prediction_id = long(row[prediction_id_field])
            score_type = row['ScoreType']
            structure_id = int(row['StructureID'])
            if (prediction_id, structure_id) not in prediction_ids_and_structs_score_count:
                prediction_ids_and_structs_score_count[(prediction_id, structure_id)] = 0
            prediction_ids_and_structs_score_count[(prediction_id, structure_id)] += 1
        structs_with_some_scores = set()
        for prediction_id, structure_id in prediction_ids_and_structs_score_count:
            if prediction_ids_and_structs_score_count[(prediction_id, structure_id)] > 0:
                structs_with_some_scores.add( (prediction_id, structure_id) )
                if prediction_ids_and_structs_score_count[(prediction_id, structure_id)] != 6:
                    print 'Missing data:', prediction_id, structure_id, prediction_ids_and_structs_score_count[(prediction_id, structure_id)]

        available_db3_files = {}
        available_db3_files_set = set()
        for task_name in job_dict:
            prediction_id, round_num, struct_type = split_task_name_path(task_name)
            if struct_type == 'wt':
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

def process_cluster_rescore_helper(ddg_output_path, prediction_id, pdb_path, chains_to_move, score_fxn, round_num, struct_type, job_data_dir, job_output_dir, rel_protocol_path, max_prediction_id):
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
    }
    return (task_name, arg_dict)

def fill_empty_score_dict(score_dict, prediction_id, structure_id, score_type, score_method_id, prediction_structure_scores_table, prediction_id_field):
    d = copy.deepcopy(score_dict)
    if prediction_id_field != None:
        d[prediction_id_field] = prediction_id
    d['ScoreMethodID'] = score_method_id
    d['ScoreType'] = score_type
    d['StructureID'] = structure_id
    return d

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
