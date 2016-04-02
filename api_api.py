from ppi_api import BindingAffinityDDGInterface, get_interface
from api_layers import *
from db_api import ddG
import db_schema as dbmodel
import numpy as np
import sys

DeclarativeBase = dbmodel.DeclarativeBase

def get_interface(passwd, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None, port = 3306, file_content_buffer_size = None):
    '''This is the function that should be used to get a DDGMonomerInterface object. It hides the private methods
       from the user so that a more traditional object-oriented API is created.'''
    return GenericUserInterface.generate(
        APIInterface,
        passwd = passwd, username = username, hostname = hostname,
        rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path,
        port = port,
        file_content_buffer_size = file_content_buffer_size,
    )

class APIInterface(BindingAffinityDDGInterface):
    def __init__(
            self,
            passwd = None, username = None, hostname = None,
            rosetta_scripts_path = None, rosetta_database_path = None,
            port = 3306,
            file_content_buffer_size = None
    ):
        super(APIInterface, self).__init__(
            passwd = passwd, username = username, hostname = hostname,
            rosetta_scripts_path = rosetta_scripts_path, rosetta_database_path = rosetta_database_path,
            port = port,
            file_content_buffer_size = file_content_buffer_size
        )

        self.run_data = None

    def set_run_data(self, run_data):
        self.run_data = run_data

    def update_prediction_id_status(self, prediction_id, verbose = True):
        if prediction_id not in self.run_data:
            return

        # Keys to use for lookup
        starting_time_key = 'Starting time'
        ending_time_key = 'Ending time'
        elapsed_time_key = 'Runtime'
        return_code_key = 'Task return code'
        memory_key = 'Max virtual memory usage'

        d = self.run_data[prediction_id]

        all_elapsed_times = [d['structIDs'][struct_id][elapsed_time_key] for struct_id in d['structIDs'] if elapsed_time_key in d['structIDs'][struct_id]]
        if len(all_elapsed_times) > 0:
            elapsed_time = np.mean( all_elapsed_times )
        else:
            elapsed_time = None

        if starting_time_key in d:
            starting_time = d[starting_time_key]
        else:
            starting_time = None

        if ending_time_key in d:
            ending_time = d[ending_time_key]
        else:
            ending_time = None

        if return_code_key in d:
            if d[return_code_key] != 0:
                status = 'failed'
                errors = '%d' % d[return_code_key]
            else:
                status = 'done'
                errors = ''
        else:
            status = 'active'
            errors = None

        if memory_key in d:
            virtual_memory_usage = d[memory_key]
        else:
            virtual_memory_usage = None

        tsession = self.importer.session
        prediction_record = tsession.query(self.PredictionTable).filter(self.PredictionTable.ID == prediction_id).one()
        prediction_record.StartDate = starting_time
        prediction_record.EndDate = ending_time
        prediction_record.Status = status
        prediction_record.maxvmem = virtual_memory_usage
        prediction_record.DDGTime = elapsed_time
        prediction_record.ERRORS = errors
        tsession.flush()
        tsession.commit()

    @job_completion
    def extract_data(self, prediction_set_id, root_directory = None, force = False, score_method_id = None):
        prediction_ids = self.get_prediction_ids(prediction_set_id)
        prediction_ids_with_data = set( self.run_data.keys() )
        missing_prediction_ids = set( prediction_ids )
        missing_prediction_ids = missing_prediction_ids.difference( prediction_ids_with_data )
        if len( missing_prediction_ids ) > 0:
            print '%d prediction ids are missing from stored data dictionary' % len( missing_prediction_ids )

        for prediction_id in self.run_data:
            self.update_prediction_id_status(prediction_id)
            self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)

    @job_completion
    def parse_prediction_scores(self, prediction_id, score_method_id = None, root_directory = None, prediction_structure_scores_table = None, prediction_id_field = None, verbose = False):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        scores = []

        for struct_id in self.run_data[prediction_id]['structIDs']:
            score_dict = self.get_score_dict(score_type='DDG', prediction_id=prediction_id, score_method_id=score_method_id, structure_id=struct_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
            subtracted_scores = subtract_wt_from_mutant(
                self.run_data[prediction_id]['structIDs'][struct_id]['WT'],
                self.run_data[prediction_id]['structIDs'][struct_id]['Mutant']
            )
            for key in score_dict:
                if key in subtracted_scores:
                    score_dict[key] = subtracted_scores[key]
            scores.append( score_dict )

        return scores

def subtract_wt_from_mutant(wt_scores, mutant_scores):
    mutant_keys = set( mutant_scores.keys() )
    wt_keys = set( wt_scores.keys() )
    common_keys = mutant_keys.intersection( wt_keys )

    return_dict = {
        key : mutant_scores[key] - wt_scores[key] for key in common_keys
    }

    if 'Sum' in return_dict:
        return_dict['total'] = return_dict['Sum']
        del return_dict['Sum']
    else:
        raise Exception("Missing 'Sum' ddG score")

    return return_dict
