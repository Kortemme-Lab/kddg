from ppi_api import BindingAffinityDDGInterface, get_interface
from api_layers import *
from db_api import ddG
import db_schema as dbmodel
import numpy as np

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

    def update_prediction_id_status(self, prediction_id, verbose = True):
        if prediction_id not in self.run_data:
            return

        # Keys to use for lookup
        runtime_key = 'Runtime'
        starting_time_key = 'starting_time'
        ending_time_key = 'ending_time'
        elapsed_time_key = 'elapsed_time'
        return_code_key = 'return_code'
        memory_key = 'memory_usage'

        task_dict = self.run_data[prediction_id]
        all_keys = set()
        for struct_id in task_dict:
            for key in task_dict[struct_id]:
                all_keys.add(key)

        if elapsed_time_key in all_keys:
            elapsed_time = np.mean( [d[struct_id][elapsed_time_key] for struct_id in d if elapsed_time_key in d[struct_id]] )
        else:
            elapsed_time = None

        if starting_time_key in all_keys:
            starting_time = min( [d[struct_id][starting_time_key] for struct_id in d if starting_time_key in d[struct_id]] )
        else:
            starting_time = None

        if ending_time_key in all_keys:
            ending_time = min( [d[struct_id][ending_time_key] for struct_id in d if ending_time_key in d[struct_id]] )
        else:
            ending_time = None

        if return_code_key in all_keys:
            status = 'done'
            errors = ''
            for return_code in [d[struct_id][return_code_key] for struct_id in d if return_code_key in d[struct_id]]:
                if return_code != 0:
                    status = 'failed'
                    errors += '%d ' % return_code
            errors = errors.strip()
        else:
            status = 'active'
            errors = None

        if memory_key in all_keys:
            virtual_memory_usage = max( [d[struct_id][memory_key] for struct_id in d if memory_key in d[struct_id]] )
        else:
            virtual_memory_usage = None

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
        prediction_ids = ppi_api.get_prediction_ids(prediction_set_id)
        prediction_ids_with_data = set( self.run_data.keys() )
        missing_prediction_ids = set( prediction_ids )
        missing_prediction_ids = missing_prediction_ids.difference( prediction_ids_with_data )
        if len( missing_prediction_ids ) > 0:
            print '%d prediction ids are missing from stored data dictionary' % len( missing_prediction_ids )

        for prediction_id in self.run_data:
            self.update_prediction_id_status(prediction_id)
            self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)

    @job_completion
    def parse_prediction_scores(self, prediction_id, score_method_id = None, prediction_structure_scores_table = None, prediction_id_field = None, verbose = False):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        scores = []

        struct_id_key = None

        raise Exception()
        for struct_id in self.run_data[struct_id_key]:
            scores_list = self.make_scores_list(self.get_score_dict(prediction_id=prediction_id, score_method_id=score_method_id, structure_id=struct_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field), wt_output_db3, mut_output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = prediction_structure_scores_table, prediction_id_field = prediction_id_field)
            scores.extend(scores_list)

        return scores

    def make_scores_list(self, empty_score_dict, wt_output_db3, mut_output_db3, prediction_id, round_num, score_method_id, prediction_structure_scores_table = None, prediction_id_field = None):
        # "left" structure has struct_id=1
        ddg_score = add_scores_from_saved_data(wt_output_db3, 1, round_num, fill_empty_score_dict(empty_score_dict, prediction_id, round_num, 'WildTypeLPartner', score_method_id, prediction_structure_scores_table, prediction_id_field))
        return [wtl_score, wtr_score, wtc_score, ml_score, mr_score, mc_score]

    def add_scores_from_saved_data(struct_id, round_num, score_dict):
        raise Exception()
        score_dict['StructureID'] = struct_id
        return score_dict
