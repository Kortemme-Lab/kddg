def fill_empty_score_dict(score_dict, prediction_id, structure_id, score_type, score_method_id, prediction_structure_scores_table, prediction_id_field):
    d = copy.deepcopy(score_dict)
    if prediction_id_field != None:
        d[prediction_id_field] = prediction_id
    d['ScoreMethodID'] = score_method_id
    d['ScoreType'] = score_type
    d['StructureID'] = structure_id
    return d
