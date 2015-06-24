#!/usr/bin/python2.4
# encoding: utf-8
"""
monomer_analysis.py

This is the code I used to generate breakdown analyses for the benchmarks website and for a group meeting in early 2015.
This code should be refactored into an analysis layer of db_api.py along with the code in analysis.py and dbstats.py.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

class AnalysisBreakdown(object):

    def __init__(self, amino_acids, pdb_details, predictions, analysis_datasets):
        self.amino_acids = amino_acids
        self.pdb_details = pdb_details

        # split the predictions over mutations to/from glycine and proline and other predictions
        GP = set(['G', 'P'])
        single_mutation_GP_predictions = {}
        single_mutation_no_GP_predictions = {}
        multiple_mutation_predictions = {}
        for p, details in predictions.iteritems():
            if not details.get('WTAA'):
                multiple_mutation_predictions[p] = details
            elif (details['WTAA'] in GP or details['MutantAA'] in GP):
                single_mutation_GP_predictions[p] = details
            else:
                single_mutation_no_GP_predictions[p] = details

        self.predictions = predictions
        self.single_mutation_GP_predictions = single_mutation_GP_predictions
        self.single_mutation_no_GP_predictions = single_mutation_no_GP_predictions
        self.multiple_mutation_predictions = multiple_mutation_predictions
        print('%d total predictions: %d Single No GP, %d Single GP, %d multiple' % (len(predictions), len(single_mutation_no_GP_predictions), len(single_mutation_GP_predictions), len(multiple_mutation_predictions)))
        self.analysis_datasets = analysis_datasets


    def analyze_subset_all(self, analysis_subset, scoring_method, prediction_details_map = {}):
        colortext.message('ANALYZING SUBSET analyze_subset_all of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.predictions, prediction_details_map)

    def analyze_subset_single_no_GP(self, analysis_subset, scoring_method):
        colortext.message('ANALYZING SUBSET single_mutation_no_GP_predictions of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.single_mutation_no_GP_predictions)

    def analyze_subset_single_GP(self, analysis_subset, scoring_method):
        colortext.message('ANALYZING SUBSET single_mutation_GP_predictions of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.single_mutation_GP_predictions)

    def analyze_subset_multiple(self, analysis_subset, scoring_method):
        colortext.message('ANALYZING SUBSET multiple_mutation_predictions of %s' % analysis_subset)
        return self._analyze_subset_sub(analysis_subset, scoring_method, self.multiple_mutation_predictions)

    def _analyze_subset_sub(self, analysis_subset, scoring_method, predictions, prediction_details_map = {}):
        'Analyzes a subset using the main datapoints.'

        analysis_dataset = self.analysis_datasets[analysis_subset]
        xvalues = []
        yvalues = []
        print('ID,Experimental,Predicted')

        # for benchmarks paper
        ddgapio = MonomericStabilityDDGInterface().DDG_db

        for prediction_id, details in sorted(predictions.iteritems()):
            if prediction_id in analysis_dataset:
                ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                predicted_score = details[scoring_method]
                if predicted_score != None:

                    # for the Benchmarks paper...
                    if False:
                        mutations = ddgapio.execute_select('''
SELECT ExperimentMutation.*
FROM Prediction
INNER JOIN ExperimentMutation ON ExperimentMutation.ExperimentID=Prediction.ExperimentID
WHERE Prediction.ID = %s''', parameters=(prediction_id,))
                        if len(mutations) == 1:
                            if abs(ExperimentalDDG - predicted_score) < 0.5:
                                if abs(ExperimentalDDG) > 2:
                                    if (mutations[0]['WildTypeAA'] == 'D'):
                                        print(mutations[0]['WildTypeAA'], mutations[0]['MutantAA'], ExperimentalDDG, predicted_score, prediction_id)


                    xvalues.append(ExperimentalDDG)
                    yvalues.append(predicted_score)
                    if prediction_details_map:
                        assert(prediction_details_map.get(int(prediction_id)))
                        #print('%f,%f,%s,%s' % (ExperimentalDDG, predicted_score, prediction_id, ','.join(map(str, prediction_details_map[int(prediction_id)]))))
                    else:
                        pass
                        #print('%f,%f,%s,%s' % (ExperimentalDDG, predicted_score, prediction_id))

        print('*' * 30)
        print(min(xvalues), max(xvalues), min(yvalues), max(yvalues))
        colortext.message('Analyzing %d values for dataset %s using scoring method %s.' % (len(xvalues), analysis_subset, scoring_method))
        stats = None
        if (len(xvalues) >= 8):
            stats = get_xy_dataset_correlations(xvalues, yvalues)
            #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
            colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
            colortext.warning('''
                              <span class="DDG_correlation_score ">%0.2f</span> /
                              <span class="DDG_stability_classification_score ">%0.2f</span> /
                              <span class="DDG_MAE_score">%0.2f</span>''' % (stats['pearsonr'][0], stats['fraction_correct'], stats['MAE']))
            pprint.pprint(stats)
        else:
            colortext.warning('Not enough data.')
        return stats


    def analyze_subset_by_specific_resolutions(self, analysis_subset, scoring_method, bins = [1.5, 2.0, 2.5]):
        ''' Analyzes a subset using specific PDB resolution bins.
            The bins argument defines the resolution bins.
            The first bin is <x where x is the smallest value in bins.
            The last bin is >x where x is the largest value in bins.
        '''

        bins = sorted(set(bins))
        assert(bins[0] > 0)
        limits = [0] + bins + [max(bins[-1], 100)]
        xvalues = {None : []}
        yvalues = {None : []}
        for x in range(0, len(limits) - 1):
            xvalues[(limits[x], limits[x + 1])] = []
            yvalues[(limits[x], limits[x + 1])] = []

        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    resolution = self.pdb_details[details['pPDB']]['Resolution']
                    if resolution:
                        for x in range(0, len(limits) - 1):
                            if limits[x] <= resolution < limits[x + 1]:
                                xvalues[(limits[x], limits[x + 1])].append(ExperimentalDDG)
                                yvalues[(limits[x], limits[x + 1])].append(predicted_score)
                                break
                    else:
                        xvalues[None].append(ExperimentalDDG)
                        yvalues[None].append(predicted_score)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for bin %s.' % (len(xvalues[k]), str(k)))
                results[k] = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                #pprint.pprint(results[k])
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results

    def analyze_subset_by_binned_resolutions(self, analysis_subset, scoring_method, num_bins = 9):
        ''' Analyzes a subset using PDB resolution bins. This function attempts to break up the result set into num_bins
            bins of somewhat equal size.
            Additionally, there is a special bin for PDB files with null resolution.'''

        assert(num_bins > 1)

        count = 0
        xyvalues = {None: []}
        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    resolution = self.pdb_details[details['pPDB']]['Resolution']
                    xyvalues[resolution] = xyvalues.get(resolution, [])
                    xyvalues[resolution].append((ExperimentalDDG, predicted_score))
                    count += 1

        # determine the ideal number per bin, ignoring the bin with null resolution
        ideal_per_bin = (count - len(xyvalues.get(None))) / num_bins

        xvalues = {None : [p[0] for p in xyvalues[None]]}
        yvalues = {None : [p[1] for p in xyvalues[None]]}
        current_bin = []
        bin_start = 0
        for bin_end in sorted(xyvalues.keys()):
            if bin_end != None:
                current_bin += xyvalues[bin_end]
                if (len(current_bin) > ideal_per_bin) or (bin_end == sorted(xyvalues.keys())[-1]):
                    xvalues[(bin_start, bin_end)] = [p[0] for p in current_bin]
                    yvalues[(bin_start, bin_end)] = [p[1] for p in current_bin]
                    bin_start = bin_end
                    current_bin = []
        assert(sum([len(xvalues[k]) for k in xvalues]) == count)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for bin %s.' % (len(xvalues[k]), str(k)))
                pprint.pprint(get_xy_dataset_correlations(xvalues[k], yvalues[k]))
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results


    def analyze_subset_by_binned_chain_length(self, analysis_subset, scoring_method, num_bins = 9):
        ''' Analyzes a subset using PDB resolution bins. This function attempts to break up the result set into num_bins
            bins of somewhat equal size.
            Additionally, there is a special bin for PDB files with null resolution.'''

        assert(num_bins > 1)

        pdb_details = self.pdb_details
        count = 0
        xyvalues = {}
        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    chain = details['Chain']
                    ppdb = details['pPDB']
                    chain_length = pdb_details[ppdb]['chains'][chain]
                    xyvalues[chain_length] = xyvalues.get(chain_length, [])
                    xyvalues[chain_length].append((ExperimentalDDG, predicted_score))
                    count += 1

        # determine the ideal number per bin, ignoring the bin with null resolution
        ideal_per_bin = (count) / num_bins

        xvalues = {}#None : [p[0] for p in xyvalues[None]]}
        yvalues = {}#sNone : [p[1] for p in xyvalues[None]]}
        current_bin = []
        bin_start = max(xyvalues.keys())
        print(sorted(xyvalues.keys()))
        print(sorted(xyvalues.keys(), reverse = True))
        for bin_end in sorted(xyvalues.keys(), reverse = True):
            current_bin += xyvalues[bin_end]
            if (len(current_bin) > ideal_per_bin) or (bin_end == sorted(xyvalues.keys(), reverse = True)[-1]):
                xvalues[(bin_start, bin_end)] = [p[0] for p in current_bin]
                yvalues[(bin_start, bin_end)] = [p[1] for p in current_bin]
                bin_start = bin_end
                current_bin = []
        print(sum([len(xvalues[k]) for k in xvalues]), count)
        assert(sum([len(xvalues[k]) for k in xvalues]) == count)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for bin %s.' % (len(xvalues[k]), str(k)))
                pprint.pprint(get_xy_dataset_correlations(xvalues[k], yvalues[k]))
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results

        # We can derive the following data:
        #   TM, Resolution, XRay per PDB
        #   for prediction_id, d in predictions.iteritems():
        #     assoc_pdb = pdb_details[d['pPDB']]
        #     d['TM'] = assoc_pdb['TM'] == 1
        #     d['XRay'] = assoc_pdb['XRay']
        #     d['Resolution'] = assoc_pdb['Resolution']
        # Derive GP (if wt or mutant is glycine or proline)
        # Derive WTPolarity, WTAromaticity, MutantPolarity, MutantAromaticity
        # Derive SL, LS, SS, LL
        # Derive ChainLength: prediction['ChainLength'] = pdbs[pc['PDBFileID']]['chains'][pc['Chain']]

        #print(self.pdb_details[details['pPDB']])
        #        print(self.pdb_details[details['pPDB']]['chains'][details['Chain']])
        #        #{u'XRay': True, u'chains': {u'A': 680}, u'Technique': u'X-RAY DIFFRACTION', u'Resolution': 2.4, u'TM': 1}


    def analyze_subset_by_exposure(self, analysis_subset, scoring_method, cut_off = 0.25):
        ''' Analyzes a subset using exposure of the wildtype residue as calculated by DSSP.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.
        '''

        assert(0 <= cut_off <= 1.0)
        xvalues = {None : [], 'B' : [], 'E' : []}
        yvalues = {None : [], 'B' : [], 'E' : []}

        main_subsets = [self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        all_results = dict.fromkeys(main_subset_names)
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {'B' : [], 'E' : []}
                yvalues = {'B' : [], 'E' : []}

                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):
                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            exposure = details.get('Exposure')
                            if exposure != None:
                                if exposure <= cut_off:
                                    xvalues['B'].append(ExperimentalDDG)
                                    yvalues['B'].append(predicted_score)
                                else:
                                    xvalues['E'].append(ExperimentalDDG)
                                    yvalues['E'].append(predicted_score)
                            else:
                                print(details)
                                xvalues[None].append(ExperimentalDDG)
                                yvalues[None].append(predicted_score)

                results = {}
                colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for exposure type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
                        #pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results
        return all_results


    def analyze_subset_by_wildtype_charge(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        xvalues = {None : [], 'C' : [], 'P' : [], 'H' : []}
        yvalues = {None : [], 'C' : [], 'P' : [], 'H' : []}

        amino_acids = self.amino_acids
        CAA = [aa for aa in amino_acids if amino_acids[aa]['Polarity'] == 'C']
        PAA = [aa for aa in amino_acids if amino_acids[aa]['Polarity'] == 'P']
        HAA = [aa for aa in amino_acids if amino_acids[aa]['Polarity'] == 'H']

        main_subsets = [self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        all_results = dict.fromkeys(main_subset_names)
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {'Change' : [], 'Polar/Charged' : [], 'Hydrophobic/Non-polar' : []}
                yvalues = {'Change' : [], 'Polar/Charged' : [], 'Hydrophobic/Non-polar' : []}

                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):

                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            wtaa = details.get('WTAA')
                            if wtaa:
                                mutaa = details.get('MutantAA')
                                if ((wtaa in CAA or wtaa in PAA) and (mutaa in HAA)) or ((mutaa in CAA or mutaa in PAA) and (wtaa in HAA)):
                                    # change in charge
                                    xvalues['Change'].append(ExperimentalDDG)
                                    yvalues['Change'].append(predicted_score)
                                elif (wtaa in CAA or wtaa in PAA) and (mutaa in CAA or mutaa in PAA):
                                    xvalues['Polar/Charged'].append(ExperimentalDDG)
                                    yvalues['Polar/Charged'].append(predicted_score)
                                elif (wtaa in HAA) and (mutaa in HAA):
                                    xvalues['Hydrophobic/Non-polar'].append(ExperimentalDDG)
                                    yvalues['Hydrophobic/Non-polar'].append(predicted_score)
                                else:
                                     raise Exception('Should not reach here.')

                results = {}
                colortext.message('Analyzing %s=%d values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), sum([len(xvalues[k]) for k in sorted(xvalues.keys())]), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for polarity type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
                        #pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results
        return all_results


    def analyze_subset_by_aromaticity(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        xvalues = {None : [], '-' : [], 'L' : [], 'R' : []}
        yvalues = {None : [], '-' : [], 'L' : [], 'R' : []}

        amino_acids = self.amino_acids
        analysis_dataset = self.analysis_datasets[analysis_subset]
        for prediction_id, details in sorted(self.single_mutation_no_GP_predictions.iteritems()):
            if prediction_id in analysis_dataset:
                predicted_score = details[scoring_method]
                if predicted_score != None:
                    ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                    wtaa = details.get('WTAA')
                    if wtaa:
                        aromaticity = amino_acids[wtaa]['Aromaticity']
                        xvalues[aromaticity].append(ExperimentalDDG)
                        yvalues[aromaticity].append(predicted_score)
                    else:
                        xvalues[None].append(ExperimentalDDG)
                        yvalues[None].append(predicted_score)

        results = {}
        colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
        for k in sorted(xvalues.keys()):
            if len(xvalues[k]) >= 8:
                colortext.warning('Analyzing %d values for aromaticity type %s.' % (len(xvalues[k]), str(k)))
                results[k] = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                #pprint.pprint(results[k])
            else:
                colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
        return results


    def analyze_subset_by_mutation_size(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        #xvalues = {None : [], 'XX' : [], 'SL' : [], 'LS' : []}
        #yvalues = {None : [], 'XX' : [], 'SL' : [], 'LS' : []}

        amino_acids = self.amino_acids
        SAA = [aa for aa in amino_acids if amino_acids[aa]['Size'] == 'small']
        LAA = [aa for aa in amino_acids if amino_acids[aa]['Size'] == 'large']

        amino_acid_volumes = {}
        for aa, details in amino_acids.iteritems():
            amino_acid_volumes[aa] = details['van_der_Waals_volume']
        assert(len(amino_acid_volumes) == 20)

        main_subsets = [self.predictions, self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['all_mutations', 'single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        main_subsets = [self.predictions]
        main_subset_names = ['all_mutations']
        all_results = dict.fromkeys(main_subset_names)
        failed_cases = 0
        non_cases = set()
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {'XX' : [], 'SL' : [], 'LS' : [], 'Failed' : []}
                yvalues = {'XX' : [], 'SL' : [], 'LS' : [], 'Failed' : []}

                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):
                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            wtaa = details.get('WTAA')
                            if wtaa:
                                mutaa = details.get('MutantAA')

                                if details.get('MutationIsReversed') != None:
                                    # todo: this is not currently an issue but it will be once we include reverse mutations
                                    if details['MutationIsReversed']:
                                        # Note: For reverse mutations, we need to switch the order since we only store the forward mutation
                                        wtaa, mutaa = mutaa, wtaa

                                if wtaa == mutaa:
                                    colortext.warning('Error in analysis: Record mutating %s to %s in Prediction #%s.' % (wtaa, mutaa, prediction_id))
                                    error = True
                                elif amino_acid_volumes[wtaa] < amino_acid_volumes[mutaa]:
                                    xvalues['SL'].append(ExperimentalDDG)
                                    yvalues['SL'].append(predicted_score)
                                elif amino_acid_volumes[wtaa] > amino_acid_volumes[mutaa]:
                                    xvalues['LS'].append(ExperimentalDDG)
                                    yvalues['LS'].append(predicted_score)
                                else:
                                    assert(amino_acid_volumes[wtaa] == amino_acid_volumes[mutaa])
                                    xvalues['XX'].append(ExperimentalDDG)
                                    yvalues['XX'].append(predicted_score)

                                if False:
                                    if wtaa in SAA and mutaa in LAA:
                                        xvalues['SL'].append(ExperimentalDDG)
                                        yvalues['SL'].append(predicted_score)
                                    elif wtaa in LAA and mutaa in SAA:
                                        xvalues['LS'].append(ExperimentalDDG)
                                        yvalues['LS'].append(predicted_score)
                                    else:
                                        non_cases.add('%s->%s' % (wtaa, mutaa))
                                        xvalues['XX'].append(ExperimentalDDG)
                                        yvalues['XX'].append(predicted_score)
                            else:
                                failed_cases += 1
                                #raise Exception('Should not reach here.')
                                xvalues['Failed'].append(ExperimentalDDG)
                                yvalues['Failed'].append(predicted_score)

                results = {}
                colortext.message('Analyzing %s=%d values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), sum([len(xvalues[k]) for k in sorted(xvalues.keys())]), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for mutation size type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE']))

                        colortext.warning('''
                                          <span class="DDG_correlation_score ">%0.2f</span> /
                                          <span class="DDG_stability_classification_score ">%0.2f</span> /
                                          <span class="DDG_MAE_score">%0.2f</span>''' % (stats['pearsonr'][0], stats['fraction_correct'], stats['MAE']))

                        pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results

        print(sorted(non_cases))
        return all_results


    def analyze_subset_by_secondary_structure(self, analysis_subset, scoring_method):
        ''' Analyzes subsets partitioned by residue charge - charged, polar, or hydrophobic.
            The cut_off argument defines the definition of burial - wildtype positions with an exposure <= cut_off are
            considered buried.'''

        xvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}
        yvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}

        main_subsets = [self.single_mutation_no_GP_predictions, self.single_mutation_GP_predictions, self.multiple_mutation_predictions]
        main_subset_names = ['single_mutation_no_GP_predictions', 'single_mutation_GP_predictions', 'multiple_mutation_predictions']
        all_results = dict.fromkeys(main_subset_names)
        for subid in range(len(main_subsets)):
            main_subset = main_subsets[subid]
            if main_subset != self.multiple_mutation_predictions:
                xvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}
                yvalues = {None : [], 'S' : [], 'H' : [], 'O' : []}
                colortext.message('ANALYZING SUBSET %s of %s' % (main_subset_names[subid], analysis_subset))
                analysis_dataset = self.analysis_datasets[analysis_subset]
                for prediction_id, details in sorted(main_subset.iteritems()):
                    if prediction_id in analysis_dataset:
                        predicted_score = details[scoring_method]
                        if predicted_score != None:
                            ExperimentalDDG = analysis_dataset[prediction_id]['ExperimentalDDG']
                            dssp = details.get('DSSP')
                            xvalues[dssp].append(ExperimentalDDG)
                            yvalues[dssp].append(predicted_score)

                results = {}
                colortext.message('Analyzing %s values for dataset %s using scoring method %s.' % ('+'.join(map(str, [len(xvalues[k]) for k in sorted(xvalues.keys())])), analysis_subset, scoring_method))
                for k in sorted(xvalues.keys()):
                    if len(xvalues[k]) >= 8:
                        colortext.warning('Analyzing %d values for secondary structure type %s.' % (len(xvalues[k]), str(k)))
                        stats = get_xy_dataset_correlations(xvalues[k], yvalues[k])
                        results[k] = stats
                        #colortext.warning('%d, %0.2f, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['gammaCC'], stats['MAE'], stats['fraction_correct']))
                        colortext.warning('%d, %0.2f, %0.2f, %0.2f' % (len(xvalues[k]), stats['pearsonr'][0], stats['fraction_correct'], stats['MAE'], ))
                        #pprint.pprint(results[k])
                    else:
                        colortext.warning('Could not analyze range %s - not enough datapoints.' % (str(k)))
                all_results[main_subset_names[subid]] = results
        return all_results
