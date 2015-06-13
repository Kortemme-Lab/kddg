#!/usr/bin/python2.4
# encoding: utf-8
"""
dbapi.py
High-level functions for interacting with the ddG database.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import os
import string
import re
import shutil
import glob
import traceback
import pickle
import md5
import random
import datetime
import zipfile
import pprint
import magic
import json

try:
    import matplotlib
    # A non-interactive backend to generate PNGs. matplotlib.use('PS') is used for PS files. If used, this command must be run before importing matplotlib.pyplot.
    matplotlib.use("AGG")
    import matplotlib.pyplot as plt
    import textwrap
except ImportError:
    plt=None

import score
import ddgdbapi
from tools.bio.pdb import PDB
from tools.bio.basics import residue_type_3to1_map as aa1, dssp_elision
from tools.bio.basics import Mutation
from tools.bio.alignment import ScaffoldModelChainMapper
#from Bio.PDB import *
from tools.fs.fsio import write_file, read_file
from tools.process import Popen
from tools.constants import rosetta_weights
from tools import colortext
from tools.stats.misc import get_xy_dataset_correlations

from tools.general.strutil import remove_trailing_line_whitespace
from tools.hash.md5 import get_hexdigest
from tools.fs.fsio import read_file, get_file_lines, write_file, write_temp_file

#import analysis
#from ddgfilters import PredictionResultSet, ExperimentResultSet, StructureResultSet

#todo: dbfields = ddgdbapi.FieldNames()

class MutationSet(object):
    def __init__(self):
        self.mutations = []

    def addMutation(self, chainID, residueID, wildtypeAA, mutantAA):
        self.mutations.append((chainID, residueID, wildtypeAA, mutantAA))

    def getChains(self):
        return sorted(list(set([m[0] for m in self.mutations])))


class ddG(object):
    '''This class is responsible for inserting prediction jobs to the database.'''

    def __init__(self, passwd = None, username = 'kortemmelab'):
        self.ddGDB = ddgdbapi.ddGDatabase(passwd = passwd, username = username)
        self.prediction_data_path = self.ddGDB.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionDataPath"')[0]['Value']

    def __del__(self):
        pass
        #self.ddGDB.close()
        #self.ddGDataDB.close()

    def _createResfile(self, pdb, mutations):
        '''The mutations here are in the original PDB numbering. pdb is assumed to use Rosetta numbering.
            We use the pdb mapping from PDB numbering to Rosetta numbering to generate the resfile.
        '''
        resfile = []
        for mutation in mutations:
            chain = mutation[0]
            resid = mutation[1]
            wt = mutation[2]
            mt = mutation[3]

            # Check that the expected wildtype exists in the PDB
            readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
            assert(wt == aa1[readwt])
            resid = resid.strip()
            resfile.append("%(resid)s %(chain)s PIKAA %(mt)s" % vars())
        if resfile:
            resfile = ["NATAA", "start"] + resfile
            return '\n'.join(resfile)
        else:
            raise Exception("An error occurred creating a resfile for the ddG job.")

    def _createMutfile(self, pdb, mutations):
        '''The mutations here are in the original PDB numbering. pdb is assumed to use Rosetta numbering.
            We use the pdb mapping from PDB numbering to Rosetta numbering to generate the mutfile.
        '''
        mutfile = []

        for mutation in mutations:
            chain = mutation.Chain
            resid = PDB.ResidueID2String(mutation.ResidueID)
            wt = mutation.WildTypeAA
            mt = mutation.MutantAA

            # Check that the expected wildtype exists in the PDB
            readwt = pdb.getAminoAcid(pdb.getAtomLine(chain, resid))
            assert(wt == aa1[readwt])
            resid = resid.strip()
            mutfile.append("%(wt)s %(resid)s %(mt)s" % vars())
        if mutfile:
            mutfile = ["total %d" % len(mutations), "%d" % len(mutations)] + mutfile
            return '\n'.join(mutfile)
        else:
            raise Exception("An error occurred creating a mutfile for the ddG job.")

    def getData(self, predictionID):
        job_data_path = os.path.join(self.prediction_data_path, '%d.zip' % predictionID)
        if os.path.exists(job_data_path):
            return read_file(job_data_path, binary = True)

    def getPublications(self, result_set):
        if result_set:
            structures = None
            experiments = None

            if result_set.isOfClass(ExperimentResultSet):
                experiments = result_set
            elif ExperimentResultSet in result_set.__class__.allowed_restrict_sets:
                experiments, experiment_map = result_set.getExperiments()

            if result_set.isOfClass(StructureResultSet):
                structures = result_set
            elif StructureResultSet in result_set.__class__.allowed_restrict_sets:
                structures, structure_map = result_set.getStructures()

            if structures:
                colortext.printf("\nRelated publications for structures:", "lightgreen")
                for id in sorted(structures.IDs):
                    pubs = self.ddGDB.callproc("GetPublications", parameters=(id,))
                    print(id)
                    for pub in pubs:
                        print("\t%s: %s" % (pub["Type"], pub["PublicationID"]))

            if experiments:
                colortext.printf("\nRelated publications for experiments:", "lightgreen")
                for id in sorted(experiments.IDs):
                    pubs = self.ddGDB.callproc("GetExperimentPublications", parameters=(id,))
                    print(id)
                    for pub in pubs:
                        print("\t%s: %s" % (pub["Type"], pub["SourceLocation.ID"]))

                experimentsets = [e[0] for e in self.ddGDB.execute_select("SELECT DISTINCT Source FROM Experiment WHERE ID IN (%s)" % ','.join(map(str, list(experiments.IDs))), cursorClass = ddgdbapi.StdCursor)]

                if experimentsets:
                    colortext.printf("\nRelated publications for experiment-set sources:", "lightgreen")
                    for id in sorted(experimentsets):
                        print(id)
                        pubs = self.ddGDB.execute_select("SELECT ID, Type FROM SourceLocation WHERE SourceID=%s", parameters = (id,))
                        for pub in pubs:
                            print("\t%s: %s" % (pub["Type"], pub["ID"]))
        else:
            raise Exception("Empty result set.")



    def dumpData(self, outfile, predictionID):
        write_file(outfile, self.getData(predictionID))

    def analyze(self, prediction_result_set, outpath = os.getcwd()):
        PredictionIDs = sorted(list(prediction_result_set.getFilteredIDs()))
        colortext.printf("Analyzing %d records:" % len(PredictionIDs), "lightgreen")
        #results = self.ddGDB.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % join(map(str, PredictionIDs), ","))

        #for r in results:
        #	r["ddG"] = pickle.loads(r["ddG"])
        #	predicted_score = r["ddG"]["data"]["ddG"]
        #	experimental_scores = [expscore["ddG"] for expscore in self.ddGDB.callproc("GetScores", parameters = r["ExperimentID"])]
        #	mean_experimental_score = float(sum(experimental_scores)) / float(len(experimental_scores))

        results = self.ddGDB.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % ','.join(map(str, PredictionIDs)))

        analysis.plot(analysis._R_mean_unsigned_error, analysis._createMAEFile, results, "my_plot1.pdf", average_fn = analysis._mean)
        analysis.plot(analysis._R_correlation_coefficient, analysis._createAveragedInputFile, results, "my_plot2.pdf", average_fn = analysis._mean)
        colortext.printf("Done", "lightgreen")



        #score.ddgTestScore

    def add_PDB_to_database(self, filepath = None, pdbID = None, contains_membrane_protein = None, protein = None, file_source = None, UniProtAC = None, UniProtID = None, testonly = False, force = False, techniques = None):
        assert(file_source)
        if filepath:
            if not os.path.exists(filepath):
                raise Exception("The file %s does not exist." % filepath)
            filename = os.path.split(filepath)[-1]
            rootname, extension = os.path.splitext(filename)
            if not extension.lower() == ".pdb":
                raise Exception("Aborting: The file does not have a .pdb extension.")
        if pdbID:
            rootname = pdbID

        try:
            dbp = ddgdbapi.PDBStructure(self.ddGDB, rootname, contains_membrane_protein = contains_membrane_protein, protein = protein, file_source = file_source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID, testonly = testonly, techniques = techniques)
            #Structure.getPDBContents(self.ddGDB)
            results = self.ddGDB.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters = (rootname,))
            if results:
                #ddgdbapi.getUniProtMapping(pdbID, storeInDatabase = True)
                #raise Exception("There is already a structure in the database with the ID %s." % rootname)
                if force:
                    dbp.commit(testonly = testonly)
                return None
            dbp.commit(testonly = testonly)
            return rootname
        except Exception, e:
            colortext.error(str(e))
            colortext.error(traceback.format_exc())
            raise Exception("An exception occurred committing %s to the database." % filepath)

    def add_pdb_file(self, filepath, pdb_id):
        #todo: use either this or add_PDB_to_database but not both
        raise Exception('deprecated in favor of add_PDB_to_database')
        existing_pdb = self.ddGDB.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
        if not existing_pdb:
            pdb_contents = read_file(filepath)
            p = PDB(pdb_contents)

            fasta = []
            for c, sequence in p.atom_sequences.iteritems():
                fasta.append('>%s:%s|PDBID|CHAIN|SEQUENCE' % (pdb_id.replace(':', '_'), c))
                fasta.append(str(sequence))
            fasta = '\n'.join(fasta)

            d = {
                'ID' : pdb_id,
                'FileSource' : 'Biosensor project',
                'Content' : read_file(filepath),
                'FASTA' : fasta,
                'Resolution' : None,
                'Techniques' : 'Rosetta model',
                'BFactors' : '',
                'Publication' : None
            }
            self.ddGDB.insertDictIfNew('PDBFile', d, ['ID'])

    def createDummyExperiment(self, pdbID, mutationset, chains, sourceID, ddG, ExperimentSetName = "DummySource"):
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        raise Exception("Out of date function.")
        Experiment = ddgdbapi.ExperimentSet(pdbID, ExperimentSetName)
        for m in mutationset.mutations:
            Experiment.addMutation(m[0], m[1], m[2], m[3])
        for c in chains:
            Experiment.addChain(c)
        Experiment.addExperimentalScore(sourceID, ddG, pdbID)
        Experiment.commit(self.ddGDB)

    def createDummyExperiment_ankyrin_repeat(self, pdbID, mutations, chain):
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        experiment = ddgdbapi.ExperimentDefinition(self.ddGDB, pdbID, interface = None)
        experiment.addChain(chain)
        for m in mutations:
            experiment.addMutation(m)
        experiment.commit(False)

    def add_mutant(self, pdb_ID, mutant_mutations):
        '''Use this function to add one set of mutations ON THE SAME CHAIN (i.e. corresponding to one mutant) to the database.
           todo: generalize this to allow different chains
        '''
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        chains = set([m.Chain for m in mutant_mutations])
        assert(len(chains) == 1)
        colortext.warning("Adding mutation: %s." % ', '.join(map(str, mutant_mutations)))
        self.createDummyExperiment_ankyrin_repeat(pdb_ID, mutant_mutations, chains.pop())


    def get_flattened_prediction_results(self, PredictionSet):
        #todo: add this as a stored procedure
        return self.ddGDB.execute_select('''
SELECT Prediction.ID AS PredictionID, Prediction.ExperimentID, Experiment.PDBFileID, ExperimentMutations.FlattenedMutations, Prediction.Scores, TIMEDIFF(Prediction.EndDate, Prediction.StartDate) AS TimeTaken FROM Prediction INNER JOIN
(
  SELECT ExperimentID, GROUP_CONCAT(Mutation SEPARATOR ', ') AS FlattenedMutations FROM
  (
    SELECT ExperimentID, CONCAT(Chain, ' ', WildTypeAA, ResidueID, MutantAA) As Mutation FROM ExperimentMutation
  ) AS FlattenedMutation
  GROUP BY ExperimentID
) AS ExperimentMutations
ON Prediction.ExperimentID=ExperimentMutations.ExperimentID
INNER JOIN Experiment ON Prediction.ExperimentID=Experiment.ID
WHERE Prediction.PredictionSet=%s AND Prediction.Scores IS NOT NULL
ORDER BY Prediction.ExperimentID''', parameters=(PredictionSet,))


    #### todo: these following functions should be refactored and renamed. In particular, the graphing functions should
    ####       be moved into the tools repository

    def create_abacus_graph_for_a_single_structure(self, PredictionSet, scoring_method, scoring_type, graph_title = None, PredictionIDs = None, graph_filename = None, cached_results = None, num_datapoints = 0):
        '''This function is meant to
            The num_datapoints is mainly for debugging - tuning the resolution/DPI to fit the number of datapoints.'''

        results = cached_results
        if not results:
            results = self.get_flattened_prediction_results(PredictionSet)

        pdb_ids = set()
        for r in results:
            pdb_ids.add(r['PDBFileID'])

        if len(pdb_ids) != 1:
            raise Exception('This function is only meant to be called when the PredictionSet or the set of results contains records for a single structure. The set of results contains %d structures.' % len(pdb_ids))

        sortable_results = {}
        for r in results:
            if (not PredictionIDs) or (r['PredictionID'] in PredictionIDs):
                sortable_results[(json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], r['ExperimentID'])] = r
        count = 0

        set_of_mutations = set()
        for k, r in sorted(sortable_results.iteritems()):
            #if r['FlattenedMutations'].find('A E141L') != -1 and r['FlattenedMutations'].find('A S142A') != -1 and r['FlattenedMutations'].find('A L78Y') != -1:
            #    print('%f, %s' % (k[0], r['FlattenedMutations']))
            #if r['FlattenedMutations'].find('A W103M') != -1 and r['FlattenedMutations'].find('A F70Y') != -1:
            #    if r['FlattenedMutations'].find('A E141L') == -1 and r['FlattenedMutations'].find('A S142A') == -1 and r['FlattenedMutations'].find('A L78Y') == -1:
            #        print('%f, %s' % (k[0], r['FlattenedMutations']))

            if r['FlattenedMutations'].find('A W103M') != -1 and r['FlattenedMutations'].find('A F70Y') != -1:
                if r['FlattenedMutations'].find('A E141L') == -1 and r['FlattenedMutations'].find('A S142A') == -1 and r['FlattenedMutations'].find('A L78Y') == -1:
                    #print('%f, %s' % (k[0], r['FlattenedMutations']))
                    count += 1
            #A E141L, A S142A

            mutations = [m for m in map(string.strip, r['FlattenedMutations'].split(',')) if m]
            for m in mutations:
                set_of_mutations.add((int(m.split()[1][1:-1]), m))
            #if r['FlattenedMutations'].find('A L78Y') == -1:
            #    print('%f, %s' % (k[0], r['FlattenedMutations']))
            #    #count += 1

        pruned_data = []
        for k, r in sorted(sortable_results.iteritems()):
            line = []
            print(json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], r['FlattenedMutations'])
            for m in sorted(set_of_mutations):
                if r['FlattenedMutations'].find(m[1]) != -1:
                    line.append(1)
                else:
                    line.append(0)
            pruned_data.append((json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], line))

        labels = [m[1].split()[1] for m in sorted(set_of_mutations)]
        graph_title = graph_title or r'$\Delta\Delta$G predictions for %s (%s.%s)' % (PredictionSet, scoring_method.replace(',0A', '.0$\AA$').replace('_', ' '), scoring_type)

        pruned_data = pruned_data[0:num_datapoints or len(pruned_data)]
        colortext.message('Creating graph with %d datapoints...' % len(pruned_data))

        number_of_non_zero_datapoints = 0
        for p in pruned_data:
            if 1 in p[1]:
                number_of_non_zero_datapoints += 1
                if number_of_non_zero_datapoints > 1:
                    break
        if number_of_non_zero_datapoints < 2:
            raise Exception('The dataset must contain at least two non-zero points.')

        if graph_filename:
            return self.write_abacus_graph(graph_filename, graph_title, labels, pruned_data, scoring_method, scoring_type)
        else:
            return self.create_abacus_graph(graph_title, labels, pruned_data, scoring_method, scoring_type)

    def write_abacus_graph(self, graph_filename, graph_title, labels, data, scoring_method, scoring_type):
        byte_stream = self.create_abacus_graph(graph_title, labels, data, scoring_method, scoring_type)
        write_file(graph_filename, byte_stream.getvalue(), 'wb')

    def create_abacus_graph(self, graph_title, labels, data, scoring_method, scoring_type):
        '''Even though this is technically a scatterplot, I call this an abacus graph because it is basically a bunch of beads on lines.'''

        from io import BytesIO
        if plt:
            assert(data)

            image_dpi = 300.0
            horizontal_margin = 400.0 # an estimate of the horizontal space not used by the graph
            horizontal_spacing = 100.0 # an estimate of the horizontal space between points on the same line
            vertical_margin = 100.0 # an estimate of the vertical space not used by the graph
            vertical_spacing = 50.0 # the rough amount of pixels between abacus lines
            point_size = 50 # the size of datapoints in points^2.
            y_offset = 1.0

            points_per_line = set([len(line[1]) for line in data])
            assert(len(points_per_line) == 1)
            points_per_line = points_per_line.pop()
            assert(len(labels) == points_per_line)

            number_of_lines = float(len(data))
            number_of_labels = float(len(labels))
            height_in_inches = max(600/image_dpi, (vertical_margin + (vertical_spacing * number_of_lines)) / image_dpi) # Use a minimum of 600 pixels in height. This avoids graphs with a small number of lines (<=10) not to become squashed.
            width_in_inches = max(700/image_dpi, (horizontal_margin + (horizontal_spacing * points_per_line)) / image_dpi) # Use a minimum of 600 pixels in width. This avoids graphs with a small number of labels (e.g. 1) not to become squashed.

            graph_color_scheme = matplotlib.cm.jet

            #y_offset = (1.75 * data_length) / 128
            #image_dpi = (400 * data_length) / 128
            #image_dpi = 400
            #point_sizes = {1 : 100, 64: 75, 128: 50, 192: 25, 256: 10}
            #index = round(data_length / 64.0) * 64
            #point_size = point_sizes.get(index, 10)


            fig = plt.figure(figsize=(width_in_inches, height_in_inches)) # figsize is specified in inches - w, h
            fig.set_dpi(image_dpi)

            # Create three identically-sized lists. Each triple is of an x-coordinate, a y-coordinate, and the DDG value
            # and corresponds to a 1 in the matrix i.e. we should draw a point/abacus bead at these coordinates.
            x_values = []
            y_values = []
            x_coordinate_skip = 3 # x-axis distance between two points
            y_coordinate_skip = 7 # y-axis distance between two points
            ddg_values = []
            y = 0
            for line in data:
                x = 0
                y += y_coordinate_skip
                w = line[0]
                #plt.text(30, y, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=9)
                for point in line[1]:
                    x += x_coordinate_skip
                    if point == 1:
                        x_values.append(x)
                        y_values.append(y)
                        ddg_values.append(line[0])

            # Draw the scatter plot
            plt.scatter(x_values, y_values, c=ddg_values, s=point_size, cmap=graph_color_scheme, edgecolors='none', zorder=99)

            # Define the limits of the cartesian coordinates. Add extra space on the right for the DDG values.
            extra_space = 1.3
            plt.axis((0, (points_per_line + 1 + extra_space) * x_coordinate_skip, -15, (y_coordinate_skip * number_of_lines) + 15))

            plt.tick_params(
                axis='both',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='off',      # ticks along the bottom edge are off
                left='off',      # ticks along the left edge are off
                labelleft='off', # labels along the bottom edge are off
                top='off',         # ticks along the top edge are off
                labelbottom='off') # labels along the bottom edge are off

            # Add the mutation labels at the botom of the diagram
            x = 1.9
            for i in range(len(labels)):
                l = labels[i]
                plt.text(x, -5 + ((i % 2) * -5), l, fontdict=None, withdash=True, fontsize=6)
                x += x_coordinate_skip

            added_zero_line = False
            last_y_value = 0
            y = 0
            for line in data:
                x = 0
                y += 7
                plt.plot([1, 25], [y, y], color='#999999', linestyle='-', linewidth=0.1)

                # Add a DDG value on every third line
                if y % 21 == 7:
                    plt.text(((points_per_line + 0.6) * x_coordinate_skip) , y-y_offset, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)
                if not added_zero_line:
                    if line[0] > 0:
                        plt.plot([1, 25], [0.5 + ((y + last_y_value) / 2), 0.5 + ((y + last_y_value) / 2)], color='k', linestyle='-', linewidth=1)
                        added_zero_line = True
                    else:
                        last_y_value = y

            plt.text(((points_per_line + 0.6) * x_coordinate_skip), y-y_offset, str('%.3f' % line[0]), fontdict=None, withdash=True, fontsize=6)

            # Set the colorbar font size and then add a colorbar
            #cbar.ax.tick_params(labelsize=6)
            #plt.colorbar(use_gridspec=True)

            #ax = fig.add_subplot(111)

            # Add a title. Note: doing this after the colorbar code below messes up the alignment.
            # Adjust the wrap length to the width of the graph
            wrap_length = 40 + max(0, (points_per_line - 3) * 14)
            graph_title = "\n".join(textwrap.wrap(graph_title, wrap_length))

            plt.title(graph_title, fontdict={'fontsize' : 6})

            from mpl_toolkits.axes_grid1 import make_axes_locatable
            ax = fig.add_subplot(111)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05) # pad specifies the padding between the right of the graph and the left of the colorbar, size seems to specify the width of the colobar relative to the graph

            CS3 = plt.contourf([[0,0],[0,0]], ddg_values, cmap=graph_color_scheme)
            #plt.colorbar(CS3)
            #cbar = fig.colorbar(CS3, format='%.2f')
            cbar = fig.colorbar(CS3, format='%.2f', cax=cax)
            cbar.set_label('$\Delta\Delta$G',size=6)
            cbar.ax.tick_params(labelsize=5)

             # Use the tight_layout command to tighten up the spaces. The pad, w_pad, and h_pad parameters are specified in fraction of fontsize.
            plt.tight_layout(pad=0.5)

            #quadmesh = ax.pcolormesh(theta,phi,data)
            #cb = fig.colorbar(quadmesh,ax=ax, shrink=.5, pad=.2, aspect=10)
            #cax = ax.imshow(ddg_values, interpolation='nearest', cmap=matplotlib.cm.coolwarm)
            #cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
            #surf = ax.contourf(X,Y,Z, 8, cmap=cm.jet)
            #cbar = fig.colorbar(surf, use_gridspec=True, shrink=0.5, aspect=20, fraction=.12,pad=.02)
            #cbar.set_label('Activation',size=18)

            byte_stream = BytesIO()
            plt.savefig(byte_stream, dpi=image_dpi, format="png")

            plt.close(fig)

            return byte_stream
        else:
            return None

    def extract_data(output_dir, PredictionID):
        assert(os.path.exists(output_dir))
        archive = self.getData(PredictionID)
        write_file(os.path.join(output_dir, '%d.zip' % PredictionID), archive, 'wb')
        p = Popen(output_dir, ['unzip', '%d.zip' % PredictionID])
        os.remove(os.path.join(output_dir, '%d.zip' % PredictionID))
        if p.errorcode != 0:
            raise colortext.Exception(p.stderr)
        else:
            colortext.warning(p.stdout)

    def test_results(output_dir, PredictionSet):
        PredictionIDs = []
        results = get_flattened_prediction_results(PredictionSet)
        mutation_lists = {}
        for r in results:
            PredictionIDs.append(r['PredictionID'])
            mutation_lists[r['PredictionID']] = r['FlattenedMutations']
        RandomPredictionIDs = [PredictionIDs[random.randint(0, len(PredictionIDs) - 1)] for k in range(10)]
        RandomPredictionIDs = [54090L, 53875L, 54085L, 54079L, 54008L, 53853L, 53952L, 54056L, 53935L, 53893L]

        # Retrieve and unzip results
        if not(os.path.exists(output_dir)):
            os.mkdir(output_dir)
        for PredictionID in PredictionIDs:#RandomPredictionIDs:
            if not(os.path.exists(os.path.join(output_dir, str(PredictionID)))):
                colortext.message('Retrieving archive for Prediction %d.' % PredictionID)
                self.extract_data(output_dir, PredictionID)

        # Get the sequences of the wildtype and mutant structures
        count = 0
        for PredictionID in PredictionIDs:#RandomPredictionIDs:
            wildtype_sequences = set()
            mutation_sequences = set()
            working_dir = os.path.join(os.path.join(output_dir, str(PredictionID)))
            for f in glob.glob(os.path.join(working_dir, '*.pdb')):
                if os.path.split(f)[1].startswith('mut_'):
                    p = PDB.from_filepath(f)
                    assert(len(p.atom_sequences) == 1)
                    sequence = str(p.atom_sequences.values()[0])
                    mutation_sequences.add(sequence)
                elif os.path.split(f)[1].startswith('repacked_wt_'):
                    p = PDB.from_filepath(f)
                    assert(len(p.atom_sequences) == 1)
                    sequence = str(p.atom_sequences.values()[0])
                    wildtype_sequences.add(sequence)

            assert(len(wildtype_sequences) == 1)
            assert(len(mutation_sequences) == 1)
            wildtype_sequence = wildtype_sequences.pop()
            mutation_sequence = mutation_sequences.pop()

            colortext.message('Prediction %d. Mutations: %s' % (PredictionID, mutation_lists[PredictionID]))
            assert(len(wildtype_sequence) == len(mutation_sequence))
            s = ''
            t = ''
            for x in range(len(wildtype_sequence)):
                if wildtype_sequence[x] != mutation_sequence[x]:
                    s += colortext.make(wildtype_sequence[x], color="green")
                    t += colortext.make(mutation_sequence[x], color="yellow")
                else:
                    s += wildtype_sequence[x]
                    t += mutation_sequence[x]
            print(s)
            print(t)

    def create_pymol_session(download_dir, PredictionID, task_number, keep_files = True):
        '''Create a PyMOL session for a pair of structures.'''

        # Retrieve and unzip results
        if not(os.path.exists(download_dir)):
            os.mkdir(download_dir)
        working_dir = os.path.join(os.path.join(download_dir, str(PredictionID)))
        if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
            self.extract_data(download_dir, PredictionID)
        if not(os.path.exists(working_dir)) or not(os.path.exists(os.path.join(working_dir, 'repacked_wt_round_%d.pdb' % task_number))):
            raise Exception('Could not extract the models for task #%d of Prediction #%d.' % (task_number, PredictionID))

        # Retrieve the two structures corresponding to the task_number
        files = sorted(glob.glob(os.path.join(working_dir, '*_round_%d.pdb' % task_number)), reverse = True)
        assert(os.path.split(files[0])[1].startswith('repacked_wt_'))
        assert(os.path.split(files[1])[1].startswith('mut_'))

        # Creator the alignment object and write the PSE file
        chain_mapper = ScaffoldModelChainMapper.from_file_paths(files[0], files[1])

        # Remove the downloaded files
        if not keep_files:
            shutil.rmtree(download_dir)
        return chain_mapper.generate_pymol_session()

    def create_pymol_session_in_memory(self, PredictionID, task_number):
        '''Create a PyMOL session for a pair of structures.
        '''
        # todo: This should replace create_pymol_session since creating a PyMOL session should be a separate task from extracting an archive. write_pymol_session should be changed to use this function

        from io import BytesIO

        # Retrieve and unzip results
        archive = self.getData(PredictionID)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            # Get the name of the files from the zip
            wildtype_filename = os.path.join(str(PredictionID), 'repacked_wt_round_%d.pdb' % task_number)
            mutant_filename = None
            for filepath in sorted(zipped_content.namelist()):
                filename = os.path.split(filepath)[1]
                if filename.startswith('mut_') and filename.endswith('_round_%d.pdb' % task_number):
                    mutant_filename = os.path.join(str(PredictionID), filename)
                    break

            PyMOL_session = None
            file_list = zipped_content.namelist()

            # If both files exist in the zip, extract their contents in memory and create a PyMOL session pair (PSE, script)
            if (mutant_filename in file_list) and (wildtype_filename in file_list):
                wildtype_pdb = zipped_content.open(wildtype_filename, 'r').read()
                mutant_pdb = zipped_content.open(mutant_filename, 'U').read()
                chain_mapper = ScaffoldModelChainMapper.from_file_contents(wildtype_pdb, mutant_pdb)
                PyMOL_session = chain_mapper.generate_pymol_session(pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol')

            zipped_content.close()
            return PyMOL_session

        except Exception, e:
            zipped_content.close()
            raise Exception(str(e))


    def extract_job_stdout_from_archive(self, PredictionID):
        '''Extract all stdout files from a Prediction. This function returns a dict mapping the type of output file to the contents of the file.'''

        from io import BytesIO

        # Retrieve and unzip results
        archive = self.getData(PredictionID)
        zipped_content = zipfile.ZipFile(BytesIO(archive), 'r', zipfile.ZIP_DEFLATED)

        try:
            stdout_file_names = {}
            stdout_file_list = [l for l in sorted(zipped_content.namelist()) if (l.find('cmd.o') != -1)]
            for f in stdout_file_list:
                tokens = os.path.split(f)
                assert(tokens[0].isdigit())
                title = tokens[1].split('_')[0]
                assert(stdout_file_names.get(title) == None)
                stdout_file_names[title] = f

            stdout_files = {}
            for stdout_type, filename in stdout_file_names.iteritems():
                stdout_files[stdout_type] = zipped_content.open(filename, 'r').read()

            zipped_content.close()
            return stdout_files

        except Exception, e:
            zipped_content.close()
            raise Exception(str(e))


    def get_ddg_monomer_scores_per_structure(self, PredictionID):
        '''Returns a dict mapping the DDG scores from a ddg_monomer run to a list of structure numbers.'''

        # Get the ddg_monomer output for the prediction
        stdout = self.extract_job_stdout_from_archive(PredictionID)['ddG']

        # Parse the stdout into two mappings (one for wildtype structures, one for mutant structures) mapping
        # structure IDs to a dict containing the score components
        wildtype_scores = {}
        mutant_scores = {}
        s1 = 'score before mutation: residue'
        s1_len = len(s1)
        s2 = 'score after mutation: residue'
        s2_len = len(s2)
        for line in stdout.split('\n'):
            idx = line.find(s1)
            if idx != -1:
                idx += s1_len
                mtchs = re.match('.*?(\d+) %s' % s1, line)
                structure_id = int(mtchs.group(1))
                assert(structure_id not in wildtype_scores)
                tokens = line[idx:].split()
                d = {'total' : float(tokens[0])}
                for x in range(1, len(tokens), 2):
                    component_name = tokens[x].replace(':', '')
                    assert(rosetta_weights.get(component_name))
                    component_value = float(tokens[x + 1])
                    d[component_name] = component_value
                wildtype_scores[structure_id] = d
            else:
                idx = line.find(s2)
                if idx != -1:
                    idx += s2_len
                    mtchs = re.match('.*?(\d+) %s' % s2, line)
                    structure_id = int(mtchs.group(1))
                    assert(structure_id not in mutant_scores)
                    tokens = line[idx:].split()
                    d = {'total' : float(tokens[1])}
                    for x in range(2, len(tokens), 2):
                        component_name = tokens[x].replace(':', '')
                        assert(rosetta_weights.get(component_name))
                        component_value = float(tokens[x + 1])
                        d[component_name] = component_value
                    mutant_scores[structure_id] = d

        # Sanity checks
        num_structures = max(wildtype_scores.keys())
        expected_keys = set(range(1, num_structures + 1))
        assert(expected_keys == set(wildtype_scores.keys()))
        assert(expected_keys == set(mutant_scores.keys()))

        # Create a list of lists - MutantScoreOrder - of structure IDs e.g. [[5,1,34], [23], [12,3], ...] which is ordered
        # by increasing energy so that each sublist contains structure IDs of equal energy and if structures have the same
        # energy then their IDs are in the same sublist
        d = {}
        for structure_id, scores in sorted(mutant_scores.iteritems()):
            d[scores['total']] = d.get(scores['total'], [])
            d[scores['total']].append(structure_id)
        MutantScoreOrder = []
        for score, structure_ids in sorted(d.iteritems()):
            MutantScoreOrder.append(structure_ids)

        # Sanity check - make sure that MutantScoreOrder is really ordered such that each set of structure IDs contains
        # structures of the same energy and of a lower energy than the following set of structure IDs in the list
        for x in range(len(MutantScoreOrder) - 1):
            s1 = set([mutant_scores[n]['total'] for n in MutantScoreOrder[x]])
            assert(len(s1) == 1)
            if x + 1 < len(MutantScoreOrder):
                s2 = set([mutant_scores[n]['total'] for n in MutantScoreOrder[x + 1]])
                assert(len(s2) == 1)
                assert(s1.pop() < s2.pop())

        return dict(
            WildType = wildtype_scores,
            Mutant = mutant_scores,
            MutantScoreOrder = MutantScoreOrder,
        )


    def determine_best_pair(self, PredictionID, ScoreMethodID = 1):
        # Iterates over the (wildtype, mutant) pairs in the PredictionStructureScore table and returns the structure ID
        # for the pair with the lowest energy mutant
        # Note: There are multiple ways to select the best pair. For example, if multiple mutants have the same minimal total
        # score, we could have multiple wildtype structures to choose from. In this case, we choose a pair where the wildtype
        # structure has the minimal total score.

        lowest_mutant_score = self.ddGDB.execute_select('SELECT total FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" ORDER BY total LIMIT 1', parameters=(PredictionID, ScoreMethodID))
        if lowest_mutant_score:
            lowest_mutant_score = lowest_mutant_score[0]['total']
            mutant_structure_ids = [r['StructureID'] for r in self.ddGDB.execute_select('SELECT StructureID FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="Mutant" AND total=%s', parameters=(PredictionID, ScoreMethodID, lowest_mutant_score))]
            if len(mutant_structure_ids) > 1:
                return self.ddGDB.execute_select(('SELECT StructureID FROM PredictionStructureScore WHERE PredictionID=%s AND ScoreMethodID=%s AND ScoreType="WildType" AND StructureID IN (' + ','.join(map(str, mutant_structure_ids)) + ') ORDER BY total LIMIT 1'), parameters=(PredictionID, ScoreMethodID ))[0]['StructureID']
            else:
                return mutant_structure_ids[0]
        return None


    def write_pymol_session(download_dir, PredictionID, task_number, keep_files = True):
        PSE_file = create_pymol_session(download_dir, PredictionID, task_number, keep_files = keep_files)
        write_file(output_filepath, PSE_file[0], 'wb')


    def get_amino_acids_for_analysis(self):
        amino_acids = {}
        polarity_map = {'polar' : 'P', 'charged' : 'C', 'hydrophobic' : 'H'}
        aromaticity_map = {'aliphatic' : 'L', 'aromatic' : 'R', 'neither' : '-'}
        results = self.ddGDB.execute_select('SELECT * FROM AminoAcid')
        for r in results:
            if r['Code'] != 'X':
                amino_acids[r['Code']] = dict(
                    LongCode = r['LongCode'],
                    Name = r['Name'],
                    Polarity = polarity_map.get(r['Polarity'], 'H'),
                    Aromaticity = aromaticity_map[r['Aromaticity']],
                    Size = r['Size'],
                    van_der_Waals_volume = r['Volume']
                )
        amino_acids['Y']['Polarity'] = 'H'
        return amino_acids


    def get_pdb_details_for_analysis(self, pdb_ids, cached_pdb_details = None):
        pdb_chain_lengths = {}
        pdb_os = {}
        pdbs = {}
        cached_pdb_ids = []
        if cached_pdb_details:
            cached_pdb_ids = cached_pdb_details.keys()
        for pdb_id in pdb_ids:
            if pdb_id in cached_pdb_ids:
                pdbs[pdb_id] = cached_pdb_details[pdb_id]
            else:
                record = self.ddGDB.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))[0]
                p = PDB(record['Content'])
                pdb_chain_lengths = {}
                for chain_id, s in p.atom_sequences.iteritems():
                    pdb_chain_lengths[chain_id] = len(s)
                pdbs[pdb_id] = dict(
                    chains = pdb_chain_lengths,
                    TM = record['Transmembrane'],
                    Technique = record['Techniques'],
                    XRay = record['Techniques'].find('X-RAY') != -1,
                    Resolution = record['Resolution'],
                )

        return pdbs


    def get_prediction_experiment_chains(self, predictionset):
        return self.ddGDB.execute_select('''
            SELECT Prediction.ID, Experiment.PDBFileID, Chain
            FROM Prediction
            INNER JOIN Experiment ON Experiment.ID=Prediction.ExperimentID
            INNER JOIN ExperimentChain ON ExperimentChain.ExperimentID=Prediction.ExperimentID
            WHERE PredictionSet=%s''', parameters=(predictionset,))







    #####
    ## Database API: start of curated functionality. Eventually all functions should be moved into here
    #####

    ##### Deprecated functions



    def create_PredictionSet(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40, allow_existing_prediction_set = False, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False): raise Exception('This function has been deprecated. Use add_prediction_set instead.')
    def charge_PredictionSet_by_number_of_residues(self, PredictionSet): raise Exception('This function has been deprecated. Use _charge_prediction_set_by_residue_count instead.')
    def createPredictionsFromUserDataSet(self, userdatasetTextID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = False, Description = {}, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False): raise Exception('This function has been deprecated. Use add_prediction_set_jobs instead.')
    def add_predictions_by_pdb_id(self, pdb_ID, PredictionSet, ProtocolID, status = 'active', priority = 5, KeepHETATMLines = False, strip_other_chains = True): raise Exception('This function has been deprecated. Use add_jobs_by_pdb_id instead.')
    def addPrediction(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False, strip_other_chains = True): raise Exception('This function has been deprecated. Use add_job instead.')



    ##### Public API: PredictionSet functions



    def add_prediction_set(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40, allow_existing_prediction_set = False, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False):
        raise Exception('This function needs to be rewritten.')

        if halted:
            Status = 'halted'
        else:
            Status = 'active'

        if allow_existing_prediction_set == False and len(self.ddGdb.execute_select('SELECT * FROM PredictionSet WHERE ID=%s', parameters=(PredictionSetID,))) > 0:
            raise Exception('The PredictionSet %s already exists.' % PredictionSetID)
        d = dict(
            ID                  = PredictionSetID,
            Status              = Status,
            Priority            = Priority,
            ProteinStability    = contains_protein_stability_predictions,
            BindingAffinity     = contains_binding_affinity_predictions,
            BatchSize           = BatchSize,
        )
        self.ddGDB.insertDictIfNew("PredictionSet", d, ['ID'])



    ##### Public API: File-related functions



    def get_file_id(self, content, hexdigest = None):
        '''Searches the database to see whether the FileContent already exists. The search uses the digest and filesize as
           heuristics to speed up the search. If a file has the same hex digest and file size then we do a straight comparison
           of the contents.
           If the FileContent exists, the value of the ID field is returned else None is returned.
           '''
        existing_filecontent_id = None
        hexdigest = hexdigest or get_hexdigest(content)
        filesize = len(content)
        for r in self.ddGDB.execute_select('SELECT * FROM FileContent WHERE MD5HexDigest=%s AND Filesize=%s', parameters=(hexdigest, filesize)):
            if r['Content'] == content:
                assert(existing_filecontent_id == None) # content uniqueness check
                existing_filecontent_id = r['ID']
        return existing_filecontent_id


    def add_file_content(self, content, rm_trailing_line_whitespace = False, forced_mime_type = None):
        '''Takes file content (and an option to remove trailing whitespace from lines e.g. to normalize PDB files), adds
           a new record if necessary, and returns the associated FileContent.ID value.'''

        if rm_trailing_line_whitespace:
            content = remove_trailing_line_whitespace(content)

        # Check to see whether the file has been uploaded before
        hexdigest = get_hexdigest(content)
        existing_filecontent_id = self.get_file_id(content, hexdigest = hexdigest)

        # Create the FileContent record if the file is a new file
        if existing_filecontent_id == None:
            mime_type = None
            if forced_mime_type:
                mime_type = forced_mime_type
            else:
                temporary_file = write_temp_file('/tmp', content, ftype = 'wb')
                m=magic.open(magic.MAGIC_MIME_TYPE) # see mime.__dict__ for more values e.g. MAGIC_MIME, MAGIC_MIME_ENCODING, MAGIC_NONE
                m.load()
                mime_type = m.file(temporary_file)
                os.remove(temporary_file)

            d = dict(
                Content = content,
                MIMEType = mime_type,
                Filesize = len(content),
                MD5HexDigest = hexdigest
            )
            self.ddGDB.insertDictIfNew('FileContent', d, ['Content'])
            existing_filecontent_id = self.get_file_id(content, hexdigest = hexdigest)
            assert(existing_filecontent_id != None)
        return existing_filecontent_id


    def add_pdb_file_content(self, pdb_content):
        return self.add_file_content(pdb_content, rm_trailing_line_whitespace = True, forced_mime_type = 'chemical/x-pdb')



    ##### Private API: File-related functions



    def _add_residue_map_to_prediction(self, prediction_table, prediction_id, residue_mapping):
        assert(type(residue_mapping) == type(self.__dict__))

        # Add the file contents to the database
        json_content = json.dumps(residue_mapping, sort_keys=True) # sorting helps to quotient the file content space over identical data
        file_content_id = self.add_file_content(json_content, forced_mime_type="application/json")

        # Link the file contents to the prediction
        d = dict(
            FileContentID = file_content_id,
            Filename = 'residue_mapping.json' % prediction_id,
            Filetype = 'RosettaPDBMapping',
            FileRole = 'Rosetta<->PDB residue mapping',
            Stage = 'Input',
        )
        if prediction_table == 'Prediction':
            d['PredictionID'] = prediction_id,
            self.ddGDB.insertDictIfNew('PredictionFile', d, ['PredictionID', 'FileRole', 'Stage'])
        elif prediction_table == 'PredictionPPI':
            d['PredictionPPIID'] = prediction_id,
            self.ddGDB.insertDictIfNew('PredictionPPIFile', d, ['PredictionPPIID', 'FileRole', 'Stage'])
        else:
            raise('Invalid table "%s" passed.' % prediction_table)







