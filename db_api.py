#!/usr/bin/python2.4
# encoding: utf-8
"""
db_api.py
High-level functions for interacting with the ddG database.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import sys
import os
import string
import re
import shutil
import glob
import traceback
import pickle
import random
import datetime
import zipfile
import StringIO
import gzip
import pprint
try:
    import magic
except ImportError:
    pass
import json

from io import BytesIO
from api_layers import *
try:
    import pandas
except ImportError:
    pass

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
from db_schema import test_schema_against_database_instance

from klab.bio.pdb import PDB
from klab.bio.basics import residue_type_3to1_map as aa1, dssp_elision
from klab.bio.basics import Mutation
from klab.bio.pdbtm import PDBTM

#from Bio.PDB import *
from klab.fs.fsio import write_file, read_file, open_temp_file
from klab.process import Popen
from klab.constants import rosetta_weights
from klab import colortext
from klab.stats.misc import get_xy_dataset_statistics

from klab.general.strutil import remove_trailing_line_whitespace
from klab.hash.md5 import get_hexdigest
from klab.fs.fsio import read_file, get_file_lines, write_file, write_temp_file


class FatalException(Exception): pass
class PartialDataException(Exception): pass
class SanityCheckException(Exception): pass

class MutationSet(object):
    '''This class is a leftover from Lin's work and should probably be folded into an API function along with the functions that call this.
       I wrote a function elsewhere to create a list of all mutations for a PDB file - maybe check out the ubiquitin complex project - which
       can probably be used to replace most of this code.'''
    def __init__(self):
        self.mutations = []

    def addMutation(self, chainID, residueID, wildtypeAA, mutantAA):
        self.mutations.append((chainID, residueID, wildtypeAA, mutantAA))

    def getChains(self):
        return sorted(list(set([m[0] for m in self.mutations])))



class ddG(object):
    '''This is the base database API class. It should not be used directly to create interface objects. Instead, use one
       of the derived classes e.g. MonomericStabilityDDGInterface or the clean user API which hides internal functionality.
       The clean API is instantiated as in the example below:

            from ddglib.monomer_api import get_interface as get_protein_stability_interface
            stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))
            stability_api.help()

       Objects of this class and derived subclasses has three main members:

          self.DDG_db - a database interface used to interact directly with the database via MySQL commands
          self.DDG_db_utf - the same interface but with UTF support. This should be used when dealing with UTF fields e.g. publication data
          self.prediction_data_path - this is the location on the file server where output form jobs of the derived class type (e.g. binding affinity jobs) should be stored.
    '''

    GET_JOB_FN_CALL_COUNTER_MAX = 10

    def __init__(self, passwd = None, username = 'kortemmelab', hostname = 'kortemmelab.ucsf.edu', rosetta_scripts_path = None, rosetta_database_path = None):
        if passwd:
            passwd = passwd.strip()
        self.DDG_db = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname)
        self.DDG_db_utf = ddgdbapi.ddGDatabase(passwd = passwd, username = username, hostname = hostname, use_utf = True)
        self.prediction_data_path = None
        self.rosetta_scripts_path = rosetta_scripts_path
        self.rosetta_database_path = rosetta_database_path

        # Before continuing, make sure that the SQLAlchemy definitions match the table definitions
        test_schema_against_database_instance(self.DDG_db)

        # This counter is used to check the number of times get_job is called and raise an exception if this exceeds a certain amount
        # If the API is misused then get_job may be called infinitely on one job - this is meant to protect against that
        self._get_job_fn_call_counter = {}
        self._get_job_fn_call_counter_max = ddG.GET_JOB_FN_CALL_COUNTER_MAX

        # Caching dictionaries
        self.cached_score_method_details = None


    def __del__(self):
        pass         #self.DDG_db.close()         #self.ddGDataDB.close()



    #########################################################################################
    ## Public API
    #########################################################################################



    #########################################################################################
    ## Broken API layer
    ##
    ## This section contains useful functions which need to be updated to work with the new
    ## schema or code
    #########################################################################################


    #== Alien functions ====================================================================
    #==
    #== These functions do not belong here.


    @alien
    def write_abacus_graph(self, graph_filename, graph_title, labels, data):
        '''NOTE: This function should be generalized and moved into the klab repository.
           This is a simple function wrapper around create_abacus_graph which writes the graph to file.'''
        byte_stream = self.create_abacus_graph(graph_title, labels, data)
        write_file(graph_filename, byte_stream.getvalue(), 'wb')


    @alien
    def create_abacus_graph(self, graph_title, labels, data):
        '''NOTE: This function should be generalized and moved into the klab repository.
           This function creates an 'abacus graph' from a set of data. Even though this is technically a scatterplot,
           I call this an abacus graph because it is looks like rows of beads on lines.
           The function takes a graph title, a set of labels (one per row of data), and an array of data where each row
           should have the same number of columns.
           A byte stream for the graph (currently PNG format but we could parameterize this) is returned. This may be
           written directly to a binary file or streamed for online display.
        '''

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


    @alien
    def get_flattened_prediction_results(self, PredictionSet):
        '''This is defined here as an API function but should be defined as a stored procedure.'''

        # @todo: this is the monomeric stability implementation - move this into that API
        #Ubiquitin scan: 1UBQ p16
        #Prediction.Scores no longer exists
        kellogg_score_id = self.get_score_method_id('global', method_type = 'protocol 16', method_authors = 'kellogg', fuzzy = True)
        noah_score_id = self.get_score_method_id('local', method_type = 'position', method_parameters = '8Å radius', method_authors = 'Noah Ollikainen', fuzzy = False)

        score_ids = {}
        score_ids['kellogg'] = kellogg_score_id
        score_ids['noah8A'] = noah_score_id

        records = self.DDG_db.execute_select('''
SELECT Prediction.ID AS PredictionID, Prediction.ExperimentID, Experiment.PDBFileID, ExperimentMutations.FlattenedMutations, TIMEDIFF(Prediction.EndDate, Prediction.StartDate) AS TimeTaken, PredictionStructureScore.ScoreMethodID, PredictionStructureScore.DDG
FROM Prediction INNER JOIN
(
  SELECT ExperimentID, GROUP_CONCAT(Mutation SEPARATOR ', ') AS FlattenedMutations FROM
  (
    SELECT ExperimentID, CONCAT(Chain, ' ', WildTypeAA, ResidueID, MutantAA) As Mutation FROM ExperimentMutation
  ) AS FlattenedMutation
  GROUP BY ExperimentID
) AS ExperimentMutations
ON Prediction.ExperimentID=ExperimentMutations.ExperimentID
INNER JOIN Experiment ON Prediction.ExperimentID=Experiment.ID
INNER JOIN PredictionStructureScore ON Prediction.ID=PredictionStructureScore.PredictionID
WHERE Prediction.PredictionSet=%s AND Prediction.Status="done" AND ScoreType="DDG" AND StructureID=-1 AND (ScoreMethodID=%s OR ScoreMethodID=%s)
ORDER BY ScoreMethodID''', parameters=(PredictionSet, kellogg_score_id, noah_score_id))

        return records, score_ids

    #== Broken functions ====================================================================


    #@informational or part of a filter API
    @brokenfn
    def get_publications_for_result_set(self, result_set):
        '''This should be fixed once the filter API has been rewritten to work with the new DB schema. It returns a list of publications associated with the filter result set.'''
        raise Exception('')
        from ddgfilters import ExperimentResultSet, StructureResultSet
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
                    pubs = self.DDG_db.callproc("GetPublications", parameters=(id,))
                    print(id)
                    for pub in pubs:
                        print("\t%s: %s" % (pub["Type"], pub["PublicationID"]))

            if experiments:
                colortext.printf("\nRelated publications for experiments:", "lightgreen")
                for id in sorted(experiments.IDs):
                    pubs = self.DDG_db.callproc("GetExperimentPublications", parameters=(id,))
                    print(id)
                    for pub in pubs:
                        print("\t%s: %s" % (pub["Type"], pub["SourceLocation.ID"]))

                experimentsets = [e[0] for e in self.DDG_db.execute_select("SELECT DISTINCT Source FROM Experiment WHERE ID IN (%s)" % ','.join(map(str, list(experiments.IDs))), cursorClass = ddgdbapi.StdCursor)]

                if experimentsets:
                    colortext.printf("\nRelated publications for experiment-set sources:", "lightgreen")
                    for id in sorted(experimentsets):
                        print(id)
                        pubs = self.DDG_db.execute_select("SELECT ID, Type FROM SourceLocation WHERE SourceID=%s", parameters = (id,))
                        for pub in pubs:
                            print("\t%s: %s" % (pub["Type"], pub["ID"]))
        else:
            raise Exception("Empty result set.")


    #@analysis_api
    @brokenfn
    def analyze(self, prediction_result_set, outpath = os.getcwd()):
        '''This function needs to be rewritten and renamed. It calls the analysis module (which creates LaTeX reports) to generate correlation and MAE graphs.'''
        raise Exception('The import of analysis was commented out - presumably some change in DB structure or API broke the import. This code probably needs to be fixed.')
        import analysis

        PredictionIDs = sorted(list(prediction_result_set.getFilteredIDs()))
        colortext.printf("Analyzing %d records:" % len(PredictionIDs), "lightgreen")
        #results = self.DDG_db.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % join(map(str, PredictionIDs), ","))

        #for r in results:
        #	r["ddG"] = pickle.loads(r["ddG"])
        #	predicted_score = r["ddG"]["data"]["ddG"]
        #	experimental_scores = [expscore["ddG"] for expscore in self.DDG_db.callproc("GetScores", parameters = r["ExperimentID"])]
        #	mean_experimental_score = float(sum(experimental_scores)) / float(len(experimental_scores))

        results = self.DDG_db.execute_select("SELECT ID, ExperimentID, ddG FROM Prediction WHERE ID IN (%s)" % ','.join(map(str, PredictionIDs)))

        analysis.plot(analysis._R_mean_unsigned_error, analysis._createMAEFile, results, "my_plot1.pdf", average_fn = analysis._mean)
        analysis.plot(analysis._R_correlation_coefficient, analysis._createAveragedInputFile, results, "my_plot2.pdf", average_fn = analysis._mean)
        colortext.printf("Done", "lightgreen")
        #score.ddgTestScore


    #== Deprecated functions =================================================================


    @deprecated
    def create_PredictionSet(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40, allow_existing_prediction_set = False, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False): raise Exception('This function has been deprecated. Use add_prediction_set instead.')

    @deprecated
    def charge_PredictionSet_by_number_of_residues(self, PredictionSet): raise Exception('This function has been deprecated. Use _charge_prediction_set_by_residue_count instead.')

    @deprecated
    def createPredictionsFromUserDataSet(self, userdatasetTextID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = False, Description = {}, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False): raise Exception('This function has been deprecated. Use add_prediction_run instead.')

    @deprecated
    def add_predictions_by_pdb_id(self, pdb_ID, PredictionSet, ProtocolID, status = 'active', priority = 5, KeepHETATMLines = False, strip_other_chains = True): raise Exception('This function has been deprecated. Use add_jobs_by_pdb_id instead.')

    @deprecated
    def addPrediction(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False, strip_other_chains = True): raise Exception('This function has been deprecated. Use add_job instead.')

    @deprecated
    def add_pdb_file(self, filepath, pdb_id): raise Exception('This function has been deprecated. Use the import_api.add_pdb_* functions instead.')

    @deprecated
    def getPublications(self, result_set): raise Exception('This function has been deprecated. Use get_publications_for_result_set instead.')

    @deprecated
    def getData(self, predictionID): raise Exception('This function has been deprecated. Use get_job_data instead.')

    @deprecated
    def dumpData(self, outfile, predictionID): raise Exception('This function has been deprecated. Use write_job_data_to_disk instead (note the change in argument order).')

    @deprecated
    def get_amino_acids_for_analysis(self): raise Exception('This function has been deprecated. Use get_amino_acid_details instead.')

    @deprecated
    def get_pdb_details_for_analysis(self, pdb_ids, cached_pdb_details = None): raise Exception('This function has been deprecated. Use get_pdb_details instead.')

    @deprecated
    def add_pdb_file_content(self, pdb_content): raise Exception('This function may never have been used and should be removed.') # return self._add_file_content(pdb_content, rm_trailing_line_whitespace = True, forced_mime_type = 'chemical/x-pdb')

    @deprecated
    def create_pymol_session(self, download_dir, prediction_id, task_number, keep_files = True): raise Exception('This function has been deprecated. Use create_pymol_session_in_memory and write_pymol_session instead.''')

    @deprecated
    def createDummyExperiment(self, pdbID, mutationset, chains, sourceID, ddG, ExperimentSetName = "DummySource"):
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        raise Exception("Out of date function.")
        Experiment = ddgdbapi.ExperimentSet(pdbID, ExperimentSetName)
        for m in mutationset.mutations:
            Experiment.addMutation(m[0], m[1], m[2], m[3])
        for c in chains:
            Experiment.addChain(c)
        Experiment.addExperimentalScore(sourceID, ddG, pdbID)
        Experiment.commit(self.DDG_db)

    @deprecated
    def createDummyExperiment_ankyrin_repeat(self, pdbID, mutations, chain):
        raise Exception("Out of date function.")
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        experiment = ddgdbapi.ExperimentDefinition(self.DDG_db, pdbID, interface = None)
        experiment.addChain(chain)
        for m in mutations:
            experiment.addMutation(m)
        experiment.commit(False)

    #@analysis_api
    @deprecated
    def test_results(self, output_dir, PredictionSet):
        PredictionIDs = []
        results = self.get_flattened_prediction_results(PredictionSet)
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
                self.write_job_data_to_disk(PredictionID, output_dir)


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

    @deprecated
    def add_mutant(self, pdb_ID, mutant_mutations):
        '''Use this function to add one set of mutations ON THE SAME CHAIN (i.e. corresponding to one mutant) to the database.
           todo: generalize this to allow different chains
        '''
        raise Exception("Out of date function.")
        #todo: elide createDummyExperiment, createDummyExperiment_ankyrin_repeat, and add_mutant
        chains = set([m.Chain for m in mutant_mutations])
        assert(len(chains) == 1)
        colortext.warning("Adding mutation: %s." % ', '.join(map(str, mutant_mutations)))
        self.createDummyExperiment_ankyrin_repeat(pdb_ID, mutant_mutations, chains.pop())


    ###########################################################################################
    ## Information layer
    ##
    ## This layer is for functions which extract data from the database.
    ###########################################################################################


    #== Information API =======================================================================


    @informational_misc
    def get_amino_acid_details(self):
        '''This function returns a dictionary of canonical amino acid details e.g. polarity, aromaticity, size etc.'''
        amino_acids = {}
        polarity_map = {'polar' : 'P', 'charged' : 'C', 'hydrophobic' : 'H'}
        aromaticity_map = {'aliphatic' : 'L', 'aromatic' : 'R', 'neither' : '-'}
        results = self.DDG_db.execute_select('SELECT * FROM AminoAcid')
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
        amino_acids['Y']['Polarity'] = 'H' # tyrosine is a special case
        return amino_acids


    @informational_misc
    def get_publication(self, ID):
        '''Returns the information (title, publication, authors etc.) for a publication.'''
        r = self.DDG_db_utf.execute_select('SELECT * FROM Publication WHERE ID=%s', parameters=(ID,))
        if not r:
            raise Exception('No publication exists with ID %s.' % str(ID))
        r = r[0]
        pubmed_id = self.DDG_db_utf.execute_select('SELECT * FROM PublicationIdentifier WHERE SourceID=%s AND Type="PMID"', parameters=(r['ID'],))
        if pubmed_id:
            pubmed_id = pubmed_id[0]['ID']
        authors = self.DDG_db_utf.execute_select('SELECT * FROM PublicationAuthor WHERE PublicationID=%s ORDER BY AuthorOrder', parameters=(r['ID'],))
        authorlist = []
        for a in authors:
            authorlist.append(dict(FirstName = a['FirstName'], MiddleNames = a['MiddleNames'], Surname = a['Surname']))
        pub_details = dict(
            Title = r['Title'],
            Publication = r['Publication'],
            Volume = r['Volume'],
            StartPage = r['StartPage'],
            EndPage = r['EndPage'],
            PublicationYear = r['PublicationYear'],
            PublicationDate = r['PublicationDate'],
            DOI = r['DOI'],
            URL = r['URL'],
            PubMedID = pubmed_id,
            Authors = authorlist,
        )
        if pub_details['PublicationDate']:
            pub_details['PublicationDate'] = pub_details['PublicationDate'].strftime('%Y-%m-%d')

        if not pub_details['URL'] and pub_details['DOI']:
            pub_details['URL'] = 'https://dx.doi.org/%s' % pub_details['DOI']
        return pub_details


    @informational_misc
    def get_publications(self):
        '''Returns the information (title, publication, authors etc.) for all publications.'''
        publications = {}
        for r in self.DDG_db.execute_select('SELECT ID FROM Publication'):
            publications[r['ID']] = self.get_publication(r['ID'])
        return publications

    def cache_all_score_method_details(self):
        '''Helper function for get_score_method_details'''
        score_methods = {}
        for r in self.DDG_db_utf.execute_select('SELECT * FROM ScoreMethod'):
            score_methods[r['ID']] = r
        self.cached_score_method_details = score_methods
    
    @informational_misc
    def get_score_method_details(self, score_method_id = None, allow_recaching = True):
        '''Returns all score method details, unless a score method id is passed, then only those details are returned'''
        if not self.cached_score_method_details:
            self.cache_all_score_method_details()

        if score_method_id:
            # Returns ScoreMethod record for specific score_method_id
            if score_method_id in self.cached_score_method_details:
                return self.cached_score_method_details[score_method_id]
            elif allow_recaching:
                # Recache and try one more time
                self.get_score_method_details(score_method_id = score_method_id, allow_recaching = False)
            else:
                # We have already tried again once, so fail
                raise Exception("score_method_id %s isn't in score methods table" % str(score_method_id))
        else:
            # Returns all defined ScoreMethod records
            return self.cached_score_method_details

    def output_score_method_information(self, score_method_id, analysis_set_id, output_directory):
        '''Outputs details about score method to a txt file in the specified output directory'''
        score_method_details = sorted([(k, v) for k, v in self.get_score_method_details(score_method_id = score_method_id).iteritems()])
        with open(os.path.join(output_directory, 'score_method.txt'), 'w') as f:
            f.write('Analysis set ID: %s\n' % str(analysis_set_id))
            f.write('Score method ID: %s\n' % str(score_method_id))
            f.write('\nScore method details\n')
            for key, value in score_method_details:
                f.write('%s: %s\n' % (str(key), str(value)))
        
    @informational_misc
    def get_score_method_id(self, method_name, method_type = None, method_parameters = None, method_authors = None, fuzzy = True):
        '''Returns the ID for the ScoreMethod with the specified parameters.
           If fuzzy is True then the string matching uses LIKE rather than equality.
           e.g. method_id = self.get_score_method_id('interface', method_authors = 'kyle')
           '''
        if fuzzy:
            match_phrase = 'LIKE %s'
            method_name = '%{0}%'.format(method_name)
            if method_type: method_type = '%{0}%'.format(method_type)
            if method_parameters: method_parameters = '%{0}%'.format(method_parameters)
            if method_authors: method_authors = '%{0}%'.format(method_authors)
        else:
            match_phrase = '=%s'

        condition_parameters = [method_name]
        conditions = ['MethodName {0}'.format(match_phrase)]
        if method_type:
            conditions.append('MethodType {0}'.format(match_phrase))
            condition_parameters.append(method_type)
        if method_parameters:
            conditions.append('Parameters {0}'.format(match_phrase))
            condition_parameters.append(method_parameters)
        if method_authors:
            conditions.append('Authors {0}'.format(match_phrase))
            condition_parameters.append(method_authors)
        conditions = ' AND '.join(conditions)
        condition_parameters = tuple(condition_parameters)

        results = self.DDG_db_utf.execute_select('SELECT ID FROM ScoreMethod WHERE {0}'.format(conditions), parameters=condition_parameters)
        if not results:
            raise Exception('Error: No ScoreMethod records were found using the criteria: {0}'.format(', '.join(map(str, [s for s in [method_name, method_type, method_parameters] if s]))))
        elif len(results) > 1:
            raise Exception('Error: Multiple ScoreMethod records were found using the criteria: {0}'.format(', '.join(map(str, [s for s in [method_name, method_type, method_parameters] if s]))))
        else:
            return results[0]['ID']


    @informational_misc
    def get_score_dict(self, prediction_id = None, score_method_id = None, score_type = None, structure_id = None, prediction_structure_scores_table = None, prediction_id_field = None):
        '''Returns a dict with keys for all fields in the Score table. The optional arguments can be used to set the
           corresponding fields of the dict. All other fields are set to None.'''

        if prediction_structure_scores_table == None:
            prediction_structure_scores_table = self._get_prediction_structure_scores_table()
        if prediction_id_field == None:
            prediction_id_field = self._get_prediction_id_field()

        # Relax the typing
        if structure_id: structure_id = int(structure_id)
        if prediction_id: prediction_id = int(prediction_id)
        if score_method_id: score_method_id = int(score_method_id)
        if score_type:
            allowed_score_types = self._get_allowed_score_types()
            if  score_type not in allowed_score_types:
                raise Exception('"{0}" is not an allowed score type. Allowed types are: "{1}".'.format(score_type, '", "'.join(sorted(allowed_score_types))))

        fieldnames = set([f for f in self.DDG_db.FieldNames.__dict__[prediction_structure_scores_table].__dict__.keys() if not(f.startswith('_'))])
        d = dict.fromkeys(fieldnames, None)
        if prediction_id_field != None:
            d[prediction_id_field] = prediction_id
        d['ScoreMethodID'] = score_method_id
        d['ScoreType'] = score_type
        d['StructureID'] = structure_id
        return d


    @informational_file
    def get_file_id(self, content, db_cursor = None, hexdigest = None):
        '''Searches the database to see whether the FileContent already exists. The search uses the digest and filesize as
           heuristics to speed up the search. If a file has the same hex digest and file size then we do a straight comparison
           of the contents.
           If the FileContent exists, the value of the ID field is returned else None is returned.
           '''
        existing_filecontent_id = None
        hexdigest = hexdigest or get_hexdigest(content)
        filesize = len(content)
        if db_cursor:
            db_cursor.execute('SELECT * FROM FileContent WHERE MD5HexDigest=%s AND Filesize=%s', (hexdigest, filesize))
            results = db_cursor.fetchall()
        else:
            results = self.DDG_db.execute_select('SELECT * FROM FileContent WHERE MD5HexDigest=%s AND Filesize=%s', parameters=(hexdigest, filesize))
        for r in results:
            if r['Content'] == content:
                assert(existing_filecontent_id == None) # content uniqueness check
                existing_filecontent_id = r['ID']
        return existing_filecontent_id


    @informational_pdb
    def get_pdb_chain_coordinates(self, pdb_id, chain_id):
        '''Read a saved dataframe.'''
        zipped_coordinates = self.DDG_db.execute_select('SELECT Coordinates FROM PDBChain WHERE PDBFileID=%s AND Chain=%s AND Coordinates IS NOT NULL', parameters=(pdb_id, chain_id))
        if zipped_coordinates:
            assert(len(zipped_coordinates) == 1)
            buf = BytesIO(zipped_coordinates[0]['Coordinates'])
            gf = gzip.GzipFile(fileobj=buf, mode="rb")

            residue_matrix = None
            try:
                store = pandas.read_hdf(gf)
                residue_matrix = store['dataframe']
                store.close()
            except NotImplementedError, e:
                # "Support for generic buffers has not been implemented"
                try:
                    nfname = None
                    f, nfname = open_temp_file('/tmp', suffix = '.hdf5')
                    f.close()
                    write_file(nfname, gf.read(), ftype = 'wb')
                    store = pandas.HDFStore(nfname)
                    residue_matrix = store['dataframe']
                    store.close()
                    os.remove(nfname)
                    print('get_pdb_chain_coordinates here')
                except:
                    if nfname: os.remove(nfname)
                    raise
            return residue_matrix
        return None


    @informational_pdb
    def get_pdb_chains_for_prediction(self, prediction_id):
        '''Returns the PDB file ID and a list of chains for the prediction.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_pdb
    def get_chain_sets_for_mutatagenesis(self, mutagenesis_id, complex_id = None):
        '''Gets a list of possibilities for the associated complex and calls get_chains_for_mutatagenesis on each.
           This function assumes that a complex structure is required i.e. that all chains in the PDB chain set are in the same PDB file.

           This is a useful method for listing the possible complexes to use in a prediction or to determine whether one
           may be missing. and we need to update the database.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_pdb
    def get_chains_for_mutatagenesis(self, mutagenesis_id, pdb_file_id, pdb_set_number, complex_id = None):
        '''Returns the PDB chains used in the mutagenesis.
           Note: At present, monomeric data e.g. protein stability does not have the notion of complex in our database
           but this abstraction is planned so that multiple choices of PDB file and chain can be easily represented.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_pdb
    def get_pdb_mutations_for_mutagenesis(self, mutagenesis_id, pdb_file_id, set_number, complex_id = None):
        '''Returns the PDB mutations for a mutagenesis experiment as well as the PDB residue information.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_pdb
    def get_pdb_details(self, pdb_ids, cached_pdb_details = None):
        '''Returns the details stored in the database about the PDB files associated with pdb_ids e.g. chains, resolution,
           technique used to determine the structure etc.'''
        raise Exception('Replace this with a call to import_api.py::DataImportInterface.get_pdb_details()')
        return self.importer.get_pdb_details(pdb_ids, cached_pdb_details = None)
        pdbs = {}
        cached_pdb_ids = []
        if cached_pdb_details:
            cached_pdb_ids = set(cached_pdb_details.keys())
        for pdb_id in pdb_ids:
            if pdb_id in cached_pdb_ids:
                pdbs[pdb_id] = cached_pdb_details[pdb_id]
            else:
                record = self.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))[0]
                p = PDB(record['Content'])
                pdb_chain_lengths = {}
                for chain_id, s in p.atom_sequences.iteritems():
                    pdb_chain_lengths[chain_id] = len(s)
                # todo: get the list of protein chains and PDB residues from the database and assert that they are the same
                #       as what were extracted from the PDB file.
                #       maybe change 'chains' below to 'protein_chains'
                pdbs[pdb_id] = dict(
                    chains = pdb_chain_lengths,
                    TM = record['Transmembrane'],
                    Technique = record['Techniques'],
                    XRay = record['Techniques'].find('X-RAY') != -1,
                    Resolution = record['Resolution'],
                )
        return pdbs


    @informational_pdb
    def get_prediction_set_pdb_chain_details(self, PredictionSet, cached_pdb_details = None, restrict_to_pdbs = set()):
        '''Used by the analysis API. This could be combined with get_pdb_details.'''

        pdb_ids = [r['PDBFileID'] for r in self.DDG_db.execute_select('SELECT DISTINCT PDBFileID FROM {0} INNER JOIN {1} ON {1}ID={1}.ID WHERE PredictionSet=%s ORDER BY PDBFileID'.format(self._get_prediction_table(), self._get_user_dataset_experiment_table()), parameters=(PredictionSet,))]

        if not pdb_ids:
            try:
                pdb_ids = [r['PDBFileID'] for r in self.DDG_db.execute_select('SELECT DISTINCT PDBFileID FROM {0} INNER JOIN {1} ON {1}ID={1}.ID WHERE PredictionSet=%s ORDER BY PDBFileID'.format(self._get_prediction_table(), 'Experiment'), parameters=(PredictionSet,))]
            except: pass

        if restrict_to_pdbs:
            pdb_ids = sorted(set(pdb_ids).intersection(restrict_to_pdbs))

        pdbs = {}
        cached_pdb_ids = []
        if cached_pdb_details:
            cached_pdb_ids = set(cached_pdb_details.keys())
        for pdb_id in pdb_ids:
            if pdb_id in cached_pdb_ids:
                pdbs[pdb_id] = cached_pdb_details[pdb_id]
            else:
                record = self.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))[0]
                p = PDB(record['Content'])

                d = {}
                chain_ids = set(p.chain_types.keys()).union(set(p.seqres_chain_order)).union(set(p.atom_sequences.keys()))
                d['Chains'] = dict.fromkeys(chain_ids)
                for chain_id in chain_ids:
                    d['Chains'][chain_id] = dict(
                        Sequence = str(p.atom_sequences.get(chain_id) or ''),
                        Type = p.chain_types.get(chain_id),
                    )
                d['Resolution'] = p.get_resolution()
                d['MethodOfDetermination'] = p.get_techniques()

                pdbs[pdb_id] = d
        return pdbs


    @informational_job
    def get_development_protocol(self, development_protocol_id):
        '''Possibly temporary function which returns a DevelopmentProtocol record from the database.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def get_complex_details(self, complex_id):
        '''Returns the database record for the given complex.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def get_job_details(self, prediction_id, include_files = True, truncate_content = None):
        '''Returns the details necessary to run the job.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @informational_job
    def get_job_files(self, prediction_id, truncate_content = None, set_pdb_occupancy_one = True):
        '''Returns a dict mapping the stages (e.g. 'input', 'output', 'analysis') of a job with the files associated with
           that stage.
           If truncate_content is set, it should be an integer specifying the amount of characters to include. This is useful
           to see if the file header is as expected.
        '''
        params = (self._get_prediction_table(),)
        qry = 'SELECT {0}File.*, FileContent.Content, FileContent.MIMEType, FileContent.Filesize, FileContent.MD5HexDigest FROM {0}File INNER JOIN FileContent ON FileContentID=FileContent.ID WHERE {0}ID=%s'.format(*params)
        results = self.DDG_db.execute_select(qry, parameters=(prediction_id,))
        job_files = {}
        for r in results:
            if truncate_content and str(truncate_content).isdigit():
                if len(r['Content']) > int(truncate_content):
                    r['Content'] = '%s...' % r['Content'][:int(truncate_content)]
            if set_pdb_occupancy_one and r['Filetype'] == 'PDB': # Set all occupancies to 1
                pdb = PDB(r["Content"].split("\n"))
                pdb.fillUnoccupied()
                r['Content'] = pdb.get_content()
            job_stage = r['Stage']
            del r['Stage']
            job_files[job_stage] = job_files.get(job_stage, [])
            job_files[job_stage].append(r)
        return job_files


    @informational_job
    def get_prediction_set_details(self, prediction_set_id):
        '''Returns the PredictionSet record from the database.'''
        results = self.DDG_db.execute_select('SELECT * FROM PredictionSet WHERE ID=%s', parameters=(prediction_set_id,))
        if len(results) == 1:
            results[0]['Job status summary'] = self._get_prediction_set_status_counts(prediction_set_id)
            return results[0]
        return None


    def _get_prediction_set_status_counts(self, prediction_set_id):
        '''Returns a summary of the prediction job statuses for the prediction set.'''
        d = {}
        for r in self.DDG_db.execute_select('SELECT Status, COUNT(Status) AS Count FROM {0} WHERE PredictionSet=%s GROUP BY Status'.format(self._get_prediction_table()), parameters=(prediction_set_id,)):
            d[r['Status']] = r['Count']
        return d


    @informational_job
    def get_prediction_ids(self, prediction_set_id):
        '''Returns the list of Prediction IDs associated with the PredictionSet.'''
        self._assert_prediction_set_is_correct_type(prediction_set_id)
        qry = 'SELECT ID FROM %s WHERE PredictionSet=%%s ORDER BY ID' % self._get_prediction_table()
        return [r['ID'] for r in self.DDG_db.execute_select(qry, parameters=(prediction_set_id,))]

    @informational_job
    def get_defined_user_datasets(self):
        '''Return a dict detailing the defined UserDataSets, their tagged subsets (if any), and the mutagenesis counts
          (i.e. the number of prediction cases) of both the user datasets and the associated tagged subsets .'''
        d = {}
        user_datasets = self.DDG_db.execute_select('SELECT * FROM UserDataSet WHERE DatasetType=%s', parameters=(self._get_prediction_dataset_type(),))
        for uds in user_datasets:
            qry = 'SELECT COUNT(ID) AS MutagenesisCount FROM %s WHERE UserDataSetID=%%s' % self._get_user_dataset_experiment_table()
            uds['MutagenesisCount'] = self.DDG_db.execute_select(qry, parameters=(uds['ID'],))[0]['MutagenesisCount']
            d[uds['TextID']] = uds

            subsets = {}
            if self._get_user_dataset_experiment_tag_table():
                qry = 'SELECT Tag, COUNT(Tag) AS MutagenesisCount FROM %s INNER JOIN %s ON %sID=%s.ID WHERE UserDataSetID=%%s GROUP BY Tag' % (self._get_user_dataset_experiment_tag_table(), self._get_user_dataset_experiment_table(), self._get_user_dataset_experiment_table(), self._get_user_dataset_experiment_table())
                for tagged_subset in self.DDG_db.execute_select(qry, parameters=(uds['ID'],)):
                    subsets[tagged_subset['Tag']] = dict(MutagenesisCount = tagged_subset['MutagenesisCount'])
            uds['Subsets'] = subsets
        return d


    @informational_job
    def get_user_dataset_experiment_ids(self, user_dataset_name, tagged_subset = None):
        '''Returns a list of UserDataSet experiment records for the given user dataset.'''
        if tagged_subset:
            qry = 'SELECT %s.* FROM %s INNER JOIN %s ON %sID=%s.ID INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%%s AND Tag=%%s' % (self._get_user_dataset_experiment_table(), self._get_user_dataset_experiment_tag_table(), self._get_user_dataset_experiment_table(), self._get_user_dataset_experiment_table(), self._get_user_dataset_experiment_table())
            return self.DDG_db.execute_select(qry, parameters=(user_dataset_name, tagged_subset))
        else:
            qry = 'SELECT %s.* FROM %s INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%%s' % (self._get_user_dataset_experiment_table(), self._get_user_dataset_experiment_table())
            return self.DDG_db.execute_select(qry, parameters=(user_dataset_name,))


    @informational_job
    def get_user_dataset_experiment_details(self, user_dataset_experiment_id, user_dataset_id = None):
        '''Returns all the data relating to a user dataset experiment.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def get_dataset_experiment_details(self, dataset_experiment_id, dataset_id = None):
        '''Returns the experimental data relating to a dataset experiment.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def export_dataset_to_json(self, dataset_id):
        '''Returns the dataset information in JSON format.'''
        return json.dumps(self._export_dataset(dataset_id))


    @informational_job
    def export_dataset_to_csv(self, dataset_id):
        '''Returns the dataset information in CSV format.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def get_predictions_experimental_details(self, prediction_id, userdatset_experiment_ids_to_subset_ddgs = None, include_files = False, reference_ids = set(), include_experimental_data = True):
        '''Returns a dict containing the experimental details for the Prediction. This is what is used by export_prediction_cases_to_json etc.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def get_experimental_ddgs_by_analysis_set(self, user_dataset_experiment_id = None, reference_ids = set()):
        '''Returns a mapping from UserPPDataSetExperimentIDs to dicts mapping analysis Subsets to a dicts containing the
           record identifier triple (subset, section, record number), the experimental DDG values, the mean of those values,
           and whether the values / one of the values are derived from other measurements e.g.

          23 : {
            'BeAtMuSiC' : {'Cases'           : set([('BeAtMuSiC', 'Main', 1408L)]),
                           'DDGs'            : [{'IsDerivedValue': 0L,
                                                 'Value': 2.9802478611},
                                                {'IsDerivedValue': 0L,
                                                 'Value': 2.1978328374}],
                           'IsDerivedValue'  : 0L,
                           'MeanDDG'         : 2.5890403492500003},
            'SKEMPI'    : {'Cases'           : set([('SKEMPI', 'Non-derivative', 1L)]),
                           'DDGs'            : [{'IsDerivedValue': 0L, 'Value': 2.9802478611},
                                              {'IsDerivedValue': 0L, 'Value': 2.1978328374}],
                           'IsDerivedValue'  : 0L,
                           'MeanDDG'         : 2.5890403492500003},
            'ZEMu'      : {'Cases'           : set([('ZEMu', 'Main', 1144L)]),
                           'DDGs'            : [{'IsDerivedValue': 0L, 'Value': 2.1978328374}],
                           'IsDerivedValue'  : 0L,
                           'MeanDDG'         : 2.1978328374}}
           ...

           This can be used to: i) generate histograms showing the spread of experimental values for a dataset; or
           ii) to add columns to an analysis dataframe so that, once created, it can be analyzed over multiple analysis sets.
        '''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_job
    def get_prediction_set_case_details(self, prediction_set_id, retrieve_references = True, include_experimental_data = True):
        '''Returns a dict containing the case information for prediction cases in the prediction set with a structure
           expected by the analysis class.'''

        # Read the Prediction details
        reference_ids = set()
        prediction_ids = self.get_prediction_ids(prediction_set_id)
        userdatset_experiment_ids_to_subset_ddgs = {}
        if include_experimental_data:
            userdatset_experiment_ids_to_subset_ddgs = self.get_experimental_ddgs_by_analysis_set(reference_ids = reference_ids)

        prediction_cases = {}
        for prediction_id in prediction_ids:
            prediction_cases[prediction_id] = self.get_predictions_experimental_details(prediction_id, userdatset_experiment_ids_to_subset_ddgs, include_experimental_data = include_experimental_data)

        references = {}
        if retrieve_references:
            for reference_id in sorted(reference_ids):
                references[reference_id] = self.get_publication(reference_id)

        return dict(
            Data = prediction_cases,
            References = references,
            PredictionSet = self.get_prediction_set_details(prediction_set_id)
            )


    @informational_job
    def export_prediction_cases_to_json(self, prediction_set_id, retrieve_references = True):
        '''A JSON wrapper to get_prediction_set_case_details.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    def _export_dataset(self, dataset_id):
        '''Returns a dict containing the dataset information.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    ###########################################################################################
    ## Prediction creation/management layer
    ##
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    @job_creator
    def add_prediction_set(self, prediction_set_id, halted = True, priority = 5, batch_size = 40,
                           allow_existing_prediction_set = False, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False,
                           series_name = None, series_color = 'ff0000', series_alpha = 1.0, description = None):
        '''Adds a new PredictionSet (a construct used to group Predictions) to the database.
           If a PredictionSet is halted then running schedulers will not kick off the jobs. Otherwise, they will be queued
           depending on the priority of the PredictionSet (higher numbers mean higher priority).
           batch_size defines the number of jobs to be added as once as an array job.
           Returns True if a PredictionSet with the same ID did not previously exist.
           The priority and batch size can be modified while a scheduler is running and will affect the next round of
           predictions to be queued.
           Raises an exception or returns False otherwise depending on the value of allow_existing_prediction_set.'''

        if halted:
            Status = 'halted'
        else:
            Status = 'active'

        existing_record = self.DDG_db.execute_select('SELECT * FROM PredictionSet WHERE ID=%s', parameters=(prediction_set_id,))
        if self.prediction_set_exists(prediction_set_id):
            if allow_existing_prediction_set == False:
                raise Exception('The PredictionSet %s already exists.' % prediction_set_id)
            else:
                return False

        d = dict(
            ID                  = prediction_set_id,
            Status              = Status,
            Priority            = priority,
            ProteinStability    = contains_protein_stability_predictions,
            BindingAffinity     = contains_binding_affinity_predictions,
            BatchSize           = batch_size,
            SeriesName          = series_name,
            SeriesColor         = series_color,
            SeriesAlpha         = series_alpha,
            Description         = description,
        )
        self.DDG_db.insertDictIfNew("PredictionSet", d, ['ID'])
        return True

    def prediction_set_exists(self, prediction_set_id):
        existing_record = self.DDG_db.execute_select('SELECT * FROM PredictionSet WHERE ID=%s', parameters=(prediction_set_id,))
        if len(existing_record) > 0:
            assert(len(existing_record) == 1)
            return True
        else:
            return False

    @job_creator
    def destroy_prediction_set(self, prediction_set_id):
        '''This function removes the PredictionSet from the database.
           THIS CANNOT BE UNDONE.
           For safety, we should only allow PredictionSets with no corresponding scores to be removed.
           It fits into the job_creator category since usually these empty PredictionSets will have been created while
           setting up a job.'''

        can_be_deleted = self.DDG_db.execute_select('SELECT CanBeDeleted FROM PredictionSet WHERE ID=%s', parameters=(prediction_set_id,))
        if len(can_be_deleted) == 0:
            raise colortext.Exception('The prediction set "%s" does not exist.' % prediction_set_id)
        elif can_be_deleted[0]['CanBeDeleted'] == 0:
            raise colortext.Exception('The prediction set "%s" is not allowed to be deleted. Change the CanBeDeleted property on its record first.' % prediction_set_id)

        params = (self._get_prediction_table(), self._get_prediction_structure_scores_table())

        qry = 'SELECT COUNT({0}.ID) AS NumRecords FROM {0} INNER JOIN {1} ON {0}.ID={1}.{0}ID WHERE PredictionSet=%s'.format(*params)
        existing_scores = self.DDG_db.execute_select(qry, parameters=(prediction_set_id,))
        if existing_scores[0]['NumRecords'] > 0:
            raise colortext.Exception('Cannot remove a prediction set with associated scores.')

        qry = 'SELECT COUNT(ID) AS NumRecords FROM {0} WHERE Status <> "queued" AND PredictionSet=%s'.format(*params)
        jobs_in_flux = self.DDG_db.execute_select(qry, parameters=(prediction_set_id,))
        if jobs_in_flux[0]['NumRecords'] > 0:
            raise colortext.Exception('Cannot remove a prediction set unless all jobs are set as "queued".')

        # Use a transaction to prevent a partial deletion
        self.DDG_db._get_connection()
        con = self.DDG_db.connection
        try:
            with con:
                cur = con.cursor()

                # Delete the associated file records
                delete_files_qry = 'DELETE {0}File FROM {0}File INNER JOIN {0} ON {0}.ID={0}File.{0}ID WHERE PredictionSet=%s'.format(*params)
                cur.execute(delete_files_qry, (prediction_set_id, ))

                # Delete the predictions
                delete_predictions_qry = 'DELETE FROM {0} WHERE PredictionSet=%s'.format(*params)
                cur.execute(delete_predictions_qry, (prediction_set_id, ))

                cur.execute('DELETE FROM PredictionSet WHERE ID=%s', (prediction_set_id, ))
        except Exception, e:
            raise colortext.Exception('An exception occurred removing the PredictionSet from the database: "%s".\n%s' % (str(e), traceback.format_exc()))


    @job_creator
    def start_prediction_set(self, PredictionSetID):
        '''Sets the Status of a PredictionSet to "active".'''
        self._set_prediction_set_status(PredictionSetID, 'active')


    @job_creator
    def stop_prediction_set(self, PredictionSetID):
        '''Sets the Status of a PredictionSet to "halted".'''
        self._set_prediction_set_status(PredictionSetID, 'halted')


    @job_creator
    def add_job(self, *args, **kwargs):
        '''Add a single prediction job to a prediction set. This should not typically be called - add_prediction_run
           is generally what should be called instead.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_creator
    def add_jobs_by_pdb_id(self, *args, **kwargs):
        ''' This function adds predictions for all Experiments corresponding to pdb_ID to the specified prediction set.
            This is useful for custom runs e.g. when we are using the DDG scheduler for design rather than for benchmarking.
            Variants of this function were used before for CypA and ubiquitin runs.
            This is currently unimplemented but ask Shane if we need this functionality again.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    def _add_prediction_run_preconditions(self, prediction_set_id, user_dataset_name, tagged_subset):
        '''Check to make sure that the prediction set, user dataset, and optional tagged subset make sense for this API.'''
        qry = 'SELECT %s FROM PredictionSet WHERE ID=%%s' % (self._get_prediction_type())
        prediction_set = self.DDG_db.execute_select(qry, parameters=(prediction_set_id,))
        if not prediction_set:
            raise colortext.Exception('The prediction set "%s" does not exist in the database.' % prediction_set_id)
        elif prediction_set[0][self._get_prediction_type()] != 1:
            raise colortext.Exception('The prediction set "%s" is not the correct type ("%s") for this API.' % (prediction_set_id, self._get_prediction_type()))

        allowed_user_datasets = self.get_defined_user_datasets()
        if user_dataset_name not in allowed_user_datasets:
            raise colortext.Exception('The user dataset "%s" does not exist in the database.' % user_dataset_name)
        if tagged_subset and tagged_subset not in allowed_user_datasets[user_dataset_name]['Subsets']:
            raise colortext.Exception('The tagged subset "%s" of user dataset "%s" does not exist in the database.' % (tagged_subset, user_dataset_name))


    @job_creator
    def add_prediction_run(self, *args, **kwargs):
        '''Adds all jobs corresponding to a user dataset e.g. add_prediction_run("my first run", "AllBindingAffinityData", tagged_subset = "ZEMu").'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_creator
    def add_job_by_user_dataset_record(self, *args, **kwargs):
        '''Uses the UserDataSet record to get most of the information needed to set up the job e.g. PDB complex, mutagenesis details.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_creator
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set, *args, **kwargs):
        '''add_prediction_run sets up a full run of dataset predictions but is slow as it needs to perform a lot of
           calculations and parsing. If you want to test the same dataset with slightly different parameters (e.g. a
           different protocol) then these calculations can be reused which reduces the overhead considerably.
           clone_prediction_run was written with this in mind. It copies the list of predictions and their setup (input
           files etc.) from an existing prediction set to an empty prediction set.'''
        raise Exception('not implemented yet')


    @job_creator
    def add_development_protocol_command_lines(self, development_protocol_id):
        '''Possibly temporary function used to add protocol command lines for methods in development.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    #== Input file generation API ===========================================================
    #
    # This part of the API is responsible for creating input files for predictions


    @job_input
    def create_resfile(self, prediction_id):
        '''This function returns the resfile content for the prediction. It is usually not called directly by the user but
           is available for convenience and debugging.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_input
    def create_mutfile(self, prediction_id):
        '''This function returns the mutfile content for the prediction. It is usually not called directly by the user but
           is available for convenience and debugging.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    #== Job execution/completion API ===========================================================
    #
    # This part of the API is responsible for starting jobs and setting them as failed or
    # completed


    @job_execution
    def get_queued_jobs(self, prediction_set_id, order_by = 'Cost', order_order_asc = False, include_files = True, truncate_content = None):
        '''An iterator to return the details of the queued prediction records in this prediction set.
           An exception is raised if the prediction set is halted.
           Assuming Cost is filled in and is representative of the expected runtime, it makes sense to request jobs ordered
           by Cost and order_order_asc = False rather than by ID as longer jobs can then be kicked off before shorter jobs.

           Usage:
               for prediction_record in ppi_api.get_queued_jobs(prediction_set_id, include_files = True, truncate_content = 30):
                   pprint.pprint(prediction_record)
        '''

        if self.get_prediction_set_details(prediction_set_id)['Status'] == 'halted':
            raise Exception('The prediction set is halted so no job details can be returned.')

        self._assert_prediction_set_exists(prediction_set_id)
        for job_id in self.get_queued_job_list(prediction_set_id, order_by = order_by, order_order_asc = order_order_asc):
            self._get_job_fn_call_counter[job_id] = self._get_job_fn_call_counter.get(job_id, 0)
            self._get_job_fn_call_counter[job_id] += 1
            if self._get_job_fn_call_counter[job_id] > self._get_job_fn_call_counter_max:
                self.DDG_db = None
                self.DDG_db_utf = None
                raise Exception('get_job was called %d times for this prediction. This is probably a bug in the calling code.' % self._get_job_fn_call_counter[job_id])
            yield(self.get_job_details(job_id, include_files = include_files, truncate_content = truncate_content))


    @job_execution
    def get_queued_job_list(self, prediction_set_id, order_by = 'Cost', order_order_asc = False):
        '''An iterator to return the list of queued prediction records in this prediction set.
           Assuming Cost is filled in and is representative of the expected runtime, it makes sense to request jobs ordered
           by Cost and order_order_asc = False rather than by ID as longer jobs can then be kicked off before shorter jobs.

           Usage:
               for prediction_id in ppi_api.get_queued_job_list(prediction_set_id):
                   print(prediction_id)
        '''

        assert((order_by in ['Cost', 'ID']) and isinstance(order_order_asc, bool))
        if order_order_asc:
            order_order_asc = 'ASC'
        else:
            order_order_asc = 'DESC'
        params = (self._get_prediction_table(), order_by, order_order_asc)
        qry = 'SELECT ID FROM {0} WHERE PredictionSet=%s AND Status="queued" ORDER BY {1} {2}'.format(*params)
        results = self.DDG_db.execute_select(qry, parameters=(prediction_set_id,))
        x = 0
        while x < len(results):
            yield results[x]['ID']
            x += 1


    @job_execution
    def start_job(self, prediction_id, prediction_set_id):
        '''Sets the job status to "active". prediction_set must be passed and is used as a sanity check.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_execution
    def get_max_number_of_cluster_jobs(self, prediction_set_id, priority):
        '''Returns the maximum number of cluster jobs that schedulers should run for this interface.'''
        return self.DDG_db.execute_select('SELECT Value FROM _DBCONSTANTS WHERE VariableName="MaxClusterJobs"')['Value']


    def _assert_prediction_set_exists(self, prediction_set_id):
        if len(self.DDG_db.execute_select('SELECT * FROM PredictionSet WHERE ID=%s', parameters=(prediction_set_id,))) != 1:
            raise Exception('The PredictionSet %s does not exist.' % prediction_set_id)


    @job_execution
    def alter_prediction_set_priority(self, prediction_set_id, priority):
        '''Modify the priority for a PredictionSet. Higher values give the PredictionSet more priority over other running PredictionSets.'''
        priority = int(priority)
        assert(priority > 0)
        self._assert_prediction_set_exists(prediction_set_id)
        self.DDG_db.execute_select('UPDATE PredictionSet SET Priority=%s WHERE ID=%s', parameters=(priority, prediction_set_id,))


    @job_execution
    def alter_prediction_set_batch_size(self, prediction_set_id, batch_size):
        '''Modify the batch size for a PredictionSet. The batch size is the number of jobs which will be submitted together
           during subsequent job submissions.'''
        batch_size = int(batch_size)
        assert(batch_size > 0)
        self._assert_prediction_set_exists(prediction_set_id)
        self.DDG_db.execute_select('UPDATE PredictionSet SET BatchSize=%s WHERE ID=%s', parameters=(batch_size, prediction_set_id,))


    @job_execution
    def set_job_temporary_protocol_field(self, prediction_id, prediction_set_id, temporary_protocol_field):
        '''Possibly temporary function which sets fields in the temporary protocol field.'''
        raise Exception('not implemented yet')


    @job_completion
    def fail_job(self, prediction_id, prediction_set, maxvmem, ddgtime, errors = None):
        '''Sets the job status to "failed". prediction_set must be passed and is used as a sanity check.'''
        self._check_prediction(prediction_id, prediction_set)
        self.DDG_db.execute('UPDATE {0} SET Status="failed", maxvmem=%s, DDGTime=%s, Errors=%s WHERE ID=%s'.format(self._get_prediction_table()), parameters=(maxvmem, ddgtime, errors, prediction_id,))


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
        prediction_ids = self.get_prediction_ids(prediction_set_id)
        for prediction_id in prediction_ids:
            job_path = os.path.join(root_directory, prediction_id)
            if not os.path.exists(job_path):
                raise Exception('The folder {0} for Prediction #{1} does not exist.'.format(job_path, prediction_id))

        for prediction_id in prediction_ids:
            self.extract_data_for_case(prediction_id, root_directory = root_directory, force = force, score_method_id = score_method_id)


    @job_completion
    def extract_data_for_case(self, prediction_id, root_directory = None, score_method_id = None, force = False):
        '''Extracts the data for the prediction case (e.g. by processing stdout) and stores it in the Prediction*StructureScore
           table.
           The scores are returned to prevent the caller having to run another query.

           If force is False and the expected number of records for the case exists in the database, these are returned.
           Otherwise, the data are extracted, stored using a database transaction to prevent partial storage, and returned.

           Note:
           We use a lot of functions here: extract_data_for_case, parse_prediction_scores, store_scores.
           This may seem like overkill but I think it could allow us to reuse a lot of the code since the tables for
           PredictionPPIStructureScore and PredictionStructureScore are very similar (but are different since the underlying
           foreign tables PredictionPPI and Prediction are at least currently separate).

           parse_prediction_scores only returns dicts for database storage so it can be useful for debugging during development.
           store_scores stores scores in the database (passed as a list of dicts) but does not care from where they came.
           extract_data_for_case calls parse_prediction_scores to get the scores and the store_scores to commit them to the database.
        '''

        root_directory = root_directory or self.prediction_data_path # defaults to /kortemmelab/shared/DDG/ppijobs
        prediction_set = self.get_job_details(prediction_id, include_files = False)['PredictionSet']

        # todo: implement force behavior

        # Create a list of dicts for the PredictionPPIStructureScore table
        scores = self.parse_prediction_scores(prediction_id, root_directory = root_directory, score_method_id = score_method_id)

        # Store the dicts as PredictionPPIStructureScore records
        if len(scores) > 0:
            self.store_scores(prediction_set, prediction_id, scores)

        return scores


    @job_completion
    def parse_prediction_scores(self, prediction_id, root_directory = None, score_method_id = None):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass. Returns a dict suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.')


    @job_completion
    def store_scores_for_many_predictions(self, prediction_set, scores, safe = True, prediction_structure_scores_table = None, prediction_id_field = None):
        '''Stores scores for many predictions.
           scores should be a list of dicts suitable for database storage e.g. PredictionStructureScore or
           PredictionPPIStructureScore records.
        '''

        prediction_id_field = prediction_id_field or self._get_prediction_id_field()
        prediction_structure_scores_table = prediction_structure_scores_table or self._get_prediction_structure_scores_table()
        
        if safe:
            # Sanity checks
            for score in scores:
                if prediction_id_field not in score:
                    raise Exception('The score record is missing a {0} field: {1}.'.format(prediction_id_field, pprint.pformat(score)))
                self._check_prediction(score[prediction_id_field], prediction_set)
        con = self.DDG_db.connection
        cursor = con.cursor()
        sql_query = None
        if safe:
            params_to_insert = set()
        else:
            params_to_insert = []
        for score in scores:
            if safe:
                sql, params, record_exists = self.DDG_db.create_insert_dict_string(prediction_structure_scores_table, score, PKfields = [prediction_id_field, 'ScoreMethodID', 'ScoreType', 'StructureID'], check_existing = True)
            else:
                sql, params, record_exists = self.DDG_db.create_insert_dict_string(prediction_structure_scores_table, score, PKfields = [prediction_id_field, 'ScoreMethodID', 'ScoreType', 'StructureID'], check_existing = False)
                
            if sql_query:
                assert( sql == sql_query )
            else:
                sql_query = sql
            if safe:
                if params in params_to_insert or record_exists:
                    print params
                    print params_list
                    raise Exception('Duplicate params')
                params_to_insert.add(params)
            else:
                params_to_insert.append(params)

        with con:
            db_cursor = con.cursor()
            if safe:
                db_cursor.executemany(sql_query, [x for x in params_to_insert])
            else:
                # print params_to_insert
                db_cursor.executemany(sql_query, params_to_insert)
           
    @job_completion
    def store_scores(self, prediction_set, prediction_id, scores, prediction_structure_scores_table = None, prediction_id_field = None):
        '''Stores scores for one prediction.
           scores should be a list of dicts suitable for database storage e.g. PredictionStructureScore or
           PredictionPPIStructureScore records.
           This function uses a transaction so if any of the insertions fail then they are all rolled back.

           The default scores table and prediction_id_field can be (evilly) overridden to put scores in the wrong table
           '''
        if prediction_set:
            # Only check prediction is in prediction set if prediction set is passed in
            self._check_prediction(prediction_id, prediction_set)
        if prediction_id_field == None:
            # Only check for self-consistency if we're not (evilly) overriding everything that is good in the world
            self._check_scores_for_main_fields(scores, prediction_id)
        if prediction_structure_scores_table == None:
            # Only check for self-consistency if we're not (evilly) overriding everything our forefathers died for
            self._check_score_fields(scores)

        prediction_structure_scores_table = prediction_structure_scores_table or self._get_prediction_structure_scores_table()
        prediction_id_field = prediction_id_field or self._get_prediction_id_field()
        
        try:
            con = self.DDG_db.connection
            with con:
                db_cursor = con.cursor()
                for score in scores:
                    sql, params, record_exists = self.DDG_db.create_insert_dict_string(prediction_structure_scores_table, score, PKfields = [prediction_id_field, 'ScoreMethodID', 'ScoreType', 'StructureID'], check_existing = True)
                    if not record_exists:
                        db_cursor.execute(sql, params)
        except Exception, e:
            raise colortext.Exception('Failed to insert scores for Prediction #{0}: "{1}".\n{2}'.format(prediction_id, str(e), traceback.format_exc()))


    @job_completion
    def complete_job(self, prediction_id, prediction_set, scores, maxvmem, ddgtime):
        '''Sets a job to 'completed' and stores scores. prediction_set must be passed and is used as a sanity check.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    ###########################################################################################
    ## Prediction results layer
    ##
    ## This part of the API for returning data about completed predictions.
    ###########################################################################################


    @job_results
    def get_ddg_scores_per_structure(self, prediction_id):
        '''Returns the list of all DDG scores for a prediction_id. NOTE: Consider allowing the score method to be passed as a parameter.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @job_results
    def get_prediction_data_path(self):
        '''Returns the file server path to the where archived prediction data is stored.'''
        return self.prediction_data_path


    @job_results
    def get_job_data(self, prediction_id):
        '''Returns (in memory) the contents of the zip file corresponding to the prediction.'''
        job_data_path = os.path.join(self.prediction_data_path, '%d.zip' % prediction_id)
        if os.path.exists(job_data_path):
            return read_file(job_data_path, binary = True)


    @job_results
    def write_job_archive_to_file(self, prediction_id, output_filename):
        '''Writes the contents of the zip file corresponding to the prediction.'''
        job_data_path = os.path.join(self.prediction_data_path, '%d.zip' % prediction_id)
        assert(output_filename != job_data_path) # do not overwrite the existing file or allow to extract in place
        write_file(output_filename, self.get_job_data(prediction_id))


    @job_results
    def write_job_data_to_disk(self, prediction_id, output_path):
        '''Saves the job output for the prediction to the specified path.'''
        assert(os.path.exists(output_path))
        assert(output_path != self.prediction_data_path) # do not overwrite the existing file or allow to extract in place
        archive = self.get_job_data(prediction_id)
        write_file(os.path.join(output_path, '%d.zip' % prediction_id), archive, 'wb')
        p = Popen(output_path, ['unzip', '%d.zip' % prediction_id])
        os.remove(os.path.join(output_path, '%d.zip' % prediction_id))
        if p.errorcode != 0:
            raise colortext.Exception(p.stderr)
        else:
            colortext.warning(p.stdout)


    @job_results
    def extract_sge_job_stdout_from_archive(self, prediction_id):
        '''Returns the stdout files created during the prediction.
           The files are returned as a dict mapping with the type of output file (e.g. ddg_monomer step) to the content
           of the stdout files.
        '''

        # Retrieve and unzip results in memory
        archive = self.get_job_data(prediction_id)
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


    ###########################################################################################
    ## Analysis layer
    ##
    ## This part of the API is responsible for running analysis on completed predictions
    ###########################################################################################


    @analysis_api
    def get_top_x_scores(self, prediction_id, score_method_id, score_type, x, component = 'total', order_by = 'ASC'):
        '''get_top_x_ddg_stability'''
        results = self.DDG_db.execute_select('SELECT * FROM {0} WHERE {1}=%s AND ScoreMethodID=%s AND ScoreType=%s ORDER BY {2} {3}'.format(self._get_prediction_structure_scores_table(), self._get_prediction_id_field(), component, order_by), parameters=(prediction_id, score_method_id, score_type))
        if len(results) < x:
            raise Exception('The top {0} best scores were requested but only {1} results are stored in the database.'.format(x, len(results)))
        results = results[:x]
        return [{
                    'PredictionID' : r[self._get_prediction_id_field()],
                    'ScoreMethodID' : score_method_id,
                    'ScoreType' : score_type,
                    'StructureID' : r['StructureID'],
                    component : r[component],
                } for r in results]


    @analysis_api
    def get_prediction_scores(self, prediction_id, expectn = None):
        '''Returns the scores for the prediction using nested dicts with the structure:
                ScoreMethodID -> StructureID -> ScoreType -> database record
        '''
        scores = {}
        for r in self.DDG_db.execute_select('SELECT * FROM {0} WHERE {1}=%s'.format(self._get_prediction_structure_scores_table(), self._get_prediction_id_field()), parameters=(prediction_id,)):
            ScoreMethodID = r['ScoreMethodID']
            ScoreType = r['ScoreType']
            StructureID = r['StructureID']
            if StructureID == -1:
                StructureID = 'None' # usually this indicates an overall or aggregate value
            scores[ScoreMethodID] = scores.get(ScoreMethodID, {})
            scores[ScoreMethodID][StructureID] = scores[ScoreMethodID].get(StructureID, {})
            scores[ScoreMethodID][StructureID][ScoreType] = r
            del scores[ScoreMethodID][StructureID][ScoreType]['ScoreMethodID']
            del scores[ScoreMethodID][StructureID][ScoreType]['StructureID']
            del scores[ScoreMethodID][StructureID][ScoreType]['ScoreType']
            del scores[ScoreMethodID][StructureID][ScoreType][self._get_prediction_id_field()]
            del scores[ScoreMethodID][StructureID][ScoreType]['ID']

        if expectn != None:
            for score_method_id, score_method_scores in scores.iteritems():
                num_cases = 0
                for k in score_method_scores.keys():
                    if isinstance(k, int) or isinstance(k, long):
                        num_cases += 1
                if num_cases < expectn:
                    raise Exception('Expected scores for at least {0} runs with score method {1}; found {2}. Prediction id: {3}.'.format(expectn, score_method_id, num_cases, prediction_id))

        return scores


    @analysis_api
    def get_top_x_ddg(self, prediction_id, score_method_id, top_x = 3, expectn = None):
        '''Returns the TopX value for the prediction. Typically, this is the mean value of the top X predictions for a
           case computed using the associated Score records in the database.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @analysis_api
    def get_analysis_dataframe(self, prediction_set_id,
            experimental_data_exists = True,
            prediction_set_series_name = None, prediction_set_description = None, prediction_set_credit = None,
            prediction_set_color = None, prediction_set_alpha = None,
            use_existing_benchmark_data = True,
            include_derived_mutations = False,
            use_single_reported_value = False,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            score_method_id = None,
            expectn = None,
            allow_failures = False,
            extract_data_for_case_if_missing = True,
            ):
        '''This function uses experimental data from the database and prediction data from the Prediction*StructureScore
           table to build a pandas dataframe and store it in the database. See .analyze for an explanation of the
           parameters.

           The dataframes mostly contain redundant data so their storage could be seen to break a key database design
           principal. However, we store the dataframe in the database as it can take a while to build it from scratch and
           pre-built dataframes can be used to run quick analysis, for rapid development of the analysis methods, or to
           plug into webservers where responsiveness is important.

           If use_existing_benchmark_data is True and the dataframe already exists then it is returned as a BenchmarkRun object.
           Otherwise, it is built from the Prediction*StructureScore records.
           If the Prediction*StructureScore records do not exist, this function falls back into extract_data_for_case
           to generate them in which case root_directory needs to be specified (this is the only use for the root_directory
           parameter).
        '''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @analysis_api
    def get_prediction_data(self, prediction_id, score_method_id, main_ddg_analysis_type, top_x = 3, expectn = None, extract_data_for_case_if_missing = True, root_directory = None):
        '''Returns a dictionary with values relevant to predictions e.g. binding affinity, monomeric stability.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    def _get_analysis_dataframe(self, benchmark_run_class,
            dataframe_type = None,
            prediction_set_id = None,
            experimental_data_exists = True,
            prediction_set_series_name = None, prediction_set_description = None, prediction_set_credit = None,
            prediction_set_color = None, prediction_set_alpha = None,
            use_existing_benchmark_data = True,
            include_derived_mutations = False,
            use_single_reported_value = False,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            report_analysis = True,
            silent = False,
            root_directory = None, # where to find the prediction data on disk
            score_method_id = None,
            expectn = None,
            allow_failures = False,
            extract_data_for_case_if_missing = True,
            debug = False
            ):
        '''This 'private' function does most of the work for get_analysis_dataframe.'''

        ddg_analysis_type = 'DDG_Top%d' % take_lowest
        assert(dataframe_type != None and prediction_set_id != None)
        hdf_store_blob = None
        if use_existing_benchmark_data:
            hdf_store_blob = self.DDG_db.execute_select('''
            SELECT PandasHDFStore FROM AnalysisDataFrame WHERE
               PredictionSet=%s AND DataFrameType=%s AND ContainsExperimentalData=%s AND ScoreMethodID=%s AND UseSingleReportedValue=%s AND TopX=%s AND BurialCutoff=%s AND
               StabilityClassicationExperimentalCutoff=%s AND StabilityClassicationPredictedCutoff=%s AND
               IncludesDerivedMutations=%s AND DDGAnalysisType=%s''', parameters=(
                    prediction_set_id, dataframe_type, experimental_data_exists, score_method_id, use_single_reported_value, take_lowest, burial_cutoff,
                    stability_classication_experimental_cutoff, stability_classication_predicted_cutoff, include_derived_mutations, ddg_analysis_type))
            if hdf_store_blob:
                assert(len(hdf_store_blob) == 1)
                mem_zip = StringIO.StringIO()
                mem_zip.write(hdf_store_blob[0]['PandasHDFStore'])
                mem_zip.seek(0)
                hdf_store_blob = gzip.GzipFile(fileobj = mem_zip, mode='rb').read()

        # This dict is similar to dataset_cases in the benchmark capture (dataset.json)
        prediction_set_case_details = None
        prediction_ids = []
        if not(use_existing_benchmark_data and hdf_store_blob):
            print('Retrieving the associated experimental data for the user dataset.')
            prediction_set_case_details = self.get_prediction_set_case_details(prediction_set_id, retrieve_references = True, include_experimental_data = experimental_data_exists)
            prediction_ids = prediction_set_case_details['Data'].keys()
            prediction_set_case_details = prediction_set_case_details['Data']

        analysis_data = {}
        top_level_dataframe_attributes = {}
        if not(use_existing_benchmark_data and hdf_store_blob):
            if extract_data_for_case_if_missing and not silent:
                print('Computing the TopX values for each prediction case, extracting data if need be.')
            elif not extract_data_for_case_if_missing and not silent:
                print('Computing the TopX values for each prediction case; skipping missing data without attempting to extract.')
            num_predictions_in_prediction_set = len(prediction_ids)
            failed_cases = set()
            for prediction_id in prediction_ids:
                try:
                    analysis_data[prediction_id] = self.get_prediction_data(prediction_id, score_method_id, ddg_analysis_type, top_x = take_lowest, expectn = expectn, extract_data_for_case_if_missing = extract_data_for_case_if_missing, root_directory = root_directory, dataframe_type = dataframe_type)
                except FatalException, e:
                    raise
                except PartialDataException, e:
                    if not allow_failures:
                        raise Exception('Prediction {0} has partial data. Skipping.'.format(prediction_id))
                    failed_cases.add(prediction_id)
                except Exception, e:
                    if not allow_failures:
                        raise Exception('An error occurred during the TopX computation: {0}.\n{1}'.format(str(e), traceback.format_exc()))
                    failed_cases.add(prediction_id)
                if debug and len(analysis_data) >= 20:
                    break

                # best_pair_id = self.determine_best_pair(prediction_id, score_method_id)

            if failed_cases:
                colortext.error('Failed to determine the TopX score for {0}/{1} predictions. Continuing with the analysis ignoring these cases.'.format(len(failed_cases), len(prediction_ids)))
            working_prediction_ids = sorted(set(prediction_ids).difference(failed_cases))
            top_level_dataframe_attributes = dict(
                num_predictions_in_prediction_set = num_predictions_in_prediction_set,
                num_predictions_in_dataframe = len(working_prediction_ids),
                dataframe_type = dataframe_type,
                contains_experimental_data = experimental_data_exists,
            )

        # Only pull PDB data for cases where we have data
        restrict_to_pdbs = set([prediction_set_case_details[k]['Structure']['PDBFileID'] for k in analysis_data])

        prediction_set_details = self.get_prediction_set_details(prediction_set_id)
        prediction_set_series_name = prediction_set_series_name or prediction_set_details['SeriesName'] or prediction_set_details['ID']
        prediction_set_description = prediction_set_description or prediction_set_details['Description']
        prediction_set_color = prediction_set_color or prediction_set_details['SeriesColor']
        prediction_set_alpha = prediction_set_alpha or prediction_set_details['SeriesAlpha']

        # Initialize the BindingAffinityBenchmarkRun object
        # Note: prediction_set_case_details, analysis_data, and top_level_dataframe_attributes will not be filled in
        benchmark_run = benchmark_run_class(
                prediction_set_series_name,
                prediction_set_case_details,
                analysis_data,
                contains_experimental_data = experimental_data_exists,
                store_data_on_disk = False,
                benchmark_run_directory = None,
                use_single_reported_value = use_single_reported_value,
                description = prediction_set_description,
                dataset_description = prediction_set_description,
                credit = prediction_set_credit,
                include_derived_mutations = include_derived_mutations,
                take_lowest = take_lowest,
                generate_plots = False,
                report_analysis = report_analysis,
                silent = silent,
                burial_cutoff = burial_cutoff,
                stability_classication_x_cutoff = stability_classication_experimental_cutoff,
                stability_classication_y_cutoff = stability_classication_predicted_cutoff,
                use_existing_benchmark_data = False,
                recreate_graphs = False,
                misc_dataframe_attributes = top_level_dataframe_attributes,
            )

        if not(use_existing_benchmark_data and hdf_store_blob):
            hdf_store_blob = benchmark_run.create_dataframe(pdb_data = self.get_prediction_set_pdb_chain_details(prediction_set_id, restrict_to_pdbs = restrict_to_pdbs))
            d = dict(
                PredictionSet                           = prediction_set_id,
                DataFrameType                           = dataframe_type,
                ContainsExperimentalData                = experimental_data_exists,
                ScoreMethodID                           = score_method_id,
                UseSingleReportedValue                  = use_single_reported_value,
                TopX                                    = take_lowest,
                BurialCutoff                            = burial_cutoff,
                StabilityClassicationExperimentalCutoff = stability_classication_experimental_cutoff,
                StabilityClassicationPredictedCutoff    = stability_classication_predicted_cutoff,
                IncludesDerivedMutations                = include_derived_mutations,
                DDGAnalysisType                         = ddg_analysis_type,
                SeriesName                              = prediction_set_series_name,
                SeriesColor                             = prediction_set_color,
                SeriesAlpha                             = prediction_set_alpha,
                Description                             = prediction_set_description,
                Credit                                  = prediction_set_credit,
                DDGAnalysisTypeDescription              = benchmark_run.ddg_analysis_type_description,
                PandasHDFStore                          = hdf_store_blob,
            )
            self.DDG_db.execute('''DELETE FROM AnalysisDataFrame WHERE PredictionSet=%s AND DataFrameType=%s AND ContainsExperimentalData=%s AND ScoreMethodID=%s AND UseSingleReportedValue=%s AND TopX=%s AND
                                    BurialCutoff=%s AND StabilityClassicationExperimentalCutoff=%s AND StabilityClassicationPredictedCutoff=%s AND
                                    IncludesDerivedMutations=%s AND DDGAnalysisType=%s''',
                                    parameters = (prediction_set_id, dataframe_type, experimental_data_exists, score_method_id, use_single_reported_value, take_lowest,
                                                  burial_cutoff, stability_classication_experimental_cutoff, stability_classication_predicted_cutoff,
                                                  include_derived_mutations, ddg_analysis_type
                                    ))
            self.DDG_db.insertDictIfNew('AnalysisDataFrame', d, ['PredictionSet', 'DataFrameType', 'ContainsExperimentalData', 'ScoreMethodID', 'UseSingleReportedValue', 'TopX', 'BurialCutoff',
                                                                 'StabilityClassicationExperimentalCutoff', 'StabilityClassicationPredictedCutoff',
                                                                 'IncludesDerivedMutations', 'DDGAnalysisType'])
        else:
            benchmark_run.read_dataframe_from_content(hdf_store_blob)

        return benchmark_run

        # if use_existing_benchmark_data and dataframe exists: return dataframe
        # else retrieve all of the Score records from the database
        #    if a record does not exist:
        #        if root_directory then call extract_data_for_case to create an analysis dataframe and store it in the database
        #    store the number of complete Score records as a column in the dataframe (to indicate whether analysis is being performed on a full set of data)
        #
        # For Shane: this extracts the dataset_description and dataset_cases data that DDGBenchmarkManager currently takes care of in the capture.
        # The analysis_data variable of DDGBenchmarkManager should be compiled via queries calls to the Prediction*StructureScore table.


    @analysis_api
    def analyze(self, prediction_set_ids,
            prediction_set_series_names = {}, prediction_set_descriptions = {}, prediction_set_credits = {}, prediction_set_colors = {}, prediction_set_alphas = {},
            use_published_data = False,
            use_existing_benchmark_data = True, recreate_graphs = False,
            include_derived_mutations = False,
            expectn = 50,
            use_single_reported_value = False,
            take_lowest = 3,
            burial_cutoff = 0.25,
            stability_classication_experimental_cutoff = 1.0,
            stability_classication_predicted_cutoff = 1.0,
            output_directory = None,
            generate_plots = True,
            report_analysis = True,
            silent = False,
            root_directory = None
            ):
        '''Runs the analyses for the specified PredictionSets and cross-analyzes the sets against each other if appropriate.

           * Analysis setup arguments *

           prediction_set_ids is a list of PredictionSet IDs. Each PredictionSet will be analyzed separately and appropriate
           pairs will be cross-analyzed.
           prediction_set_series_names, prediction_set_descriptions, and prediction_set_credits are mappings from PredictionSet IDs
           to series names (in plots), descriptions, and credits respectively. These details are stored in PredictionSet so
           they are optional arguments. If passed, these mappings will override the PredictionSet values in the database
           which allows the user to customize the analysis reports. Likewise, prediction_set_colors and prediction_set_alphas
           are mappings to series colors and transparency values for use in the plots.

           use_published_data. todo: implement later. This should include any published data e.g. the Kellogg et al. data for protein stability.
           use_existing_benchmark_data and recreate_graphs are data creation arguments i.e. "should we use existing data or create it from scratch?"
           include_derived_mutations is used to filter out dataset cases with derived mutations.
           expectn declares how many predictions we expect to see per dataset case. If the actual number is less than expectn
           then a warning will be included in the analysis.

           * Dataframe arguments *

           use_single_reported_value is specific to ddg_monomer. If this is True then the DDG value reported by the application is used and take_lowest is ignored. This is inadvisable - take_lowest = 3 is a better default.
           take_lowest AKA Top_X. Specifies how many of the best-scoring groups of structures to consider when calculating the predicted DDG value.
           burial_cutoff defines what should be considered buried (DSSPExposure field). Values around 1.0 are fully exposed, values of 0.0 are fully buried. For technical reasons, the DSSP value can exceed 1.0 but usually not by much.
           stability_classication_experimental_cutoff AKA x_cutoff. This defines the neutral mutation range for experimental values in kcal/mol i.e. values between -1.0 and 1.0 kcal/mol are considered neutral by default.
           stability_classication_predicted_cutoff AKA y_cutoff. This defines the neutral mutation range for predicted values in energy units.

           * Reporting arguments *

           output_directory : The directory in which to save plots and reports.
           generate_plots   : if plots are not needed, setting this to False can shorten the analysis time.
           report_analysis  : Whether or not to print analysis to stdout.
           silent = False   : Whether or not anything should be printed to stdout (True is useful for webserver interaction).
        '''

        raise Exception('Abstract method. This needs to be overridden by a subclass.')

        # colors, alpha, and default series name and descriptions are taken from PredictionSet records
        # The order (if p1 before p2 then p1 will be on the X-axis in comparative plots) in comparative analysis plots is determined by the order in PredictionSets
        assert(take_lowest > 0 and (int(take_lowest) == take_lowest))
        assert(0 <= burial_cutoff <= 2.0)
        assert(stability_classication_experimental_cutoff > 0)
        assert(stability_classication_predicted_cutoff > 0)
        # assert PredictionSet for PredictionSet in PredictionSets is in the database

        # calls get_analysis_dataframe(options) over all PredictionSets
        # if output_directory is set, save files
        # think about how to handle this in-memory. Maybe return a dict like:
            #"run_analyis" -> benchmark_name -> {analysis_type -> object}
            #"comparative_analysis" -> (benchmark_name_1, benchmark_name_2) -> {analysis_type -> object}
        # comparative analysis
        #   only compare dataframes with the exact same points
        #   allow cutoffs, take_lowest to differ but report if they do so


    @analysis_api
    def determine_best_pair(self, prediction_id, score_method_id = 1):
        '''This returns the best wildtype/mutant pair for a prediction given a scoring method. NOTE: Consider generalising this to the n best pairs.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @analysis_api
    def create_abacus_graph_for_a_single_structure(self, PredictionSet, scoring_method, scoring_type, graph_title = None, PredictionIDs = None, graph_filename = None, cached_results = None, num_datapoints = 0):
        '''This function creates an abacus graph for one PDB file. It is useful when scanning all mutations at all positions
           on small proteins e.g. ubiquitin to show which mutations at which positions are likely to improve the stability or
           binding affinity.
           The num_datapoints variable is mainly for debugging - I was tuning the resolution/DPI to fit the number of datapoints.'''
        raise Exception('This should work or nearly work. Test it again when we have real data. Does it assume single point mutations?')

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
            #print(json.loads(r['Scores'])['data'][scoring_method][scoring_type]['ddG'], r['FlattenedMutations'])
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


    ################################################################################################
    ## Application layer
    ## These functions combine the database and prediction data with useful klab
    ################################################################################################


    #== PyMOL API ===========================================================


    @app_pymol
    def create_pymol_session_in_memory(self, prediction_id, task_number, pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol'):
        '''Returns (in memory) a PyMOL session for a pair of structures.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @app_pymol
    def write_pymol_session(self, prediction_id, task_number, output_filepath, pymol_executable = '/var/www/tg2/tg2env/designdb/pymol/pymol/pymol'):
        '''Writes the PyMOL session for a pair of structures to disk.'''
        PSE_file_contents = self.create_pymol_session_in_memory(prediction_id, task_number, pymol_executable = pymol_executable)
        write_file(output_filepath, PSE_file_contents, 'wb')


    ################################################################################################
    ## Private API layer
    ## These are helper functions used internally by the class but which are not intended for export
    ################################################################################################


    ###########################################################################################
    ## Subclass layer
    ##
    ## These functions need to be implemented by subclasses
    ###########################################################################################


    def _get_prediction_table(self): return None
    def _get_prediction_structure_scores_table(self): return None
    def _get_prediction_id_field(self): return self._get_prediction_table() + 'ID'
    def _get_prediction_type(self): return None
    def _get_prediction_dataset_type(self): return None
    def _get_prediction_type_description(self): return None
    def _get_user_dataset_experiment_table(self): return None
    def _get_user_dataset_experiment_tag_table(self): return None
    def _get_allowed_score_types(self): return None


    ###########################################################################################
    ## Assertion layer
    ##
    ## These functions check pre- and post-conditions
    ###########################################################################################


    def _check_prediction(self, prediction_id, prediction_set):
        '''Sanity check: Asserts that a Prediction belongs in the expected PredictionSet.'''
        prediction_table = self._get_prediction_table()
        if not self.DDG_db.execute_select('SELECT * FROM {0} WHERE ID=%s AND PredictionSet=%s'.format(prediction_table), parameters=(prediction_id, prediction_set)):
            raise Exception('{0} record #{1} does not belong to PredictionSet {2}.'.format(prediction_table, prediction_id, prediction_set))


    def _check_scores_for_main_fields(self, scores, prediction_id):
        '''Sanity check: Asserts that the identifying fields for the scores make sense for this interface.'''
        prediction_id_field = self._get_prediction_id_field()
        score_method_details = self.get_score_method_details()
        allowed_score_types = self._get_allowed_score_types()
        int_type = type(1)
        for score in scores:
            assert(prediction_id_field in score and score[prediction_id_field] == prediction_id)
            assert('ScoreMethodID' in score and score['ScoreMethodID'] in score_method_details)
            assert('ScoreType' in score and score['ScoreType'] in allowed_score_types)
            assert('StructureID' in score and type(score['StructureID']) == int_type)


    def _check_score_fields(self, scores):
        '''Sanity check: Asserts that the fields for the scores are represented in the database table.'''
        fieldnames = set([f for f in self.DDG_db.FieldNames.__dict__[self._get_prediction_structure_scores_table()].__dict__.keys() if not(f.startswith('_'))])
        for score in scores:
            score_keys = score.keys()
            if sorted(fieldnames.intersection(score_keys)) != sorted(score_keys):
                print score_keys
                print fieldnames
                raise Exception('These score table fieldnames were not recognized: %s.'.format(', '.join(sorted(set(score_keys).difference(fieldnames)))))


    ###########################################################################################
    ## File management layer
    ##
    ## This part of the API is responsible for file content abstraction
    ###########################################################################################


    def _add_file_content(self, content, db_cursor = None, rm_trailing_line_whitespace = False, forced_mime_type = None):
        '''Takes file content (and an option to remove trailing whitespace from lines e.g. to normalize PDB files), adds
           a new record if necessary, and returns the associated FileContent.ID value.'''

        if rm_trailing_line_whitespace:
            content = remove_trailing_line_whitespace(content)

        # Check to see whether the file has been uploaded before
        hexdigest = get_hexdigest(content)
        existing_filecontent_id = self.get_file_id(content, db_cursor = db_cursor, hexdigest = hexdigest)

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

            if db_cursor:
                sql, params, record_exists = self.DDG_db.create_insert_dict_string('FileContent', d, ['Content'])
                db_cursor.execute(sql, params)
            else:
                self.DDG_db.insertDictIfNew('FileContent', d, ['Content'])
            existing_filecontent_id = self.get_file_id(content, db_cursor = db_cursor, hexdigest = hexdigest)
            assert(existing_filecontent_id != None)
        return existing_filecontent_id


    ###########################################################################################
    ## Prediction layer
    ##
    ## This part of the API is responsible for inserting prediction jobs in the database via
    ## the trickle-down proteomics paradigm.
    ###########################################################################################


    #== Job creation/management API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via the
    # trickle-down proteomics paradigm.


    # PredictionSet interface


    def _assert_prediction_set_is_correct_type(self, PredictionSetID):
        '''Returns the list of Prediction IDs associated with the PredictionSet.'''
        assert(self._get_prediction_type() and self._get_prediction_type_description())
        if (self.get_prediction_set_details(PredictionSetID) or {}).get(self._get_prediction_type()) != 1:
            raise Exception('This PredictionSet either does not exist or else contains no %s predictions.' % self._get_prediction_type_description())


    def _set_prediction_set_status(self, PredictionSetID, status):
        '''Sets the Status of a PredictionSet.'''
        assert(status == 'halted' or status == 'active')
        assert(self.get_prediction_set_details(PredictionSetID))
        self.DDG_db.execute('UPDATE PredictionSet SET Status=%s WHERE ID=%s', parameters=(status, PredictionSetID))


    # Prediction setup interface


    def _add_prediction_file(self, prediction_id, file_content, filename, filetype, filerole, stage, db_cursor = None, rm_trailing_line_whitespace = False, forced_mime_type = None):
        '''This function adds file content to the database and then creates a record associating that content with a prediction.
           If db_cursor is passed then we call that directly. This is crucial for transactions as many of the database functions
           create new cursors which commit changes so transactions do not work properly.'''
        prediction_table = self._get_prediction_table()

        # Add the file contents to the database
        if filetype == 'PDB':
            forced_mime_type = forced_mime_type or 'chemical/x-pdb'

        file_content_id = self._add_file_content(file_content, db_cursor = db_cursor, rm_trailing_line_whitespace = rm_trailing_line_whitespace, forced_mime_type = forced_mime_type)

        # Link the file contents to the prediction
        d = dict(
            FileContentID = file_content_id,
            Filename = filename,
            Filetype = filetype,
            FileRole = filerole,
            Stage = stage,
        )
        if prediction_table == 'Prediction':
            d['PredictionID'] = prediction_id
            if db_cursor:
                # Note: When less tired, add select statement here to see if info already in database
                sql, params, record_exists = self.DDG_db.create_insert_dict_string('PredictionFile', d, ['PredictionID', 'FileRole', 'Stage'])
                db_cursor.execute(sql, params)
            else:
                self.DDG_db.insertDictIfNew('PredictionFile', d, ['PredictionID', 'FileRole', 'Stage'])
        elif prediction_table == 'PredictionPPI':
            d['PredictionPPIID'] = prediction_id
            if db_cursor:
                sql, params, record_exists = self.DDG_db.create_insert_dict_string('PredictionPPIFile', d, ['PredictionPPIID', 'FileRole', 'Stage'])
                db_cursor.execute(sql, params)
            else:
                self.DDG_db.insertDictIfNew('PredictionPPIFile', d, ['PredictionPPIID', 'FileRole', 'Stage'])
        else:
            raise('Invalid table "%s" passed.' % prediction_table)


    def _add_residue_map_to_prediction(self, prediction_id, residue_mapping):
        assert(type(residue_mapping) == type(self.__dict__))
        json_content = json.dumps(residue_mapping, sort_keys=True) # sorting helps to quotient the file content space over identical data
        self._add_prediction_file(prediction_id, json_content, 'residue_mapping.json', 'RosettaPDBMapping', 'Rosetta<->PDB residue mapping', 'Input', forced_mime_type = "application/json")


    def _strip_pdb(self, pdb_file_id, chains):
        raise Exception('assert that chains exist in PDBChain table. reads PDB content from the database. call PDB class functions to strip to chains.')


    def _add_stripped_pdb_to_prediction(self, prediction_id):
        pdb_file_id, chains = self.get_pdb_chains_for_prediction(prediction_id)
        pdb_content = self._strip_pdb(pdb_file_id, chains)
        filename = '%s_%s' % (pdb_file_id, ''.join(sorted(chains)))
        self._add_prediction_file(prediction_id, pdb_content, filename, 'PDB', 'StrippedPDB', 'Input', rm_trailing_line_whitespace = True, forced_mime_type = 'chemical/x-pdb')


    def _add_resfile_to_prediction(self, prediction_id):
        resfile_content = self.create_resfile(prediction_id)
        self._add_prediction_file(prediction_id, resfile_content, 'mutations.resfile', 'Resfile', 'Resfile', 'Input', rm_trailing_line_whitespace = True)


    def _add_mutfile_to_prediction(self, prediction_id):
        mutfile_content = self.create_mutfile(prediction_id)
        self._add_prediction_file(prediction_id, mutfile_content, 'mutations.mutfile', 'Mutfile', 'Mutfile', 'Input', rm_trailing_line_whitespace = True)


    def _create_resfile_from_pdb_mutations(self, stripped_pdb, pdb_mutations):
        '''This function takes a PDB object to be used in a DDG job (i.e. usually stripped to certain chains but with the
           original PDB numbering) and a list of mutations using the original PDB numbering. Resfiles use PDB numbering
           so no mapping needs to be done.'''

        if not pdb_mutations:
            raise Exception("There needs to be at least one mutation.")

        try:
            resfile = []
            for mutation in pdb_mutations:
                # Check that the expected wildtype exists in the PDB
                stripped_pdb.assert_wildtype_matches(mutation)
                chain, resid, mt = mutation.Chain, mutation.ResidueID.strip(), mutation.MutantAA
                #resfile.append("%(resid)s %(chain)s PIKAA %(mt)s" % vars())
                resfile.append("%(resid)s %(chain)s PIKAA %(mt)s" % vars())
            assert(resfile)
            return '\n'.join(["NATAA", "start"] + resfile)
        except:
            raise Exception("An error occurred creating a resfile for the ddG job.")


    def _create_mutfile_from_pdb_mutations(self, stripped_pdb, pdb_mutations):
        '''This function takes a PDB object to be used in a DDG job (i.e. usually stripped to certain chains but with the
           original PDB numbering)) and a list of mutations using the original PDB numbering. Since mutfiles use Rosetta
           numbering, we need to map the residue IDs from PDB numbering to Rosetta numbering.'''

        if not pdb_mutations:
            raise Exception("There needs to be at least one mutation.")

        try:
            # Map the mutations from PDB numbering to Rosetta numbering
            rosetta_mutations = stripped_pdb.map_pdb_residues_to_rosetta_residues(pdb_mutations)
            assert(len(rosetta_mutations) == len(pdb_mutations))

            mutfile = []
            for x in len(pdb_mutations):
                pdb_mutation = pdb_mutations[x]
                rosetta_mutation = pdb_mutations[x]

                # Check that the expected wildtype exists in the PDB
                stripped_pdb.assert_wildtype_matches(pdb_mutation)

                wt, resid, mt = rosetta_mutation.WildTypeAA, rosetta_mutation.ResidueID, rosetta_mutation.MutantAA
                mutfile.append("%(wt)s %(resid)s %(mt)s" % vars())

            assert(mutfile)
            return '\n'.join(["total %d" % len(rosetta_mutations), "%d" % len(rosetta_mutations)] + mutfile)
        except:
            raise Exception("An error occurred creating a mutfile for the ddG job.")
