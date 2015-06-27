#!/usr/bin/python2.4
# encoding: utf-8
"""
db_api.py
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

from io import BytesIO
from api_layers import *

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


    def __init__(self, passwd = None, username = 'kortemmelab', rosetta_scripts_path = None, rosetta_database_path = None):
        if passwd:
            passwd = passwd.strip()
        self.DDG_db = ddgdbapi.ddGDatabase(passwd = passwd, username = username)
        self.DDG_db_utf = ddgdbapi.ddGDatabase(passwd = passwd, username = username, use_utf = True)
        self.prediction_data_path = None
        self.rosetta_scripts_path = rosetta_scripts_path
        self.rosetta_database_path = rosetta_database_path


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
        '''NOTE: This function should be generalized and moved into the tools repository.
           This is a simple function wrapper around create_abacus_graph which writes the graph to file.'''
        byte_stream = self.create_abacus_graph(graph_title, labels, data)
        write_file(graph_filename, byte_stream.getvalue(), 'wb')


    @alien
    def create_abacus_graph(self, graph_title, labels, data):
        '''NOTE: This function should be generalized and moved into the tools repository.
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
    def add_PDB_to_database(self, filepath = None, pdbID = None, contains_membrane_protein = None, protein = None, file_source = None, UniProtAC = None, UniProtID = None, testonly = False, force = False, techniques = None):
        '''NOTE: This API is used to create and analysis predictions or retrieve information from the database.
                 This function adds new raw data to the database and does not seem to belong here. It should be moved into
                 an admin API instead.
           This function adds imports a PDB into the database, creating the associated molecule, chain and residue etc. records.'''

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
            dbp = ddgdbapi.PDBStructure(self.DDG_db, rootname, contains_membrane_protein = contains_membrane_protein, protein = protein, file_source = file_source, filepath = filepath, UniProtAC = UniProtAC, UniProtID = UniProtID, testonly = testonly, techniques = techniques)
            #Structure.getPDBContents(self.DDG_db)
            results = self.DDG_db.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters = (rootname,))
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


    @alien
    def get_flattened_prediction_results(self, PredictionSet):
        '''This is defined here as an API function but should be defined as a stored procedure.'''
        return self.DDG_db.execute_select('''
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
    def add_pdb_file(self, filepath, pdb_id): raise Exception('This function has been deprecated. Use add_PDB_to_database instead.')

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


    @informational_file
    def get_file_id(self, content, hexdigest = None):
        '''Searches the database to see whether the FileContent already exists. The search uses the digest and filesize as
           heuristics to speed up the search. If a file has the same hex digest and file size then we do a straight comparison
           of the contents.
           If the FileContent exists, the value of the ID field is returned else None is returned.
           '''
        existing_filecontent_id = None
        hexdigest = hexdigest or get_hexdigest(content)
        filesize = len(content)
        for r in self.DDG_db.execute_select('SELECT * FROM FileContent WHERE MD5HexDigest=%s AND Filesize=%s', parameters=(hexdigest, filesize)):
            if r['Content'] == content:
                assert(existing_filecontent_id == None) # content uniqueness check
                existing_filecontent_id = r['ID']
        return existing_filecontent_id


    @informational_pdb
    def get_pdb_chains_for_prediction(self, prediction_id):
        '''Returns the PDB file ID and a list of chains for the prediction.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_pdb
    def get_chains_for_mutatagenesis(self, mutagenesis_id, pdb_file_id, pdb_set_number, complex_id = None):
        '''Returns the PDB chains used in the mutagenesis.
           Note: At present, monomeric data e.g. protein stability does not have the notion of complex in our database
           but this abstraction is planned so that multiple choices of PDB file and chain can be easily represented.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @informational_pdb
    def get_pdb_details(self, pdb_ids, cached_pdb_details = None):
        '''Returns the details stored in the database about the PDB files associated with pdb_ids e.g. chains, resolution,
           technique used to determine the structure etc.'''
        pdbs = {}
        cached_pdb_ids = []
        if cached_pdb_details:
            cached_pdb_ids = cached_pdb_details.keys()
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


    @informational_job
    def get_prediction_set_details(self, PredictionSetID):
        '''Returns the PredictionSet record from the database.'''
        results = self.DDG_db.execute_select('SELECT * FROM PredictionSet WHERE ID=%s', parameters=(PredictionSetID,))
        if len(results) == 1:
            return results[0]
        return None


    @informational_job
    def get_prediction_ids(self, PredictionSetID):
        '''Returns the list of Prediction IDs associated with the PredictionSet.'''
        self._assert_prediction_set_is_correct_type(PredictionSetID)
        qry = 'SELECT ID FROM %s WHERE PredictionSet=%%s ORDER BY ID' % self._get_prediction_table()
        return [r['ID'] for r in self.DDG_db.execute_select(qry, parameters=(PredictionSetID,))]


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


    ###########################################################################################
    ## Prediction creation/management layer
    ##
    ###########################################################################################


    #== Job creation API ===========================================================
    #
    # This part of the API is responsible for inserting prediction jobs in the database via
    # the trickle-down proteomics paradigm.


    @job_creator
    def add_prediction_set(self, prediction_set_id, halted = True, priority = 5, batch_size = 40, allow_existing_prediction_set = False, contains_protein_stability_predictions = True, contains_binding_affinity_predictions = False):
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
        if len(existing_record) > 0:
            assert(len(existing_record) == 1)
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
        )
        self.DDG_db.insertDictIfNew("PredictionSet", d, ['ID'])
        return True


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
    def clone_prediction_run(self, existing_prediction_set, new_prediction_set, *args, **kwargs):
        '''add_prediction_run sets up a full run of dataset predictions but is slow as it needs to perform a lot of
           calculations and parsing. If you want to test the same dataset with slightly different parameters (e.g. a
           different protocol) then these calculations can be reused which reduces the overhead considerably.
           clone_prediction_run was written with this in mind. It copies the list of predictions and their setup (input
           files etc.) from an existing prediction set to an empty prediction set.'''
        raise Exception('not implemented yet')


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
    def get_job(self, prediction_set):
        '''Returns None if no queued jobs exist or if the PredictionSet is halted otherwise return details necessary to run the job.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')
        #returns None if no queued jobs exist or if the PredictionSet is halted otherwise return details necessary to run the job


    @job_execution
    def start_job(self, prediction_id, prediction_set):
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


    @job_completion
    def fail_job(self, prediction_id, prediction_set, maxvmem, ddgtime):
        '''Sets the job status to "failed". prediction_set must be passed and is used as a sanity check.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')
        # sets a job to 'failed'.


    @job_completion
    def parse_prediction_scores(self, *args, **kwargs):
        '''Returns a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass. Returns a dict suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.')


    @job_completion
    def store_scores(self, scores, prediction_set, prediction_id):
        '''Stores a list of dicts suitable for database storage e.g. PredictionStructureScore or PredictionPPIStructureScore records.'''
        raise Exception('Abstract method. This needs to be overridden by a subclass.')


    @job_completion
    def complete_job(self, prediction_id, prediction_set, scores, maxvmem, ddgtime):
        '''Sets a job to 'completed' and stores scores. prediction_set must be passed and is used as a sanity check.'''
        raise Exception('This function needs to be implemented by subclasses of the API.')


    @job_results
    def add_ddg_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores):
        '''Add the DDG scores for a given score method for a prediction.'''
        raise Exception('not implemented yet. This should only need to be implemented for the base class')


    @job_results
    def _add_structure_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores, is_wildtype):
        '''Add the scores for a given score method for one structure of a prediction.'''
        raise Exception('not implemented yet. This should only need to be implemented for the base class')
        #if ScoreMethodID is None, raise an exception but report the score method id for the default score method
        #assert(scores fields match db fields)


    @job_results
    def add_wildtype_structure_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores):
        '''Add the scores for a given score method for one wildtype structure of a prediction.'''
        return self._add_structure_score(prediction_set, prediction_id, structure_id, ScoreMethodID, scores, True)
        raise Exception('not implemented yet. This should only need to be implemented for the base class')


    @job_results
    def add_mutant_structure_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores):
        '''Add the scores for a given score method for one mutant structure of a prediction.'''
        return self._add_structure_score(prediction_set, prediction_id, structure_id, ScoreMethodID, scores, False)
        raise Exception('not implemented yet. This should only need to be implemented for the base class')


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


    ################################################################################################
    ## Application layer
    ## These functions combine the database and prediction data with useful tools
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
    def _get_prediction_type(self): return None
    def _get_prediction_dataset_type(self): return None
    def _get_prediction_type_description(self): return None
    def _get_user_dataset_experiment_table(self): return None
    def _get_user_dataset_experiment_tag_table(self): return None


    ###########################################################################################
    ## File management layer
    ##
    ## This part of the API is responsible for file content abstraction
    ###########################################################################################


    def _add_file_content(self, content, rm_trailing_line_whitespace = False, forced_mime_type = None):
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
            self.DDG_db.insertDictIfNew('FileContent', d, ['Content'])
            existing_filecontent_id = self.get_file_id(content, hexdigest = hexdigest)
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
            raise Exception('This PredictionSet either does not exist or else contains no %s predictions.' % self.self._get_prediction_type_description())


    def _set_prediction_set_status(self, PredictionSetID, status):
        '''Sets the Status of a PredictionSet.'''
        assert(status == 'halted' or status == 'active')
        assert(self.get_prediction_set_details(PredictionSetID))
        self.DDG_db.execute('UPDATE PredictionSet SET Status=%s WHERE ID=%s', parameters=(status, PredictionSetID))


    # Prediction setup interface


    def _add_prediction_file(self, prediction_id, file_content, filename, filetype, filerole, stage, rm_trailing_line_whitespace = False, forced_mime_type = None):
        '''This function adds file content to the database and then creates a record associating that content with a prediction.'''
        prediction_table = self._get_prediction_table()

        # Add the file contents to the database
        if filetype == 'PDB':
            forced_mime_type = forced_mime_type or 'chemical/x-pdb'

        file_content_id = self._add_file_content(file_content, rm_trailing_line_whitespace = rm_trailing_line_whitespace, forced_mime_type = forced_mime_type)

        # Link the file contents to the prediction
        d = dict(
            FileContentID = file_content_id,
            Filename = filename,
            Filetype = filetype,
            FileRole = filerole,
            Stage = stage,
        )
        if prediction_table == 'Prediction':
            d['PredictionID'] = prediction_id,
            self.DDG_db.insertDictIfNew('PredictionFile', d, ['PredictionID', 'FileRole', 'Stage'])
        elif prediction_table == 'PredictionPPI':
            d['PredictionPPIID'] = prediction_id,
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


