#!/usr/bin/python2.4
# encoding: utf-8
"""
ppi_api.py
High-level functions for interacting with the protein-protein interaction sections of the ddG database.
The API for the database is overgrown and all over the place. It should be refactored properly but for now I am
separating out the protein-protein functionality so that we can at least design that from a clean slate.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

from tools.fs.fsio import read_file
from tools import colortext
import ddgdbapi
import pprint
from dbapi import ddG, jobcreator, inputfiles, analysisfn, pymolapi, deprecated
import inspect
import functools

def bind_object_function(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs): return fn(*args, **kwargs)
    return wrapper


class BindingAffinityDDGUserInterface(object):
    '''This is the class that should be used to interface with the database. It hides functions that should only be called
       within this other API functions.

       The class contains a private copy of the internal API and wraps the public functions of that API so that the
       functions of BindingAffinityDDGUserInterface contain only the public functions of the internal API. Private functions
       are denoted as such by a leading underscore in the function name.
       '''

    def __init__(self, passwd = None, username = 'kortemmelab'):

        self._ddg_interface = BindingAffinityDDGInterface(passwd = passwd, username = username)
        self._api_functions = []
        self._api_function_args = {}
        self.DDG_db = self._ddg_interface.DDG_db
        self.DDG_db_utf = self._ddg_interface.DDG_db_utf

        for m in inspect.getmembers(BindingAffinityDDGInterface, predicate=inspect.ismethod):
            if m[0][0] != '_':
                fn_name = m[0]
                self._api_function_args[fn_name] = getattr(self._ddg_interface, fn_name).func_code.co_varnames
                self._api_functions.append(fn_name)
                self.__dict__[fn_name] = bind_object_function(getattr(self._ddg_interface, fn_name))


    def help(self):
        helpstr = []
        helpstr.append(colortext.mcyan('\n*****\n*** %s API\n*****\n' % self.__class__.__name__))

        doc_strings = {}
        for fn_name in sorted(self._api_functions):
            fn = self.__dict__[fn_name]
            function_class = 'Miscellanous'
            try:
                function_class = fn._helptype
            except: pass
            doc_strings[function_class] = doc_strings.get(function_class, {})
            doc_strings[function_class][fn_name] = fn.__doc__

        for function_class, fn_names in sorted(doc_strings.iteritems()):
            helpstr.append(colortext.mlightpurple('* %s *\n' % function_class))
            for fn_name, docstr in sorted(fn_names.iteritems()):
                helpstr.append(colortext.mgreen('  %s(%s)' % (fn_name, ', '.join(self._api_function_args[fn_name]))))
                if docstr:
                    helpstr.append(colortext.myellow('    %s' % ('\n    '.join([s.strip() for s in docstr.split('\n') if s.strip()]))))
                else:
                    helpstr.append(colortext.mred('    <not documented>'))
                helpstr.append('')
        return '\n'.join(helpstr)


class BindingAffinityDDGInterface(ddG):
    '''This is the internal API class that should be NOT used to interface with the database.'''


    def __init__(self, passwd = None, username = 'kortemmelab'):
        super(BindingAffinityDDGInterface, self).__init__(passwd = passwd, username = username)

    def get_prediction_table(self):
        return 'PredictionPPI'


    ##### Public API: Rosetta-related functions


    @inputfiles
    def create_resfile(self, prediction_id):
        raise Exception('This needs to be implemented.')


    @inputfiles
    def create_mutfile(self, prediction_id):
        raise Exception('This needs to be implemented.')



    ##### Public API: PDB-related functions



    def get_pdb_chains_for_prediction(self, prediction_id):
        '''<!Job insertion>Returns the PDB file ID and a list of chains for the prediction.'''
        raise Exception('This needs to be implemented.')



    ##### Public API: PredictionSet functions



    def add_prediction_set(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40, allow_existing_prediction_set = False):
        raise Exception('need to implement this')
        #call super function with contains_protein_stability_predictions = False, contains_binding_affinity_predictions = True


    ##### Private API: Job insertion helper functions



    def _charge_prediction_set_by_residue_count(self, PredictionSet):
        '''This function assigns a cost for a prediction equal to the number of residues in the chains.'''
        raise Exception('This function needs to be rewritten.')
        from tools.bio.rcsb import parseFASTAs

        DDG_db = self.DDG_db
        predictions = DDG_db.execute_select("SELECT ID, ExperimentID FROM Prediction WHERE PredictionSet=%s", parameters=(PredictionSet,))

        PDB_chain_lengths ={}
        for prediction in predictions:
            chain_records = DDG_db.execute_select('SELECT PDBFileID, Chain FROM Experiment INNER JOIN ExperimentChain ON ExperimentID=Experiment.ID WHERE ExperimentID=%s', parameters=(prediction['ExperimentID']))
            num_residues = 0
            for chain_record in chain_records:
                key = (chain_record['PDBFileID'], chain_record['Chain'])

                if PDB_chain_lengths.get(key) == None:
                    fasta = DDG_db.execute_select("SELECT FASTA FROM PDBFile WHERE ID=%s", parameters = (chain_record['PDBFileID'],))
                    assert(len(fasta) == 1)
                    fasta = fasta[0]['FASTA']
                    f = parseFASTAs(fasta)
                    PDB_chain_lengths[key] = len(f[chain_record['PDBFileID']][chain_record['Chain']])
                chain_length = PDB_chain_lengths[key]
                num_residues += chain_length

            print("UPDATE Prediction SET Cost=%0.2f WHERE ID=%d" % (num_residues, prediction['ID']))

            predictions = DDG_db.execute("UPDATE Prediction SET Cost=%s WHERE ID=%s", parameters=(num_residues, prediction['ID'],))



    ##### Public API: Job insertion functions


    @jobcreator
    def add_prediction_set_jobs(self, userdatasetTextID, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = False, Description = {}, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False):
        raise Exception('This function needs to be rewritten.')

        assert(self.DDG_db.execute_select("SELECT ID FROM PredictionSet WHERE ID=%s", parameters=(PredictionSet,)))

        #results = self.DDG_db.execute_select("SELECT * FROM UserDataSet WHERE TextID=%s", parameters=(userdatasetTextID,))
        results = self.DDG_db.execute_select("SELECT UserDataSetExperiment.* FROM UserDataSetExperiment INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%s", parameters=(userdatasetTextID,))
        if not results:
            return False

        if not(quiet):
            colortext.message("Creating predictions for UserDataSet %s using protocol %s" % (userdatasetTextID, ProtocolID))
            colortext.message("%d records found in the UserDataSet" % len(results))

        count = 0
        showprogress = not(quiet) and len(results) > 300
        if showprogress:
            print("|" + ("*" * (int(len(results)/100)-2)) + "|")
        for r in results:

            existing_results = self.DDG_db.execute_select("SELECT * FROM Prediction WHERE PredictionSet=%s AND UserDataSetExperimentID=%s", parameters=(PredictionSet, r["ID"]))
            if len(existing_results) > 0:
                #colortext.warning('There already exist records for this UserDataSetExperimentID. You probably do not want to proceed. Skipping this entry.')
                continue

            PredictionID = self.addPrediction(r["ExperimentID"], r["ID"], PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = r["PDBFileID"], StoreOutput = StoreOutput, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = testonly)
            count += 1
            if showprogress:
                if count > 100:
                    colortext.write(".", "cyan", flush = True)
                    count = 0
            if shortrun and count > 4:
                break
        print("")
        return(True)


    @jobcreator
    def add_jobs_by_pdb_id(self, pdb_ID, PredictionSet, ProtocolID, status = 'active', priority = 5, KeepHETATMLines = False, strip_other_chains = True):
        ''' This function adds predictions for all Experiments corresponding to pdb_ID to the specified prediction set.
            This is useful for custom runs e.g. when we are using the DDG scheduler for design rather than for benchmarking.
        '''
        raise Exception('This function needs to be rewritten.')
        colortext.printf("\nAdding any mutations for this structure which have not been queued/run in the %s prediction set." % PredictionSet, "lightgreen")

        d = {
            'ID' : PredictionSet,
            'Status' : status,
            'Priority' : priority,
            'BatchSize' : 40,
            'EntryDate' : datetime.datetime.now(),
        }
        self.DDG_db.insertDictIfNew('PredictionSet', d, ['ID'])

        # Update the priority and activity if necessary
        self.DDG_db.execute('UPDATE PredictionSet SET Status=%s, Priority=%s WHERE ID=%s', parameters = (status, priority, PredictionSet))

        # Determine the set of experiments to add
        ExperimentIDs = set([r['ID'] for r in self.DDG_db.execute_select('SELECT ID FROM Experiment WHERE PDBFileID=%s', parameters=(pdb_ID,))])
        ExperimentIDsInPredictionSet = set([r['ExperimentID'] for r in self.DDG_db.execute_select('SELECT ExperimentID FROM Prediction WHERE PredictionSet=%s', parameters=(PredictionSet,))])
        experiment_IDs_to_add = sorted(ExperimentIDs.difference(ExperimentIDsInPredictionSet))

        if experiment_IDs_to_add:
            colortext.printf("\nAdding %d jobs to the prediction set." % len(experiment_IDs_to_add), "lightgreen")
            count = 0
            for experiment_ID in experiment_IDs_to_add:
                colortext.write('.', "lightgreen")
                self.addPrediction(experiment_ID, None, PredictionSet, ProtocolID, KeepHETATMLines, StoreOutput = True, strip_other_chains = strip_other_chains)
                count +=1
        else:
            colortext.printf("\nAll jobs are already in the queue or have been run.", "lightgreen")
        print('')


    @jobcreator
    def add_job(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False, strip_other_chains = True):
        '''This function inserts a prediction into the database.
            The parameters define:
                the experiment we are running the prediction for;
                the name of the set of predictions for later grouping;
                the short description of the Command to be used for prediction;
                whether HETATM lines are to be kept or not.
            We strip the PDB based on the chains used for the experiment and KeepHETATMLines.
            We then add the prediction record, including the stripped PDB and the inverse mapping
            from Rosetta residue numbering to PDB residue numbering.'''
        raise Exception('This function needs to be rewritten.')

        parameters = (experimentID,)
        assert(ReverseMutation == False) # todo: allow this later
        try:
            predictionPDB_ID = None

            sql = "SELECT PDBFileID, Content FROM Experiment INNER JOIN PDBFile WHERE Experiment.PDBFileID=PDBFile.ID AND Experiment.ID=%s"
            results = self.DDG_db.execute_select(sql, parameters = parameters)
            if len(results) != 1:
                raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
            experimentPDB_ID = results[0]["PDBFileID"]
            pdbID = results[0]["PDBFileID"]

            if PDB_ID:
                #sql = "SELECT ID, Content FROM PDBFile WHERE ID=%s"
                results = self.DDG_db.execute_select("SELECT ID, Content FROM PDBFile WHERE ID=%s", parameters=(PDB_ID))
                if len(results) != 1:
                    raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
                predictionPDB_ID = results[0]["ID"]
                pdbID = results[0]["ID"]
            else:
                predictionPDB_ID = experimentPDB_ID

            # Get the related PDB ID and file
            assert(len(results) == 1)
            result = results[0]
            contents = result["Content"]

            pdb = PDB(contents.split("\n"))

            # Check that the mutated positions exist and that the wild-type matches the PDB
            mutations = self.DDG_db.call_select_proc("GetMutations", parameters = parameters)

            # todo: Hack. This should be removed when PDB homologs are dealt with properly.
            mutation_objects = []
            for mutation in mutations:
                if experimentPDB_ID == "1AJ3" and predictionPDB_ID == "1U5P":
                    assert(int(mutation['ResidueID']) < 1000)
                    mutation['ResidueID'] = str(int(mutation['ResidueID']) + 1762)
                mutation_objects.append(Mutation(mutation['WildTypeAA'], mutation['ResidueID'], mutation['MutantAA'], mutation['Chain']))

            #todo: a
            #checkPDBAgainstMutations(pdbID, pdb, mutations)
            pdb.validate_mutations(mutation_objects)

            #for mutation in mutations:
            #    if experimentPDB_ID == "ub_OTU":
            #        mutation['ResidueID'] = str(int(mutation['ResidueID']) + 172)

            # Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
            chains = [result['Chain'] for result in self.DDG_db.call_select_proc("GetChains", parameters = parameters)]
            if strip_other_chains:
                pdb.stripForDDG(chains, KeepHETATMLines, numberOfModels = 1)
            else:
                pdb.stripForDDG(True, KeepHETATMLines, numberOfModels = 1)

            #print('\n'.join(pdb.lines))
            # - Post stripping checks -
            # Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
            # then check again that the mutated positions exist and that the wild-type matches the PDB
            colortext.warning('mutations %s' % (str(mutations)))

            remappedMutations = pdb.remapMutations(mutations, pdbID)

            #resfile = self._createResfile(pdb, remappedMutations)
            mutfile = self._createMutfile(pdb, remappedMutations)

            # Check to make sure that we haven't stripped all the ATOM lines
            if not pdb.GetAllATOMLines():
                raise colortext.Exception("No ATOM lines remain in the stripped PDB file of %s." % pdbID)

            # Check to make sure that CSE and MSE are not present in the PDB
            badresidues = pdb.CheckForPresenceOf(["CSE", "MSE"])
            if badresidues:
                raise colortext.Exception("Found residues [%s] in the stripped PDB file of %s. These should be changed to run this job under Rosetta." % (', '.join(badresidues), pdbID))

            # Turn the lines array back into a valid PDB file
            strippedPDB = '\n'.join(pdb.lines)
        except Exception, e:
            colortext.error("Error in %s, %s: .\n%s" % (experimentID, UserDataSetExperimentID, traceback.format_exc()))
            colortext.warning(str(e))
            return
            colortext.error("\nError: '%s'.\n" % (str(e)))
            colortext.error(traceback.format_exc())
            raise colortext.Exception("An exception occurred retrieving the experimental data for Experiment ID #%s." % experimentID)

        #InputFiles["RESFILE"] = resfile
        InputFiles["MUTFILE"] = mutfile

        ExtraParameters = {}
        InputFiles = pickle.dumps(InputFiles)
        Description = pickle.dumps(Description)
        ExtraParameters = pickle.dumps(ExtraParameters)

        PredictionFieldNames = self.DDG_db.FieldNames.Prediction
        params = {
            PredictionFieldNames.ExperimentID		: experimentID,
            PredictionFieldNames.UserDataSetExperimentID : UserDataSetExperimentID,
            PredictionFieldNames.PredictionSet		: PredictionSet,
            PredictionFieldNames.ProtocolID			: ProtocolID,
            PredictionFieldNames.KeptHETATMLines	: KeepHETATMLines,
            PredictionFieldNames.StrippedPDB		: strippedPDB,
            PredictionFieldNames.ResidueMapping		: pickle.dumps(pdb.get_ddGInverseResmap()),
            PredictionFieldNames.InputFiles			: InputFiles,
            PredictionFieldNames.Description		: Description,
            PredictionFieldNames.Status 			: "queued",
            PredictionFieldNames.ExtraParameters	: ExtraParameters,
            PredictionFieldNames.StoreOutput		: StoreOutput,
        }
        if not testonly:
            self.DDG_db.insertDict('Prediction', params)

            # Add cryptID string
            predictionID = self.DDG_db.getLastRowID()
            entryDate = self.DDG_db.execute_select("SELECT EntryDate FROM Prediction WHERE ID=%s", parameters = (predictionID,))[0]["EntryDate"]
            rdmstring = ''.join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16))
            cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
            cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
            entryDate = self.DDG_db.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))
            return predictionID


    @jobcreator
    def clone_prediction_set(self, existing_prediction_set, new_prediction_set):
        raise Exception('not implemented yet')
        #assert(existing_prediction_set exists and has records)
        #assert(new_prediction_set is empty)
        #for each prediction record, add the record and all associated predictionfile records,



    ##### Private API: Job completion functions



    def _add_structure_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores, is_wildtype):
        raise Exception('not implemented yet')
        #if ScoreMethodID is None, raise an exception but report the score method id for the default score method
        #assert(scores fields match db fields)



    ##### Public API: Job completion functions



    def add_ddg_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores):
        raise Exception('not implemented yet')
        #if ScoreMethodID is None, raise an exception but report the score method id for the default score method
        #assert(scores fields match db fields)


    def add_wildtype_structure_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores):
        raise Exception('not implemented yet')


    def add_mutant_structure_score(self, prediction_set, prediction_id, structure_id, ScoreMethodID, scores):
        raise Exception('not implemented yet')




class PPIJobCreator(ddG):
    '''This class is responsible for inserting prediction jobs in the database via the trickle-down proteomics paradigm.'''

    def __init__(self, passwd = None, username = 'kortemmelab'):
        self.ddGdb = ddgdbapi.ddGDatabase(passwd = passwd, username = username)
        self.ddGdb_utf = ddgdbapi.ddGDatabase(passwd = passwd, username = username, use_utf=True)
        self.prediction_data_path = self.ddGDB.execute('SELECT Value FROM _DBCONSTANTS WHERE VariableName="PredictionDataPath"')[0]['Value']


    def create_prediction_set(self, PredictionSetID, halted = True, Priority = 5, BatchSize = 40, allow_existing_prediction_set = False, contains_protein_stability_predictions = False, contains_binding_affinity_predictions = True):
        # Re-use functions from the ddG object
        return self.create_PredictionSet(PredictionSetID, halted = halted, Priority = Priority, BatchSize = BatchSize, allow_existing_prediction_set = allow_existing_prediction_set, contains_protein_stability_predictions = contains_protein_stability_predictions, contains_binding_affinity_predictions = contains_binding_affinity_predictions)


    def add_prediction(self, experimentID, UserDataSetExperimentID, PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = None, StoreOutput = False, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = False, strip_other_chains = True):
        '''This function inserts a prediction into the database.
            The parameters define:
                the experiment we are running the prediction for;
                the name of the set of predictions for later grouping;
                the short description of the Command to be used for prediction;
                whether HETATM lines are to be kept or not.
            We strip the PDB based on the chains used for the experiment and KeepHETATMLines.
            We then add the prediction record, including the stripped PDB and the inverse mapping
            from Rosetta residue numbering to PDB residue numbering.'''

        parameters = (experimentID,)
        assert(ReverseMutation == False) # todo: allow this later
        try:
            predictionPDB_ID = None

            sql = "SELECT PDBFileID, Content FROM Experiment INNER JOIN PDBFile WHERE Experiment.PDBFileID=PDBFile.ID AND Experiment.ID=%s"
            results = self.ddGDB.execute_select(sql, parameters = parameters)
            if len(results) != 1:
                raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
            experimentPDB_ID = results[0]["PDBFileID"]
            pdbID = results[0]["PDBFileID"]

            if PDB_ID:
                #sql = "SELECT ID, Content FROM PDBFile WHERE ID=%s"
                results = self.ddGDB.execute_select("SELECT ID, Content FROM PDBFile WHERE ID=%s", parameters=(PDB_ID))
                if len(results) != 1:
                    raise colortext.Exception("The SQL query '%s' returned %d results where 1 result was expected." % (sql, len(results)))
                predictionPDB_ID = results[0]["ID"]
                pdbID = results[0]["ID"]
            else:
                predictionPDB_ID = experimentPDB_ID

            # Get the related PDB ID and file
            assert(len(results) == 1)
            result = results[0]
            contents = result["Content"]

            pdb = PDB(contents.split("\n"))

            # Check that the mutated positions exist and that the wild-type matches the PDB
            mutations = self.ddGDB.call_select_proc("GetMutations", parameters = parameters)

            # todo: Hack. This should be removed when PDB homologs are dealt with properly.
            mutation_objects = []
            for mutation in mutations:
                if experimentPDB_ID == "1AJ3" and predictionPDB_ID == "1U5P":
                    assert(int(mutation['ResidueID']) < 1000)
                    mutation['ResidueID'] = str(int(mutation['ResidueID']) + 1762)
                mutation_objects.append(Mutation(mutation['WildTypeAA'], mutation['ResidueID'], mutation['MutantAA'], mutation['Chain']))

            #todo: a
            #checkPDBAgainstMutations(pdbID, pdb, mutations)
            pdb.validate_mutations(mutation_objects)

            #for mutation in mutations:
            #    if experimentPDB_ID == "ub_OTU":
            #        mutation['ResidueID'] = str(int(mutation['ResidueID']) + 172)

            # Strip the PDB to the list of chains. This also renumbers residues in the PDB for Rosetta.
            chains = [result['Chain'] for result in self.ddGDB.call_select_proc("GetChains", parameters = parameters)]
            if strip_other_chains:
                pdb.stripForDDG(chains, KeepHETATMLines, numberOfModels = 1)
            else:
                pdb.stripForDDG(True, KeepHETATMLines, numberOfModels = 1)

            #print('\n'.join(pdb.lines))
            # - Post stripping checks -
            # Get the 'Chain ResidueID' PDB-formatted identifier for each mutation mapped to Rosetta numbering
            # then check again that the mutated positions exist and that the wild-type matches the PDB
            colortext.warning('mutations %s' % (str(mutations)))

            remappedMutations = pdb.remapMutations(mutations, pdbID)

            #resfile = self._createResfile(pdb, remappedMutations)
            mutfile = self._createMutfile(pdb, remappedMutations)

            # Check to make sure that we haven't stripped all the ATOM lines
            if not pdb.GetAllATOMLines():
                raise colortext.Exception("No ATOM lines remain in the stripped PDB file of %s." % pdbID)

            # Check to make sure that CSE and MSE are not present in the PDB
            badresidues = pdb.CheckForPresenceOf(["CSE", "MSE"])
            if badresidues:
                raise colortext.Exception("Found residues [%s] in the stripped PDB file of %s. These should be changed to run this job under Rosetta." % (', '.join(badresidues), pdbID))

            # Turn the lines array back into a valid PDB file
            strippedPDB = '\n'.join(pdb.lines)
        except Exception, e:
            colortext.error("Error in %s, %s: .\n%s" % (experimentID, UserDataSetExperimentID, traceback.format_exc()))
            colortext.warning(str(e))
            return
            colortext.error("\nError: '%s'.\n" % (str(e)))
            colortext.error(traceback.format_exc())
            raise colortext.Exception("An exception occurred retrieving the experimental data for Experiment ID #%s." % experimentID)

        #InputFiles["RESFILE"] = resfile
        InputFiles["MUTFILE"] = mutfile

        ExtraParameters = {}
        InputFiles = pickle.dumps(InputFiles)
        Description = pickle.dumps(Description)
        ExtraParameters = pickle.dumps(ExtraParameters)

        PredictionFieldNames = self.ddGDB.FieldNames.Prediction
        params = {
            PredictionFieldNames.ExperimentID		: experimentID,
            PredictionFieldNames.UserDataSetExperimentID : UserDataSetExperimentID,
            PredictionFieldNames.PredictionSet		: PredictionSet,
            PredictionFieldNames.ProtocolID			: ProtocolID,
            PredictionFieldNames.KeptHETATMLines	: KeepHETATMLines,
            PredictionFieldNames.StrippedPDB		: strippedPDB,
            PredictionFieldNames.ResidueMapping		: pickle.dumps(pdb.get_ddGInverseResmap()),
            PredictionFieldNames.InputFiles			: InputFiles,
            PredictionFieldNames.Description		: Description,
            PredictionFieldNames.Status 			: "queued",
            PredictionFieldNames.ExtraParameters	: ExtraParameters,
            PredictionFieldNames.StoreOutput		: StoreOutput,
        }
        if not testonly:
            self.ddGDB.insertDict('Prediction', params)

            # Add cryptID string
            predictionID = self.ddGDB.getLastRowID()
            entryDate = self.ddGDB.execute_select("SELECT EntryDate FROM Prediction WHERE ID=%s", parameters = (predictionID,))[0]["EntryDate"]
            rdmstring = ''.join(random.sample('0123456789abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', 16))
            cryptID = "%(predictionID)s%(experimentID)s%(PredictionSet)s%(ProtocolID)s%(entryDate)s%(rdmstring)s" % vars()
            cryptID = md5.new(cryptID.encode('utf-8')).hexdigest()
            entryDate = self.ddGDB.execute("UPDATE Prediction SET cryptID=%s WHERE ID=%s", parameters = (cryptID, predictionID))
            return predictionID



    def create_predictions_from_userdataset(self, PredictionSet, userdataset_textid, ProtocolID = None, tagged_subset = None, KeepHETATMLinesStoreOutput = False, Description = {}, InputFiles = {}, quiet = False, testonly = False, only_single_mutations = False, shortrun = False):

        assert(ProtocolID == None) # this should hold until we add Kyle's protocols into the database properly

        assert(self.ddGDB.execute_select("SELECT ID FROM PredictionSet WHERE ID=%s", parameters=(PredictionSet,)))

        #results = self.ddGDB.execute_select("SELECT * FROM UserDataSet WHERE TextID=%s", parameters=(userdataset_textid,))
        results = self.ddGDB.execute_select("SELECT UserDataSetExperiment.* FROM UserDataSetExperiment INNER JOIN UserDataSet ON UserDataSetID=UserDataSet.ID WHERE UserDataSet.TextID=%s", parameters=(userdataset_textid,))
        if not results:
            return False

        if not(quiet):
            colortext.message("Creating predictions for UserDataSet %s using protocol %s" % (userdataset_textid, ProtocolID))
            colortext.message("%d records found in the UserDataSet" % len(results))

        count = 0
        showprogress = not(quiet) and len(results) > 300
        if showprogress:
            print("|" + ("*" * (int(len(results)/100)-2)) + "|")
        for r in results:

            existing_results = self.ddGDB.execute_select("SELECT * FROM Prediction WHERE PredictionSet=%s AND UserDataSetExperimentID=%s", parameters=(PredictionSet, r["ID"]))
            if len(existing_results) > 0:
                #colortext.warning('There already exist records for this UserDataSetExperimentID. You probably do not want to proceed. Skipping this entry.')
                continue

            PredictionID = self.addPrediction(r["ExperimentID"], r["ID"], PredictionSet, ProtocolID, KeepHETATMLines, PDB_ID = r["PDBFileID"], StoreOutput = StoreOutput, ReverseMutation = False, Description = {}, InputFiles = {}, testonly = testonly)
            count += 1
            if showprogress:
                if count > 100:
                    colortext.write(".", "cyan", flush = True)
                    count = 0
            if shortrun and count > 4:
                break
        print("")
        return(True)



if __name__ == '__main__':

    a='''
    PPIJobCreator = PPIJobCreator(passwd = read_file('ddgdb.pw').strip())

    1.
    PPIJobCreator.create_prediction_set("Shane's test prediction set")

    2.
    PPIJobCreator.create_predictions_from_userdataset("Shane's test prediction set", 'AllBindingAffinity', tagged_subset = 'ZEMu')

    # Create a list of PredictionPPI records with:
    - PDB file with stripped chains
       - PPMutagenesis.ID
       - UserPPDataSetExperimentID (specifies PPMutagenesisID and PDB complex definition (PDB ID, PPComplexID, SetNumber))
       - ProtocolID none at present
       - Cost (num residues in stripped PDB)
       - KeptHETATMLines?
       - ResidueMapping (JSON from Rosetta numbering to PDB numbering)
       - InputFiles - mutfile/resfile?
       - Description
       - ScoreVersion
       - ddG (NULL from 23505 to 76632) - 1860 records
       - Scores (NULL from 23505 to 76632)
       - StructureScores (only non-NULL on records 55808, 55809 )

       Each prediction has a set of PredictionStructureScores
         - per prediction, score method (Global p16, Local 8A Noah, ...), score type (DDG, Mutant, Wildtype), run number e.g. 1-50, we store:
           - score components
           - DDG

    3.
    Kyle's runner populates these records with command lines?
    Kyle's runner runs the jobs and saves the results back into the database


    '''