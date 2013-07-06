#!/usr/bin/python2.4
# -*- coding: iso-8859-15 -*-

import sys, os
import MySQLdb
import MySQLdb.cursors
import traceback
import pickle
import time
from datetime import datetime, date
from string import join, letters
import math
import getpass
import itertools
if __name__ == "__main__":
    sys.path.insert(0, "../common")
from tools.db.mysql import DatabaseInterface
from tools.fs.io import read_file, write_file
from tools.bio import rcsb
from tools import colortext
from tools.pdb import PDB, relaxed_amino_acid_codes
from tools.bio import uniprot
from tools.hash import CRC64
from ddgobjects import DBObject

sqrt = math.sqrt
DictCursor = MySQLdb.cursors.DictCursor
StdCursor = MySQLdb.cursors.Cursor

UniProtACToID = {}
PDBToUniProt = {}
UniProtKBACToPDB = {}
uniprotmapping = os.path.join("rawdata", "uniprotmapping.csv")
UniqueIDs = {}

PDBChains = {}

aas = [
    ["A", "ALA", "Alanine",			"non-polar",	"small"],
    ["C", "CYS", "Cysteine",		"non-polar",	"small"],
    ["D", "ASP", "Aspartic acid",	"polar",		"small"],
    ["E", "GLU", "Glutamic acid",	"polar",		"large"],
    ["F", "PHE", "Phenylalanine",	"non-polar",	"large"],
    ["G", "GLY", "Glycine",			"non-polar",	"small"],
    ["H", "HIS", "Histidine",		"polar",		"large"],
    ["I", "ILE", "Isoleucine",		"non-polar",	"large"],
    ["K", "LYS", "Lysine",			"polar",		"large"],
    ["L", "LEU", "Leucine",			"non-polar",	"large"],
    ["M", "MET", "Methionine",		"non-polar",	"large"],
    ["N", "ASN", "Asparagine",		"polar",		"small"],
    ["P", "PRO", "Proline",			"non-polar",	"small"],
    ["Q", "GLN", "Glutamine",		"polar",		"large"],
    ["R", "ARG", "Arginine",		"polar",		"large"],
    ["S", "SER", "Serine",			"polar",		"small"],
    ["T", "THR", "Threonine",		"polar",		"small"],
    ["V", "VAL", "Valine",			"non-polar",	"small"],
    ["W", "TRP", "Tryptophan",		"non-polar",	"large"],
    ["Y", "TYR", "Tyrosine",		"polar",		"large"]
]
# These lists are used to verify record imports from the different databases
# We scan these lists in order to find matches so reorder the chain letters and insertion codes according to what (I'm guessing here) is typical usage in PDB files.
AllowedAminoAcids = [aa[0] for aa in aas]
CommonChainLetters = ['_', '-'] + list(itertools.chain(*[[letters[i+26], letters[i]] for i in range(26)])) + ['4', '1']
AllowedChainLetters = [chr(i) for i in range(32, 127)]
AllowedChainLetters = [c for c in AllowedChainLetters if c not in CommonChainLetters]
AllowedChainLetters = CommonChainLetters + AllowedChainLetters
AllowedInsertionCodes = list(itertools.chain(*[[letters[i+26], letters[i]] for i in range(26)]))

def getTransitiveClosureOfUniProtandPDBIDs(startingIDs, fromtype, totype, numIterations = 3):
    '''This does not resolve very well. I ran it with numIterations = 17 on one PDB ID ('107L') and
    it never resolved. In the last iteration, it had 11188 PDB IDs and 1759 ProTherm IDs.'''
    allowedArgs = ('PDB_ID', 'ACC')
    assert(fromtype in allowedArgs)
    assert(totype in allowedArgs)
    assert(fromtype != totype)

    if type(startingIDs) == type(""):
        startingIDs = set([startingIDs])
    if type(startingIDs) == type([]):
        startingIDs = set(startingIDs)
    assert(type(startingIDs) == type(set()))

    # Start with a set of PDB IDs
    if fromtype == 'ACC':
        startingIDs, MappedUniProtACs, ACCsToPDBs = mapBetweenUniProtandPDB(startingIDs, 'ACC', 'PDB_ID')

    for i in range(numIterations):
        # Starting with a set of PDB IDs
        MappedPDBIDs, MappedUniProtACs1, PDBsToACCs = mapBetweenUniProtandPDB(startingIDs, 'PDB_ID', 'ACC')
        MappedPDBIDs1, MappedUniProtACs, ACCsToPDBs = mapBetweenUniProtandPDB(MappedUniProtACs1, 'ACC', 'PDB_ID')
        assert(len(startingIDs.difference(MappedPDBIDs)) == 0 and len(MappedPDBIDs.difference(startingIDs)) >= 0)
        MappedPDBIDs, MappedUniProtACs2, PDBsToACCs = mapBetweenUniProtandPDB(MappedPDBIDs1, 'PDB_ID', 'ACC')
        MappedPDBIDs2, MappedUniProtACs, ACCsToPDBs = mapBetweenUniProtandPDB(MappedUniProtACs1, 'ACC', 'PDB_ID')
        if (MappedPDBIDs1 == MappedPDBIDs2) and (MappedUniProtACs1 == MappedUniProtACs2):
            print("Mapping from PDBs to ACCs resolved after %d iterations." % numIterations)
            ACCsToIDs = mapFromUniProtACs2IDs(MappedUniProtACs)
            print(PDBsToACCs)
            PDBsToACCs, ACCsToPDBs, ACCsToIDs
        else:
            print(len(MappedPDBIDs1), len(MappedPDBIDs2), len(MappedUniProtACs1), len(MappedUniProtACs2))
            startingIDs = MappedPDBIDs2
    raise Exception("Mapping from PDBs to ACCs did not resolve after %d iterations." % numIterations)

def mapFromUniProtACs2IDs(ACs):
    if type(ACs) == type(""):
        ACs = set([ACs])
    if type(ACs) == type([]):
        ACs = set(ACs)
    if type(ACs) == type(set([])):
        ACIDMapping = {}
        import urllib,urllib2
        url = 'http://www.uniprot.org/mapping/'
        params = {
            'from': 'ACC',
            'to': 'ID',
            'format':'tab',
            'query':join(ACs, " ")
        }
        data = urllib.urlencode(params)
        request = urllib2.Request(url, data)
        response = urllib2.urlopen(request)
        lines = [l for l in response.read(200000).split("\n") if l]
        assert(lines[0]=="From\tTo")
        lines = [line.split("\t") for line in lines[1:]]

        for line in lines:
            assert(len(line) == 2)
            AC = line[0]
            ID = line[1]
            ACIDMapping[AC] = ACIDMapping.get(AC) or []
            ACIDMapping[AC].append(ID)
        assert(set(ACIDMapping.keys()) == ACs)
        return ACIDMapping

def mapBetweenUniProtandPDB(IDs, fromtype, totype):
    '''Takes in  a list of PDB IDs or one single ID string and returns a tuple of two values:
        - a table mapping PDB IDs to UniProt ACCs (not necessarily a bijection) and
        - a table mapping the above UniProt ACCs to UniProt IDs (necessarily a bijection).'''
    allowedArgs = ('PDB_ID', 'ACC')
    assert(fromtype in allowedArgs)
    assert(totype in allowedArgs)
    assert(fromtype != totype)
    if fromtype == 'PDB_ID':
        pdbIndex, ACIndex = (0, 1)
    else:
        pdbIndex, ACIndex = (1, 0)
    if type(IDs) == type(""):
        IDs = set([IDs])
    if type(IDs) == type([]):
        IDs = set(IDs)
    if type(IDs) == type(set([])):
        MappingBetweenIDs = {}
        import urllib,urllib2
        url = 'http://www.uniprot.org/mapping/'
        params = {
            'from': fromtype,
            'to': totype,
            'format':'tab',
            'query':join(IDs, " ")
        }
        data = urllib.urlencode(params)
        request = urllib2.Request(url, data)
        response = urllib2.urlopen(request)
        lines = [l for l in response.read(200000).split("\n") if l]
        if not lines[0]=="From\tTo":
            print(join(lines, "\n"))
        assert(lines[0]=="From\tTo")
        lines = [line.split("\t") for line in lines[1:]]
        ACCs = set()
        MappedPDBIDs = set([line[pdbIndex] for line in lines])
        MappedUniProtACs = set([line[ACIndex] for line in lines])
        if fromtype == "PDB_ID":
            for line in lines:
                assert(len(line) == 2)
                MappingBetweenIDs[line[pdbIndex]] = MappingBetweenIDs.get(line[pdbIndex]) or []
                MappingBetweenIDs[line[pdbIndex]].append(line[ACIndex])
            assert(MappedPDBIDs == IDs)
        else:
            for line in lines:
                assert(len(line) == 2)
                MappingBetweenIDs[line[ACIndex]] = MappingBetweenIDs.get(line[ACIndex]) or []
                MappingBetweenIDs[line[ACIndex]].append(line[pdbIndex])
            assert(MappedUniProtACs == IDs)
        return MappedPDBIDs, MappedUniProtACs, MappingBetweenIDs

def commitUniProtMapping(ddGdb, ACtoID_mapping, PDBtoAC_mapping):
    for uACC, uID in ACtoID_mapping.iteritems():
        results = ddGdb.locked_execute(("SELECT * FROM UniProtKB WHERE %s" % ddGdb.FlatFieldNames.UniProtKB_AC)+"=%s", parameters=(uACC,))
        if results:
            if results[0][ddGdb.FlatFieldNames.UniProtKB_ID] != uID:
                raise Exception("Existing UniProt mapping (%s->%s) does not agree with the passed-in parameters (%s->%s)." % (results[0][ddGdb.FlatFieldNames.UniProtKB_AC],results[0][ddGdb.FlatFieldNames.UniProtKB_ID],uACC,uID))
        else:
            UniProtMapping = {
                ddGdb.FlatFieldNames.UniProtKB_AC : uACC,
                ddGdb.FlatFieldNames.UniProtKB_ID : uID,
            }
            ddGdb.insertDict('UniProtKB', UniProtMapping)
    for updbID, uACCs in PDBtoAC_mapping.iteritems():
        for uACC in uACCs:
            results = ddGdb.locked_execute(("SELECT * FROM UniProtKBMapping WHERE %s" % ddGdb.FlatFieldNames.UniProtKB_AC)+"=%s", parameters=(uACC,))
            associatedPDBsInDB = []
            if results:
                associatedPDBsInDB = [r[ddGdb.FlatFieldNames.PDB_ID] for r in results]
            if updbID not in associatedPDBsInDB:
                UniProtPDBMapping = {
                    ddGdb.FlatFieldNames.UniProtKB_AC : uACC,
                    ddGdb.FlatFieldNames.PDB_ID : updbID,
                }
                ddGdb.insertDict('UniProtKBMapping', UniProtPDBMapping)

class UPFatalException(Exception): pass
def getUniProtMapping(pdbIDs, storeInDatabase = False):
    '''Takes in  a list of PDB IDs or one single ID string and returns a tuple of two values:
        - a table mapping PDB IDs to UniProt ACCs (not necessarily a bijection) and
        - a table mapping the above UniProt ACCs to UniProt IDs (necessarily a bijection).'''
    numtries = 1
    maxtries = 3
    if type(pdbIDs) == type(""):
        pdbIDs = [pdbIDs]
    if type(pdbIDs) == type([]):
        db = ddGDatabase()
        for numtries in range(1, maxtries + 1):
            try:
                import urllib,urllib2
                url = 'http://www.uniprot.org/mapping/'
                print("Querying UniProt: ")
                PDBtoAC_mapping = {}
                params = {
                    'from':'PDB_ID',
                    'to':'ACC',
                    'format':'tab',
                    'query':join(pdbIDs, " ")
                }
                data = urllib.urlencode(params)
                request = urllib2.Request(url, data)
                response = urllib2.urlopen(request)
                lines = [l for l in response.read(200000).split("\n") if l]
                if len(lines) <= 1:
                    raise UPFatalException("No PDB->ACC mapping returned for PDB IDs: %s." % join(pdbIDs, ", "))
                assert(lines[0]=="From\tTo")
                ACCs = set()
                for line in lines[1:]:
                    line = line.split("\t")
                    assert(len(line) == 2)
                    pdbID = line[0]
                    uniprotACC = line[1]
                    PDBtoAC_mapping[pdbID] = PDBtoAC_mapping.get(pdbID) or []
                    PDBtoAC_mapping[pdbID].append(uniprotACC)
                    ACCs.add(uniprotACC)
                ACCs = sorted(list(ACCs))

                ACtoID_mapping = {}
                params = {
                    'from':'ACC',
                    'to':'ID',
                    'format':'tab',
                    'query':join(ACCs, " ")
                }
                data = urllib.urlencode(params)
                request = urllib2.Request(url, data)
                response = urllib2.urlopen(request)
                lines = [l for l in response.read(200000).split("\n") if l]
                if len(lines) <= 1:
                    raise UPFatalException("No ACC->ID mapping returned for ACCs: %s." % join(pdbIDs, ", "))
                assert(lines[0]=="From\tTo")
                assert(len(lines) == len(ACCs) + 1)
                for line in lines[1:]:
                    line = line.split("\t")
                    assert(len(line) == 2)
                    uniprotACC = line[0]
                    uniprotID = line[1]
                    assert(not(ACtoID_mapping.get(uniprotACC)))
                    ACtoID_mapping[uniprotACC] = uniprotID

                if storeInDatabase:
                    commitUniProtMapping(db, ACtoID_mapping, PDBtoAC_mapping)
                colortext.message("success")
                return (ACtoID_mapping, PDBtoAC_mapping)
            except UPFatalException, e:
                raise(str(e).strip())
            except Exception, e:
                emsg = str(e).strip()
                if emsg and emsg.startswith("HTTP Error 500"):
                    colortext.error("failed (HTTP Error 500)")
                else:
                    colortext.error("failed")
                    if emsg:
                        colortext.error(emsg)
                    colortext.error(traceback.format_exc())
                    #print("Lines:\n%s" % lines)
                    #print("ACCs:\n%s" % ACCs)
                time.sleep(1)
    else:
        raise Exception("Expected a list of PDB IDs or a string with a single ID.")

    db.close()
    raise Exception("The request to UniProt failed %d times." % (numtries))


def computeStandardDeviation(values):
    sum = 0
    n = len(values)

    for v in values:
        sum += v

    mean = sum / n
    sumsqdiff = 0

    for v in values:
        t = (v - mean)
        sumsqdiff += t * t

    variance = sumsqdiff / n
    stddev = math.sqrt(variance)

    return stddev, variance

def readUniProtMap(ddGdb):
    if not (UniProtACToID) or not(PDBToUniProt) or not(UniProtKBACToPDB):
        results = ddGdb.locked_execute("SELECT * FROM UniProtKB")
        for r in results:
            assert(not(UniProtACToID.get(r[ddGdb.FlatFieldNames.UniProtKB_AC])))
            UniProtACToID[r[ddGdb.FlatFieldNames.UniProtKB_AC]] = r[ddGdb.FlatFieldNames.UniProtKB_ID]

        results = ddGdb.locked_execute("SELECT * FROM UniProtKBMapping")
        for r in results:
            pdbID = r[ddGdb.FlatFieldNames.PDB_ID]
            UPAC = r[ddGdb.FlatFieldNames.UniProtKB_AC]
            UPID = UniProtACToID[UPAC]
            if PDBToUniProt.get(pdbID):
                PDBToUniProt[pdbID].append((UPAC, UPID))
            else:
                PDBToUniProt[pdbID] = [(UPAC, UPID)]
            if UniProtKBACToPDB.get(UPAC):
                UniProtKBACToPDB[UPAC].append(pdbID)
            else:
                UniProtKBACToPDB[UPAC] = [pdbID]

class PDBStructure(DBObject):

    # At the time of writing, these PDB IDs had no UniProt entries
    NoUniProtIDs = set(['1GTX', '1UOX', '1WSY', '1YYJ', '2IMM', '2MBP'])

    # At the time of writing, these PDB IDs had no JRNL lines
    NoPublicationData = ['2FX5']

    def __init__(self, ddGdb, pdbID, content = None, protein = None, source = None, filepath = None, UniProtAC = None, UniProtID = None, testonly = False):
        '''UniProtACs have forms like 'P62937' whereas UniProtIDs have forms like 'PPIA_HUMAN.'''

        super(PDBStructure, self).__init__(ddGdb)
        FieldNames_ = ddGdb.FlatFieldNames
        self.dict = {
            FieldNames_.PDB_ID : pdbID,
            FieldNames_.Content : content,
            FieldNames_.FASTA : None,
            FieldNames_.Resolution : None,
            FieldNames_.NumberOfMolecules : None,
            FieldNames_.Length : None,
            FieldNames_.Protein : protein,
            FieldNames_.Source : source,
            FieldNames_.Techniques : None,
            FieldNames_.BFactors : None,
            FieldNames_.Publication : None,
        }
        self.testonly = testonly
        self.filepath = filepath
        self.UniProtAC = UniProtAC
        self.UniProtID = UniProtID
        self.ACtoID_mapping = None
        self.PDBtoAC_mapping = None
        self.chains = None

    def getPDBContents(self):
        ddGdb = self.ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        d = self.dict
        id = d[FieldNames_.PDB_ID]
        if len(id) != 4:
            print(id)
        assert(len(id) <= 10)

        if self.filepath:
            filename = self.filepath
        else:
            filename = os.path.join("../pdbs", id + ".pdb")
        contents = None
        chains = {}

        if not os.path.exists(filename):
            sys.stdout.write("The file for %s is missing. Retrieving it now from RCSB: " % (id))
            sys.stdout.flush()
            try:
                contents = rcsb.getPDB(id)
                write_file(os.path.join("../pdbs", id + ".pdb"), contents)
                colortext.message("success")
            except:
                colortext.error("failure")
                raise Exception("Error retrieving %s." % filename)

        else:
            contents = read_file(filename)

        resolution = p.get_resolution()
        print(p.get_techniques())
        print(p.get_ATOM_and_HETATM_chains())

        resolution = None
        lines = contents.split("\n")
        pdb_o = PDB(lines)
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chains[line[21]] = True
            elif line.startswith("EXPDTA"):
                techniques = line[10:71].split(";")
                for k in range(len(techniques)):
                    techniques[k] = techniques[k].strip()
                techniques = join(techniques, ";")
            elif line[0:6] == "REMARK" and line[9] == "2" and line[11:22] == "RESOLUTION.":
                #if id == :
                #	line = "REMARK   2 RESOLUTION. 3.00 ANGSTROMS.

                                # This code SHOULD work but there are badly formatted PDBs in the RCSB database.
                # e.g. "1GTX"
                #if line[31:41] == "ANGSTROMS.":
                #	try:
                #		resolution = float(line[23:30])
                #	except:
                #		raise Exception("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard. Expected data for diffraction experiments." % line )
                #if line[23:38] == "NOT APPLICABLE.":
                #	resolution = "N/A"
                #else:
                #	raise Exception("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard." % line )
                #
                # Instead, we use the code below:
                strippedline = line[22:].strip()
                Aindex = strippedline.find("ANGSTROMS.")
                NA = strippedline == "NOT APPLICABLE."
                if NA:
                    resolution = "N/A"
                elif Aindex != -1 and strippedline.endswith("ANGSTROMS."):
                    if strippedline[:Aindex].strip() == "NULL":
                        resolution = "N/A" # Yes, yes, yes, I know. Look at 1WSY.pdb.
                    else:
                        try:
                            resolution = float(strippedline[:Aindex].strip())
                        except:
                            raise Exception("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard. Expected data for diffraction experiments." % line )
                else:
                    raise Exception("Error parsing PDB file to determine resolution. The resolution line\n  '%s'\ndoes not match the PDB standard." % line )

        self.chains = chains.keys()
        PDBChains[d[FieldNames_.PDB_ID]] = self.chains

        if not resolution:
            raise Exception("Could not determine resolution for %s." % filename)
        if resolution == "N/A":
            resolution = None

        UniqueIDs[id] = True

        if id not in self.NoUniProtIDs:
            readUniProtMap(self.ddGdb)
            if not PDBToUniProt.get(id):
                if not (self.UniProtAC and self.UniProtID):
                    ACtoID_mapping, PDBtoAC_mapping = getUniProtMapping(id, storeInDatabase = False)
                    if not (PDBtoAC_mapping and ACtoID_mapping):
                        raise Exception("Could not find a UniProt mapping for %s in %s." % (id, uniprotmapping))
                    self.ACtoID_mapping = ACtoID_mapping
                    self.PDBtoAC_mapping = PDBtoAC_mapping

        d[FieldNames_.Content] = contents
        d[FieldNames_.Resolution] = resolution
        d[FieldNames_.Techniques] = techniques

        if id not in self.NoPublicationData:
            self.getPublication()
        self.getFASTA()


        pdbID = d[FieldNames_.PDB_ID]
        foundRes = pdb_o.CheckForPresenceOf(["CSE", "MSE"])
        if foundRes:
            colortext.error("The PDB %s contains residues which could affect computation (%s)." % (pdbID, join(foundRes, ", ")))
            if "CSE" in foundRes:
                colortext.error("The PDB %s contains CSE. Check." % pdbID)
            if "MSE" in foundRes:
                colortext.error("The PDB %s contains MSE. Check." % pdbID)

        d[FieldNames_.BFactors] = pickle.dumps(pdb_o.ComputeBFactors())
        return contents

    def getPublication(self):
        '''Extracts the PDB source information.'''
        ddGdb = self.ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        d = self.dict

        PUBTYPES = ['ISSN', 'ESSN']

        p = PDB(d[FieldNames_.Content].split("\n"))
        j = p.getJournal()
        pdbID = d[FieldNames_.PDB_ID].strip()

        if False:
            if j["published"]:
                print(j["REFN"]["type"])
                print(j["REFN"]["ID"])
            print("doi: %s" % j["DOI"])
            print("--------")
            if j["DOI"]:
                pass
            print(join(j["lines"],"\n"))
            for k, v in j.iteritems():
                if k != "lines":
                    print("%s: %s" % (k, v))
            print("********\n")

        # We identify the sources for a PDB identifier with that identifier
        SourceID = "PDB:%s" % pdbID
        sourceExists = self.ddGdb.locked_execute("SELECT ID FROM Source WHERE ID=%s", parameters=(SourceID,))
        if not sourceExists:
            if not self.testonly:
                self.ddGdb.insertDict(FieldNames_.Source, {FieldNames_.ID : SourceID})

        d[FieldNames_.Publication] = SourceID

        locations = self.ddGdb.locked_execute("SELECT * FROM SourceIdentifier WHERE SourceID=%s", parameters=(SourceID,))
        publocations = [location for location in locations if location[FieldNames_.Type] in PUBTYPES]
        doilocations = [location for location in locations if location[FieldNames_.Type] == "DOI"]
        assert(len(publocations) <= 1)
        assert(len(doilocations) <= 1)
        if j["published"]:
            skip = False
            if publocations:
                location = publocations[0]
                if j["REFN"]["type"] == location[FieldNames_.Type]:
                    if j["REFN"]["ID"] != location[FieldNames_.ID]:
                        colortext.warning("REFN: Check that the SourceIdentifier data ('%s') matches the PDB REFN data ('%s')." % (str(location), j["REFN"]))
            else:
                assert(j["REFN"]["type"] in PUBTYPES)
                source_location_dict  =  {
                    FieldNames_.SourceID	: SourceID,
                    FieldNames_.ID			: j["REFN"]["ID"],
                    FieldNames_.Type		: j["REFN"]["type"],
                }
                if not self.testonly:
                    print(source_location_dict)
                    self.ddGdb.insertDict(ddGdb.FieldNames.SourceIdentifier._name, source_location_dict)
                else:
                    print(source_location_dict)
        if j["DOI"]:
            if doilocations:
                location = doilocations[0]
                if j["DOI"] != location[FieldNames_.ID]:
                    colortext.warning("DOI: Check that the SourceIdentifier data ('%s') matches the PDB DOI data ('%s')." % (str(doilocations), j["DOI"]))
            else:
                source_location_dict = {
                    FieldNames_.SourceID	: SourceID,
                    FieldNames_.ID			: j["DOI"],
                    FieldNames_.Type		: "DOI",
                }
                if not self.testonly:
                    self.ddGdb.insertDict(ddGdb.FieldNames.SourceIdentifier._name, source_location_dict)
                else:
                    print(source_location_dict)

    def getFASTA(self):
        ddGdb = self.ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        pdbID = self.dict[FieldNames_.PDB_ID]
        if len(pdbID) == 4:
            try:
                results = ddGdb.execute_select("SELECT FASTA FROM Structure WHERE PDB_ID=%s", parameters=(pdbID,))
                if results:
                    assert(len(results) == 1)
                    fastafile = results[0]["FASTA"]
                else:
                    fastafile = rcsb.getFASTA(pdbID)
                self.dict[FieldNames_.FASTA] = fastafile
                return
            except:
                pass
        raise Exception("No FASTA file could be found for %s." % pdbID)

    def parseFASTA(self):
        pdbID = self.dict[self.ddGdb.FlatFieldNames.PDB_ID]
        results = self.ddGdb.execute_select("SELECT FASTA FROM Structure WHERE PDB_ID = %s", parameters = (pdbID,))
        if results:
            assert(len(results) == 1)
            return rcsb.parseFASTAs(results[0]["FASTA"])
        else:
            return rcsb.parseFASTAs(pdbID)

    def find(self):
        results = self.ddGdb.locked_execute('''SELECT PDB_ID FROM Structure WHERE PDB_ID=%s''', parameters = (self.PDB_ID, ))
        if results:
            return self.PDB_ID, self.PDB_ID
        return None, None

    def test(self):
        # #todo: Implement this later
        pass

    def commit(self, testonly = False):
        '''Returns the database record ID if an insert occurs but will typically return None if the PDB is already in the database.'''
        ddGdb = self.ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        d = self.dict
        testonly = testonly or self.testonly

        self.getPDBContents()
        assert(PDBChains.get(d[FieldNames_.PDB_ID]))

        if self.UniProtAC and self.UniProtID:
            # todo: Append to uniprotmapping.csv file
            results = self.ddGdb.locked_execute(("SELECT * FROM %s WHERE %s" % (FieldNames_.UniProtKB, FieldNames_.UniProtKB_AC))+"=%s", parameters=(self.UniProtAC,))
            if results:
                if results[0][FieldNames_.UniProtKB_ID] != self.UniProtID:
                    raise Exception("Existing UniProt mapping (%s->%s) does not agree with the passed-in parameters (%s->%s)." % (results[0][FieldNames_.UniProtKB_AC],results[0][FieldNames_.UniProtKB_ID],self.UniProtAC,self.UniProtID))
            else:
                UniProtMapping = {
                    FieldNames_.UniProtKB_AC : self.UniProtAC,
                    FieldNames_.UniProtKB_ID : self.UniProtID,
                }
                if not testonly:
                    self.ddGdb.insertDict(FieldNames_.UniProtKB, UniProtMapping)

        results = self.ddGdb.locked_execute("SELECT * FROM Structure WHERE PDB_ID=%s", parameters = (d[FieldNames_.PDB_ID]))

        if results:
            assert(len(results) == 1)
            result = results[0]
            pdbID = result[FieldNames_.PDB_ID]
            for k, v in d.iteritems():
                if k != FieldNames_.PDB_ID:
                    if k == FieldNames_.Techniques and result[k] == "":
                        SQL = "UPDATE Structure SET %s" % k
                        SQL += "=%s WHERE PDB_ID=%s"
                        if not testonly:
                            results = self.ddGdb.locked_execute(SQL, parameters = (v, pdbID))
                    if d[k] and not(result[k]):
                        SQL = "UPDATE Structure SET %s" % k
                        SQL += "=%s WHERE PDB_ID=%s"
                        if not testonly:
                            results = self.ddGdb.locked_execute(SQL, parameters = (v, pdbID))
        else:
            if not testonly:
                pdb_dict = {
                    FieldNames_.ID	    	: d[FieldNames_.PDB_ID],
                    FieldNames_.Content		: d[FieldNames_.Content],
                    FieldNames_.FASTA		: d[FieldNames_.FASTA],
                    FieldNames_.Resolution	: d[FieldNames_.Resolution],

                    FieldNames_.NumberOfMolecules	: None, # todo
                    FieldNames_.Length			 	: None, # todo

                    FieldNames_.Protein		: d[FieldNames_.Protein],
                    FieldNames_.Source		: d[FieldNames_.Source],
                    FieldNames_.Techniques	: d[FieldNames_.Techniques],
                    FieldNames_.BFactors	: d[FieldNames_.BFactors],
                    FieldNames_.Publication	: d[FieldNames_.Publication],
                }
                self.ddGdb.insertDict(FieldNames_.PDBFile, pdb_dict)
                self.databaseID = self.ddGdb.getLastRowID()

        if self.UniProtAC and self.UniProtID:
            results = self.ddGdb.locked_execute(("SELECT * FROM UniProtKBMapping WHERE %s" % FieldNames_.UniProtKB_AC)+"=%s", parameters=(self.UniProtAC,))
            if results:
                if results[0][FieldNames_.PDB_ID] != d[FieldNames_.PDB_ID]:
                    raise Exception("Existing UniProt mapping (%s->%s) does not agree with the passed-in parameters (%s->%s)." % (results[0][FieldNames_.UniProtKB_AC],results[0][FieldNames_.PDB_ID],self.UniProtAC,d[FieldNames_.PDB_ID]))
            else:
                UniProtPDBMapping = {
                    FieldNames_.UniProtKB_AC : self.UniProtAC,
                    FieldNames_.PDB_ID : d[FieldNames_.PDB_ID],
                }
                if not testonly:
                    self.ddGdb.insertDict('UniProtKBMapping', UniProtPDBMapping)

        # Store the UniProt mapping in the database
        if d[FieldNames_.PDB_ID] not in self.NoUniProtIDs:
            if not (self.ACtoID_mapping and self.PDBtoAC_mapping):
                self.ACtoID_mapping, self.PDBtoAC_mapping = getUniProtMapping(d[FieldNames_.PDB_ID], storeInDatabase = False)
            assert(self.ACtoID_mapping and self.PDBtoAC_mapping)
            commitUniProtMapping(self.ddGdb, self.ACtoID_mapping, self.PDBtoAC_mapping)

        if not testonly:
            return self.databaseID
        else:
            return None

        return self.databaseID

    def remove(self):
        # We do not usually want to remove PDBs from the database
        pass

    def __repr__(self):
        ddGdb = self.ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        d = self.dict
        str = []
        str.append("%s: %s" % (FieldNames_.PDB_ID, d[FieldNames_.PDB_ID]))
        str.append("%s: %s" % (FieldNames_.Protein, d[FieldNames_.Protein]))
        return join(str, "\n")

MutantMatchFudges = {
    # This dict is used to help match up mutants where the length does not match.
    # This could be due to:
    #   tagging of the wildtype or mutant sequence on one or both ends;
    #   missing coordinates;
    #   a bad mutant PDB ID;
    #   or something else.
    # The values are slice co-ordinates for the sequence strings.
    ('2CI2', 'I', '1YPC', 'I') : ((  19, None),(None, None),'Wildtype has a longer head (SSVEKKPEGVNTGAGDRHN) with an L20M mutation. The mutant also has an E78Q mutation.'),
    ('2CI2', 'I', '1YPB', 'I') : ((  19, None),(None, None),'Wildtype has a longer head (SSVEKKPEGVNTGAGDRHN) with an L20M mutation. The mutant also has an E78Q mutation.'),
    ('2CI2', 'I', '1YPA', 'I') : ((  19, None),(None, None),'Wildtype has a longer head (SSVEKKPEGVNTGAGDRHN) with an L20M mutation. The mutant also has an E78Q mutation.'),
    ('2CI2', 'I', '1COA', 'I') : ((  19, None),(None, None),'Wildtype has a longer head (SSVEKKPEGVNTGAGDRHN) with an L20M mutation.'),
    ('1BNI', 'A', '1BRG', 'A') : ((   2, None),(None, None),'Wildtype has a longer head (AQ).'),
    ('1CEY', 'A', '1E6K', 'A') : ((None, None),(   2, None),'Mutant has a longer head (MR) followed by an A->S mutation.'),
    ('1CEY', 'A', '1E6L', 'A') : ((   1, None),(None, None),'Wildtype has information about an A residue at the start.'),
    ('1CEY', 'A', '1E6M', 'A') : ((None, None),(None, None),'Mutant has an A->S mutation at the start.'),
    ('1N0J', 'A', '1VAR', 'A') : ((   1, None),(None, None),'Wildtype has information about an M at the start.'),
    ('1IHB', 'A', '1MX2', 'A') : ((None, None),(None,   -6),'Mutant has information about the last six residues (GATNLQ).'),
    ('1IHB', 'A', '1MX4', 'A') : ((None, None),(None,   -6),'Mutant has information about the last six residues (GATNLQ).'),
    ('1IHB', 'A', '1MX6', 'A') : ((None, None),(None,   -6),'Mutant has information about the last six residues (GATNLQ).'),
}

class MutantChainMap(object):
    def __init__(self, ddGdb, wildtype_PDB_ID, mutant_PDB_ID, wildtype_chain, mutant_chain):
        self.ddGdb = ddGdb
        self.wildtype_PDB_ID = wildtype_PDB_ID
        self.mutant_PDB_ID = mutant_PDB_ID
        self.wildtype_chain = wildtype_chain
        self.mutant_chain = mutant_chain
        self.notes = None
        fudge = MutantMatchFudges.get((wildtype_PDB_ID, wildtype_chain, mutant_PDB_ID, mutant_chain))
        if fudge:
            self.notes = fudge[2]

    def __repr__(self):
        return "%(wildtype_PDB_ID)s chain %(wildtype_chain)s -> %(mutant_PDB_ID)s chain %(mutant_chain)s" % self.__dict__

    def test(self, experiment_mutations):
        fasta = rcsb.parseFASTAs([self.wildtype_PDB_ID, self.mutant_PDB_ID], silent = True, database_ref = self.ddGdb, database_table = "Structure", database_field = "FASTA", database_IDfield = "PDB_ID")

        wildtype_PDB_ID = self.wildtype_PDB_ID
        wildtype_chain = self.wildtype_chain
        if not fasta[wildtype_PDB_ID].get(wildtype_chain):
            raise colortext.Exception("Chain %s does not exist in wildtype %s.\n%s" % (wildtype_chain, wildtype_PDB_ID, fasta[wildtype_PDB_ID]))
        wildtype_sequence = fasta[self.wildtype_PDB_ID][wildtype_chain]

        mutant_PDB_ID = self.mutant_PDB_ID
        mutant_chain = self.mutant_chain
        if not fasta[mutant_PDB_ID].get(mutant_chain):
            raise colortext.Exception("Chain %s does not exist in mutant %s.\n%s" % (mutant_chain, mutant_PDB_ID, fasta[mutant_PDB_ID]))
        mutant_sequence = fasta[mutant_PDB_ID][mutant_chain]

        fudge = MutantMatchFudges.get((wildtype_PDB_ID, wildtype_chain, mutant_PDB_ID, mutant_chain))
        if fudge:
            wildtype_sequence = wildtype_sequence[fudge[0][0]:fudge[0][1]]
            mutant_sequence = mutant_sequence[fudge[1][0]:fudge[1][1]]

        mutations_by_diff = []
        mutation_pairs = []

        wildtype_leftpadding = 0
        mutant_leftpadding = 0

        if len(wildtype_sequence) == len(mutant_sequence):
            for i in range(len(wildtype_sequence)):
                if wildtype_sequence[i] != mutant_sequence[i]:
                    if fudge:
                        mutations_by_diff.append("%s->%s at position %d of the wt FASTA" % (wildtype_sequence[i], mutant_sequence[i], i+1+(fudge[0][0] or 0)))
                    else:
                        mutations_by_diff.append("%s->%s at position %d of the wt FASTA" % (wildtype_sequence[i], mutant_sequence[i], i+1))
                    mutation_pairs.append((wildtype_sequence[i], mutant_sequence[i]))
        else:
            # Simple try to match up the sequences
            if len(wildtype_sequence) < len(mutant_sequence) and len(wildtype_sequence) > 15:
                idx = mutant_sequence.find(wildtype_sequence[3:10])
                if idx != -1:
                    wildtype_leftpadding = idx - 3
                else:
                    idx = mutant_sequence.find(wildtype_sequence[6:13])
                    if idx != -1:
                        wildtype_leftpadding = idx - 6
            elif len(mutant_sequence) < len(wildtype_sequence) and len(mutant_sequence) > 15:
                idx = wildtype_sequence.find(mutant_sequence[3:10])
                if idx != -1:
                    mutant_leftpadding = idx - 3
                else:
                    idx = wildtype_sequence.find(mutant_sequence[6:13])
                    if idx != -1:
                        mutant_leftpadding = idx - 6
            colortext.printf("%s%s" % (" " * wildtype_leftpadding, wildtype_sequence), "lightpurple")
            colortext.printf("%s%s" % (" " * mutant_leftpadding, mutant_sequence), "cyan")
            raise Exception("Failed matching. wildtype_leftpadding = %d, mutant_leftpadding = %d." % (wildtype_leftpadding, mutant_leftpadding))


        if len(mutations_by_diff) != len(experiment_mutations):
            colortext.error("The number of experiment mutations and the diff-list from the sequences disagree in length:")
            print(wildtype_PDB_ID, mutant_PDB_ID)
            colortext.printf(wildtype_sequence, "lightpurple")
            colortext.printf(mutant_sequence, "cyan")
            print(mutations_by_diff)
            print(experiment_mutations)
            raise colortext.Exception("The number of experiment mutations and the diff-list from the sequences disagree in length:")
        else:
            if sorted(mutation_pairs) != sorted([(m.WildTypeAA, m.MutantAA) for m in experiment_mutations]):
                colortext.error("The wildtype/mutant amino acid pairs between the experiment mutations and the diff-list from the sequences disagree disagree:")
                print(wildtype_PDB_ID, mutant_PDB_ID)
                print(mutations_by_diff)
                print(experiment_mutations)
                raise colortext.Exception("The wildtype/mutant amino acid pairs between the experiment mutations and the diff-list from the sequences disagree disagree:")


class ExperimentDefinition(DBObject):

    AllowedSecondaryStructurePositions = ['Coil', 'Helix', 'Sheet', 'Turn', '3_10-Helix']

    def __init__(self, ddGdb, pdbid, interface = None):
        super(ExperimentDefinition, self).__init__(ddGdb)
        self.PDB_ID = pdbid
        self.interface = interface
        self.mutantmaps = []
        self.mutations = []
        self.chains = []

    def __repr__(self):
        s = []
        s.append("PDB ID: %s" % self.PDB_ID)
        if self.interface:
            s.append("Interface: %s" % self.interface)
        if self.mutantmaps:
            s.append("Mutants: %s" % join(map(str, self.mutantmaps), ","))
        if self.chains:
            s.append("Chains: %s" % join(self.chains, ","))
        if self.mutations:
            s.append("Mutations:\n  %s" % join(map(str, self.mutations), "\n  "))
        return join(s, "\n")

    def addMutant(self, mutant_PDB_ID, wildtype_chain, mutant_chain):
        if wildtype_chain in [m.wildtype_chain for m in self.mutantmaps if m.mutant_PDB_ID == mutant_PDB_ID]:
            raise colortext.Exception("Error: Adding wildtype chain %s twice to the mutant %s of %s." % (wildtype_chain, mutant_PDB_ID, self.PDB_ID))
        if mutant_chain in [m.mutant_chain for m in self.mutantmaps if m.mutant_PDB_ID == mutant_PDB_ID]:
            raise colortext.Exception("Error: Adding mutant chain %s twice to the mutant %s of %s." % (mutant_chain, mutant_PDB_ID, self.PDB_ID))

        self.mutantmaps.append(MutantChainMap(self.ddGdb, self.PDB_ID, mutant_PDB_ID, wildtype_chain, mutant_chain))

    def addChain(self, chainID, ID = ""):
        if not chainID in AllowedChainLetters:
            raise Exception("An exception occurred processing a chain (ID='%s'):%s.\n\tThe chain '%s' is invalid." % (ID, errors, chainID))
        if chainID in self.chains:
            raise Exception("Error: Adding chain %s twice to %s." % (chainID, self.PDB_ID))
        self.chains.append(chainID)

    def getChains(self):
        return self.chains

    def addMutation(self, mutation, ID = None):
        errors = []
        residueID = ("%s" % mutation.ResidueID).strip()
        if not mutation.Chain in AllowedChainLetters:
            errors.append("The chain '%s' is invalid." % mutation.Chain)
        if not mutation.WildTypeAA in AllowedAminoAcids:
            errors.append("The wildtype amino acid '%s' is invalid." % mutation.WildTypeAA)
        if not mutation.MutantAA in AllowedAminoAcids:
            errors.append("The mutant amino acid '%s' is invalid." % mutation.MutantAA)
        residueID = mutation.ResidueID
        if not residueID.isdigit():
            if not residueID[0:-1].isdigit():
                errors.append("The residue '%s' is invalid." % residueID)
            elif residueID[-1] not in AllowedInsertionCodes:
                errors.append("The insertion code '%s' of residue '%s' is invalid." % (residue[-1], residueID))
        if mutation.SecondaryStructurePosition and not(mutation.SecondaryStructurePosition in self.AllowedSecondaryStructurePositions):
            errors.append("The secondary structure location '%s' is invalid." % (mutation.SecondaryStructurePosition))
        if mutation.AccessibleSurfaceArea:
            assert(type(mutation.AccessibleSurfaceArea) == type(0.1))

        if errors:
            ID = ID or ""
            if ID:
                ID = ", ID %s" % ID
            errors = join(['\t%s\n' % e for e in errors], "")
            raise Exception("An exception occurred processing a mutation (ID='%s'):\n%s" % (ID, errors))

        self.mutations.append(mutation)

    def remove(self):
        ddGdb = self.ddGdb
        if not self.databaseID:
            raise Exception("Cannot remove a record with no corresponding database ID.")
        ddGdb.locked_execute("DELETE FROM ExperimentChain WHERE ExperimentID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM ExperimentInterface WHERE ExperimentID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM ExperimentMutant WHERE ExperimentID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM ExperimentMutation WHERE ExperimentID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM Experiment WHERE ID=%s", parameters = (self.databaseID,))

    def find(self):

        FieldNames = self.ddGdb.FieldNames
        results = self.ddGdb.locked_execute('''SELECT ID, Chain, ResidueID, WildTypeAA, MutantAA FROM Experiment INNER JOIN ExperimentMutation ON Experiment.ID = ExperimentMutation.ExperimentID WHERE Experiment.Structure=%s''', parameters = (self.PDB_ID, ))
        groupedresults = {}

        ExperimentIDsWithMatchingMutations = {}
        mymutations = []
        for mutation in self.mutations:
            mymutations.append((mutation.Chain.upper(), mutation.ResidueID.upper(), mutation.WildTypeAA.upper(), mutation.MutantAA.upper()))
        mymutations = sorted(mymutations)
        for r in results:
            ExperimentID = r[FieldNames.Experiment.ID]
            groupedresults[ExperimentID] = groupedresults.get(ExperimentID, [])
            groupedresults[ExperimentID].append((r['Chain'].upper(), r['ResidueID'].upper(), r['WildTypeAA'].upper(), r['MutantAA'].upper()))

        for ExperimentID, dbmutations in groupedresults.iteritems():
            if sorted(dbmutations) == mymutations:
                ExperimentIDsWithMatchingMutations[ExperimentID] = sorted(dbmutations)

        matchedAllInformation = True
        if ExperimentIDsWithMatchingMutations:
            #print("Database record found with matched mutations: %s" % groupedresults[ExperimentID])
            for ExperimentID, foundMutations in sorted(ExperimentIDsWithMatchingMutations.iteritems()):
                results = self.ddGdb.locked_execute('''SELECT Chain FROM Experiment INNER JOIN ExperimentChain ON Experiment.ID = ExperimentChain.ExperimentID WHERE Experiment.ID=%s''', parameters = (ExperimentID, ), cursorClass = StdCursor)
                foundChains = sorted([r[0] for r in results])
                if foundChains != sorted(self.chains):
                    matchedAllInformation = False
                else:
                    pass
                    #print("Database record found with matched chains: %s" % foundChains)

                myinterfaces = []
                if self.interface:
                    myinterfaces = [self.interface]
                results = self.ddGdb.locked_execute('''SELECT Interface FROM Experiment INNER JOIN ExperimentInterface ON Experiment.ID = ExperimentInterface.ExperimentID WHERE Experiment.ID=%s''', parameters = (ExperimentID, ), cursorClass = StdCursor)
                foundInterfaces = sorted([r[0] for r in results])
                if foundInterfaces != myinterfaces:
                    matchedAllInformation = False
                else:
                    pass
                    #print("Database record found with matched interface: %s" % foundInterfaces)

                results = self.ddGdb.locked_execute('''SELECT Mutant, WildTypeChainID, MutantChainID FROM Experiment INNER JOIN ExperimentMutant ON Experiment.ID = ExperimentMutant.ExperimentID WHERE Experiment.ID=%s''', parameters = (ExperimentID, ))
                storedMutants = []
                for r in results:
                    storedMutants.append((r["Mutant"], r["Mutant"], r["Mutant"]))
                    foundMutant = False
                    for mm in self.mutantmaps:
                        if mm.mutant_PDB_ID == r["Mutant"] and mm.wildtype_chain == r["WildTypeChainID"] and mm.mutant_chain == r["MutantChainID"]:
                            foundMutant = True
                            break
                    if not foundMutant:
                        matchedAllInformation = False
                storedMutants = sorted(storedMutants)

                if not matchedAllInformation:
                    raise Exception('''Database record %d with PDB ID %s was found with:\n  chains %s,\n  mutations %s,\n  interfaces %s, and\n  mutant maps %s\nwhich did not match the record being added which has:\n  chains %s,\n  mutations %s,\n  interfaces %s, and\n  mutant maps %s.''' % (ExperimentID, self.PDB_ID, foundChains, foundMutations, foundInterfaces, storedMutants, self.chains, mymutations, myinterfaces, self.mutantmaps))

                return ExperimentID, groupedresults[ExperimentID]

        return None, None

    def test(self, testonly = False, pdbPath = None, mutationAllowedToBeStoredDespiteMissingCoordinates = False):

        # Assertions
        wildtype_chain_pairs = [(m.mutant_PDB_ID, m.wildtype_chain) for m in self.mutantmaps]
        mutant_chains_pairs = [(m.mutant_PDB_ID, m.mutant_chain) for m in self.mutantmaps]
        assert(len(wildtype_chain_pairs) == len(set(wildtype_chain_pairs)))
        assert(len(mutant_chains_pairs) == len(set(mutant_chains_pairs)))
        assert(len(self.chains) == len(set(self.chains)))
        if self.mutantmaps:
            assert(sorted(list(set([m.wildtype_chain for m in self.mutantmaps]))) == sorted(self.chains))

        # Check that the mutant structures have the correct sequences
        for mutant_chain_map in self.mutantmaps:
            chain_mutations = [m for m in self.mutations if m.Chain == mutant_chain_map.wildtype_chain]
            mutant_chain_map.test(chain_mutations)

        # Look for an existing record
        ExperimentID = super(ExperimentDefinition, self).test()
        if ExperimentID:
            return ExperimentID

        failed = False

        # Make sure all mutant structures are homologous
        mutantsequences = None
        MutantStructures = []
        mutant_pdbs = []

        mutants = [mm.mutant_PDB_ID for mm in self.mutantmaps]
        for mutant in mutants:
            # Create a PDBStructure object for the mutant structure
            MutantStructure = PDBStructure(self.ddGdb, mutant)

            fastasequences = MutantStructure.parseFASTA().sequences
            rawsequences = sorted(list(set([fs[2] for fs in fastasequences])))

            if mutantsequences:
                if mutantsequences != rawsequences:
                    # todo: This could be relaxed to make sure that the mutant sequence (original sequence + mutation) is present as a chain in all mutant structures
                    raise Exception("The experiment has multiple mutants but they are not homologous:\n%s\n%s." % (mutantsequences, rawsequences))
            else:
                mutantsequences = rawsequences
            MutantStructures.append(MutantStructure)
            contents = MutantStructure.getPDBContents()
            mutant_pdbs.append(PDB(contents.split("\n")))

        # Add all mutant structures to the database
        for MutantStructure in MutantStructures:
            results = self.ddGdb.execute_select("SELECT PDB_ID FROM Structure WHERE PDB_ID = %s", parameters = (MutantStructure.dict[self.ddGdb.FlatFieldNames.PDB_ID]))
            if not results:
                MutantStructure.commit()

        # Add the wildtype structure to the database
        if pdbPath:
            WildTypeStructure = PDBStructure(self.ddGdb, self.PDB_ID, filepath = os.path.join(pdbPath, "%s.pdb" % self.PDB_ID))
        else:
            WildTypeStructure = PDBStructure(self.ddGdb, self.PDB_ID)
        results = self.ddGdb.execute_select("SELECT PDB_ID FROM Structure WHERE PDB_ID = %s", parameters = (WildTypeStructure.dict[self.ddGdb.FlatFieldNames.PDB_ID]))
        if not results:
            WildTypeStructure.commit()

        # Create a PDB object from the wildtype
        contents = WildTypeStructure.getPDBContents()
        wildtype_pdb = PDB(contents.split("\n"))

        # Check for the existance of CSE or MSE residues which can cause problems for Rosetta
        badResidues = ["CSE", "MSE"]
        foundRes = wildtype_pdb.CheckForPresenceOf(badResidues)
        if foundRes:
            colortext.warning("The PDB %s contains residues which could affect computation (%s)." % (pdbID, join(foundRes, ", ")))
            failed = True
            for res in foundRes:
                colortext.warning("The PDB %s contains %s. Check." % (pdbID, res))

        # Sanity check that the wildtypes of all mutations are correct
        FieldNames = self.ddGdb.FieldNames.ExperimentMutation
        for mutation in self.mutations:
            foundMatch = False
            foundAlternativeWildType = None
            for resid, wtaa in sorted(wildtype_pdb.ProperResidueIDToAAMap().iteritems()):
                c = resid[0]
                resnum = resid[1:].strip()
                if mutation.Chain == c and mutation.ResidueID == resnum:
                    if mutation.WildTypeAA == wtaa:
                        foundMatch = True
                    else:
                        foundAlternativeWildType = wtaa
            if not foundMatch and not(mutationAllowedToBeStoredDespiteMissingCoordinates):
                #raise colortext.Exception("%s: Could not find a match for mutation %s %s:%s -> %s in %s." % (pdbID, mutation[FieldNames.Chain], mutation[FieldNames.ResidueID], mutation[FieldNames.WildTypeAA], mutation[FieldNames.MutantAA], associatedRecordsStr ))
                if foundAlternativeWildType:
                    colortext.error("%s: Could not find a match for mutation %s. A different wildtype (%s) seems to be at that position." % (self.PDB_ID, mutation, foundAlternativeWildType))
                else:
                    colortext.error("%s: Could not find a match for mutation %s. No residue information seems available at that position." % (self.PDB_ID, mutation))
                failed = True

            # Sanity check that the mutant structure has the mutation.
            # Note: This assumes the residue numbering is identical which may not be the case. Otherwise, mappings must be provided in the database.
            for mutant_pdb in mutant_pdbs:
                for resid, mutantaa in sorted(mutant_pdb.ProperResidueIDToAAMap().iteritems()):
                    c = resid[0]
                    resnum = resid[1:].strip()
                    if mutation.Chain == c and mutation.ResidueID == resnum:
                        if mutation.MutantAA == mutantaa:
                            foundMatch = True
                        else:
                            foundAlternativeMutantType = mutantaa
                if not foundMatch and not(mutationAllowedToBeStoredDespiteMissingCoordinates):
                    #raise colortext.Exception("%s: Could not find a match for mutation %s %s:%s -> %s in %s." % (pdbID, mutation[FieldNames.Chain], mutation[FieldNames.ResidueID], mutation[FieldNames.WildTypeAA], mutation[FieldNames.MutantAA], associatedRecordsStr ))
                    if foundAlternativeMutantType:
                        colortext.error("%s: Could not find a match for mutation %s. A different mutant type (%s) seems to be at that position in the mutant %s." % (self.PDB_ID, mutation, foundAlternativeMutantType, mutant_pdb["PDB_ID"]))
                    else:
                        colortext.error("%s: Could not find a match for mutation %s. No residue information seems available at that position in the mutant %s." % (self.PDB_ID, mutation, mutant_pdb["PDB_ID"]))
                    failed = True

        # Sanity check that the chain information is correct (ProTherm has issues)
        chainsInPDB = WildTypeStructure.chains
        if not chainsInPDB:
            raise Exception("The chains for %s were not read in properly." % associatedRecordsStr)
        for c in self.chains:
            if not c in chainsInPDB:
                if len(chainsInPDB) == 1 and len(self.chains) == 1:
                    colortext.warning("%s: Chain '%s' does not exist in the PDB %s. PDB has one chain - %s. Use that chain instead." % (self.PDB_ID, c, self.PDB_ID, chainsInPDB[0]))
                    self.ddGdb.addChainWarning(pdbID, associatedRecords, c)
                    failed = True
                else:
                    self.ddGdb.addChainError(self.PDB_ID, c)
                    raise colortext.Exception("Error committing experiment:\n%s: Chain '%s' does not exist in the PDB %s. Chains %s exist." % (self.PDB_ID, c, self.PDB_ID, join(chainsInPDB, ", ")))
        if failed:
            raise colortext.Exception("Failed adding the record to the database.")

    def commit(self, mutationAllowedToBeStoredDespiteMissingCoordinates, testonly = False, pdbPath = None, quiet = False):
        '''Commits the set of experiments associated with the mutation to the database. Returns the unique ID of the associated Experiment record.'''
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames
        failed = False

        # Look for an existing record
        ExperimentID = super(ExperimentDefinition, self).commit()
        if ExperimentID:
            return ExperimentID

        #ExperimentID, dbMutations = self.find()
        #if ExperimentID != None:
        #	if not quiet:
        #		colortext.error("\nExperiment already exists in the database with Experiment ID=%s." % ExperimentID)
        #		colortext.warning("*** This record ***%s\n" % self)
        #		colortext.warning("*** Database record ***\n%s\n" % dbMutations)
        #	return ExperimentID

        # Test mutant homology, whether mutations exist etc.
        self.test(mutationAllowedToBeStoredDespiteMissingCoordinates = mutationAllowedToBeStoredDespiteMissingCoordinates)

        try:
            d = {
                FieldNames.Experiment.Structure : self.PDB_ID,
            }
            if not testonly:
                ddGdb.insertDict(FieldNames.Experiment._name, d)
                ExperimentID = ddGdb.getLastRowID()
                #print("ExperimentID", ExperimentID)
            self.databaseID = ExperimentID
        except Exception, e:
            raise Exception("\nError committing Experiment to database.\n***\n%s\n%s\n***" % (self, str(e)))

        if self.interface:
            try:
                d = {
                    FieldNames.ExperimentInterface.ExperimentID	: ExperimentID,
                    FieldNames.ExperimentInterface.Interface	: self.interface,
                }
                if not testonly:
                    ddGdb.insertDict(FieldNames.ExperimentInterface._name, d)
            except Exception, e:
                self.remove()
                raise Exception("\nError committing Experiment to database.\n***\n%s\n%s\n***" % (self, str(e)))

        for mutantmap in self.mutantmaps:
            try:
                d = {
                    FieldNames.ExperimentMutant.ExperimentID	: ExperimentID,
                    FieldNames.ExperimentMutant.Mutant			: mutantmap.mutant_PDB_ID,
                    FieldNames.ExperimentMutant.WildTypeChainID	: mutantmap.wildtype_chain,
                    FieldNames.ExperimentMutant.MutantChainID	: mutantmap.mutant_chain,
                    FieldNames.ExperimentMutant.Notes			: mutantmap.notes,
                }
                if not testonly:
                    ddGdb.insertDict(FieldNames.ExperimentMutant._name, d)
            except Exception, e:
                self.remove()
                raise Exception("\nError committing Experiment to database.\n***\n%s\n%s\n***" % (self, str(e)))

        for chain in self.chains:
            try:
                d = {
                    FieldNames.ExperimentChain.ExperimentID : ExperimentID,
                    FieldNames.ExperimentChain.Chain		: chain,
                }
                if not testonly:
                    ddGdb.insertDict(FieldNames.ExperimentChain._name, d)
            except Exception, e:
                self.remove()
                raise Exception("\nError committing Experiment to database.\n***\n%s\n%s\n***" % (self, str(e)))

        emFieldNames = self.ddGdb.FieldNames.ExperimentMutation
        for mutation in self.mutations:
            if not testonly:
                try:
                    d = {
                        emFieldNames.ExperimentID				: ExperimentID,
                        emFieldNames.Chain 						: mutation.Chain,
                        emFieldNames.ResidueID					: mutation.ResidueID,
                        emFieldNames.WildTypeAA					: mutation.WildTypeAA,
                        emFieldNames.MutantAA					: mutation.MutantAA,
                        emFieldNames.SecondaryStructurePosition	: mutation.SecondaryStructurePosition,
                        emFieldNames.AccessibleSurfaceArea		: mutation.AccessibleSurfaceArea,
                    }
                    ddGdb.insertDict(FieldNames.ExperimentMutation._name, d)
                except Exception, e:
                    self.remove()
                    raise Exception("\nError committing Experiment to database.\n***\n%s\n%s\n***" % (self, str(e)))

        return self.databaseID


class ExperimentAssayDefinition(DBObject):

    ExperimentalConditionsFields = [
        'Temperature', 'pH', 'Buffer', 'BufferConcentration',
        'Ion1', 'Ion1Concentration', 'Ion2', 'Ion2Concentration', 'Ion3', 'Ion3Concentration',
        'Additives', 'ProteinConcentration',
        'Measure1', 'Measure2', 'Measure3',
        'MethodOfDenaturation1','MethodOfDenaturation2'
    ]
    ThermodynamicFields = set(['DG', 'DG_H2O', 'Tm', 'dTm', 'dHvH', 'dHcal', 'm', 'Cm', 'dCp', 'activity', 'activity_Km', 'activity_Kcat', 'activity_Kd'])

    def __init__(self, ddGdb, SecondaryID, ExperimentID, Username, Publication = None, DuplicateOf = None):
        self.ddGdb = ddGdb
        self.DDGs = {}
        self.Thermodynamics = {}
        FieldNames = self.ddGdb.FieldNames.ExperimentAssay
        d = {}
        for k, v in FieldNames.__dict__.iteritems():
            if k != "_name" and k != "ID":
                d[k] = None
        self.d = d

        time_now = datetime.now()
        self.SecondaryID = SecondaryID
        self.ExperimentID = ExperimentID
        self.Publication = Publication
        self.DuplicateOf = DuplicateOf
        self.AddedBy = Username
        self.AddedDate = time_now
        self.LastModifiedBy = Username
        self.LastModifiedDate = time_now

    def __setattr__(self, k, v):
        if k == "ddGdb" or k == "d" or k == "DDGs" or k == "Thermodynamics" or k == "databaseID":
            super(ExperimentAssayDefinition, self).__setattr__(k, v)
        else:
            assert(k in self.d.keys())
            self.d[k] = v

    def __getattr__(self, k):
        if k == "ddGdb" or k == "d":
            super(ExperimentAssayDefinition, self).__getattr__(k)
        else:
            return self.d[k]

    def remove(self):
        ddGdb = self.ddGdb
        if not self.databaseID:
            raise Exception("Cannot remove a record with no corresponding database ID.")
        ddGdb.locked_execute("DELETE FROM ExperimentAssayDDG WHERE ExperimentAssayID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM ExperimentAssayThermodynamic WHERE ExperimentAssayID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM ExperimentAssay WHERE ID=%s", parameters = (self.databaseID,))

    def find(self):
        results = self.ddGdb.locked_execute('''SELECT ID FROM ExperimentAssay WHERE SecondaryID=%s''', parameters = (self.SecondaryID, ))
        if results:
            assert(len(results) == 1)
            return results[0]["ID"], self
        return None, None

    def __repr__(self):
        s = []
        s.append("SecondaryID: %s" % self.SecondaryID)
        s.append("ExperimentID: %s" % self.ExperimentID)
        s.append("Publication: %s" % self.Publication)
        s.append("DuplicateOf: %s" % self.DuplicateOf)
        if self.DDGs.get("DDG_H2O"):
            DDG_H2O = self.DDGs["DDG_H2O"]
            s.append("DDG_H2O: %s kcal/mol, Rosetta convention (published as %s %s)" % (DDG_H2O["Value"], DDG_H2O["PublishedValue"], DDG_H2O["PublishedUnit"]))
        if self.DDGs.get("DDG"):
            DDG = self.DDGs["DDG"]
            s.append("DDG: %s kcal/mol, Rosetta convention (published as %s %s)" % (DDG["Value"], DDG["PublishedValue"], DDG["PublishedUnit"]))
        if self.DDGs.get("Unknown"):
            someDDG = self.DDGs["Unknown"]
            s.append("DDG_H2O or DDG: %s kcal/mol, Rosetta convention (published as %s %s)" % (someDDG["Value"], someDDG["PublishedValue"], someDDG["PublishedUnit"]))
        for fieldname in ExperimentAssayDefinition.ExperimentalConditionsFields:
            if self.__getattr__(fieldname):
                s.append("%s: %s" % (fieldname, self.__getattr__(fieldname)))
        if len(self.Thermodynamics) > 0:
            s.append("Thermodynamics:")
            for fieldname in ExperimentAssayDefinition.ThermodynamicFields:
                if self.Thermodynamics.get(fieldname):
                    s.append("  %s: %s" % (fieldname, self.Thermodynamics.get(fieldname)["PublishedValue"]))
        s.append("Number of transition states: %s" % self.NumberOfTransitionStates)
        s.append("Reversibility: %s" % self.Reversibility)
        s.append("ReversibilityLevel: %s" % self.ReversibilityLevel)
        s.append("AddedBy: %s on %s" % (self.AddedBy, self.AddedDate))
        s.append("LastModifiedBy: %s on %s" % (self.LastModifiedBy, self.LastModifiedDate))
        return join(s, "\n")

    def addDDG(self, ddg_type, LocationOfValueInPublication, DDGValue, PublishedValue, PublishedUnit, PublishedError, NumberOfMeasurements = None, Remarks = None):
        '''
        DDGValue - Value is in kcal/mol using the Rosetta convention i.e. the lower the value, the more stabilizing the mutation.
        PublishedValue - Use the exact published value here even if wrong. Use the original sign convention and values in the original unit without conversion (unless conversion is necessary to convert back to the original unit).
        '''
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames.ExperimentAssayDDG

        assert(ddg_type == "DDG" or ddg_type == "DDG_H2O" or ddg_type == "Unknown")
        assert(not(self.DDGs.get(ddg_type)))
        assert(len(LocationOfValueInPublication) <= 96)
        assert(type(DDGValue) == type(0.1))
        assert(type(PublishedValue) == type(0.1))
        assert(PublishedUnit == "kcal/mol" or PublishedUnit == "kJ/mol" or PublishedUnit == "cal/mol")
        if NumberOfMeasurements:
            assert(type(NumberOfMeasurements) == type(1))

        self.DDGs[ddg_type] = {
            FieldNames.ExperimentAssayID			: None,
            FieldNames.Type							: ddg_type,
            FieldNames.LocationOfValueInPublication	: LocationOfValueInPublication,
            FieldNames.Value						: DDGValue,
            FieldNames.PublishedValue				: PublishedValue,
            FieldNames.PublishedUnit				: PublishedUnit,
            FieldNames.PublishedError				: PublishedError,
            FieldNames.NumberOfMeasurements			: NumberOfMeasurements,
            FieldNames.Remarks						: Remarks,
        }

    def addThermodynamic(self, t_type, PublishedValue, PublishedUnit = 'Unknown', PublishedError = None):
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames.ExperimentAssayThermodynamic

        assert(t_type in ExperimentAssayDefinition.ThermodynamicFields)
        assert(not(self.Thermodynamics.get(t_type)))
        assert(len(str(PublishedValue)) <= 256)
        # todo: Change the assertions below when we actually store the proper information
        assert(PublishedUnit == "Unknown")
        assert(PublishedError == None)

        self.Thermodynamics[t_type] = {
            FieldNames.ExperimentAssayID			: None,
            FieldNames.Type							: t_type,
            FieldNames.PublishedValue				: str(PublishedValue),
            FieldNames.PublishedUnit				: PublishedUnit,
            FieldNames.PublishedError				: PublishedError,
        }

    def test(self):
        pass

    def commit(self, UserID, testonly = False, pdbPath = None, quiet = False):
        '''Commits the experimental assay to the database. Returns the unique ID of the associated ExperimentAssay record.'''
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames
        failed = False

        # Look for an existing record
        ExperimentAssayID = super(ExperimentAssayDefinition, self).commit()
        if ExperimentAssayID:
            return ExperimentAssayID

        # Look for an existing record
        ExperimentAssayID, dbRecord = self.find()
        if ExperimentAssayID != None:
            if not quiet:
                colortext.error("\nExperimentAssay already exists in the database with ID=%s." % ExperimentAssayID)
                colortext.warning("*** This record ***%s\n" % self)
                colortext.warning("*** Database record ***\n%s\n" % dbRecord)
            return ExperimentAssayID

        self.test()

        ExperimentAssayID = None
        try:
            if not testonly:
                ddGdb.insertDict(FieldNames.ExperimentAssay._name, self.d)
                ExperimentAssayID = ddGdb.getLastRowID()
                #print("ExperimentAssayID", ExperimentAssayID)
            self.databaseID = ExperimentAssayID
        except Exception, e:
            raise Exception("\nError committing ExperimentAssay to database.\n***\n%s\n%s\n%s\n***" % (str(e), traceback.format_exc(), self))

        for DDGtype, ddGdict in self.DDGs.iteritems():
            try:
                if not testonly:
                    ddGdict["ExperimentAssayID"] = ExperimentAssayID
                    ddGdb.insertDict(FieldNames.ExperimentAssayDDG._name, ddGdict)
            except Exception, e:
                self.remove()
                raise Exception("\nError committing ExperimentAssayDDG to database.\n***\n%s\n%s\n***" % (self, str(e)))

        for thermoType, thermodict in self.Thermodynamics.iteritems():
            try:
                if not testonly:
                    thermodict["ExperimentAssayID"] = ExperimentAssayID
                    ddGdb.insertDict(FieldNames.ExperimentAssayThermodynamic._name, thermodict)
            except Exception, e:
                self.remove()
                raise Exception("\nError committing ExperimentAssayThermodynamic to database.\n***\n%s\n%s\n***" % (self, str(e)))

        return self.databaseID


class DataSetDDG(DBObject):

    def __init__(self, ddGdb, DataSetID, Section, RecordNumber, PublishedValue, PDB_ID, PublishedPDB_ID = None, AggregateType = 'SingleValue', MutationIsReversed = False, PossibleError = False, Remark = None, CorrectionRemark = None):
        self.ddGdb = ddGdb
        self.DDGSources = {}

        FieldNames = self.ddGdb.FieldNames.DataSetDDG
        d = {}
        for k, v in FieldNames.__dict__.iteritems():
            if k != "_name" and k != "ID":
                d[k] = None
        self.d = d

        assert(len(DataSetID) <= 128)
        assert(len(Section) <= 64)
        assert(type(RecordNumber) == type(1))
        assert(AggregateType == "SingleValue" or AggregateType == "MeanValue")
        assert(type(PublishedValue) == type(0.1))
        assert(type(MutationIsReversed) == type(True))
        assert(len(PDB_ID) <= 10)
        if not PublishedPDB_ID:
            PublishedPDB_ID = PDB_ID
        assert(len(PublishedPDB_ID) <= 10)
        assert(type(PossibleError) == type(True))

        self.DataSetID = DataSetID
        self.Section = Section
        self.RecordNumber = RecordNumber
        self.AggregateType = AggregateType
        self.PublishedValue = PublishedValue
        self.MutationIsReversed = MutationIsReversed
        self.PDB_ID = PDB_ID
        self.PublishedPDB_ID = PublishedPDB_ID
        self.PossibleError = PossibleError
        self.Remark = Remark
        self.CorrectionRemark = CorrectionRemark

    def __setattr__(self, k, v):
        if k == "ddGdb" or k == "d" or k == "DDGSources" or k == "databaseID":
            super(DataSetDDG, self).__setattr__(k, v)
        else:
            assert(k in self.d.keys())
            self.d[k] = v

    def __getattr__(self, k):
        if k == "ddGdb" or k == "d":
            super(DataSetDDG, self).__getattr__(k)
        else:
            return self.d[k]

    def remove(self):
        ddGdb = self.ddGdb
        if not self.databaseID:
            raise Exception("Cannot remove a record with no corresponding database ID.")
        ddGdb.locked_execute("DELETE FROM DataSetDDGSource WHERE DataSetDDGID=%s", parameters = (self.databaseID,))
        ddGdb.locked_execute("DELETE FROM DataSetDDG WHERE ID=%s", parameters = (self.databaseID,))

    def find(self):
        results = self.ddGdb.locked_execute('''SELECT ID FROM DataSetDDG WHERE DataSetID=%s AND Section=%s and RecordNumber=%s''', parameters = (self.DataSetID, self.Section, self.RecordNumber))
        if results:
            assert(len(results) == 1)
            return results[0]["ID"], self
        return None, None

    def __repr__(self):
        s = []
        s.append("DataSetID: %s" % self.DataSetID)
        s.append("Section: %s" % self.Section)
        s.append("RecordNumber: %s" % self.RecordNumber)
        if self.AggregateType != 'SingleValue':
            assert(len(self.DDGSources) > 1)
            s.append("AggregateType: %s relating to the experimental assays with the IDs %s" % (self.AggregateType, join(map(str, sorted(self.DDGSources)), ", ")))
        else:
            assert(len(self.DDGSources) == 1)
            s.append("Relates to the experimental assays with ID %s" % list(self.DDGSources.keys())[0])
            for ExperimentAssayID, e_type in sorted(self.DDGSources.iteritems()):
                s.append("%s: %s" % (ExperimentAssayID, e_type))
        s.append("PublishedValue: %s" % self.PublishedValue)
        s.append("PDB_ID: %s" % self.PDB_ID)
        if self.PublishedPDB_ID != self.PDB_ID:
            s.append("PublishedPDB_ID: %s" % self.PublishedPDB_ID)
        if self.Remark:
            s.append(colortext.make("Remark: %s" % self.Remark, color="cyan"))
        if self.CorrectionRemark:
            s.append(colortext.make("CorrectionRemark: %s" % self.CorrectionRemark, color="lightblue"))
        if self.MutationIsReversed:
            s.append(colortext.make("[*The mutation is reversed from the published experimental assay*]", color = "orange"))
        if self.PossibleError:
            s.append(colortext.make("[*This record is possibly erroneous*]", color = "pink"))
        return join(s, "\n")

    def addDDGSource(self, ExperimentAssayID, e_type):
        ddGdb = self.ddGdb
        assert(not(self.DDGSources.get(ExperimentAssayID)))
        assert(e_type == 'Unknown' or e_type == 'DDG' or e_type == 'DDG_H2O')
        self.DDGSources[ExperimentAssayID] = e_type
        if self.AggregateType == 'SingleValue':
            assert(len(self.DDGSources) == 1)

    def test(self):
        pass

    def commit(self, testonly = False, quiet = False):
        '''Commits the DataSet DDG record to the database. Returns the unique ID of the associated DataSetDDG record.'''
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames
        failed = False

        # Look for an existing record
        DataSetDDGID = super(DataSetDDG, self).commit()
        if DataSetDDGID:
            return DataSetDDGID

        # Look for an existing record
        DataSetDDGID, dbRecord = self.find()
        if DataSetDDGID != None:
            if not quiet:
                colortext.error("\nDataSetDDG already exists in the database with ID=%s." % DataSetDDGID)
                colortext.warning("*** This record ***%s\n" % self)
                colortext.warning("*** Database record ***\n%s\n" % dbRecord)
            return DataSetDDGID

        self.test()

        DataSetDDGID = None
        try:
            if not testonly:
                ddGdb.insertDict(FieldNames.DataSetDDG._name, self.d)
                DataSetDDGID = ddGdb.getLastRowID()
            else:
                print(self.d)
            self.databaseID = DataSetDDGID
        except Exception, e:
            raise Exception("\nError committing DataSetDDG to database.\n***\n%s\n%s\n%s\n***" % (str(e), traceback.format_exc(), self))

        for DDGSource, e_type in self.DDGSources.iteritems():
            try:
                t = {
                    FieldNames.DataSetDDGSource.DataSetDDGID		: DataSetDDGID,
                    FieldNames.DataSetDDGSource.ExperimentAssayID	: DDGSource,
                    FieldNames.DataSetDDGSource.Type				: e_type,
                }
                if not testonly:
                    ddGdb.insertDict(FieldNames.DataSetDDGSource._name, t)
                else:
                    print(t)
            except Exception, e:
                self.remove()
                raise Exception("\nError committing DataSetDDGSource to database.\n***\n%s\n%s\n***" % (self, str(e)))

        return self.databaseID





class Prediction(DBObject):

    def __init__(self, ddGdb, ExperimentID, PredictionSet, ProtocolID, ddG, status, NumberOfMeasurements = 1):
        self.ddGdb = ddGdb
        FieldNames_ = ddGdb.FlatFieldNames
        self.dict = {
            FieldNames_.ExperimentID		: ExperimentID,
            FieldNames_.PredictionSet		: PredictionSet,
            FieldNames_.ProtocolID			: ProtocolID,
            FieldNames_.KeptHETATMLines		: None,
            FieldNames_.StrippedPDB			: None,
            FieldNames_.ResidueMapping		: None,
            FieldNames_.InputFiles			: {},
            FieldNames_.Description			: {},
            FieldNames_.ddG					: ddG,
            FieldNames_.NumberOfMeasurements: NumberOfMeasurements,
            FieldNames_.Status				: status,
            FieldNames_.ExtraParameters		: pickle.dumps({}),
        }
        if ExperimentID == None:
            raise Exception("Cannot create the following Prediction - Missing ExperimentID:\n***\n%s\n***" % self)

    def setOptional(self, KeptHETATMLines = None, StrippedPDB = None, ResidueMapping = None, InputFiles = None, Description = None):
        d = self.dict
        if KeptHETATMLines:
            d[FieldNames_.KeptHETATMLines] = KeptHETATMLines
        if StrippedPDB:
            d[FieldNames_.StrippedPDB] = StrippedPDB
        if ResidueMapping:
            d[FieldNames_.ResidueMapping] = ResidueMapping
        if InputFiles:
            d[FieldNames_.InputFiles] = InputFiles
        if Description:
            d[FieldNames_.Description] = Description

    def commit(self):
        d = self.dict
        d[FieldNames_.InputFiles] = pickle.dumps(d[FieldNames_.InputFiles])
        d[FieldNames_.Description] = pickle.dumps(d[FieldNames_.Description])
        fields = [FieldNames_.ExperimentID, FieldNames_.PredictionSet, FieldNames_.ProtocolID, FieldNames_.KeptHETATMLines,
                FieldNames_.StrippedPDB, FieldNames_.ResidueMapping, FieldNames_.InputFiles, FieldNames_.Description,
                FieldNames_.ddG, FieldNames_.NumberOfMeasurements, FieldNames_.Status, FieldNames_.ExtraParameters]
        try:
            self.ddGdb.insertDict('Prediction', d, fields)
        except Exception, e:
            raise Exception("\nError committing prediction to database.\n***\n%s\n%s\n***" % (self, str(e)))
        self.databaseID = self.ddGdb.getLastRowID()
        return self.databaseID

    def __repr__(self):
        raise Exception('''This is unlikely to work as I have not tested it in a while.''')
        d = self.dict
        str = []
        str.append("%s: %s" % (FieldNames_.ExperimentID, d[FieldNames_.ExperimentID]))
        str.append("%s: %s" % (FieldNames_.PredictionSet, d[FieldNames_.PredictionSet]))
        str.append("%s: %d" % (FieldNames_.ProtocolID, d[FieldNames_.ProtocolID]))
        if d[FieldNames_.KeptHETATMLines] == None:
            str.append("%s: NULL" % (FieldNames_.KeptHETATMLines))
        else:
            str.append("%s: %d" % (FieldNames_.KeptHETATMLines, d[FieldNames_.KeptHETATMLines]))
        n = d[FieldNames_.NumberOfMeasurements]
        if n > 1:
            str.append("%s: %0.2f (%d measurements)" % (FieldNames_.ddG, d[FieldNames_.ddG], n))
        else:
            str.append("%s: %0.2f" % (FieldNames_.ddG, d[FieldNames_.ddG]))

        str.append("%s:" % (FieldNames_.InputFiles))
        if d[FieldNames_.InputFiles]:
            ifiles = d[FieldNames_.InputFiles]
            if type(ifiles) == type(""):
                ifiles = pickle.loads(ifiles)
            for k,v in ifiles.iteritems():
                str.append("\t%s" % k)
        else:
            str.append("\tEmpty")
        str.append("%s:" % (FieldNames_.Description))
        if d[FieldNames_.Description]:
            idesc = d[FieldNames_.Description]
            if type(idesc) == type(""):
                idesc = pickle.loads(idesc)
            for k,v in idesc.iteritems():
                str.append("\t%s: %s" % (k, v))
        else:
            str.append("\tEmpty")

        return join(str, "\n")

    def __getitem__(self, key):
        return dict_[key]

class DataSet(DBObject):

    def __init__(self, ddGdb, ID, ShortID, Description, CreationDateStart, CreationDateEnd, DDGConvention, UserID = None):
        self.ddGdb = ddGdb
        FieldNames = ddGdb.FieldNames.DataSet
        self.dict = {
            FieldNames.ID					: ID,
            FieldNames.ShortID				: ShortID,
            FieldNames.UserID				: UserID,
            FieldNames.Description			: Description,
            FieldNames.CreationDateStart	: CreationDateStart,
            FieldNames.CreationDateEnd		: CreationDateEnd,
            FieldNames.DDGConvention		: DDGConvention,
        }
        self.References = []

    def addReference(self, SourceID):
        self.References.append(SourceID)

    def find(self):
        results = self.ddGdb.locked_execute('''SELECT ID FROM DataSet WHERE ID=%s''', parameters = (self.dict[self.ddGdb.FieldNames.DataSet.ID]))
        if results:
            assert(len(results) == 1)
            return results[0]["ID"], self
        return None, None

    def commit(self):
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames
        failed = False

        # Look for an existing record
        DataSetID = super(DataSet, self).commit()
        if DataSetID:
            return DataSetID

        # Look for an existing record
        DataSetID, dbRecord = self.find()
        if DataSetID != None:
            if not quiet:
                colortext.error("\nDataSet already exists in the database with ID=%s." % DataSetID)
                colortext.warning("*** This record ***%s\n" % self)
                colortext.warning("*** Database record ***\n%s\n" % dbRecord)
            return DataSetID

        d = self.dict
        IDfield = FieldNames.DataSet.ID
        DataSetID = None
        try:
            DataSetID = self.ddGdb.insertDictIfNew(FieldNames.DataSet._name, d, IDfield)[1][IDfield]
        except Exception, e:
            raise Exception("\nError committing DataSet to database.\n***\n%s\n%s\n***" % (self, str(e)))

        self.databaseID = self.ddGdb.getLastRowID()
        for reference in self.References:
            try:
                d = {
                    FieldNames.DataSetReference.DataSetID	: DataSetID,
                    FieldNames.DataSetReference.Publication	: reference
                }
                self.ddGdb.insertDictIfNew(FieldNames.DataSetReference._name, d, [FieldNames.DataSetReference.DataSetID, FieldNames.DataSetReference.Publication])
            except Exception, e:
                raise Exception("\nError committing DataSetReference to database.\n***\n%s\n%s\n***" % (self, str(e)))

        return self.databaseID

    def __repr__(self):
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames.DataSet

        s = []
        s.append("ID: %s" % self.dict[FieldNames.ID])
        s.append("ShortID: %s" % self.dict[FieldNames.ShortID])
        s.append("UserID: %s" % self.dict[FieldNames.UserID])
        s.append("Description: %s" % self.dict[FieldNames.Description])
        s.append("CreationDateStart: %s" % self.dict[FieldNames.CreationDateStart])
        s.append("CreationDateEnd: %s" % self.dict[FieldNames.CreationDateEnd])
        s.append("DDGConvention: %s" % self.dict[FieldNames.DDGConvention])
        if self.References:
            s.append("References:")
            for reference in self.References:
                s.append("  %s" % reference)
        return join(s, "\n")

    def __getitem__(self, key):
        return dict_[key]


class Publication(DBObject):

    def __init__(self, ddGdb, ID, DDGUnit = None, DDGConvention = None, Notes = None, DGNotes = None, DGUnitUsedInProTherm = None, DDGProThermSignNotes = None, DDGValuesNeedToBeChecked = None, RIS = None):
        self.ddGdb = ddGdb
        FieldNames = ddGdb.FieldNames.Source
        self.dict = {
            FieldNames.ID							: ID,
            FieldNames.DGUnit						: DDGUnit,
            FieldNames.DDGConvention				: DDGConvention,
            FieldNames.Notes						: Notes,
            FieldNames.DGNotes						: DGNotes,
            FieldNames.DGUnitUsedInProTherm			: DGUnitUsedInProTherm,
            FieldNames.DDGProThermSignNotes			: DDGProThermSignNotes,
            FieldNames.DDGValuesNeedToBeChecked		: DDGValuesNeedToBeChecked,
            FieldNames.RIS							: RIS,
        }
        self.IDs = []
        self.DDGLocations = []

    def addIdentifier(self, IDType, ID):
        self.IDs.append((ID, IDType))

    def addDDGValueLocation(self, Location, Notes):
        self.IDs.append((Location, Notes))

    def commit(self):
        d = self.dict
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames
        IDfield = FieldNames.Source.ID
        SourceID = None
        try:
            SourceID = self.ddGdb.insertDictIfNew(FieldNames.Source._name, d, [IDfield])[1][IDfield]
        except Exception, e:
            raise Exception("\nError committing Source to database.\n***\n%s\n%s\n***" % (self, str(e)))

        self.databaseID = self.ddGdb.getLastRowID()

        for ID in self.IDs:
            try:
                d = {
                    FieldNames.SourceIdentifier.SourceID	: SourceID,
                    FieldNames.SourceIdentifier.ID 			: ID[0],
                    FieldNames.SourceIdentifier.Type		: ID[1],
                }
                self.ddGdb.insertDictIfNew(FieldNames.SourceIdentifier._name, d, [FieldNames.SourceIdentifier.SourceID, FieldNames.SourceIdentifier.ID])
            except Exception, e:
                raise Exception("\nError committing SourceIdentifier to database.\n***\n%s\n%s\n***" % (self, str(e)))

        for Location in self.DDGLocations:
            try:
                d = {
                    FieldNames.SourceDDGValueLocation.SourceID	: SourceID,
                    FieldNames.SourceDDGValueLocation.Location 	: Location[0],
                    FieldNames.SourceDDGValueLocation.Notes		: Location[1],
                }
                self.ddGdb.insertDictIfNew(FieldNames.SourceDDGValueLocation._name, d, [FieldNames.SourceDDGValueLocation.SourceID, FieldNames.SourceDDGValueLocation.Location])
            except Exception, e:
                raise Exception("\nError committing SourceDDGValueLocation to database.\n***\n%s\n%s\n***" % (self, str(e)))

        return self.databaseID

    def __repr__(self):
        ddGdb = self.ddGdb
        FieldNames = ddGdb.FieldNames.Source

        s = []
        s.append("ID: %s" % self.dict[FieldNames.ID])
        if self.dict[FieldNames.DGUnit]:
            s.append("DGUnit: %s" % self.dict[FieldNames.DGUnit])
        if self.dict[FieldNames.DDGConvention]:
            s.append("DDGConvention: %s" % self.dict[FieldNames.DDGConvention])
        if self.dict[FieldNames.Notes]:
            s.append("Notes: %s" % self.dict[FieldNames.Notes])
        if self.dict[FieldNames.DGNotes]:
            s.append("DGNotes: %s" % self.dict[FieldNames.DGNotes])
        if self.dict[FieldNames.DGUnitUsedInProTherm]:
            s.append("DGUnitUsedInProTherm: %s" % self.dict[FieldNames.DGUnitUsedInProTherm])
        if self.dict[FieldNames.DDGProThermSignNotes]:
            s.append("DDGProThermSignNotes: %s" % self.dict[FieldNames.DDGProThermSignNotes])
        if self.dict[FieldNames.DDGValuesNeedToBeChecked]:
            s.append("DDGValuesNeedToBeChecked: Yes")
        s.append("RIS: %s" % self.dict[FieldNames.RIS])
        if self.IDs:
            s.append("IDs:")
            for ID in self.IDs:
                s.append("  %s: %s" % (ID[1], ID[0]))
        if self.DDGLocations:
            s.append("Locations:")
            for Location in self.DDGLocations:
                s.append("  %s: %s" % (Location))
        return join(s, "\n")

    def __getitem__(self, key):
        return dict_[key]

class ddGPredictionDataDatabase(object):

    def __init__(self, passwd = None):
        if not passwd:
            if os.path.exists("pw"):
                F = open("pw")
                passwd = F.read().strip()
                F.close()
            else:
                passwd = getpass.getpass("Enter password to connect to MySQL database:")

        self.passwd = passwd
        self.connectToServer()
        self.numTries = 32
        self.lastrowid = None

    def close(self):
        self.connection.close()

    def connectToServer(self):
        print("[CONNECTING TO SQL SERVER]")
        self.connection = MySQLdb.Connection(host = "kortemmelab.ucsf.edu", db = "ddGPredictionData", user = "kortemmelab", passwd = self.passwd, port = 3306, unix_socket = "/var/lib/mysql/mysql.sock")

    def execute(self, sql, parameters = None, cursorClass = MySQLdb.cursors.DictCursor, quiet = False):
        """Execute SQL query. This uses DictCursor by default."""
        i = 0
        errcode = 0
        caughte = None
        while i < self.numTries:
            i += 1
            try:
                cursor = self.connection.cursor(cursorClass)
                if parameters:
                    errcode = cursor.execute(sql, parameters)
                else:
                    errcode = cursor.execute(sql)
                self.connection.commit()
                results = cursor.fetchall()
                self.lastrowid = int(cursor.lastrowid)
                cursor.close()
                return results
            except MySQLdb.OperationalError, e:
                caughte = str(e)
                errcode = e[0]
                if errcode == 2006 or errcode == 2013:
                    self.connectToServer()
                self.connection.ping()
                continue
            except Exception, e:
                caughte = str(e)
                traceback.print_exc()
                break

        if not quiet:
            sys.stderr.write("\nSQL execution error in query %s at %s:" % (sql, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            sys.stderr.write("\nErrorcode/Error: %d - '%s'.\n" % (errcode, str(caughte)))
            sys.stderr.flush()
        raise MySQLdb.OperationalError(caughte)


class DatabaseMissingKeyException(Exception):
    def __init__(self, key, tbl):
        self.missing_key = key
        self.table = tbl
    def __str__(self):
        return ("Missing database key '%s' for table %s." % (str(self.missing_key), self.table))

class DatabaseBadDataException(Exception):
    def __init__(self, tbl, msg = ""):
        self.table = tbl
        self.msg = msg
    def __str__(self):
        if self.msg:
            return "A consistency check on table %s failed:\n%s" % (self.table, self.msg)
        else:
            return "A consistency check on table %s failed." % self.table

class DatabaseCannotDeleteRecordException(Exception):
    def __init__(self, id, tbl):
        self.id = id
        self.table = tbl
    def __str__(self):
        return "The records associated with ID '%s' in table %s could not be deleted." % (self.id, self.table)

class ddGDatabase(DatabaseInterface):

    chainErrors = {}
    chainWarnings= {}

    def __init__(self, passwd = None):
        if not passwd:
            if os.path.exists("pw"):
                F = open("pw")
                passwd = F.read().strip()
                F.close()
            else:
                passwd = getpass.getpass("Enter password to connect to MySQL database:")
        self.passwd = passwd

        super(ddGDatabase, self).__init__({},
            isInnoDB = True,
            numTries = 32,
            host = "kortemmelab.ucsf.edu",
            db = "ddG",
            user = "kortemmelab",
            passwd = passwd,
            port = 3306,
            unix_socket = "/var/lib/mysql/mysql.sock")

    def _create_protein_deletion_stored_procedure(self):
        '''This stored procedure returns 1 on error, -1 when there was no associated Protein record, and 0 on success.'''

        self.execute('DROP PROCEDURE IF EXISTS _DELETE_PROTEIN')
        if '_DELETE_PROTEIN' not in self.list_stored_procedures():
            self.execute("""
CREATE PROCEDURE _DELETE_PROTEIN(IN protein_ID VARCHAR(18), OUT error_code TINYINT)
BEGIN

DECLARE number_of_deleted_rows INT;

DECLARE EXIT HANDLER FOR SQLEXCEPTION ROLLBACK;
DECLARE EXIT HANDLER FOR SQLWARNING ROLLBACK;

START TRANSACTION;
SET error_code = 1;

DELETE FROM ProteinResidue WHERE ProteinID=protein_ID;
DELETE FROM ProteinName WHERE ProteinID=protein_ID;
DELETE FROM ProteinOrganism WHERE ProteinID=protein_ID;
DELETE FROM ProteinSegment WHERE ProteinID=protein_ID;
DELETE FROM Protein WHERE ID=protein_ID;

SET number_of_deleted_rows = ROW_COUNT();

IF number_of_deleted_rows > 0 THEN
    SET error_code = 0;
ELSE
    SET error_code = -1;
END IF;

COMMIT;

END
""", allow_unsafe_query = True)

    def _delete_protein(self, ProteinID):
        '''Returns True if deletion was successful and False if no records were deleted.
        An exception is raised if deletion was not possible.'''
        if '_DELETE_PROTEIN' not in self.list_stored_procedures():
            self._create_protein_deletion_stored_procedure()
        results = self.callproc('_DELETE_PROTEIN', parameters=(ProteinID, '@return_value'))
        assert(len(results) == 1)
        if int(results[0]['return_value']) > 0:
            raise DatabaseCannotDeleteRecordException(ProteinID, 'Protein')
        return int(results[0]['return_value']) == 0

    def _add_protein_residues(self, ProteinID, update = False):
        results = self.execute_select("SELECT * FROM Protein WHERE ID=%s", parameters=(ProteinID,))
        if len(results) != 1:
            raise DatabaseMissingKeyException(ProteinID, 'Protein')
        sequence = results[0]['Sequence']

        # Sanity check stored data
        results = self.execute_select("SELECT COUNT(ResidueID) AS ResidueCount FROM ProteinResidue WHERE ProteinID=%s", parameters=(ProteinID,))
        assert(len(results) == 1)
        if results[0]['ResidueCount'] > 0:
            if results[0]['ResidueCount'] == len(sequence):
                results = self.execute_select("SELECT ResidueAA FROM ProteinResidue WHERE ProteinID=%s ORDER BY ResidueID", parameters=(ProteinID,))
                stored_sequence = "".join([r['ResidueAA'] for r in results])
                assert(stored_sequence == sequence)
                if not(update):
                    return
            else:
                raise DatabaseBadDataException('ProteinResidue', 'Expected %d results, got %d results' % (len(sequence), results[0]['ResidueCount']))

        # Delete if updating
        if update:
            results = self.execute("DELETE FROM ProteinResidue WHERE ProteinID=%s", parameters=(ProteinID,))

        for c in range(len(sequence)):
            x = sequence[c]
            if x not in relaxed_amino_acid_codes: # Allow X for some proteins e.g. UPI000012EE21 / P00346
                raise Exception("Unknown amino acid '%s' at position %d of the sequence of protein %s." % (x, c + 1, ProteinID))
            assert(x in relaxed_amino_acid_codes)

        for c in range(len(sequence)):
            d = {
                'ProteinID': ProteinID,
                'ResidueID': c + 1,
                'ResidueAA': sequence[c],
                }
            self.insertDictIfNew('ProteinResidue', d, ['ProteinID', 'ResidueID'])

        # Sanity check again
        stored_sequence = "".join([r['ResidueAA'] for r in self.execute_select("SELECT ResidueAA FROM ProteinResidue WHERE ProteinID=%s ORDER BY ResidueID", parameters=(ProteinID,))])
        assert(stored_sequence == sequence)

    def look_for_protein_sequence(self, sequence):
        digest = CRC64.CRC64digest(sequence)
        results = self.execute_select("SELECT * FROM Protein WHERE CRC_64_ISO_Digest=%s", parameters=(digest,))
        for r in results:
            if r['Sequence'] == sequence:
                return r
        return None

    def add_raw_protein_sequence(self, sequence, IDScheme, IsASegmentOf = None, AddProteinResidues = True):

        # Determine the ID prefix
        ID_prefix = None
        if IDScheme == 'Kortemme Lab':
            ID_prefix = 'KOR'
        elif IDScheme == 'UniParcSegment':
            ID_prefix = 'SEG'
        else:
            raise Exception("The ID scheme %s is not recognized." % IDScheme)

        # Sanity check the sequence, allowing unknown residues 'X'
        for x in sequence:
            assert(x in relaxed_amino_acid_codes)

        # Sanity check IsASegmentOf
        if IsASegmentOf:
            results = self.execute_select("SELECT Sequence FROM Protein WHERE ID=%s", parameters=(IsASegmentOf,))
            if not results:
                raise DatabaseMissingKeyException(IsASegmentOf, 'Protein')
            else:
                assert(len(results) == 1)
                if results[0]['Sequence'].find(sequence) == -1:
                    raise Exception("The protein with sequence\n%s\ndoes not seem to be a segment of protein %s with sequence:\n%s" % (sequence, IsASegmentOf, results[0]['Sequence']))

        # Sanity check existing records
        ProteinID = None
        existing_sequence = self.look_for_protein_sequence(sequence)
        if existing_sequence:
            # We already have this sequence stored in the database
            ProteinID = existing_sequence['ID']
            if existing_sequence['IDScheme'] == 'UniParc':
                raise Exception("You are trying to add a raw protein sequence but an existing UniParc record with the same sequence exists with ID %s." % ProteinID)
            assert(existing_sequence['UniParcID'] == None)
            assert(existing_sequence['Sequence'] == sequence)
            # We do not need to check the CRC64 digest
            assert(existing_sequence['Mass'] == None)
            assert(existing_sequence['IsASegmentOf'] == IsASegmentOf)
        else:
            # Create a ProteinID if none exist. These are increasing integers just like the UniParc ID scheme.
            if not ProteinID:
                ProteinID = 1
                results = self.execute_select("SELECT ID FROM Protein WHERE ID LIKE '%s%%' ORDER BY ID DESC" % ID_prefix)
                if results:
                    ProteinID = int(results[0]['ID'][3:], base=16) + 1
                ProteinID = "%s%015x" % (ID_prefix, ProteinID)
                assert(len(ProteinID) == 18)

            # Protein record
            d = {
                'ID': ProteinID,
                'IDScheme': IDScheme,
                'UniParcID': None,
                'Sequence': sequence,
                'CRC_64_ISO_Digest': CRC64.CRC64digest(sequence),
                'Mass': None,
                'IsASegmentOf': IsASegmentOf,
                }
            self.insertDictIfNew('Protein', d, ['ID'])

        # ProteinResidue records
        if AddProteinResidues:
            self._add_protein_residues(ProteinID, False)
        return ProteinID

    def add_protein_from_UniParc_ID(self, UniParcID, cache_dir = None):
        '''Adds a Protein record and related records for the UniParc sequence with ID UniParcID. Returns the ID of the protein.'''

        uniparco = uniprot.UniParcEntry(UniParcID, cache_dir = cache_dir)

        # Sanity check the sequence, allowing unknown residues 'X'
        for x in uniparco.sequence:
            assert(x in relaxed_amino_acid_codes)

        # Sanity check existing records
        existing_sequence = self.look_for_protein_sequence(uniparco.sequence)
        if existing_sequence:
            # We already have this sequence stored in the database
            assert(existing_sequence['ID'] == UniParcID)
            assert(existing_sequence['IDScheme'] == 'UniParc')
            assert(existing_sequence['UniParcID'] == int(UniParcID[3:], 16))
            assert(existing_sequence['Sequence'] == uniparco.sequence)
            assert(existing_sequence['CRC_64_ISO_Digest'] == uniparco.CRC64Digest)
            assert(existing_sequence['Mass'] == uniparco.atomic_mass)
            assert(existing_sequence['IsASegmentOf'] == None)
        else:
            # Protein record
            d = {
                'ID': UniParcID,
                'IDScheme': 'UniParc',
                'UniParcID': int(UniParcID[3:], 16),
                'Sequence': uniparco.sequence,
                'CRC_64_ISO_Digest': uniparco.CRC64Digest,
                'Mass': uniparco.atomic_mass,
                'IsASegmentOf': None,
                }
            self.insertDictIfNew('Protein', d, ['ID'])

        for protein_segment in uniparco.subsections.sections:
            d = protein_segment.to_db()

            SegmentProteinID = None
            if d['StartResidue'] == 1 and d['EndResidue'] == len(uniparco.sequence):
                subsequence = uniparco.sequence
                SegmentProteinID = UniParcID
            else:
                subsequence = uniparco.sequence[d['StartResidue'] - 1:d['EndResidue']]
                SegmentProteinID = self.add_raw_protein_sequence(subsequence, 'UniParcSegment', IsASegmentOf = None, AddProteinResidues = False)

            d['ProteinID'] = UniParcID
            d['DefinedBy'] = 'UniProt'
            d['SegmentProteinID'] = SegmentProteinID
            self.insertDictIfNew('ProteinSegment', d, ['ProteinID', 'StartResidue', 'EndResidue'])

        # Fill in the Occurrence fields for each segment
        db_segments = self.execute("SELECT * FROM ProteinSegment WHERE ProteinID=%s ORDER BY StartResidue", parameters=(UniParcID,))
        occurrences = {}
        for db_segment in db_segments:
            occurrences[db_segment['SegmentProteinID']] = occurrences.get(db_segment['SegmentProteinID'], 0) + 1
            occurrence = occurrences[db_segment['SegmentProteinID']]
            self.execute("UPDATE ProteinSegment SET Occurrence=%s WHERE ProteinID=%s AND StartResidue=%s AND EndResidue=%s", parameters=(occurrence, UniParcID, db_segment['StartResidue'], db_segment['EndResidue']))

        for AC, names in uniparco.organisms.iteritems():
            d = names
            d['UniProt_ACC'] = AC
            self.insertDictIfNew('ProteinOrganism', d, ['UniProt_ACC', 'scientific'])

        # ProteinName record
        assert(uniparco.recommended_name)
        ECNumber = None
        if len(uniparco.recommended_name['EC numbers']) == 1:
            ECNumber = uniparco.recommended_name['EC numbers'][0]
        ProteinName = uniparco.recommended_name['Name']
        d = {
            'ProteinID' : UniParcID,
            'NameOrder' : 0,
            'Name' : ProteinName,
            'ECNumber' : ECNumber,
            'Validity' : None,
        }
        self.insertDictIfNew("_ProteinName", {'ProteinName' : ProteinName}, ['ProteinName'])
        self.insertDictIfNew("ProteinName", d, ['ProteinID', 'NameOrder'])

        # ProteinDatabaseIdentifier records
        for UniProtAC in uniparco.UniProtACs:
            d = {
                'ProteinID' : UniParcID,
                'Scheme'    : 'UniProt_ACC',
                'SchemeID'  : UniProtAC,
            }
            self.insertDictIfNew('ProteinDatabaseIdentifier', d, ['ProteinID', 'Scheme', 'SchemeID'])

        for UniProtID in uniparco.UniProtIDs:
            d = {
                'ProteinID' : UniParcID,
                'Scheme'    : 'UniProt_ID',
                'SchemeID'  : UniProtID,
            }
            self.insertDictIfNew('ProteinDatabaseIdentifier', d, ['ProteinID', 'Scheme', 'SchemeID'])

        # ProteinResidue records
        self._add_protein_residues(UniParcID, False)


    def addChainWarning(self, pdbID, associatedRecords, c):
        chainWarnings = self.chainWarnings
        chainWarnings[pdbID] = chainWarnings.get(pdbID) or []
        chainWarnings[pdbID].append((associatedRecords, c))

    def addChainError(self, pdbID, c):
        chainErrors = self.chainErrors
        chainErrors[pdbID] = chainErrors.get(pdbID) or []
        chainErrors[pdbID].append(c)

    def addTechniquesFields(self):
        '''Used to update missing Techniques fields as this field was added after the initial PDB import.'''
        return
        results = self.locked_execute("SELECT * FROM Structure")
        for result in results:
            pdbID = result[FieldNames_.PDB_ID]
            contents = result[FieldNames_.Content]
            lines = contents.split("\n")
            for line in lines:
                if line.startswith("EXPDTA"):
                    techniques = line[10:71].split(";")
                    for k in range(len(techniques)):
                        techniques[k] = techniques[k].strip()
                    techniques = join(techniques, ";")
                    break
            if not result[FieldNames_.Techniques]:
                SQL = "UPDATE Structure SET %s" % FieldNames_.Techniques
                SQL += "=%s WHERE PDB_ID=%s"
                self.locked_execute(SQL, parameters = (techniques, pdbID))

    def REMOVE_THIS_FUNCTION_insert(self, table, fieldnames, values):
        try:
            sql = None
            valuestring = join(["%s" for field in fieldnames], ", ")
            sql = "INSERT INTO %s (%s) VALUES (%s)" % (table, join(fieldnames, ", "), valuestring)
            #print(sql, values)
            if not len(fieldnames) == len(values):
                raise Exception("Fieldnames and values lists are not of equal size.")
            return
            self.locked_execute(sql, parameters)
        except Exception, e:
            if sql:
                sys.stderr.write("\nSQL execution error in query '%s' %% %s at %s:" % (sql, values, datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            sys.stderr.write("\nError: '%s'.\n" % (str(e)))
            sys.stderr.flush()
            raise Exception("Error occurred during database insertion.")

    def getStandardDeviation(self, ID):
        raise Exception("getStandardDeviation needs to be updated.")
        results = self.callproc("GetScores", ID)
        scores = []
        if len(results) == 1:
            return 0
        else:
            for result in results:
                if result["NumberOfMeasurements"] != 1:
                    raise Exception("Need to add logic for calculating standard deviation.")
            scores = [result["ddG"] for result in results]
            stddev, variance = computeStandardDeviation(scores)
            return stddev

class DatabasePrimer(object):
    '''This class fills in initial values for Tool, AminoAcid, UniProtKB, and UniProtKBMapping. The last will print errors if the corresponding PDB is not in the database.'''

    def __init__(self, ddGdb):
        self.ddGdb = ddGdb
        if False:
            self.insertTools()
            self.insertAminoAcids()
            self.insertUniProtKB()

    def computeBFactors(self):
        SQL = 'SELECT PDB_ID, Content FROM Structure'
        results = self.ddGdb.locked_execute(SQL)
        for result in results:
            pdbID = result["PDB_ID"]
            colortext.message(pdbID)
            pdb = PDB(result["Content"].split("\n"))
            BF = pickle.dumps(pdb.ComputeBFactors())
            SQL = ('UPDATE Structure SET %s=' % FieldNames_.BFactors) + '%s WHERE PDB_ID = %s'
            results = self.ddGdb.locked_execute(SQL, parameters = (BF, pdbID))

    def checkForCSEandMSE(self):
        SQL = 'SELECT PDB_ID, Content FROM Structure'
        results = self.ddGdb.locked_execute(SQL)
        for result in results:
            pdbID = result["PDB_ID"]
            pdb = PDB(result["Content"].split("\n"))
            foundRes = pdb.CheckForPresenceOf(["CSE", "MSE"])
            if foundRes:
                colortext.warning("The PDB %s contains residues which could affect computation (%s)." % (pdbID, join(foundRes, ", ")))
                if "CSE" in foundRes:
                    colortext.printf(pdbID, color = 'lightpurple')
                    colortext.warning("The PDB %s contains CSE. Check." % pdbID)
                if "MSE" in foundRes:
                    colortext.printf(pdbID, color = 'lightpurple')
                    colortext.warning("The PDB %s contains MSE. Check." % pdbID)

    def insertTools(self):
        emptydict = pickle.dumps({})

        # Tool name, Version, SVN information
        Tools = [
            ('CC/PBSA',		'Unknown', 	0, emptydict),
            ('EGAD', 		'Unknown', 	0, emptydict),
            ('FoldX',		'3.0', 		0, emptydict),
            ('Hunter',		'Unknown', 	0, emptydict),
            ('IMutant2',	'2.0', 		0, emptydict),
        ]

        Tools.append(('Rosetta', '2.1', 8075,
            pickle.dumps({
                "FirstBranchRevision"			:	9888,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1",
                "SourceSVNRevision"				:	8075,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/rosetta++",
                "DatabaseSVNRevisionInTrunk"	:	7966,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/rosetta_database",
            })))

        Tools.append(('Rosetta', '2.1.1', 13074,
            pickle.dumps({
                "FirstBranchRevision"			:	13894 ,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1.1",
                "SourceSVNRevision"				:	13074,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1/rosetta++",
                "DatabaseSVNRevisionInTrunk"	:	13074,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1/rosetta_database",
            })))

        Tools.append(('Rosetta', '2.1.2', 15393 ,
            pickle.dumps({
                "FirstBranchRevision"			:	15394,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1.2",
                "SourceSVNRevision"				:	15393,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1.1/rosetta++",
                "DatabaseSVNRevisionInTrunk"	:	15393,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.1.1/rosetta_database",
            })))

        Tools.append(('Rosetta', '2.2.0', 16310,
            pickle.dumps({
                "FirstBranchRevision"			:	16427,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.2.0",
                "SourceSVNRevision"				:	16310,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/rosetta++",
                "DatabaseSVNRevisionInTrunk"	:	15843,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/rosetta_database",
            })))

        Tools.append(('Rosetta', '2.3.0', 20729,
            pickle.dumps({
                "FirstBranchRevision"			:	20798,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.0",
                "SourceSVNRevision"				:	20729,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/rosetta++",
                "DatabaseSVNRevisionInTrunk"	:	20479,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/rosetta_database",
            })))

        Tools.append(('Rosetta', '2.3.1', 0,
            pickle.dumps({
                "FirstBranchRevision"			:	36012,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.1",
                "SourceSVNRevision"				:	None,
                "Source_SVN_URL"				:	None,
                "DatabaseSVNRevisionInTrunk"	:	None,
                "Database_SVN_URL"				:	None,
            })))

        Tools.append(('Rosetta', '3.0', 26316,
            pickle.dumps({
                "FirstBranchRevision"			:	26323 ,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.0",
                "SourceSVNRevision"				:	26316,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/mini",
                "DatabaseSVNRevisionInTrunk"	:	26298,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/minirosetta_database",
            })))

        Tools.append(('Rosetta', 'r32231', 32231,
            pickle.dumps({
                "FirstBranchRevision"			:	None,
                "Branch_SVN_URL"				:	None,
                "SourceSVNRevision"				:	32231,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/mini",
                "DatabaseSVNRevisionInTrunk"	:	32231,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/minirosetta_database",
            })))

        Tools.append(('Rosetta', 'r32257', 32257,
            pickle.dumps({
                "FirstBranchRevision"			:	None,
                "Branch_SVN_URL"				:	None,
                "SourceSVNRevision"				:	32257,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/mini",
                "DatabaseSVNRevisionInTrunk"	:	32257,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/minirosetta_database",
            })))

        Tools.append(('Rosetta', '3.1', 32528,
            pickle.dumps({
                "FirstBranchRevision"			:	30467,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.1",
                "SourceSVNRevision"				:	32528,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/mini",
                "DatabaseSVNRevisionInTrunk"	:	32509,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/minirosetta_database",
            })))


        Tools.append(('Rosetta', '3.2', 39284,
            pickle.dumps({
                "FirstBranchRevision"			:	39352,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.2",
                "SourceSVNRevision"				:	39284,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/mini",
                "DatabaseSVNRevisionInTrunk"	:	39117,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/minirosetta_database",
            })))

        Tools.append(('Rosetta', '3.2.1', 40878,
            pickle.dumps({
                "FirstBranchRevision"			:	40885,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.2.1",
                "SourceSVNRevision"				:	40878,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.2/rosetta_source",
                "DatabaseSVNRevisionInTrunk"	:	40878,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.2/rosetta_database",
            })))

        Tools.append(('Rosetta', '3.3', 42941,
            pickle.dumps({
                "FirstBranchRevision"			:	42943,
                "Branch_SVN_URL"				:	"https://svn.rosettacommons.org/source/branches/releases/rosetta-3.3",
                "SourceSVNRevision"				:	42941,
                "Source_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/mini",
                "DatabaseSVNRevisionInTrunk"	:	42940,
                "Database_SVN_URL"				:	"https://svn.rosettacommons.org/source/trunk/minirosetta_database",
            })))

        for t in Tools:
            SQL = 'SELECT * FROM Tool WHERE Name=%s AND Version=%s AND SVNRevision=%s'
            numresults = len(self.ddGdb.locked_execute(SQL, parameters = (t[0], t[1], t[2])))
            assert(numresults == 0 or numresults == 1)
            if numresults == 0:
                SQL = 'INSERT INTO Tool (Name, Version, SVNRevision, SVNRevisionInfo) VALUES (%s, %s, %s, %s)'
                self.ddGdb.locked_execute(SQL, parameters = t)

    def deleteAllExperimentalData(self):
        '''THIS WILL REMOVE *ALL* EXPERIMENTAL DATA FROM THE DATABASE. USE AT GREAT RISK!
           This function runs much quicker than the selective data removal function removeExperimentalData.
           It should fail when there are associated Predictions as this breaks a foreign key constraint.
           This is by design; Prediction data should not be deleted lightly.
           To avoid deleting the other records associated with the Experiment, we raise an exception first.
          '''

        predictions = self.ddGdb.locked_execute('SELECT * FROM Prediction')
        if not predictions:
            results = self.ddGdb.locked_execute('DELETE FROM ExperimentInterface')
            results = self.ddGdb.locked_execute('DELETE FROM ExperimentChain')
            results = self.ddGdb.locked_execute('DELETE FROM ExperimentMutation')
            results = self.ddGdb.locked_execute('DELETE FROM ExperimentMutant')
            results = self.ddGdb.locked_execute('DELETE FROM ExperimentScore')
            results = self.ddGdb.locked_execute('DELETE FROM Experiment')
        else:
            raise Exception("Database integrity failure: Cannot delete an Experiment (ID = %s) with an associated Prediction (ID = %s)." % (ID, predictions[0]['ID']))

    def deleteExperimentalDataByDataset(self):
        '''THIS WILL REMOVE ALL EXPERIMENTAL DATA FROM THE removethese ARRAY FROM THE DATABASE. USE AT GREAT RISK!
           It should fail when there are associated Predictions as this breaks a foreign key constraint.
           This is by design; Prediction data should not be deleted lightly.
           To avoid deleting the other records associated with the Experiment, we raise an exception first.
          '''

        removethese = ["Potapov-2009", "SenLiu-ComplexExperimentalDataset", "ProTherm-2008-09-08-23581"]
        experimentIDs = []
        for dataset in removethese:
            SQL = 'SELECT ID FROM Experiment WHERE Source=%s'
            results = self.ddGdb.locked_execute(SQL, parameters = (dataset,))
            for result in results:
                experimentIDs.append(result['ID'])

        for ID in experimentIDs:
            predictions = self.ddGdb.locked_execute('SELECT * FROM Prediction WHERE ExperimentID=%s', parameters = (ID,))
            if not predictions:
                results = self.ddGdb.locked_execute('DELETE FROM ExperimentInterface WHERE ExperimentID=%s', parameters = (ID,))
                results = self.ddGdb.locked_execute('DELETE FROM ExperimentChain WHERE ExperimentID=%s', parameters = (ID,))
                results = self.ddGdb.locked_execute('DELETE FROM ExperimentMutation WHERE ExperimentID=%s', parameters = (ID,))
                results = self.ddGdb.locked_execute('DELETE FROM ExperimentMutant WHERE ExperimentID=%s', parameters = (ID,))
                results = self.ddGdb.locked_execute('DELETE FROM ExperimentScore WHERE ExperimentID=%s', parameters = (ID,))
                results = self.ddGdb.locked_execute('DELETE FROM Experiment WHERE ID=%s', parameters = (ID,))
            else:
                raise Exception("Database integrity failure: Cannot delete an Experiment (ID = %s) with an associated Prediction (ID = %s)." % (ID, predictions[0]['ID']))

    def insertUniProtKB(self):
        uniprot = os.path.join("..", "rawdata", "uniprotmapping.csv")
        F = open(uniprot)
        lines = F.read().split("\n")[1:]
        F.close()
        UniProtKB = {}
        UniProtKBMapping = {}

        for line in lines:
            data = line.split("\t")
            if len(data) == 3:
                PDBID, AC, ID = data
                PDBID = PDBID.strip()
                AC = AC.strip()
                ID = ID.strip()
                UniProtKB[AC] = ID
                UniProtKBMapping[AC] = UniProtKBMapping.get(AC, []) or []
                UniProtKBMapping[AC].append(PDBID)

        for AC, ID in sorted(UniProtKB.iteritems()):
            if not self.ddGdb.locked_execute("SELECT * FROM UniProtKB WHERE UniProtKB_AC=%s", parameters = (AC,)):
                SQL = 'INSERT INTO UniProtKB (UniProtKB_AC, UniProtKB_ID) VALUES (%s, %s);'
                self.ddGdb.locked_execute(SQL, parameters = (AC, ID))

        for AC, pdbIDs in sorted(UniProtKBMapping.iteritems()):
            for pdbID in pdbIDs:

                if not self.ddGdb.locked_execute("SELECT * FROM UniProtKBMapping WHERE UniProtKB_AC=%s AND PDB_ID=%s", parameters = (AC, pdbID)):
                    SQL = 'INSERT INTO UniProtKBMapping (UniProtKB_AC, PDB_ID) VALUES (%s, %s);'
                    try:
                        self.ddGdb.locked_execute(SQL, parameters = (AC, pdbID), quiet = True)
                    except:
                        print("Error inserting UniProt record AC %s for PDB ID %s." % (AC, pdbID))


    def insertAminoAcids(self):
        global aas
        for aa in aas:
            SQL = 'INSERT INTO AminoAcids (Code, LongCode, Name, Polarity, Size) VALUES (%s, %s, %s, %s, %s);'
            self.ddGdb.locked_execute(SQL, parameters = tuple(aa))

    def insertKelloggLeaverFayBakerProtocols(self):

        protocols = [{} for i in range(0,21)]
                #
        commonstr = [
            '-in:file:s', '%(in:file:s)s',
            '-resfile', '%(resfile)s',
            '-database', '%(DATABASE_DIR)s',
            '-ignore_unrecognized_res',
            '-in:file:fullatom',
            '-constraints::cst_file', '%(constraints::cst_file)s'
        ]

        softrep = ['-score:weights', 'soft_rep_design']
        hardrep = ['-score:weights standard', '-score:patch score12']
        minnohardrep = ['-ddg::minimization_scorefunction', 'standard', '-ddg::minimization_patch', 'score12']

        protocols1617 = [
            '-ddg::weight_file', 'soft_rep_design',
            '-ddg::iterations', '50',
            '-ddg::local_opt_only', 'false',
            '-ddg::min_cst', 'true',
            '-ddg::mean', 'false',
            '-ddg::min', 'true',
            '-ddg::sc_min_only', 'false', # Backbone and sidechain minimization
            '-ddg::ramp_repulsive', 'true',
            '-ddg::minimization_scorefunction', 'standard',
            '-ddg::minimization_patch', 'score12'
        ]

        # Command for protocol 16 preminimization
        preminCmd = {
            FieldNames_.Type : "CommandLine",
            FieldNames_.Command : pickle.dumps([
                '%(BIN_DIR)s/minimize_with_cst.static.linuxgccrelease',
                '-in:file:l', '%(-in:file:l)s',
                '-in:file:fullatom',
                '-ignore_unrecognized_res',
                '-fa_max_dis', '9.0',
                '-database', '%(DATABASE_DIR)s',
                '-ddg::harmonic_ca_tether', '0.5',
                '-score:weights', 'standard',
                '-ddg::constraint_weight','1.0',
                '-ddg::out_pdb_prefix', 'min_cst_0.5',
                '-ddg::sc_min_only', 'false',
                '-score:patch', 'score12']),
            FieldNames_.Description : "Preminimization for Kellogg:10.1002/prot.22921:protocol16:32231",
        }
        alreadyExists = self.ddGdb.locked_execute("SELECT ID FROM Command WHERE Type=%s AND Command=%s", parameters = (preminCmd[FieldNames_.Type], preminCmd[FieldNames_.Command]))
        if not alreadyExists:
            self.ddGdb.insertDict('Command', preminCmd)
            preminCmdID = self.ddGdb.getLastRowID()
        else:
            preminCmdID = alreadyExists[0]["ID"]

        # Command for protocol 16 ddG
        ddGCmd = {
            FieldNames_.Type : "CommandLine",
            FieldNames_.Command : pickle.dumps(['%(BIN_DIR)s/fix_bb_monomer_ddg.linuxgccrelease'] + commonstr + softrep +  protocols1617),
            FieldNames_.Description : "ddG for Kellogg:10.1002/prot.22921:protocol16:32231",
        }
        alreadyExists = self.ddGdb.locked_execute("SELECT ID FROM Command WHERE Type=%s AND Command=%s", parameters = (ddGCmd[FieldNames_.Type], ddGCmd[FieldNames_.Command]))
        if not alreadyExists:
            self.ddGdb.insertDict('Command', ddGCmd)
            ddGCmdID = self.ddGdb.getLastRowID()
        else:
            ddGCmdID = alreadyExists[0]["ID"]

        # Protocol 16
        name = "Kellogg:10.1002/prot.22921:protocol16:32231"
        alreadyExists = self.ddGdb.locked_execute("SELECT ID FROM Protocol WHERE ID=%s", parameters = (name,))
        if not alreadyExists:
            PreMinTool = self.ddGdb.locked_execute("SELECT ID FROM Tool WHERE Name=%s and Version=%s", parameters = ("Rosetta", 3.3))
            ddGTool = self.ddGdb.locked_execute("SELECT ID FROM Tool WHERE Name=%s and SVNRevision=%s", parameters = ("Rosetta", 32231))
            ddGDatabaseToolID = self.ddGdb.locked_execute("SELECT ID FROM Tool WHERE Name=%s and SVNRevision=%s", parameters = ("Rosetta", 32257))
            if PreMinTool and ddGTool and ddGDatabaseToolID:
                PreMinTool = PreMinTool[0]["ID"]
                ddGTool = ddGTool[0]["ID"]
                ddGDatabaseToolID = ddGDatabaseToolID[0]["ID"]
            else:
                raise Exception("Cannot add protocol %s." % name)
            print("Inserting %s." % name)
            proto = {
                FieldNames_.ID : name,
                FieldNames_.Description : "Protocol 16 from Kellogg, Leaver-Fay, and Baker",
            }
            self.ddGdb.insertDict('Protocol', proto)
            pstep = {
                FieldNames_.ProtocolID : name,
                FieldNames_.StepID : "preminimization",
                FieldNames_.ToolID : PreMinTool,
                FieldNames_.CommandID : preminCmdID,
                FieldNames_.DatabaseToolID : PreMinTool,
                FieldNames_.DirectoryName : "",
                FieldNames_.ClassName : None,
                FieldNames_.Description : "Preminimization step",
            }
            self.ddGdb.insertDict('ProtocolStep', pstep)
            pstep = {
                FieldNames_.ProtocolID : name,
                FieldNames_.StepID : "ddG",
                FieldNames_.ToolID : ddGTool,
                FieldNames_.CommandID : ddGCmdID,
                FieldNames_.DatabaseToolID : ddGDatabaseToolID,
                FieldNames_.DirectoryName : "",
                FieldNames_.ClassName : None,
                FieldNames_.Description : "ddG step",
            }
            self.ddGdb.insertDict('ProtocolStep', pstep)
            pedge = {
                FieldNames_.ProtocolID : name,
                FieldNames_.FromStep : 1,
                FieldNames_.ToStep : 2,
            }
            self.ddGdb.insertDict('ProtocolGraphEdge', pedge)

    def updateCommand(self):
        # Command for protocol 16 ddG
        commonstr = [
            '-in:file:s', '%(in:file:s)s',
            '-resfile', '%(resfile)s',
            '-database', '%(DATABASE_DIR)s',
            '-ignore_unrecognized_res',
            '-in:file:fullatom',
            '-constraints::cst_file', '%(constraints::cst_file)s'
        ]

        softrep = ['-score:weights', 'soft_rep_design']
        hardrep = ['-score:weights standard', '-score:patch score12']
        minnohardrep = ['-ddg::minimization_scorefunction', 'standard', '-ddg::minimization_patch', 'score12']

        protocols1617 = [
            '-ddg::weight_file', 'soft_rep_design',
            '-ddg::iterations', '50',
            '-ddg::local_opt_only', 'false',
            '-ddg::min_cst', 'true',
            '-ddg::mean', 'false',
            '-ddg::min', 'true',
            '-ddg::sc_min_only', 'false', # Backbone and sidechain minimization
            '-ddg::ramp_repulsive', 'true',
            '-ddg::minimization_scorefunction', 'standard',
            '-ddg::minimization_patch', 'score12'
        ]

        newcmd = pickle.dumps(['%(BIN_DIR)s/fix_bb_monomer_ddg.linuxgccrelease'] + commonstr + softrep +  protocols1617)
        ddGdb.locked_execute("UPDATE Command SET Command=%s WHERE ID=5;", parameters = (newcmd,))
        newcmd = pickle.dumps([
                '%(BIN_DIR)s/minimize_with_cst.static.linuxgccrelease',
                '-in:file:l', '%(in:file:l)s',
                '-in:file:fullatom',
                '-ignore_unrecognized_res',
                '-fa_max_dis', '9.0',
                '-database', '%(DATABASE_DIR)s',
                '-ddg::harmonic_ca_tether', '0.5',
                '-score:weights', 'standard',
                '-ddg::constraint_weight','1.0',
                '-ddg::out_pdb_prefix', 'min_cst_0.5',
                '-ddg::sc_min_only', 'false',
                '-score:patch', 'score12'])
        ddGdb.locked_execute("UPDATE Command SET Command=%s WHERE ID=4;", parameters = (newcmd,))

        newcmd = "%(BIN_DIR)s/minimize_with_cst.static.linuxgccrelease -in:file:l %(in:file:l)s -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -database %(DATABASE_DIR)s -ddg::harmonic_ca_tether 0.5 -score:weights standard -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false -score:patch score12"
        ddGdb.locked_execute("UPDATE Command SET Command=%s WHERE ID=4;", parameters = (newcmd,))

        newcmd = "%(BIN_DIR)s/fix_bb_monomer_ddg.linuxgccrelease -in:file:s %(in:file:s)s -resfile %(resfile)s -database %(DATABASE_DIR)s -ignore_unrecognized_res -in:file:fullatom -constraints::cst_file %(constraints::cst_file)s -score:weights soft_rep_design -ddg::weight_file soft_rep_design -ddg::iterations 50 -ddg::local_opt_only false -ddg::min_cst true -ddg::mean false -ddg::min true -ddg::sc_min_only false -ddg::ramp_repulsive true -ddg::minimization_scorefunction standard -ddg::minimization_patch score12"
        ddGdb.locked_execute("UPDATE Command SET Command=%s WHERE ID=5;", parameters = (newcmd,))


if __name__ == "__main__":
    ddGdb = ddGDatabase()

    #primer = DatabasePrimer(ddGdb)
    #primer.insertUniProtKB()
    #primer.checkForCSEandMSE()
    #primer.computeBFactors()
    #print("Removing all data")
    #primer.deleteAllExperimentalData()
    #primer.insertKelloggLeaverFayBakerProtocols()
    #primer.insertTools()
    #primer.addPDBSources()
    #primer.updateCommand()
