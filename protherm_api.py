#!/usr/bin/python2.4
# -*- coding: utf-8 -*-

import sys
import os
import re
import string
sys.path.insert(0, "..")
from string import join
import common.colortext as colortext
import ddgdbapi
from common.rosettahelper import kJtokcal
from common.rosettahelper import NUMBER_KELVIN_AT_ZERO_CELSIUS

sometimesFields = ["ION_NAME_2", "ION_CONC_2", "ION_NAME_3", "ION_CONC_3", "ION_CON_1C"] # ION_CON_1C seems to be a weird typo in version 23581

field_order = [
	["Protein information",
		"PROTEIN"				, 
		"SOURCE"				, 
		"LENGTH"				,
		"MOL-WEIGHT"			, 
		"PIR_ID"				, 
		"SWISSPROT_ID"			, 
		"E.C.NUMBER"			, 
		"PMD.NO"				, 
		"PDB_wild"				, 
		"PDB_mutant"			, 
		"MUTATION"				, 
		"MUTATED_CHAIN"			, 
		"NO_MOLECULE"			, 
		"SEC.STR."				, 
		"ASA"					,
	], 
	["Experimental condition",
		"T"						, 
		"pH"					, 
		"BUFFER_NAME"			, 
		"BUFFER_CONC"			, 
		"ION_NAME_1"			, 
		"ION_NAME_1C"			, 
		"ION_CONC_1"			, 
		"ION_NAME_2"			, 
		"ION_CONC_2"			, 
		"ION_NAME_3"			, 
		"ION_CONC_3"			, 
		"PROTEIN_CONC"			, 
		"MEASURE"				, 
		"METHOD"				,
	], 
	["Thermodynamic data", 
		["Denaturant denaturation",
			"dG_H2O"			, 
			"ddG_H2O"			, 
			"Cm"				, 
			"m"					,
		],
		["Thermal denaturation",
			"dG"				, 
			"ddG"				, 
			"Tm"				, 
			"dTm"				, 
			"dHvH"				, 
			"dHcal"				, 
			"dCp"				,
		],
		["Other",
			"STATE"				, 
			"REVERSIBILITY"		, 
			"ACTIVITY"			, 
			"ACTIVITY_Km"		, 
			"ACTIVITY_Kcat"		, 
			"ACTIVITY_Kd"		,
		],
	], 
	["Literature",
		"REFERENCE"				, 
		"AUTHOR"				, 
		"KEY_WORDS"				, 
		"REMARKS"				, 
		"RELATED_ENTRIES"		, 
	],
]

field_descriptions = {
	# Protein information
	"PROTEIN"				: "Protein name",
	"SOURCE"				: "Protein source",
	"LENGTH"				: "Protein length",
	"MOL-WEIGHT"			: "",
	"PIR_ID"				: "PIR ID",
	"SWISSPROT_ID"			: "SWISSPROT ID",
	"E.C.NUMBER"			: "E.C. number",
	"PMD.NO"				: "PMD number",
	"PDB_wild"				: "PDB code of wildtype",
	"PDB_mutant"			: "PDB code of mutant",
	"MUTATION"				: "Mutation detail (wildtype, position, mutation)",
	"MUTATED_CHAIN"			: "Mutation detail (chain)",
	"NO_MOLECULE"			: "Number of molecules?",
	"SEC.STR."				: "Mutation detail? (secondary structure)",
	"ASA"					: "Accessible surface area (ASA)",
	# Experimental condition:
	"T"						: "Temperature",
	"pH"					: "pH",
	"BUFFER_NAME"			: "Buffer",
	"BUFFER_CONC"			: "Buffer concentration",
	"ION_NAME_1"			: "Ion?",
	"ION_CONC_1"			: "Ion concentration?",
	"ION_NAME_2"			: "Ion?",
	"ION_CONC_2"			: "Ion concentration?",
	"ION_NAME_3"			: "Ion?",
	"ION_CONC_3"			: "Ion concentration?",
	"ION_NAME_1C"			: "Ion? (typo)",
	"PROTEIN_CONC"			: "Protein concentration",
	"MEASURE"				: "Measure (DSC, CD and so on)",
	"METHOD"				: "Method of denaturation",
	# Thermodynamic data 
		# Denaturant denaturation",
		"dG_H2O"			: "Free energy of unfolding: DGH2O",
		"ddG_H2O"			: "Difference in DGH2O : DDGH2O",
		"Cm"				: "Denaturation concentration: Cm",
		"m"					: "Slope of denaturation curve: m",
		# Thermal denaturation:",
		"dG"				: "Free energy of unfolding: DG",
		"ddG"				: "Difference in DG: DDG",
		"Tm"				: "Transition temperature: Tm",
		"dTm"				: "Change in Tm: DTm",
		"dHvH"				: "Enthalpy change: DHvH",
		"dHcal"				: "Enthalpy change: DHcal",
		"dCp"				: "Heat capacity change: DCp",
		# Unknown
		"STATE"				: "?",
		"REVERSIBILITY"		: "?",
		"ACTIVITY"			: "?",
		"ACTIVITY_Km"		: "?",
		"ACTIVITY_Kcat"		: "?",
		"ACTIVITY_Kd"		: "?",
	# Literature:
	"REFERENCE"				: "Reference",
	"AUTHOR"				: "Author",
	"KEY_WORDS"				: "Keywords",
	"REMARKS"				: "Remarks",
	"RELATED_ENTRIES"		: "Related entries",
}

fields_of_interest = {
	"PROTEIN" 			: True,
	"SOURCE" 			: True,
	"LENGTH" 			: True,
	"PDB_wild" 			: True,
	"PDB_mutant" 		: True,
	"MUTATION" 			: True,
	"MUTATED_CHAIN" 	: True,
	"T" 				: True,
	"pH" 				: True,
	"BUFFER_NAME" 		: True,
	"BUFFER_CONC" 		: True,
	"ION_NAME_1" 		: True,
	"ION_CONC_1" 		: True,
	"ION_NAME_2" 		: True,
	"ION_CONC_2" 		: True,
	"ION_NAME_3" 		: True,
	"ION_CONC_3" 		: True,
	"ION_NAME_1C" 		: True,
	"PROTEIN_CONC" 		: True,
	"MEASURE" 			: True,
	"METHOD" 			: True,
	"ddG" 				: True,
	"REFERENCE" 		: True,
}

fields_requiring_qualification = ["ddG"]  
		
secondary_structure_values = {
	"c" 		: "Coil",
	"Coil"		: "Coil",
	"s"			: "Sheet",
	"Sheet"		: "Sheet",
	"h"			: "Helix",
	"Helix"		: "Helix",
	"Turn"		: "Turn",
}

MeasureMapping = {
	'Abs' : 'Absorbance',
	'Absorbance' : 'Absorbance',
	'Absorption' : 'Absorbance',
	'Activity' : 'Activity',
	'activity' : 'Activity',
	'Activity assay' : 'Activity',
	'Anisotropy' : 'Anisotropy',
	'ANS binding' : 'ANS binding',
	'Capillary electrophoresis' : 'Capillary electrophoresis',
	'CD' : 'CD',
	'CD (far-UV)' : 'CD (far-UV)',
	'CD(far-UV)' : 'CD (far-UV)',
	'CD (near-UV)' : 'CD (near-UV)',
	'CD(near-UV)' : 'CD (near-UV)',
	'Chromatography' : 'Chromatography',
	'DSC' : 'DSC',
	'DSMC' : 'DSMC',
	'Emission' : 'Emission',
	'Enzyme activity' : 'Enzyme activity',
	'enzyme activity' : 'Enzyme activity',
	'Enzyme assay' : 'Enzyme assay',
	'EPR' : 'EPR',
	'ESR' : 'ESR',
	'far-UV CD' : 'CD (far-UV)',
	'Fluorescence' : 'Fluorescence',
	'Fluorescence (ANS)' : 'Fluorescence (ANS)',
	'Fluorescence (Trp)' : 'Fluorescence (Trp)',
	'FTIR' : 'FTIR',
	'Gel electrophoresis' : 'Gel electrophoresis',
	'HPLC' : 'HPLC',
	'Hydrogen exchange' : 'Hydrogen exchange',
	'IATC' : 'IATC',
	'IR spectroscopy' : 'IR spectroscopy',
	'Isothermal denaturation' : 'Isothermal denaturation',
	'Light scattering' : 'Light-scattering',
	'Light-scattering' : 'Light-scattering',
	'Magnetic Relaxation Dispersion' : 'Magnetic relaxation dispersion',
	'near-UV CD' : 'CD (near-UV)',
	'NMR' : 'NMR',
	'NMR amide hydrogen exchange' : 'NMR amide hydrogen exchange',
	'NMR hydrogen exchange' : 'NMR hydrogen exchange',
	'NMR Hydrogen exchange' : 'NMR hydrogen exchange',
	'optical' : 'Optical',
	'Optical Density' : 'Optical density',
	'Pulse Protolysis' : 'Pulse protolysis',
	'Quantitative cysteine reactivity' : 'Quantitative cysteine reactivity',
	'Refraction' : 'Refraction',
	'SAXS' : 'SAXS',
	'SEC' : 'SEC',
	'SUPREX' : 'SUPREX',
	'Thiol reactivity' : 'Thiol reactivity',
	'UV spectroscopy' : 'UV spectroscopy',
}

MethodMapping = {
	'2-Propanol' : '2-Propanol',
	'Acid' : 'Acid',
	'Activity' : 'Activity',
	'DimethylUrea' : 'Dimethylurea',
	'Dynamic fluctuation' : 'Dynamic fluctuation',
	'GSH' : 'GSH',
	'GSSG' : 'GSSG',
	'GdnHCl' : 'GdnHCl',
	'GdnHSCN' : 'GdnHSCN',
	'GdnSCN' : 'GdnSCN',
	'HClO4' : 'HClO4',
	'Heat treatment' : 'Heat treatment',
	'heat treatment' : 'Heat treatment',
	'KSCN' : 'KSCN',
	'LiCl' : 'LiCl',
	'NaCl' : 'NaCl',
	'NaClO4 titration' : 'NaClO4 titration',
	'Pressure' : 'Pressure',
	'pressure' : 'Pressure',
	'Pressure denaturation' : 'Pressure denaturation',
	'SDS' : 'SDS',
	'TFE' : 'TFE',
	'Thermal' : 'Thermal',
	'Urea' : 'Urea',
	'Urea/GdnHCl(0.5 M)' : 'Urea/GdnHCl(0.5 M)',
	'Urea/GdnHCl(0.9 M)' : 'Urea/GdnHCl(0.9 M)',
	'Urea/GdnHCl(1.35 M)' : 'Urea/GdnHCl(1.35 M)',
}

# These are the records where the mutations include insertion codes
# It turns out that none of these are eligible for inclusion in the database
# as they are all missing ddG values.
# As a precaution, we cast all residue IDs to int which will raise an exception
# if an insertion code is parsed.
iCodeRecords = [5438, 5439, 5440, 5441, 8060, 13083, 13084]

# NOTE: THESE ARE PATCHES FOR MISSING DATA IN ProTherm
# This is used to print out information to the admin on what needs patching
patchfields = {
	"MUTATED_CHAIN" : ("PDB_wild", "MUTATION"), #"%(PDB_wild)s-%(MUTATION)s",
	"LENGTH" 		: ("PDB_wild"),
	"PDB_wild" 		: ("SWISSPROT_ID"),
}
	
patch = {
	2396 : {'PDB_wild' : None}, # -> 2405. P08505 No related PDB entry. 
	2397 : {'PDB_wild' : None}, #P08505
	2398 : {'PDB_wild' : None}, #P08505
	2400 : {'PDB_wild' : None}, #P08505
	2401 : {'PDB_wild' : None}, #P08505
	2403 : {'PDB_wild' : None}, #P08505
	2404 : {'PDB_wild' : None}, #P08505
	2405 : {'PDB_wild' : None}, #P08505
	4216 : {'PDB_wild' : None}, #P00912 No related PDB entry.
	8588 : {'LENGTH' : 104}, #1ONC
	8589 : {'LENGTH' : 104}, #1ONC
	8590 : {'LENGTH' : 104}, #1ONC
	8591 : {'LENGTH' : 104}, #1ONC
	8592 : {'LENGTH' : 104}, #1ONC
	8593 : {'LENGTH' : 104}, #1ONC
	8594 : {'LENGTH' : 104}, #1ONC
	8595 : {'LENGTH' : 104}, #1ONC
	8596 : {'LENGTH' : 104}, #1ONC
	14229 : {'PDB_wild' : None}, # -> 14233. P08821 No related PDB entry.
	14230 : {'PDB_wild' : None}, #P08821
	14231 : {'PDB_wild' : None}, #P08821
	14232 : {'PDB_wild' : None}, #P08821
	14233 : {'PDB_wild' : None}, #P08821
	14978 : {'LENGTH' : 238}, #1CHK
	14979 : {'LENGTH' : 238}, #1CHK
	14980 : {'LENGTH' : 238}, #1CHK
	14981 : {'LENGTH' : 238}, #1CHK
	14987 : {'LENGTH' : 238}, #1CHK
	14988 : {'LENGTH' : 238}, #1CHK
	14989 : {'LENGTH' : 238}, #1CHK
	14990 : {'LENGTH' : 238}, #1CHK
	14996 : {'LENGTH' : 238}, #1CHK
	14997 : {'LENGTH' : 238}, #1CHK
	14998 : {'LENGTH' : 238}, #1CHK
	14999 : {'LENGTH' : 238}, #1CHK
	16597 : {'LENGTH' : 435}, #1KFW
	16598 : {'LENGTH' : 435}, #1KFW
	16599 : {'LENGTH' : 435}, #1KFW
	16600 : {'LENGTH' : 435}, #1KFW
	19423 : {'MUTATED_CHAIN' : None},# -> 19538. 1OTR A33 - I cannot determine what chain this is
	19424 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19425 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19426 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19427 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19428 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19429 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19430 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19431 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19432 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19433 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19434 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19449 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19450 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19451 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19452 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19453 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19454 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19455 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19456 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19457 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19458 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19459 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19460 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19475 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19476 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19477 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19478 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19479 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19480 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19481 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19482 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19483 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19484 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19485 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19486 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19501 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19502 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19503 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19504 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19505 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19506 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19507 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19508 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19509 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19510 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19511 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19512 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19527 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19528 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19529 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19530 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19531 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19532 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19533 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19534 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19535 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19536 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19537 : {'MUTATED_CHAIN' : None},# 1OTR A33
	19538 : {'MUTATED_CHAIN' : None},# 1OTR A33
	21040 : {'MUTATED_CHAIN' : None},# -> 21332. 1CSP Cannot determine what the mutation is
	21041 : {'MUTATED_CHAIN' : None},# 1CSP
	21097 : {'MUTATED_CHAIN' : None},# 1CSP
	21098 : {'MUTATED_CHAIN' : None},# 1CSP
	21157 : {'MUTATED_CHAIN' : None},# 1CSP
	21158 : {'MUTATED_CHAIN' : None},# 1CSP
	21215 : {'MUTATED_CHAIN' : None},# 1CSP
	21216 : {'MUTATED_CHAIN' : None},# 1CSP
	21273 : {'MUTATED_CHAIN' : None},# 1CSP
	21274 : {'MUTATED_CHAIN' : None},# 1CSP
	21331 : {'MUTATED_CHAIN' : None},# 1CSP
	21332 : {'MUTATED_CHAIN' : None},# 1CSP
}

# These PDB files have exactly one chain but ProTherm lists the wrong chain e.g. '-' rather than 'A'
singleChainPDBs = {
	'1A23' : {'MUTATED_CHAIN' : 'A'},
	'1A2I' : {'MUTATED_CHAIN' : 'A'},
	'1A43' : {'MUTATED_CHAIN' : 'A'},
	'1A53' : {'MUTATED_CHAIN' : 'A'},
	'1A5E' : {'MUTATED_CHAIN' : 'A'},
	'1A70' : {'MUTATED_CHAIN' : 'A'},
	'1ABE' : {'MUTATED_CHAIN' : 'A'},
	'1AG2' : {'MUTATED_CHAIN' : 'A'},
	'1AG4' : {'MUTATED_CHAIN' : 'A'},
	'1AG6' : {'MUTATED_CHAIN' : 'A'},
	'1AIN' : {'MUTATED_CHAIN' : 'A'},
	'1AJ3' : {'MUTATED_CHAIN' : 'A'},
	'1AKK' : {'MUTATED_CHAIN' : 'A'},
	'1AMQ' : {'MUTATED_CHAIN' : 'A'},
	'1APC' : {'MUTATED_CHAIN' : 'A'},
	'1APS' : {'MUTATED_CHAIN' : 'A'},
	'1AQH' : {'MUTATED_CHAIN' : 'A'},
	'1AVR' : {'MUTATED_CHAIN' : 'A'},
	'1AX1' : {'MUTATED_CHAIN' : 'A'},
	'1AXB' : {'MUTATED_CHAIN' : 'A'},
	'1AYE' : {'MUTATED_CHAIN' : 'A'},
	'1B0O' : {'MUTATED_CHAIN' : 'A'},
	'1B5M' : {'MUTATED_CHAIN' : 'A'},
	'1BAH' : {'MUTATED_CHAIN' : 'A'},
	'1BCX' : {'MUTATED_CHAIN' : 'A'},
	'1BD8' : {'MUTATED_CHAIN' : 'A'},
	'1BET' : {'MUTATED_CHAIN' : 'A'},
	'1BKE' : {'MUTATED_CHAIN' : 'A'},
	'1BLC' : {'MUTATED_CHAIN' : 'A'},
	'1BMC' : {'MUTATED_CHAIN' : 'A'},
	'1BOY' : {'MUTATED_CHAIN' : 'A'},
	'1BP2' : {'MUTATED_CHAIN' : 'A'},
	'1BPI' : {'MUTATED_CHAIN' : 'A'},
	'1BPR' : {'MUTATED_CHAIN' : 'A'},
	'1BTA' : {'MUTATED_CHAIN' : 'A'},
	'1BVC' : {'MUTATED_CHAIN' : 'A'},
	'1BZO' : {'MUTATED_CHAIN' : 'A'},	
	'1C52' : {'MUTATED_CHAIN' : 'A'},
	'1C53' : {'MUTATED_CHAIN' : 'A'},
	'1C5G' : {'MUTATED_CHAIN' : 'A'},
	'1CAH' : {'MUTATED_CHAIN' : 'A'},
	'1CEY' : {'MUTATED_CHAIN' : 'A'},
	'1CLW' : {'MUTATED_CHAIN' : 'A'},
	'1CMS' : {'MUTATED_CHAIN' : 'A'},
	'1COK' : {'MUTATED_CHAIN' : 'A'},
	'1CPM' : {'MUTATED_CHAIN' : 'A'},
	'1CSP' : {'MUTATED_CHAIN' : 'A'},
	'1CTS' : {'MUTATED_CHAIN' : 'A'},
	'1CUS' : {'MUTATED_CHAIN' : 'A'},
	'1CYO' : {'MUTATED_CHAIN' : 'A'},
	'1D0X' : {'MUTATED_CHAIN' : 'A'},
	'1DE3' : {'MUTATED_CHAIN' : 'A'},
	'1DEC' : {'MUTATED_CHAIN' : 'A'},
	'1DFX' : {'MUTATED_CHAIN' : 'A'},
	'1DHN' : {'MUTATED_CHAIN' : 'A'},
	'1DIL' : {'MUTATED_CHAIN' : 'A'},
	'1DIV' : {'MUTATED_CHAIN' : 'A'},
	'1DLC' : {'MUTATED_CHAIN' : 'A'},
	'1DO9' : {'MUTATED_CHAIN' : 'A'},
	'1DTO' : {'MUTATED_CHAIN' : 'A'},
	'1DVC' : {'MUTATED_CHAIN' : 'A'},
	'1EKG' : {'MUTATED_CHAIN' : 'A'},
	'1ELV' : {'MUTATED_CHAIN' : 'A'},
	'1EQ1' : {'MUTATED_CHAIN' : 'A'},
	'1ERU' : {'MUTATED_CHAIN' : 'A'},
	'1EVQ' : {'MUTATED_CHAIN' : 'A'},
	'1EW4' : {'MUTATED_CHAIN' : 'A'},
	'1EXG' : {'MUTATED_CHAIN' : 'A'},
	'1EZA' : {'MUTATED_CHAIN' : 'A'},
	'1FAJ' : {'MUTATED_CHAIN' : 'A'},
	'1FEP' : {'MUTATED_CHAIN' : 'A'},
	'1FGA' : {'MUTATED_CHAIN' : 'A'},
	'1FKJ' : {'MUTATED_CHAIN' : 'A'},
	'1FLV' : {'MUTATED_CHAIN' : 'A'},
	'1FMM' : {'MUTATED_CHAIN' : 'S'},
	'1FNF' : {'MUTATED_CHAIN' : 'A'},
	'1FRD' : {'MUTATED_CHAIN' : 'A'},
	'1FTG' : {'MUTATED_CHAIN' : 'A'},
	'1FTT' : {'MUTATED_CHAIN' : 'A'},
	'1GAL' : {'MUTATED_CHAIN' : 'A'},
	'1GKG' : {'MUTATED_CHAIN' : 'A'},
	'1GLH' : {'MUTATED_CHAIN' : 'A'},
	'1GLM' : {'MUTATED_CHAIN' : 'A'},
	'1GPC' : {'MUTATED_CHAIN' : 'A'},
	'1GRX' : {'MUTATED_CHAIN' : 'A'},
	'1H09' : {'MUTATED_CHAIN' : 'A'},
	'1H7M' : {'MUTATED_CHAIN' : 'A'},
	'1HCD' : {'MUTATED_CHAIN' : 'A'},
	'1HEV' : {'MUTATED_CHAIN' : 'A'},
	'1HIC' : {'MUTATED_CHAIN' : 'A'},
	'1HK0' : {'MUTATED_CHAIN' : 'X'},
	'1HME' : {'MUTATED_CHAIN' : 'A'},
	'1HML' : {'MUTATED_CHAIN' : 'A'},
	'1HXN' : {'MUTATED_CHAIN' : 'A'},
	'1HYW' : {'MUTATED_CHAIN' : 'A'},
	'1IFB' : {'MUTATED_CHAIN' : 'A'},
	'1IGS' : {'MUTATED_CHAIN' : 'A'},
	'1IMQ' : {'MUTATED_CHAIN' : 'A'},
	'1IOB' : {'MUTATED_CHAIN' : 'A'},
	'1IOJ' : {'MUTATED_CHAIN' : 'A'},
	'1IRL' : {'MUTATED_CHAIN' : 'A'},
	'1IRO' : {'MUTATED_CHAIN' : 'A'},
	'1JAE' : {'MUTATED_CHAIN' : 'A'},
	'1JBK' : {'MUTATED_CHAIN' : 'A'},
	'1JHN' : {'MUTATED_CHAIN' : 'A'},
	'1JNK' : {'MUTATED_CHAIN' : 'A'},
	'1KDU' : {'MUTATED_CHAIN' : 'A'},
	'1KFD' : {'MUTATED_CHAIN' : 'A'},
	'1KKJ' : {'MUTATED_CHAIN' : 'A'},
	'1KTQ' : {'MUTATED_CHAIN' : 'A'},
	'1KUM' : {'MUTATED_CHAIN' : 'A'},
	'1L63' : {'MUTATED_CHAIN' : 'A'},
	'1LFO' : {'MUTATED_CHAIN' : 'A'},
	'1LPS' : {'MUTATED_CHAIN' : 'A'},
	'1LRE' : {'MUTATED_CHAIN' : 'A'},
	'1LS4' : {'MUTATED_CHAIN' : 'A'},
	'1LVE' : {'MUTATED_CHAIN' : 'A'},
	'1LZ1' : {'MUTATED_CHAIN' : 'A'},
	'1MAX' : {'MUTATED_CHAIN' : 'A'},
	'1MBD' : {'MUTATED_CHAIN' : 'A'},
	'1MBG' : {'MUTATED_CHAIN' : 'A'},
	'1MGR' : {'MUTATED_CHAIN' : 'A'},
	'1MJC' : {'MUTATED_CHAIN' : 'A'},
	'1MSI' : {'MUTATED_CHAIN' : 'A'},
	'1MUL' : {'MUTATED_CHAIN' : 'A'},
	'1N02' : {'MUTATED_CHAIN' : 'A'},
	'1OLR' : {'MUTATED_CHAIN' : 'A'},
	'1OMU' : {'MUTATED_CHAIN' : 'A'},
	'1ONC' : {'MUTATED_CHAIN' : 'A'},
	'1ORC' : {'MUTATED_CHAIN' : 'A'},
	'1OSA' : {'MUTATED_CHAIN' : 'A'},
	'1P2P' : {'MUTATED_CHAIN' : 'A'},
	'1PAH' : {'MUTATED_CHAIN' : 'A'},
	'1PBA' : {'MUTATED_CHAIN' : 'A'},
	'1PCA' : {'MUTATED_CHAIN' : 'A'},
	'1PDO' : {'MUTATED_CHAIN' : 'A'},
	'1PGA' : {'MUTATED_CHAIN' : 'A'},
	'1PHP' : {'MUTATED_CHAIN' : 'A'},
	'1PII' : {'MUTATED_CHAIN' : 'A'},
	'1PK2' : {'MUTATED_CHAIN' : 'A'},
	'1PMC' : {'MUTATED_CHAIN' : 'A'},
	'1POH' : {'MUTATED_CHAIN' : 'A'},
	'1PPI' : {'MUTATED_CHAIN' : 'A'},
	'1PPN' : {'MUTATED_CHAIN' : 'A'},
	'1PPP' : {'MUTATED_CHAIN' : 'A'},
	'1PRR' : {'MUTATED_CHAIN' : 'A'},
	'1QLX' : {'MUTATED_CHAIN' : 'A'},
	'1QQV' : {'MUTATED_CHAIN' : 'A'},
	'1RBP' : {'MUTATED_CHAIN' : 'A'},
	'1RCB' : {'MUTATED_CHAIN' : 'A'},
	'1RH1' : {'MUTATED_CHAIN' : 'A'},
	'1RHD' : {'MUTATED_CHAIN' : 'A'},
	'1RIL' : {'MUTATED_CHAIN' : 'A'},
	'1RIS' : {'MUTATED_CHAIN' : 'A'},
	'1RRO' : {'MUTATED_CHAIN' : 'A'},
	'1RTB' : {'MUTATED_CHAIN' : 'A'},
	'1RX4' : {'MUTATED_CHAIN' : 'A'},
	'1SAP' : {'MUTATED_CHAIN' : 'A'},
	'1SEE' : {'MUTATED_CHAIN' : 'A'},
	'1SFP' : {'MUTATED_CHAIN' : 'A'},
	'1SHG' : {'MUTATED_CHAIN' : 'A'},
	'1SMD' : {'MUTATED_CHAIN' : 'A'},
	'1SSO' : {'MUTATED_CHAIN' : 'A'},
	'1STN' : {'MUTATED_CHAIN' : 'A'},
	'1SUP' : {'MUTATED_CHAIN' : 'A'},
	'1TCA' : {'MUTATED_CHAIN' : 'A'},
	'1TEN' : {'MUTATED_CHAIN' : 'A'},
	'1TFE' : {'MUTATED_CHAIN' : 'A'},
	'1TGN' : {'MUTATED_CHAIN' : 'A'},
	'1TIN' : {'MUTATED_CHAIN' : 'A'},
	'1TIT' : {'MUTATED_CHAIN' : 'A'},
	'1TML' : {'MUTATED_CHAIN' : 'A'},
	'1TMY' : {'MUTATED_CHAIN' : 'A'},
	'1TOF' : {'MUTATED_CHAIN' : 'A'},
	'1TPE' : {'MUTATED_CHAIN' : 'A'},
	'1TTG' : {'MUTATED_CHAIN' : 'A'},
	'1TUR' : {'MUTATED_CHAIN' : 'A'},
	'1UBQ' : {'MUTATED_CHAIN' : 'A'},
	'1UCU' : {'MUTATED_CHAIN' : 'A'},
	'1UOX' : {'MUTATED_CHAIN' : 'A'},
	'1URK' : {'MUTATED_CHAIN' : 'A'},
	'1UW3' : {'MUTATED_CHAIN' : 'A'},
	'1VIE' : {'MUTATED_CHAIN' : 'A'},
	'1VQB' : {'MUTATED_CHAIN' : 'A'},
	'1W3D' : {'MUTATED_CHAIN' : 'A'},
	'1W4H' : {'MUTATED_CHAIN' : 'A'},
	'1WIT' : {'MUTATED_CHAIN' : 'A'},
	'1WRP' : {'MUTATED_CHAIN' : 'R'},
	'1XAS' : {'MUTATED_CHAIN' : 'A'},
	'1YAL' : {'MUTATED_CHAIN' : 'A'},
	'1YCC' : {'MUTATED_CHAIN' : 'A'},
	'1YEA' : {'MUTATED_CHAIN' : 'A'},
	'1YMB' : {'MUTATED_CHAIN' : 'A'},
	'2A36' : {'MUTATED_CHAIN' : 'A'},
	'2ABD' : {'MUTATED_CHAIN' : 'A'},
	'2ACE' : {'MUTATED_CHAIN' : 'A'},
	'2ACY' : {'MUTATED_CHAIN' : 'A'},
	'2ADA' : {'MUTATED_CHAIN' : 'A'},
	'2AIT' : {'MUTATED_CHAIN' : 'A'},
	'2AKY' : {'MUTATED_CHAIN' : 'A'},
	'2ASI' : {'MUTATED_CHAIN' : 'A'},
	'2B4Z' : {'MUTATED_CHAIN' : 'A'},
	'2BRD' : {'MUTATED_CHAIN' : 'A'},
	'2CBR' : {'MUTATED_CHAIN' : 'A'},
	'2CPP' : {'MUTATED_CHAIN' : 'A'},
	'2CRK' : {'MUTATED_CHAIN' : 'A'},
	'2CRO' : {'MUTATED_CHAIN' : 'A'},
	'2DRI' : {'MUTATED_CHAIN' : 'A'},
	'2EQL' : {'MUTATED_CHAIN' : 'A'},
	'2FAL' : {'MUTATED_CHAIN' : 'A'},
	'2FHA' : {'MUTATED_CHAIN' : 'A'},
	'2GA5' : {'MUTATED_CHAIN' : 'A'},
	'2HMB' : {'MUTATED_CHAIN' : 'A'},
	'2HPR' : {'MUTATED_CHAIN' : 'A'},
	'2IFB' : {'MUTATED_CHAIN' : 'A'},
	'2IMM' : {'MUTATED_CHAIN' : 'A'},
	'2LZM' : {'MUTATED_CHAIN' : 'A'},
	'2NUL' : {'MUTATED_CHAIN' : 'A'},
	'2PDD' : {'MUTATED_CHAIN' : 'A'},
	'2PEC' : {'MUTATED_CHAIN' : 'A'},
	'2PGK' : {'MUTATED_CHAIN' : 'A'},
	'2PRD' : {'MUTATED_CHAIN' : 'A'},
	'2RN2' : {'MUTATED_CHAIN' : 'A'},
	'2TRT' : {'MUTATED_CHAIN' : 'A'},
	'2TS1' : {'MUTATED_CHAIN' : 'A'},
	'3MBP' : {'MUTATED_CHAIN' : 'A'},
	'3PGK' : {'MUTATED_CHAIN' : 'A'},
	'3PSG' : {'MUTATED_CHAIN' : 'A'},
	'3SSI' : {'MUTATED_CHAIN' : 'A'},
	'3VUB' : {'MUTATED_CHAIN' : 'A'},
	'451C' : {'MUTATED_CHAIN' : 'A'},
	'4GCR' : {'MUTATED_CHAIN' : 'A'},
	'4LYZ' : {'MUTATED_CHAIN' : 'A'},
	'4TLN' : {'MUTATED_CHAIN' : 'A'},
	'4TMS' : {'MUTATED_CHAIN' : 'A'},
	'5CPV' : {'MUTATED_CHAIN' : 'A'},
	'5PEP' : {'MUTATED_CHAIN' : 'A'},
	'6TAA' : {'MUTATED_CHAIN' : 'A'},
	'9PCY' : {'MUTATED_CHAIN' : 'A'},
}

# These PDB files have more than one chain but are sequence-identical and ProTherm sometimes lists the wrong chain e.g. '-' rather than 'A'
identicalChainPDBs = {
	'1AAR' : {'MUTATED_CHAIN' : 'A'},
	'1ARR' : {'MUTATED_CHAIN' : 'A'},
	'1AZP' : {'MUTATED_CHAIN' : 'A'},
	'1B26' : {'MUTATED_CHAIN' : 'A'},
	'1BNI' : {'MUTATED_CHAIN' : 'A'},
	'1FC1' : {'MUTATED_CHAIN' : 'A'},
	'1G6N' : {'MUTATED_CHAIN' : 'A'},
	'1HFY' : {'MUTATED_CHAIN' : 'A'},
	'1HFZ' : {'MUTATED_CHAIN' : 'A'},
	'1LRP' : {'MUTATED_CHAIN' : 'A'},
	'1N0J' : {'MUTATED_CHAIN' : 'A'},
	'1RGG' : {'MUTATED_CHAIN' : 'A'},
	#'1RN1' : {'MUTATED_CHAIN' : 'B'}, # Three identical chains A, B, C. This case is an odd one - we choose chain B since while the PDB file contains identical chains, residue 45 is missing in chain A but required for records 10057 and 10058.
	'1WQ5' : {'MUTATED_CHAIN' : 'A'},
	'2AFG' : {'MUTATED_CHAIN' : 'A'},
	'2TRX' : {'MUTATED_CHAIN' : 'A'},
	'2ZTA' : {'MUTATED_CHAIN' : 'A'},
	'1HTI' : {'MUTATED_CHAIN' : 'A'},
}

# In these cases, the data in ProTherm is incorrect according to the publication
overridden = {
	# These cases have ambiguous entries
	
	3047  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1TUP'}, # Two identical chains A, B, and C (and two DNA chains E and F) but '-' specified in ProTherm
	3048  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1TUP'},
	3049  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1TUP'},
	3050  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1TUP'},
	3051  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1TUP'},
	3469  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1G6N'}, # Two identical chains A and B but '-' specified in ProTherm
	3470  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1G6N'},
	2418  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'}, # Two identical chains A and B but '-' specified in ProTherm
	5979  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5980  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5981  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5982  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5983  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5984  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5985  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5986  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},
	5987  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'},	
	3629  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'}, # Two identical chains A and B but '-' specified in ProTherm
	3630  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3631  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3632  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3633  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3634  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3635  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3636  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3637  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3638  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3639  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3640  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3641  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3642  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3643  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3644  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'},
	3604  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'}, # Three identical chains A, B, and C but '-' specified in ProTherm
	3605  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	3606  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	3607  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	3608  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	3609  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	3610  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	3611  : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13412 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13413 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13414 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13415 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13416 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13417 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13418 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13419 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13420 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13985 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13986 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	13421 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	14253 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	14254 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	14255 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'},
	8302  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'}, # Four identical chains A, B, C, and D but '-' specified in ProTherm
	8303  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	8304  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	8305  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	8306  : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	11864 : {'dCp' 			 : None,'PDB' : '2FHA'},
	11865 : {'dCp' 			 : None,'PDB' : '2FHA'},
	11866 : {'dCp' 			 : None,'PDB' : '2FHA'},
	11867 : {'dCp' 			 : None,'PDB' : '2FHA'},
	24388 : {'dCp' 			 : None,'PDB' : ''},
	16900 : {'m'			 : '6.86 kJ/mol/M', 'PDB' : '1C9O'},
	
	889   : {'Tm'			 : '53.4 C',	'PDB' : '1ARR'},
	890   : {'Tm'			 : '67.3 C',	'PDB' : '1ARR'},
	5303  : {'Tm'			 : '<= 10.0 C',	'PDB' : '1YCC'},
	23589  : {'Tm'			 : '> 80 C', 	'PDB' : ''},
	25269  : {'Tm'			 : '52-54 C',	'PDB' : ''},
	
	17877  : {'dTm'			 : None,	'PDB' : ''},
	23676  : {'dTm'			 : None,	'PDB' : ''},
	
	14474 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	14475 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	14476 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	24298 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'}, # As above but ProTherm has no entry at all
	24299 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	24300 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	24301 : {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'},
	14153 : {'MUTATED_CHAIN' : 'A', 'PDB' : '1N0J'}, # Two identical chains A and B but '-' specified in ProTherm
	17394 : {'MUTATED_CHAIN' : 'B', 'PDB' : '1OTR'}, # Two distinct chains. There's no mutation here so there's no harm in specifying B as the chain.
	17395 : {'MUTATED_CHAIN' : 'B', 'PDB' : '1OTR'}, # 
	17396 : {'MUTATED_CHAIN' : 'B', 'PDB' : '1OTR'}, # 
	8498  : {'ddG_H2O' 		 : '-2.82', 		'PDB' : '1TEN'}, # Typo 
	14192 : {'ddG' 		 	 : '-0.6453154876', 'PDB' : '1LZ1'}, # Bad computation
	15807 : {'ddG_H2O'			 : '2.72', 'PDB' : '1FKJ'}, # Wrong sign
	15808 : {'ddG_H2O'			 : '2.35', 'PDB' : '1FKJ'}, # Wrong sign
	17873 : {'ddG_H2O'		 : None, 'PDB' : '1RN1'},# I would need to check the reference
	
	# PMID: 11513583 - dG values entered as ddG
	11745 : {'ddG_H2O' :  '3.1', 'dG_H2O' : '12.1', 'PDB' : '2TRX'}, # D26I is a stabilizing mutation
	11746 : {'ddG_H2O' : '-3.7', 'dG_H2O' : '5.3',  'PDB' : '2TRX'}, # Destabilizing mutation
	11747 : {'ddG_H2O' : '-3.1', 'dG_H2O' : '5.9',  'PDB' : '2TRX'}, # Destabilizing mutation
	11748 : {'ddG_H2O' : '-1.4', 'dG_H2O' : '7.6',  'PDB' : '2TRX'}, # Destabilizing mutation
	11749 : {'ddG_H2O' : '-1.1', 'dG_H2O' : '7.9',  'PDB' : '2TRX'}, # Destabilizing mutation
	
	# PMID:11695900 - ProTherm appears to use the wrong sign
	12235 : {'ddG_H2O' :  '-11.9', 'PDB' : '1OH0'}, # Y14S  is a destabilizing mutation
	12236 : {'ddG_H2O' :  '-13.7', 'PDB' : '1OH0'}, # Destabilizing mutation
	12237 : {'ddG_H2O' :   '-9.5', 'PDB' : '1OH0'}, # Destabilizing mutation

	
	# PMID:8652517 - ProTherm appears to use the wrong sign
	5038 : {'ddG_H2O' :  '1.2', 'PDB' : '1YCC'}, # N52I is a stabilizing mutation
	
	# PMID:11964251 - ProTherm uses kJ/mol but does not specify them
	13086 : {'ddG_H2O' :   '-9.9 kJ/mol', 'PDB' : '5AZU'},
	13087 : {'ddG_H2O' :   '-10.9 kJ/mol', 'PDB' : '5AZU'},
	
	# PMID:18189416 - ProTherm uses the wrong sign # todo: put this in a separate dict so we can check that the values are the same as ProTherm, just the opposite sign i.e. sanity-check against bad data entry here
	23154 : {'ddG_H2O' :  '-9.2', 'PDB' : '1CAH'},
	23155 : {'ddG_H2O' :  '-8.3', 'PDB' : '1CAH'},
	23156 : {'ddG_H2O' :  '-8.6', 'PDB' : '1CAH'},
	23157 : {'ddG_H2O' :  '-5.4', 'PDB' : '1CAH'},
	23158 : {'ddG_H2O' :  '-4.7', 'PDB' : '1CAH'},
	23159 : {'ddG_H2O' :  '-3.0', 'PDB' : '1CAH'},
	23160 : {'ddG_H2O' :  '-3.2', 'PDB' : '1CAH'},
	23161 : {'ddG_H2O' :  '-1.4', 'PDB' : '1CAH'},
	23162 : {'ddG_H2O' :  '-2.8', 'PDB' : '1CAH'},
	23163 : {'ddG_H2O' :  '-0.8', 'PDB' : '1CAH'},
	23164 : {'ddG_H2O' :  '-1.8', 'PDB' : '1CAH'},
	
	# PMID:8771183 - ProTherm uses the wrong sign # todo: put this in a separate dict so we can check that the values are the same as ProTherm, just the opposite sign i.e. sanity-check against bad data entry here
	11246 : {'ddG_H2O' :  '-1.92', 'PDB' : '4LYZ'},
	11247 : {'ddG_H2O' :  '-1.27', 'PDB' : '4LYZ'},
	11248 : {'ddG_H2O' :  '-2.00', 'PDB' : '4LYZ'},
	
	
}

RoundingErrors = {
	# PMID: 1404369. Rounding errors
	180 : {'ddG' 		 	 : '-0.71', 'PDB'	:	'1BNI'},
	179 : {'ddG' 		 	 : '-0.69', 'PDB'	:	'1BNI'},
	186 : {'ddG' 		 	 : '-0.91', 'PDB'	:	'1BNI'},
	181 : {'ddG' 		 	 : '-0.78', 'PDB'	:	'1BNI'},
	183 : {'ddG' 		 	 : '-0.81', 'PDB'	:	'1BNI'},
	172 : {'ddG' 		 	 : '-0.19', 'PDB'	:	'1BNI'},
	174 : {'ddG' 		 	 : '-0.35', 'PDB'	:	'1BNI'},
	173 : {'ddG' 		 	 : '-0.31', 'PDB'	:	'1BNI'},
	178 : {'ddG' 		 	 : '-0.66', 'PDB'	:	'1BNI'},
	189 : {'ddG' 		 	 : '-4.08', 'PDB'	:	'1BNI'},
	176 : {'ddG' 		 	 : '-0.48', 'PDB'	:	'1BNI'},
	171 : {'ddG' 		 	 : '-0.14', 'PDB'	:	'1BNI'},
	175 : {'ddG' 		 	 : '-0.41', 'PDB'	:	'1BNI'},
	182 : {'ddG' 		 	 : '-0.79', 'PDB'	:	'1BNI'},
	185 : {'ddG' 		 	 : '-0.88', 'PDB'	:	'1BNI'},
	187 : {'ddG' 		 	 : '-0.98', 'PDB'	:	'1BNI'},
	184 : {'ddG' 		 	 : '-0.82', 'PDB'	:	'1BNI'},
	
	# PMID: 8358293. Rounding errors
	13328 : {'ddG'	: '1.05', 'PDB' : '1BVC'},
	13329 : {'ddG'	: '0.75', 'PDB' : '1BVC'},
	13330 : {'ddG'	: '0.59', 'PDB' : '1BVC'},
	13331 : {'ddG'	: '0.16', 'PDB' : '1BVC'},
	13332 : {'ddG'	: '-0.67', 'PDB' : '1BVC'},
	13333 : {'ddG'	: '-0.26', 'PDB' : '1BVC'},
	13334 : {'ddG'	: '-0.26', 'PDB' : '1BVC'},
	13335 : {'ddG'	: '-0.44', 'PDB' : '1BVC'},
	13336 : {'ddG'	: '-0.41', 'PDB' : '1BVC'},
	13337 : {'ddG'	: '-1.12', 'PDB' : '1BVC'},
	13338 : {'ddG'	: '-1.78', 'PDB' : '1BVC'},
	13339 : {'ddG'	: '-1.45', 'PDB' : '1BVC'},
	13340 : {'ddG'	: '-1.41', 'PDB' : '1BVC'},
	13341 : {'ddG'	: '-1.60', 'PDB' : '1BVC'},
	13342 : {'ddG'	: '-1.92', 'PDB' : '1BVC'},

	# More from this reference needs to be updated
	# PMID: 1569557. Precision lost during data entry.
	96 : {'ddG'	: '-2.07', 'PDB' : '1BNI'},
	40 : {'ddG'	: '-1.35', 'PDB' : '1BNI'},
	41 : {'ddG'	: '-1.85', 'PDB' : '1BNI'},
	42 : {'ddG'	: '-1.24', 'PDB' : '1BNI'},
	43 : {'ddG'	: '-2.15', 'PDB' : '1BNI'},
	44 : {'ddG'	: '-0.89', 'PDB' : '1BNI'},
	45 : {'ddG'	: '-2.48', 'PDB' : '1BNI'},
	46 : {'ddG'	: '-3.39', 'PDB' : '1BNI'},
	47 : {'ddG'	: '-0.31', 'PDB' : '1BNI'},
	48 : {'ddG'	: '-3.34', 'PDB' : '1BNI'},
	49 : {'ddG'	: '-4.32', 'PDB' : '1BNI'},
	50 : {'ddG'	: '-1.68', 'PDB' : '1BNI'},
	51 : {'ddG'	:  '0.54', 'PDB' : '1BNI'},
	52 : {'ddG'	: '-2.03', 'PDB' : '1BNI'},
	53 : {'ddG'	: '-2.25', 'PDB' : '1BNI'},
	54 : {'ddG'	: '-0.02', 'PDB' : '1BNI'},
	55 : {'ddG'	: '-1.12', 'PDB' : '1BNI'},
	56 : {'ddG'	: '-3.52', 'PDB' : '1BNI'},
	57 : {'ddG'	: '-1.46', 'PDB' : '1BNI'},
	58 : {'ddG'	: '-1.94', 'PDB' : '1BNI'},
	59 : {'ddG'	: '-0.44', 'PDB' : '1BNI'},
	60 : {'ddG'	: '-1.76', 'PDB' : '1BNI'},
	# Missing
	65 : {'ddG'	: '-1.15', 'PDB' : '1BNI'},
	# Missing
	67 : {'ddG'	: '-1.75', 'PDB' : '1BNI'},
	68 : {'ddG'	: '-2.44', 'PDB' : '1BNI'},
	# Missing
	71 : {'ddG'	: '-2.97', 'PDB' : '1BNI'},
	72 : {'ddG'	: '-2.42', 'PDB' : '1BNI'},
	73 : {'ddG'	: '-0.27', 'PDB' : '1BNI'},
	74 : {'ddG'	: '-1.15', 'PDB' : '1BNI'},
	# Missing
	77 : {'ddG'	:  '0.47', 'PDB' : '1BNI'},
	# Missing
	79 : {'ddG'	: '-0.82', 'PDB' : '1BNI'},
	80 : {'ddG'	: '-1.89', 'PDB' : '1BNI'},
	81 : {'ddG'	: '-1.65', 'PDB' : '1BNI'},
	82 : {'ddG'	: '-1.35', 'PDB' : '1BNI'},
	83 : {'ddG'	: '-2.02', 'PDB' : '1BNI'},
	84 : {'ddG'	: '-1.34', 'PDB' : '1BNI'},
	# Missing
	88 : {'ddG'	: '-1.93', 'PDB' : '1BNI'},
	# Missing
	90 : {'ddG'	: '-0.88', 'PDB' : '1BNI'},
	91 : {'ddG'	: '-3.17', 'PDB' : '1BNI'},
	92 : {'ddG'	: '-2.67', 'PDB' : '1BNI'},
	# Missing
	94 : {'ddG'	: '-2.24', 'PDB' : '1BNI'},
	95 : {'ddG'	: '-0.76', 'PDB' : '1BNI'},
	# Missing
	
	# More from this reference needs to be updated
	# PMID: 1404369. Precision lost during data entry.
	171 : {'ddG'	: '-0.14', 'PDB' : '1BNI'}, # A 32 R ARG
	172 : {'ddG'	: '-0.19', 'PDB' : '1BNI'}, # A 32 K LYS
	173 : {'ddG'	: '-0.31', 'PDB' : '1BNI'}, # A 32 M MET
	174 : {'ddG'	: '-0.35', 'PDB' : '1BNI'}, # A 32 L LEU
	175 : {'ddG'	: '-0.41', 'PDB' : '1BNI'}, # A 32 S SER
	176 : {'ddG'	: '-0.48', 'PDB' : '1BNI'}, # A 32 Q GLN
	177 : {'ddG'	: '-0.55', 'PDB' : '1BNI'}, # A 32 E GLU
	178 : {'ddG'	: '-0.66', 'PDB' : '1BNI'}, # A 32 N ASN
	179 : {'ddG'	: '-0.69', 'PDB' : '1BNI'}, # A 32 F PHE
	180 : {'ddG'	: '-0.71', 'PDB' : '1BNI'}, # A 32 D ASP
	181 : {'ddG'	: '-0.78', 'PDB' : '1BNI'}, # A 32 H HIS
	182 : {'ddG'	: '-0.79', 'PDB' : '1BNI'}, # A 32 T THR
	183 : {'ddG'	: '-0.81', 'PDB' : '1BNI'}, # A 32 I ILE
	184 : {'ddG'	: '-0.82', 'PDB' : '1BNI'}, # A 32 Y TYR
	185 : {'ddG'	: '-0.88', 'PDB' : '1BNI'}, # A 32 V VAL
	186 : {'ddG'	: '-0.91', 'PDB' : '1BNI'}, # A 32 G GLY
	187 : {'ddG'	: '-0.98', 'PDB' : '1BNI'}, # A 32 W TRP
	188 : {'ddG'	: '-1.00', 'PDB' : '1BNI'}, # A 32 C CYS
	189 : {'ddG'	: '-4.08', 'PDB' : '1BNI'}, # A 32 P PRO
	
	
	# More from these references need to be updated
	19236 : {'ddG_H2O'	: str(-6.0/4.184), 'PDB' : '5AZU'},
	21991 : {'ddG_H2O'	: str(-10.0/4.184), 'PDB' : '5AZU'},
}
PMIDReferencesInWhichToFixDDGPrecision = [
	1569557, "Protherm 19236", "ProTherm 21991"

]

for ID, data in RoundingErrors.iteritems():
	assert(ID not in overridden.keys())
	overridden[ID] = data
	

# Records with the wrong PMID
badPublicationReferences = {}
for recordID in range(13376, 13381 + 1):
	badPublicationReferences[recordID] = 8390295
for recordID in range(15714, 15733 + 1):
	badPublicationReferences[recordID] = 12473461

badSecondaryStructure = dict.fromkeys(
	[2747, 4611, 12310, 12701, 12702, 12979, 12980, 12982, 12983, 16895, 19886, 19887, 19888, 19889, 19893, 22231, 24335]
	+ range(15529, 15534 + 1)
	+ range(24921, 24929 + 1)
	+ range(24931, 24939 + 1)
	+ range(24962, 24964 + 1) + [24966]
	+ range(24968, 24983 + 1)
	+ range(24985, 25000 + 1)
	,True)

# In these cases, the protein is elongated at the 67th position. This is more than a mutation so I ignore it. 	
skipTheseCases = [12156, 12159, 12161, 12188, 12191, 12193, 14468]

# Note: The pair of records 12193 and 14468 are one example of duplicated data

#These cases fail parsing the mutation line - need to write a new regex
#skipTheseCases.extend([19893,19894,19895])

# In these cases, the structural information needed for the mutations (residues A57, A62) is missing in the PDB file
# Some of the information is present in chain B (A57's corresponding residue) but I decided against remapping the residues. 
skipTheseCases.extend([13451, 13452])

# In this case, I couldn't map the paper's structure well enough onto 1YCC. The mutations mentioned e.g. A57 N->I do not correspond with the PDB structure (attributed to a later paper).
skipTheseCases.append(11817)

# The DDG value here is a limit rather than a value (the paper states that the mutant (G107V) "was too unstable to yield an accurate determination of DDG_H20".
skipTheseCases.append(2043)

skipTheseCases = set(skipTheseCases)

# Mutations involving cysteine which have a different format involving bridges (S-H or S-S)
CysteineMutationCases = [13663, 13664, 13677, 13678]

# Mutations with different parsing requirements and their regexes
multimapCases1 = [16597, 16598, 16599, 16600, 19893, 19894, 22383, 22384]
mmapCases1regex = re.compile("PDB:(.+[,]{0,1})+\)")

multimapCases2Ranges = ((17681, 17687), (17692, 17698), (18104, 18105), (18108, 18136), (18138, 18175))
multimapCases2 = []
for m2r in multimapCases2Ranges:
	multimapCases2.extend(range(m2r[0], m2r[1] + 1))
#multimapCases2 = range(17681, 17687 + 1) + range(17692, 17698 + 1) + range(18104, 18105 + 1) +
#				  range(18108, 18136 + 1) + range(18138, 18175 + 1)
mmapCases2regex = re.compile("^.*PDB:(.*);PIR.*$")	

		
multimapCases3Ranges = ((3352, 3383), (14215, 14223), (6366, 6368), (8470, 8504), (12235, 12237),
					(12308, 12310), (15402, 15408), (16251, 16253), (16255, 16257), 
					(16259, 16261), (16263, 16265))
multimapCases3 = []
for m3r in multimapCases3Ranges:
	multimapCases3.extend(range(m3r[0], m3r[1] + 1))
multimapCases3.extend([10384, 16991, 17678, 17679, 17680, 17689, 17690, 17691])
mmapCases3regex = re.compile("PDB:(\w.*?\w)[;)]")

# 18125 - The PIR mutation has the wrong mutant (it should be A, not I)

# A mapping for missing ProTherm references to their PubMed IDs
missingRefMap = {
	"BIOCHEMISTRY 34, 7094-7102 (1995)" 		: ("PMID", 7766619),
	"BIOCHEMISTRY 34, 7103-7110 (1995)" 		: ("PMID", 7766620), # This is probably the correct record. Pubmed has a different end page (two extra pages)
	"BIOCHEMISTRY 37, 2477-2487 (1998)" 		: ("PMID", 9485396),
	"BIOCHIM BIOPHYS ACTA 1429, 365-376 (1999)" : ("PMID", 9989221),
	"PROTEIN SCI 6, 2196-2202 (1997)" 			: ("PMID", 9336842),
	"J MOL BIOL 224, 819-835 (1992)"			: ("PMID", 1569559),
	"STRUCTURE 14, 1401-1410 (2006)"			: ("PMID", 16962971),
	"J AM CHEM SOC 115, 8523-8526 (1993) "		: ("PMID", 0), # No PMID for this article. "Phospholipase A2 Engineering.  10.  The Aspartate...Histidine Catalytic Diad Also Plays an Important Structural Role."  Y. Li and M.-D. Tsai, J. Am. Chem. Soc. 115, 8523-8526 (1993).
	"J MOL BIOL 351, 402-416 (2005)"			: ("PMID", 16002092),
	"PROTEIN SCI 16, 227-238 (2007)"			: ("PMID", 17189482),
} 

def getDDGUnitsUsedInDB(ddGDB):
	results = ddGDB.locked_execute('SELECT SourceID, DGUnitUsedInProTherm FROM ProThermUnits')
	unitsUsed = {}
	for r in results:
		unitsUsed[r["SourceID"]] = r["DGUnitUsedInProTherm"]
	return unitsUsed
	
def getIDsInDB(ddGDB = None, source = "ProTherm-2008-09-08-23581"):
	records_in_database = set([int(r[0]) for r in ddGDB.execute('SELECT SourceID FROM ExperimentScore INNER JOIN Experiment on ExperimentID=Experiment.ID WHERE Source=%s', parameters =(source,), cursorClass = ddgdbapi.StdCursor)])
	return records_in_database
	
def summarizeList(l):
	s = ""
	if len(l) > 1:
		oldn = l[0]
		s = "%d" % oldn
		for n in l[1:]:
			if n > oldn + 1:
				s += "-%d,%d" % (oldn, n)
			oldn = n
		s += "-%d" % n
	elif len(l) == 1:
		s = l[0]
	return s
		
	
class ProThermReader(object):

	updated_dates = {
		23581 : "2008-09-08",
		25616 : "2011-12-21",
	}
	
	def __init__(self, infilepath, ddGDB = None, quiet = False, skipIndexStore = False):
		
		if not ddGDB:
			ddGDB = ddgdbapi.ddGDatabase()
		self.ddGDB = ddGDB
		
		mtchs = re.match(".*(ProTherm)(\d+)[.]dat$", infilepath, re.IGNORECASE)
		if mtchs:
			lastrecord = int(mtchs.group(2))
		if lastrecord in [23581, 25616]:
			# These fields of ProTherm records cannot be empty for our purposes
			self.requiredFields = ["NO.", "PDB_wild", "LENGTH", "ddG", "MUTATION", "MUTATED_CHAIN"]
			self.quiet = quiet
			# For bad data
			self.iCodeRecords = iCodeRecords
			self.patch = patch
			self.patchthis = {}
			self.singleChainPDBs = singleChainPDBs
			self.identicalChainPDBs = identicalChainPDBs
			self.overridden = overridden
			self.skipTheseCases = skipTheseCases
			self.CysteineMutationCases = CysteineMutationCases
			self.multimapCases1 = multimapCases1
			self.multimapCases2 = multimapCases2
			self.multimapCases3 = multimapCases3
			self.mmapCases1regex = mmapCases1regex
			self.mmapCases2regex = mmapCases2regex
			self.mmapCases3regex = mmapCases3regex
			#self.mutationregex = re.compile("^(\w\d+\w,)+(\w\d+\w)[,]{0,1}$")
			self.singlemutationregex = re.compile("^(\w)(\d+\w?)(\w)$")
			self.missingRefMap = missingRefMap
			self.updated_date = ProThermReader.updated_dates[lastrecord]
			self.missingddGUnits = {}
			self.secondary_structure_values = secondary_structure_values
			# Experimental data
			self.missingExpData = {}
			self.maxDBfieldlengths = {}
			self.kelvin_regex = re.compile("^(\d*[.]?\d*) K$")
			self.celsius_regex = re.compile("^(\d*[.]?\d*) C$")
			# References
			self.missingReferences = 0
			self.noPMIDs = {}
			self.ExistingDBIDs = {}
			self.ddGUnitsUsed = getDDGUnitsUsedInDB(self.ddGDB)
			
			self.ExistingScores = {}		
		else:
			raise Exception("No patch data is available for %s. Run a diff against the most recent database to determine if any changes need to be made." % infilepath)
		
		self.fieldnames = None
		self.infilepath = infilepath
		self.indices = None
		self.fhandle = None
		
		self.open()
		self.readFieldnames()
		self.singleErrors = dict.fromkeys(self.fieldnames, 0)
		if not skipIndexStore:
			if not self.quiet:
				colortext.write("[Storing indices for %s: " % self.infilepath, "green")
			colortext.flush()
			self.storeRecordIndices()
			if not quiet:
				colortext.printf("done]", "green")
			colortext.flush()
		self.close()
		
		if not quiet:
			self.printSummary()
		
		self.setUpExperimentalFieldNames()
		
	def setUpExperimentalFieldNames(self):
		dbAssayFieldNames = self.ddGDB.FieldNames.ExperimentAssay
		self.exp2DBfield = {
			"T"				: dbAssayFieldNames.Temperature,
			"pH"			: dbAssayFieldNames.pH,
			"BUFFER_NAME"	: dbAssayFieldNames.Buffer,
			"BUFFER_CONC"	: dbAssayFieldNames.BufferConcentration,
			"ION_NAME_1"	: dbAssayFieldNames.Ion1,
			"ION_NAME_1C"	: dbAssayFieldNames.Ion1,
			"ION_CONC_1"	: dbAssayFieldNames.Ion1Concentration,
			"ION_NAME_2"	: dbAssayFieldNames.Ion2, 
			"ION_CONC_2"	: dbAssayFieldNames.Ion2Concentration,
			"ION_NAME_3"	: dbAssayFieldNames.Ion3,
			"ION_CONC_3"	: dbAssayFieldNames.Ion3Concentration, 
			"ADDITIVES"		: dbAssayFieldNames.Additives,
			"PROTEIN_CONC"	: dbAssayFieldNames.ProteinConcentration,
			"MEASURE"		: "Measure",
			"METHOD"		: "MethodOfDenaturation",
		}
		self.expfields = self.exp2DBfield.keys()
		self.numericexpfields = ["T", "pH"]  
	
	def open(self):
		if not self.fhandle:
			self.fhandle = open(self.infilepath, "r")
		else:
			raise Exception("Trying to reopen an open file.")
		
	def close(self):
		if self.fhandle:
			self.fhandle.close()
			self.fhandle = None
		else:
			raise Exception("Trying to close a null file handle.")
	
	def test(self):
		success = True
		expected_results = {
		# PLAIN
			# L 121 A, A 129 M, F 153 L
			1163  : [('L',  '121', 'A'), ('A',  '129', 'M'), ('F',  '153', 'L')],
			# S 31 G, E 33 A, E 34 A
			1909  : [('S',   '31', 'G'), ('E',   '33', 'A'), ('E',   '34', 'A')],
			#Q 15 I, T 16 R, K 19 R, G 65 S, K 66 A, K 108 R
			2227  : [('Q',   '15', 'I'), ('T',   '16', 'R'), ('K',   '19', 'R'), ('G',   '65', 'S'), ('K',   '66', 'A'), ('K',  '108', 'R')],
			# I 6 R, T 53 E, T 44 A
			12848 : [('I',    '6', 'R'), ('T',   '53', 'E'), ('T',   '44', 'A')],
			# E 128 A, V 131 A, N 132 A
			17608 : [('E',  '128', 'A'), ('V',  '131', 'A'), ('N',  '132', 'A')],
		# PLAIN BUT BADLY FORMATTED
			11146 : [('G',   '23', 'A'), ('G',   '25', 'A')],
			21156 : [('E',    '3', 'R'), ('F',   '15', 'A'), ('D',   '25', 'K')],
		# MUTATION (PDB: MUTATION; PIR MUTATION), repeat
			#S 15 R (PDB: S 14 R; PIR: S 262 R), H 19 E (PDB: H 18 E; PIR: H 266 E), N 22 R ( PDB: N 21 R; PIR: N 269 R)
			14215 : [('S',   '14', 'R'), ('H',   '18', 'E'), ('N',   '21', 'R')],
		# MUTATION_LIST (PDB: MUTATION_LIST; PIR MUTATION_LIST)
			# N 51 H, D 55 H (PDB: N 47 H, D 51 H; PIR: N 135 H, D 139 H)
			17681 : [('N',   '47', 'H'), ('D',   '51', 'H')],
		# MUTATION_LIST (PDB: MUTATION_LIST)
			# A 2 K, L 33 I (PDB: A 1 K, L 32 I)
			22383 : [('A',    '1', 'K'), ('I',   '23', 'V')],
			16597 : [('G',   '92', 'P')],
			19893 : [('Q',  '11A', 'E'), ('H',  '53A', 'E'), ('S',  '76A', 'K'), ('M',  '78A', 'K'), ('D',  '81A', 'K')],
		}
		
		errors = []
		colortext.write("Testing %d mutation formats: " % len(expected_results), "green")	
		for ID, expected in sorted(expected_results.iteritems()):
			numerrormsgs = len(errors)
			for i in range(len(expected)):
				expected[i] = {
					"WildTypeAA"	: expected[i][0],
					"ResidueID"		: expected[i][1],
					"MutantAA"		: expected[i][2]
				}
			record = self._getRecord(ID, None)
			mutations, secondary_structure_positions = self.getMutations(ID, record)
			if expected and (expected != mutations):
				errors.append("Record %d: Expected %s, got %s." % (ID, expected, mutations))
			if (not badSecondaryStructure.get(ID)):
				if len(mutations) != len(secondary_structure_positions):
					errors.append("%s: \n" % record["MUTATION"])
					errors.append("Record %d: Read %d mutations but %d secondary structure positions." % (ID, len(mutations), len(secondary_structure_positions)))
					errors.append(secondary_structure_positions)
			if numerrormsgs < len(errors):
				colortext.write(".", "red")
			else:
				colortext.write(".", "green")
		if errors:
			colortext.error(" [failed]")
			colortext.error("Failed mutation parsing test.")
			for e in errors:
				colortext.error("\t%s" % e)
			success = False
		else:
			colortext.message(" [passed]")
			
		expected_results = {
			897   : -0.08, # -0.08
			5980  : 1.1 / 4.184, # 1.1 kJ/mol
			12144 : -3.9 / 4.184, #-3.9 kJ/mol
			13535 : -0.14, #-0.14
			14451 : -20.7 / 4.184, # -20.7 kJ/mol
			14751 : 3.35, # 3.35
			17849 : 12.7, # 12.7 kcal/mol
			17858 : 0.06, # 0.06 kcal/mol
			20129 : 0.08, # 0.08 
			21100 : 9.5 / 4.184, # 9.5 kJ/mol,
			22261 : -0.4 / 4.184, # -0.4 kJ/mol
			23706 : -2.8, # -2.8
			25089 : -3.3, # -3.3
		}
		errors = []
		colortext.write("Testing %d ddG conversions: " % len(expected_results), "green")	
		for ID, expectedDDG in sorted(expected_results.iteritems()):
			numerrormsgs = len(errors)
			record = self._getRecord(ID, None)
			referenceID = self.getReference(ID, record)
			record["dbReferencePK"] = "PMID:%s" % referenceID
			ddG = self.getDDGInKcal(ID, record)
			if expectedDDG != ddG:
				print("Record %d: Expected %s, got %s." % (ID, str(expectedDDG), str(ddG)))
				errors.append("Record %d: Expected %f, got %f." % (ID, expectedDDG, ddG))
			if numerrormsgs < len(errors):
				colortext.write(".", "red")
			else:
				colortext.write(".", "green")
		if errors:
			colortext.error(" [failed]")
			colortext.error("Failed ddG parsing test.")
			for e in errors:
				colortext.error("\t%s" % e)
			success = False
		else:
			colortext.message(" [passed]")
		return success
					
	def printSummary(self):
		colortext.printf("File %s: " % self.infilepath, "green")
		colortext.printf("Number of records: %d" % len(self.indices), "yellow")
		colortext.printf("First record ID: %d" % sorted(self.indices.keys())[0], "yellow")
		colortext.printf("Last record ID: %d" % sorted(self.indices.keys())[-1], "yellow")
		i = 1
		count = 0
		while i < self.list_of_available_keys[-1]:
			if i not in self.set_of_available_keys:
				count += 1
				#colortext.write("%d " % i)
				#colortext.flush()
			i += 1
		colortext.printf("Missing record IDs in this range: %d" % count, "yellow")
		colortext.write("\n")
				
			
	def readFieldnames(self):
		'''Reads fieldnames from the first record in the file without affecting the current position.''' 
		fhandle = self.fhandle
		oldpos = fhandle.tell()
		fhandle.seek(0)
		fieldnames = []
		line = fhandle.readline()
		while line and not(line.startswith("//")):
			if line[0] != "*" and line[0] != " ":
				fieldnames.append(line.split()[0])
			line = fhandle.readline()
		fhandle.seek(oldpos)
		self.fieldnames = fieldnames
		
	def saveIndicesToFile(self, filename):
		import pickle
		F = open(filename, "w")
		F.write(pickle.dumps(self.indices))
		F.close()

	def loadIndicesFromFile(self, filename):
		import pickle
		F = open(filename, "r")
		self.indices = pickle.loads(F.read())
		F.close()
		ikeys = self.indices.keys()
		self.list_of_available_keys = sorted(ikeys)
		self.set_of_available_keys = set(ikeys)
		
	def storeRecordIndices(self):
		'''Reads record indices from the file without affecting the current position.'''
		mustclose = False
		if not self.fhandle:
			self.open()
			mustclose = True
		fhandle = self.fhandle
		oldpos = fhandle.tell()
		fhandle.seek(0)
		self.indices = {}
		
		startofline = 0
		line = fhandle.readline()
		while line:
			if line.startswith("NO.  "):
				n = re.match("^NO[.]\s*(\d+)\s*$", line).group(1)
				assert(not(self.indices.get(n)))
				self.indices[int(n)] = startofline 
			startofline = fhandle.tell()
			line = fhandle.readline()
		ikeys = self.indices.keys()
		self.list_of_available_keys = sorted(ikeys)
		self.set_of_available_keys = set(ikeys)
		
		fhandle.seek(oldpos)
		if mustclose:
			self.close()
		
	def readRecordsConsecutively(self):
		'''Reads a record from the current position in the file.'''
		
		fhandle = self.fhandle
		record = {}
		lastkey = None
		line = fhandle.readline()
		while line and not(line.startswith("//")):
			if line[0] != "*":
				tokens = line.split()
				if len(tokens) > 1:
					v = join(tokens[1:], " ")
					if line[0] == " ":
						assert(lastkey)
						assert(lastkey in record.keys())
						record[lastkey] = "%s %s" % (record[lastkey], join(tokens[0:], " ").strip())
					else:
						k = tokens[0]
						try:
							assert(k in self.fieldnames)
						except:
							# This throws an error in ProTherm->23581 at record 19609 in the REMARKS field. I manually edited the entry, adding spaces to the record conformed to the rest of the database
							# In ProTherm->23581, records 21571 and 21739 may have an error (errant field called "kJ/mol") or that may have been an accidental local modification
							# In ProTherm->23581, record 22054 may have an error (errant field called "PDB_mutan") or that may have been an accidental local modification
							if k not in sometimesFields:
								colortext.error("Could not find data for field '%s' in record %s of %s." % (k, str(record.get("NO.")), self.infilepath))
								raise
						record[k] = v.strip()
				else:
					if line[0] == " ":
						assert(lastkey)
						assert(k in record.keys())
					else:
				 		k = tokens[0]
						try:
							assert(k in self.fieldnames)
						except:
							if k not in sometimesFields:
								colortext.error("Could not find data for field '%s' in record %s of %s." % (k, str(record.get("NO.")), self.infilepath))
								raise
						record[k] = None
				lastkey = k or lastkey
			line = fhandle.readline()
		return record
	
	def _getRecord(self, ID, record):
		if not record:
			if ID:
				record = self.readRecord(ID)
		if not record:
			raise Exception("Could not read a record with number #%s" % str(ID))
		return record	
	
	def searchFor(self, pdbID = None, chain = None, startIndex = None, endIndex = None, numMutations = None, hasddG = None, mutations = None, residueIDs = None):
		matching_records = []
		
		list_of_available_keys = sorted(list(self.set_of_available_keys.difference(set(skipTheseCases))))
		for ID in list_of_available_keys:
			record = self.readRecord(ID)
			self.fixRecord(ID, record)
			mts = None
			mtlocs = None
			try:
				mts, mtlocs = self.getMutations(ID, record)
			except:
				pass
			if mts and len(mts) == 1:
				#if str(mts[0]["ResidueID"]) == "32" and  m["WildTypeAA"] == 'A':
				#	print(ID, mts, record["ddG"], record["PDB_wild"])
				if str(mts[0]["ResidueID"]) == "45" and mts[0]["WildTypeAA"] == 'A' and mts[0]["MutantAA"] == 'G':
					print(ID, mts, record["ddG"], record["PDB_wild"])
			
			if (pdbID == None) or (pdbID != None and record["PDB_wild"] == pdbID):
				#if mts:
				#	print(ID, mts)
				if not(hasddG) or (hasddG == True and record["ddG"] != None):
					if (numMutations == None) or (mts and numMutations == len(mts)):
						foundResidueID = False
						if residueIDs:
							for m in mts:
								if m["ResidueID"] in residueIDs and m["WildTypeAA"] == 'A':
									foundResidueID = True 
									print(ID, mts)
						if not(residueIDs) or foundResidueID: 
							matching_records.append(ID)
		return matching_records

	def _printRecordSection(self, field, record, level = -1):
		if type(field) == type(""):
			if record.get(field):
				colortext.write("%s%s: " % ("   " * level, field), "yellow")
				print(record[field])
		else:
			if type(field[0]) == type(""):
				for f in field[1:]:
					if type(f) != type("") or record.get(f):
						colortext.printf("%s%s" % ("   " * level, field[0]), "cyan")
						break
			else:
				self._printRecordSection(field[0], record, level + 1)
			for field in field[1:]:
				self._printRecordSection(field, record, level + 1)
	
	def printRecord(self, ID, record = None):
		record = self._getRecord(ID, record)
		colortext.message("Record: %d" % ID)
		self._printRecordSection(field_order, record)
	
	def _getRecordHTMLSection(self, html, field, record, maxlevel, level = -1):
		#print(field)
		if type(field) == type(""):
			if record.get(field):
				extra=""
				if level < maxlevel:
					extra='colspan="%d"' % (maxlevel - level + 1)
				html.append('''<tr>%s<td style="background-color:#bbbbbb; border:1px solid black;" %s>%s</td>''' % ('''<td style="border:0"></td>''' * level, extra, field))
				# Wrap long strings
				splitstr = ""
				for i in range(0, len(record[field]), 80):
					splitstr = splitstr + record[field][i:i+80] + " "
				if field =="REFERENCE":
					refnum = self.getReference(0, record)
					if refnum:
						splitstr = '''<a target="ddGpubmed" href="http://www.ncbi.nlm.nih.gov/pubmed?term=%s[uid]">%s</a>''' % (refnum, splitstr)
				html.append('''<td style="background-color:#bbbbee; border:1px solid black;">%s</td></tr>''' % splitstr)
		else:
			if type(field[0]) == type(""):
				for f in field[1:]:
					if type(f) != type("") or record.get(f):
						html.append('''<tr>%s<td style="background-color:#993333;border:1px solid black;"><b>%s</b></td></tr>''' % ('''<td style="border:0"></td>''' * level, field[0]))
						
						if type(f) == type(""):
							extra=""
							if level < maxlevel:
								extra='colspan="%d"' % (maxlevel - level)
							html.append('''<tr>%s<td style="background-color:#bb8888; border:1px solid black;" %s><b>Field</b></td><td style="background-color:#bb8888; border:1px solid black;"><b>Value</b></td></tr>''' % ('''<td style="border:0"></td>''' * (level + 1), extra))
						break
			else:
				self._getRecordHTMLSection(html, field[0], record, maxlevel, level + 1)
			for field in field[1:]:
				self._getRecordHTMLSection(html, field, record, maxlevel, level + 1)
	
	def getRecordHTML(self, ID, record = None):
		record = self._getRecord(ID, record)
		html = []
		html.append('''<br><br><table width="900px" style="text-align:left"><tr><td align="center" colspan="4"><b>Record #%d</b></td></tr>''' % ID)
		self._getRecordHTMLSection(html, field_order, record, 2)
		html.append("</table>")
		return html
	
	def fixThermodynamic_dG(self, ID, record):
		# Leave these as are for the moment since we only really care about the DDG values
		pass
	
	def fixThermodynamic_dG_H2O(self, ID, record):
		# Leave these as are for the moment since we only really care about the DDG_H2O values
		pass
	
	def fixThermodynamic_Tm(self, ID, record):
		# Values are in degrees Celsius 
		Tm = record["Tm"]
		if Tm:
			if overridden.get(ID) and overridden[ID].get("Tm"):
				return
			if Tm.endswith(" K"):
				record["Tm"] = "%g C" % (float(Tm[:-2]) - NUMBER_KELVIN_AT_ZERO_CELSIUS)
			elif Tm.endswith("K"):
				record["Tm"] = "%g C" % (float(Tm[:-1]) - NUMBER_KELVIN_AT_ZERO_CELSIUS)
			elif Tm.endswith(" C"):
				TmInCelsius = float(Tm[:-2])
			else:
				try:
					TmInCelsius = float(Tm)
					if TmInCelsius > 200 and ID != 8530:
						colortext.warning("Record %d has a high Tm (%s) which may be in Celsius or Kelvin. Check." % (ID, TmInCelsius))
					record["Tm"] = "%s C" % record["Tm"]
				except:
					colortext.error("Error parsing Tm %s of record %d." % (Tm, ID))

	def fixThermodynamic_dTm(self, ID, record):
		# Values are in degrees Celsius 
		dTm = record["dTm"]
		if dTm:
			if dTm.endswith(" K") or dTm.endswith(" k"):
				dTm = dTm[:-2]
			elif dTm.endswith("K"):
				dTm = dTm[:-1]
			elif dTm.endswith(" C"):
				dTm = dTm[:-2]
			try:
				dTmInCelsius = float(dTm)
				record["dTm"] = dTm
			except:
				colortext.error("Error parsing dTm %s of record %d." % (record["dTm"], ID))
	
	def fixThermodynamic_dHvH(self, ID, record):
		# Probably in kcal/mol or kJ/mol unless otherwise specified but ProTherm does not explicitly specify the units. Unqualified values appear to be in both kJ/mol and kcal/mol. 
		pass
	
	def fixThermodynamic_dHcal(self, ID, record):
		# Probably in kcal/mol or kJ/mol unless otherwise specified but ProTherm does not explicitly specify the units. Unqualified values appear to be in both kJ/mol and kcal/mol. 
		pass
		
	def fixThermodynamic_m(self, ID, record):
		m = record["m"]
		m_unit = None
		if m:
			m = m.strip()
			tokens = [t for t in record["m"].split(' ') if t]
			tokens[0] = tokens[0].replace(",", "")
			assert(len(tokens) <= 3)
			if len(tokens) == 1:
				m = float(tokens[0])
				m_unit = 'kcal/mol/M'
			else:
				m = float(tokens[0])
				m_unit = join(map(string.strip, tokens[1:]), " ")
			record['m'] = m
			record['m_UNIT'] = m_unit

	def fixThermodynamic_Cm(self, ID, record):
		Cm = record["Cm"]
		if Cm:
			if Cm.endswith(" M"):
				Cm = Cm[:-2]
			if ID != 4919:
				# We usually expect a float value so we hack the odd cases
				floatvalue = float(Cm)
			record["Cm"] = Cm
	
	def fixThermodynamic_dCp(self, ID, record):
		dCp = record["dCp"]
		if dCp:
			kjIndex = dCp.find(" kJ/mol/K")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/Kmol")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/K/mol")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/K mol")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/K.mol")
			if kjIndex == -1:
				kjIndex = dCp.find("kJ/K/mol")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/mol K")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/m/K")
			if kjIndex == -1:
				kjIndex = dCp.find(" kJ/mol")
				
			if kjIndex != -1:
				dCp = kJtokcal(float(dCp[:kjIndex]))
			else:
				kcalIndex = dCp.find(" kcal/mol/K")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" kcal/mole/K")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" kcal/mol/deg")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" kcal/mol deg")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" Kcal/(mol.deg)")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" kcla/molK")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" kcal/K mol")
				if kcalIndex == -1:
					kcalIndex = dCp.find(" Kcal/K mol")
				if kcalIndex != -1:
					dCp = float(dCp[:kcalIndex])
				else:
					calIndex = dCp.find(" cal/K/mol")
					if calIndex != -1:
						dCp = float(dCp[:calIndex])/1000.0
					else:
						# We expect a float value here
						dCp = float(dCp)
			record["dCp"] = dCp
	
	def fixThermodynamic_state(self, ID, record):
		if record["STATE"]:
			record["STATE"] = int(record["STATE"])
			
	def fixThermodynamic_reversibility(self, ID, record):
		reversibility = record["REVERSIBILITY"]
		level = None
		if reversibility:
			reversibility = reversibility.lower()
			if reversibility == 'unknown' or reversibility == 'unknownnouwn':
				reversibility = 'Unknown'
			elif reversibility == 'yes, 90%':
				# Hack for special case
				reversibility = 'Yes'
				level = "90%"
			elif reversibility.startswith('yes'):
				level = reversibility[3:].strip() or None
				if level:
					assert(level[0] == '(' and level[-1] == ')')
					level = level[1:-1]
				reversibility = 'Yes'
			elif reversibility.startswith('no'):
				level = reversibility[2:].strip() or None
				if level:
					assert(level[0] == '(' and level[-1] == ')')
					level = level[1:-1]
				reversibility = 'No'
		else:
			reversibility = 'Unknown'
		record["REVERSIBILITY"] = reversibility
		record["REVERSIBILITY_LEVEL"] = level
	
	def fixThermodynamic_activity(self, ID, record):
		pass
	
	def fixThermodynamic_activity_Km(self, ID, record):
		pass
	
	def fixThermodynamic_activity_Kcat(self, ID, record):
		pass
	
	def fixThermodynamic_activity_Kd(self, ID, record):
		pass
	
	def fixRecord(self, ID, record = None):
		'''Override bad data in ProTherm.'''
		record = self._getRecord(ID, record)
		
		# Patch for typos
		if 24390 <= ID <= 24420 and record["ddG_H2O"] and record["ddG_H2O"].find("kal/mol") != -1:
			record["ddG_H2O"] = record["ddG_H2O"].replace("kal/mol", "kcal/mol")
		
		passed = True
		
		overridden = self.overridden
		singleChainPDBs = self.singleChainPDBs
		identicalChainPDBs = self.identicalChainPDBs
		patch = self.patch
		if overridden.get(ID):
			if record.get('PDB_wild'):
				if record["PDB_wild"] != overridden[ID]["PDB"]:
					raise colortext.Exception("Error in overridden table: Record %d. Read '%s' for PDB_wild, expected '%s'." % (ID, record["PDB_wild"], overridden[ID]["PDB"]))
			for k, v in overridden[ID].iteritems():
				if k != "PDB":
					record[k] = v
		if record["PDB_wild"]:
			pdbID = record["PDB_wild"].upper()
			if not(overridden.get(ID)) or not(overridden[ID].get("MUTATED_CHAIN")):
				# Update PDB IDs so long as we don't have them in the override dict 
				if singleChainPDBs.get(pdbID):
					for k,v in singleChainPDBs[pdbID].iteritems():
						record[k] = v
				elif identicalChainPDBs.get(pdbID):
					for k,v in identicalChainPDBs[pdbID].iteritems():
						record[k] = v
		
		# ** Experimental conditions **
		
		# Normalize measures values
		measures = []
		msrs = record['MEASURE']
		if msrs.find('+') != -1:
			msrs = msrs.split('+')
		elif msrs.find(',') != -1:
			msrs = msrs.split(',')
		else:
			msrs = [msrs]
		for m in msrs:
			m = m.strip()
			if not MeasureMapping.get(m):
				passed = False
				raise colortext.Exception("Cannot find match for measure '%s' in record %d." % (m, ID))
			else:
				measures.append(MeasureMapping[m])
		if len(measures) > 3:
				raise colortext.Exception("More than three measures (%s) were found for record %d. Our DDG database does not have enough fields to store all measures." % (str(measures), ID))
		record['MEASURES'] = measures 
		
		# Normalize methods values
		methods = []
		mthds = record['METHOD']
		if mthds != None:
			if mthds.find('+') != -1:
				raise colortext.Exception("Cannot find match for method '%s' in record %d." % (m, ID))
			elif mthds.find(',') != -1:
				mthds = mthds.split(',')
			else:
				mthds = [mthds]
			for m in mthds:
				m = m.strip()
				if not MethodMapping.get(m):
					passed = False
					raise colortext.Exception("Cannot find match for method '%s' in record %d." % (m, ID))
				else:
					methods.append(MethodMapping[m])
		if len(methods) > 2:
				raise colortext.Exception("More than two methods (%s) were found for record %d. Our DDG database does not have enough fields to store all methods." % (str(methods), ID))
		record['METHODS'] = methods 
		
		# ** Thermodynamic data **
		self.fixThermodynamic_dG(ID, record)
		self.fixThermodynamic_dG_H2O(ID, record)
		self.fixThermodynamic_Tm(ID, record)
		self.fixThermodynamic_dTm(ID, record)
		self.fixThermodynamic_dHvH(ID, record)
		self.fixThermodynamic_dHcal(ID, record)
		self.fixThermodynamic_m(ID, record)
		self.fixThermodynamic_Cm(ID, record)
		self.fixThermodynamic_dCp(ID, record)
		self.fixThermodynamic_state(ID, record)
		self.fixThermodynamic_reversibility(ID, record)
		self.fixThermodynamic_activity(ID, record)
		self.fixThermodynamic_activity_Km(ID, record)
		self.fixThermodynamic_activity_Kcat(ID, record)
		self.fixThermodynamic_activity_Kd(ID, record)
		
		missingFields = []
		for field in self.requiredFields:
			if not record[field]:
				passed = False
				missingFields.append(field)
		
		if not passed:
			# Recover when one field is missing
			if len(missingFields) == 1:
				# Recover when there is a ddG value and the mutation is specified
				# The patch dict lets us either skip entries with missing information or else correct them
				foundpatch = False
				if record["MUTATION"] and record["MUTATION"] != "wild" and ("ddG" not in missingFields):
					for patchfield, patchfield_info in patchfields.iteritems():
						if not record[patchfield]:
							foundpatch = patch.get(ID) and (patchfield in patch[ID].keys())
							if not patch.get(ID):
								colortext.error("Error processing record ID %d; no %s" %  (ID, patchfield))
								self.singleErrors[patchfield] += 1
								self.patchthis[ID] = "%s %s" % (patchfield, join(["%s" % record.get(pfi) for pfi in patchfield_info],"-"))
							elif patch[ID][patchfield]:
								record[patchfield] = patch[ID][patchfield]
								passed = True
							break
					if not foundpatch:
						#colortext.error("Error processing structure: ID %d, no %s " % (ID, missingFields[0]))
						self.singleErrors[missingFields[0]] += 1
		return passed
		
	def getMutations(self, ID, record = None):
		'''Either returns a list of triples (WildTypeAA, ResidueID, MutantAA) and a corresponding list of secondary structure locations or else None.'''
		record = self._getRecord(ID, record)
			
		mutations = []
		if ID == 24283:
			record["MUTATION"] = record["MUTATION"].replace("Y236 F", "Y 236 F") # Hack for this special case
			
		mutationline = record["MUTATION"].split()
		#print(len(mutationline))
		cline = join(mutationline, "")
		if ID in self.CysteineMutationCases: # Hack for input which includes information on the bonds generated from mutation to cysteine
			if not mutationline[1].isdigit():
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID)) # This will raise an exception if there is an insertion code - this is fine. I want to examine these cases manually when they occur.
			mutations = [{"WildTypeAA" : mutationline[0], "ResidueID" : mutationline[1], "MutantAA" : mutationline[2]}]  
		elif ID in self.multimapCases1:
			mtchslst = self.mmapCases1regex.findall(cline)
			if mtchslst:
				assert(len(mtchslst) == 1)
				mtchslst = mtchslst[0].split(',')
				for mtch in mtchslst:
					assert(mtch)
					residueID = mtch[1:-1]
					if not residueID.isdigit():
						if (not residueID[:-1].isdigit()) or (not residueID[-1].isalpha()):
							raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
					mutations.append({"WildTypeAA" : mtch[0], "ResidueID" : mtch[1:-1], "MutantAA" : mtch[-1]}) 
			else:
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
		elif ID in self.multimapCases2:
			mtchslst = self.mmapCases2regex.match(cline)
			if mtchslst:
				mtchslst = mtchslst.group(1).split(",")
				for mtch in mtchslst:
					assert(mtch)
					if not mtch[1:-1].isdigit():
						raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID)) # This will raise an exception if there is an insertion code - this is fine. I want to examine these cases manually when they occur.
					mutations.append({"WildTypeAA" : mtch[0], "ResidueID" : mtch[1:-1], "MutantAA" : mtch[-1]})
			else:
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
		elif ID in self.multimapCases3:
			mtchslst = self.mmapCases3regex.findall(cline)
			if mtchslst:
				for mtch in mtchslst:
					assert(mtch)
					if not mtch[1:-1].isdigit():
						raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID)) # This will raise an exception if there is an insertion code - this is fine. I want to examine these cases manually when they occur.
					mutations.append({"WildTypeAA" : mtch[0], "ResidueID" : mtch[1:-1], "MutantAA" : mtch[-1]}) 
			else:
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
		elif len(mutationline) == 3:
			if (not mutationline[1].isdigit()) and (not mutationline[1][:-1].isdigit()):
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
			mutations = [{"WildTypeAA" : mutationline[0], "ResidueID" : mutationline[1], "MutantAA" : mutationline[2]}] 
		elif len(mutationline) == 1:
			mline = mutationline[0]
			if mline == "wild" or mline == "wild*" or mline == "wild**":
				return [], []
			else:
				raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
		elif len(mutationline) % 3 == 0 or len(mutationline) == 5: # 2nd case is a hack for 11146
			cline = [entry for entry in cline.split(",") if entry] # Hack for extra comma in 21156
			for mutation in cline:
				m = self.singlemutationregex.match(mutation)
				if m:
					if (not m.group(2).isdigit()) and (not m.group(2)[:-1].isdigit()):
						raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
					mutations.append({"WildTypeAA" : m.group(1), "ResidueID" : m.group(2), "MutantAA" : m.group(3)})
				else:
					# todo: I realized after the fact that the regexes below do not deal properly with insertion codes in
					# the mutation e.g. "MUTATION        Y 27D D, S 29 D" in record 5438. However, none of these records
					# have ddG values for ProTherm on 2008-09-08 (23581 entries) so we can ignore this issue unless the
					# database is updated.
					# In this case, maybe singlemutationregex = re.compile("^(\w)(\d+\w?)(\w)$") would work
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
				 
		else:
			raise Exception("We need to add a case to handle this mutation string: '%s'." % cline)
		
		mutation_locations = []
		if (not badSecondaryStructure.get(ID)) and record["SEC.STR."]:
			mutation_locations = record["SEC.STR."].split(",")
			numlocations = len(mutation_locations)
			if not numlocations == len(mutations):
				#todo raise Exception
				colortext.error("The mutations '%s' do not have corresponding locations '%s' in record %d." % (mutations, mutation_locations, ID))
			for i in range(numlocations):
				n_location = self.secondary_structure_values.get(mutation_locations[i].strip())
				if not n_location:
					raise Exception("Bad secondary structure information '%s' for the mutation in record %d." % (record["SEC.STR."], ID))
				mutation_locations[i] = n_location
		
		return mutations, mutation_locations
	
	def getDDGInKcal(self, ID, record = None, useRosettaConvention = False, getDDGH2OInstead = False):
		record = self._getRecord(ID, record)
		dbReferencePK = record.get("dbReferencePK", "publication undetermined")
		
		ddGKey = "ddG"
		if getDDGH2OInstead:
			ddGKey = "ddG_H2O"
		
		ddGline = record[ddGKey]
		if getDDGH2OInstead and not ddGline:
			# todo: I should handle this better
			return None
		if ddGline.find("kJ/mol") == -1 and ddGline.find("kcal/mol") == -1:
			try:
				x = float(ddGline)
			except:
				colortext.error("Error processing %s: ID %d, %s" % (ddGKey, ID, record["ddG"]))
			if self.ddGUnitsUsed.get(dbReferencePK):
				unitsUsed = self.ddGUnitsUsed[dbReferencePK]
				# todo: These cases should be checked
				if unitsUsed[-1] == '?':
					unitsUsed = unitsUsed[:-1].strip()
				ddGline = "%s %s" % (ddGline, unitsUsed)
					
		idx = ddGline.find("kJ/mol")
		if idx != -1:
			try:
				val = kJtokcal(float(ddGline[0:idx]))
				if useRosettaConvention:
					return -val
				return val
			except:
				colortext.error("Error processing %s in kJ/mol: ID %d, %s" % (ddGKey, ID, record["ddG"]))
				return None
	
		idx = ddGline.find("kcal/mol")
		if idx != -1:
			try:
				val = float(ddGline[0:idx])
				if useRosettaConvention:
					return -val
				return val
			except:
				colortext.error("Error processing %s in kcal/mol: ID %d, %s" % (ddGKey, ID, record["ddG"]))
				return None
		
		mutationline = record["MUTATION"]
		if not(mutationline == "wild" or mutationline == "wild*" or mutationline == "wild**"):
			if self.ExistingDBIDs.get(ID):
				colortext.printf("No %s unit specified for existing record: ID %d, %s, publication ID='%s'" % (ddGKey, ID, ddGline, dbReferencePK), "pink")
			else:
				pass
				colortext.printf("No %s unit specified for new record: ID %d, %s, publication ID='%s'" % (ddGKey, ID, ddGline, dbReferencePK), "cyan")
			mutations, mutation_locations = self.getMutations(ID, record)
			mutations = join([join(map(str, m),"") for m in mutations], ",")
			#sys.exit(0)
			self.missingddGUnits[dbReferencePK] = self.missingddGUnits.get(dbReferencePK) or []
			self.missingddGUnits[dbReferencePK].append((ID, "*" + mutations + "=" + ddGline + "*"))
		return None

	def getReference(self, ID, record = None):
		record = self._getRecord(ID, record)
		
		if badPublicationReferences.get(ID):
			return badPublicationReferences[ID]
		
		# Parse reference
		if not record["REFERENCE"]:
			self.missingReferences += 1
		referenceData = record["REFERENCE"] or ""
		referenceData = referenceData.strip()
		idx = referenceData.find("PMID:")
		referenceID = None
		if idx != -1:
			refPMID = referenceData[idx+5:].strip()
			if refPMID:
				if refPMID[-1] == ".":
					refPMID = refPMID[0:-1] # Hack for record 25600
				if refPMID.isdigit():
					referenceID = int(refPMID)
				else:
					colortext.error("Check reference for record %(ID)d: %(referenceData)s. '%(refPMID)s' is not numeric." % vars())				
		if not referenceID:
			refdetails = None
			if idx == -1:
				refdetails = referenceData
			else:
				refdetails = referenceData[:idx]
			if self.missingRefMap.get(refdetails) and self.missingRefMap[refdetails][0] == "PMID":
				referenceID = self.missingRefMap[refdetails][1]
			else:
				self.missingReferences += 1
				authors = record["AUTHOR"]
				colortext.warning("No PMID reference for record %(ID)d: '%(referenceData)s', %(authors)s." % vars())
				self.noPMIDs[referenceData.strip()] = authors
		return referenceID
			
	def getExperimentalConditions(self, ID, record):
		record = self._getRecord(ID, record)
		
		exp2DBfield = self.exp2DBfield
		expfields = self.expfields
		
		#if record['MEASURE'] and
		if not record.get('MEASURES'):
			raise Exception("No measures stored for record %d." % ID)
		#if record['METHOD'] and
		if not record.get('METHODS'):
			raise Exception("No methods stored for record %d." % ID)

		ExperimentalConditions = {}
		for h in expfields:
			if not record[h]:
				self.missingExpData[h] = self.missingExpData.get(h, 0) + 1
			fielddata = record[h]
			if record[h]:
				self.maxDBfieldlengths[h] = max(self.maxDBfieldlengths.get(h, 0), len(fielddata))
			if fielddata == "Unknown":
				fielddata = None
			if fielddata and h == "METHOD":
				assert(record['METHODS'])
				for i in range(len(record['METHODS'])):
					ExperimentalConditions["%s%d" % (exp2DBfield[h], i + 1)] = record['METHODS'][i]
			elif fielddata and h == "MEASURE":
				assert(record['MEASURES'])
				for i in range(len(record['MEASURES'])):
					ExperimentalConditions["%s%d" % (exp2DBfield[h], i + 1)] = record['MEASURES'][i]
			elif fielddata and h in self.numericexpfields:
				try:
					fielddata = float(fielddata)
				except:
					if h == "T":
						K = self.kelvin_regex.match(fielddata)
						C = self.celsius_regex.match(fielddata)
						if K:
							fielddata = float(K.group(1)) - NUMBER_KELVIN_AT_ZERO_CELSIUS 
						elif C:
							fielddata = float(C.group(1))
						else:
							fielddata = None
				if fielddata == None:
					colortext.error("Failed to convert field %s's number %s for record %d." % (h, fielddata, ID))
				ExperimentalConditions[exp2DBfield[h]] = fielddata or None
			else:
				ExperimentalConditions[exp2DBfield[h]] = fielddata or None
		
		return ExperimentalConditions
	
	def readRecord(self, recordnumber):
		'''Reads a record from the current position in the file.'''
		if not recordnumber in self.indices.keys():
			if not self.quiet:
				colortext.error("Record %d not found" % recordnumber)
			return None
		openhandlehere = False
		if not self.fhandle:
			openhandlehere = True
			self.open()
		self.fhandle.seek(self.indices[recordnumber])
		record = self.readRecordsConsecutively()
		if openhandlehere:
			self.close()
			self.fhandle = False
		return record
	
	def oldPrintRecord(self, record):
		assert(sorted(record.keys()) == sorted(self.fieldnames))
		for fieldname in self.fieldnames:
			print("%s: %s" % (fieldname, record[fieldname]))

	def diff(self, secondDB, fields_to_ignore = [], start_at_id = 0, end_at_id = None, showall = False):
		mykeys = self.set_of_available_keys
		theirkeys = secondDB.set_of_available_keys
		
		fields_to_ignore = set(fields_to_ignore)
		
		removedKeys = sorted(list(mykeys.difference(theirkeys)))
		addedKeys = sorted(list(theirkeys.difference(mykeys)))
		if removedKeys:
			colortext.printf("Records in %s but not in %s:" % (self.infilepath, secondDB.infilepath), "green")
			colortext.printf(summarizeList(removedKeys))
		if addedKeys:
			colortext.printf("Records in %s but not in %s:" % (secondDB.infilepath, self.infilepath), "green")
			colortext.printf(summarizeList(addedKeys))
		
		removedFields = set(self.fieldnames) - set(secondDB.fieldnames) 
		addedFields = set(secondDB.fieldnames) - set(self.fieldnames) 
		if removedFields:
			colortext.printf("Fields in %s but not in %s:" % (self.infilepath, secondDB.infilepath), "green")
			colortext.printf(join(sorted(removedFields), ","))
		if addedFields:
			colortext.printf("Fields in %s but not in %s:" % (secondDB.infilepath, self.infilepath), "green")
			colortext.printf(join(sorted(addedFields), ","))
			
		colortext.printf("\nChanges from %s to %s:" % (self.infilepath, secondDB.infilepath), "green")
		colortext.printf("(Ignoring fields %s)" % join(sorted(list(fields_to_ignore.union(removedFields))), ", "), "green")
		IDs = list(self.set_of_available_keys - mykeys.difference(theirkeys))
		commonFieldsOfInterest = [f for f in self.fieldnames if (f not in removedFields) and (f not in fields_to_ignore)]
		foundDifference = False
		count = 0
		differing_fields = {}
		
		added_info = []
		changed_info = []
		deleted_info = []
		records_with_ddG_values = {}
		
		openedhere = [False, False]
		if not self.fhandle:
			openedhere[0] = True
			self.open()
		if not secondDB.fhandle:
			openedhere[1] = True
			secondDB.open()
			
		if IDs:
			if start_at_id:
				startID = start_at_id
			else:
				startID = IDs[0]
			if end_at_id == None:
				end_at_id = IDs[-1]
			for id in IDs:
				if start_at_id <= id <= end_at_id:
					myrecord = self.readRecord(id)
					theirrecord = secondDB.readRecord(id)
					for fieldname in commonFieldsOfInterest:
						if myrecord["ddG"] or theirrecord["ddG"]:
							records_with_ddG_values[id] = True
						if not fieldname in myrecord.keys():
							colortext.error("Record %s is missing field %s in %s." % (str(id), fieldname, self.infilepath))
						elif not fieldname in theirrecord.keys():
							colortext.error("Record %s is missing field %s in %s." % (str(id), fieldname, secondDB.infilepath))
						elif myrecord[fieldname] != theirrecord[fieldname]:
							differing_fields[fieldname] = True
							if theirrecord[fieldname] == None:
								deleted_info.append((id, fieldname, myrecord[fieldname]))
							elif myrecord[fieldname] == None:
								added_info.append((id, fieldname, theirrecord[fieldname]))
							else:
								changed_info.append((id, fieldname, myrecord[fieldname], theirrecord[fieldname]))
					count += 1
					if count % 1000 == 0:
						colortext.printf("[%d records read (#%d-#%d)]" % (count, startID, id), "yellow")
						startID = id
		if openedhere[0]:
			self.close()
		if openedhere[1]:
			secondDB.close()
					
		existingIDs = getIDsInDB(self.ddGDB, source = "ProTherm-2008-09-08-23581")
		for i in added_info:
			if showall:
				colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
				if i[0] not in existingIDs:
					colortext.printf("\tIgnored as record %d not in our database." % i[0], "cyan")
				if not records_with_ddG_values.get(i[0]):
					colortext.printf("\tIgnored as record %d has no ddG value." % i[0], "cyan")
				if not fields_of_interest.get(i[1]):
					colortext.printf("\tIgnored as field '%s' is not of interest." % i[1], "cyan")
			else:
				if i[0] in existingIDs and records_with_ddG_values.get(i[0]) and fields_of_interest.get(i[1]):
					colortext.message("Record %10.d: %s (%s). Value '%s' added." % (i[0], i[1], field_descriptions[i[1]], i[2]))			
		for i in deleted_info:
			if showall:
				colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
				if i[0] not in existingIDs:
					colortext.printf("\tIgnored as record %d not in our database." % i[0], "cyan")
				if not records_with_ddG_values.get(i[0]):
					colortext.printf("\tIgnored as record %d has no ddG value." % i[0], "cyan")
				if not fields_of_interest.get(i[1]):
					colortext.printf("\tIgnored as field '%s' is not of interest." % i[1], "cyan")
			else:
				if i[0] in existingIDs and records_with_ddG_values.get(i[0]) and fields_of_interest.get(i[1]):
					colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
		for i in changed_info:
			if showall:
				colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
				if i[0] not in existingIDs:
					colortext.printf("\tIgnored as record %d not in our database." % i[0], "cyan")
				if not records_with_ddG_values.get(i[0]):
					colortext.printf("\tIgnored as record %d has no ddG value." % i[0], "cyan")
				if not fields_of_interest.get(i[1]):
					colortext.printf("\tIgnored as field '%s' is not of interest." % i[1], "cyan")
			else:
				if i[0] in existingIDs and records_with_ddG_values.get(i[0]) and fields_of_interest.get(i[1]):
					colortext.message("Record %10.d: %s (%s)" % (i[0], i[1], field_descriptions[i[1]]))
					print("Original value: '%s'" % i[2])
					print("New value: '%s'" % i[3])
		if not (changed_info or added_info or deleted_info):
			print("None found.")
		else:
			print("Changed fields: %s" % join(sorted(differing_fields.keys()), ", "))
		
		#changed_records = set([i[0] for i in added_info]).union(set([i[0] for i in deleted_info])).union(set([i[0] for i in changed_info]))
		#records_in_database = getIDsInDB()
		#print(changed_records.intersection(records_in_database))

if __name__ == "__main__":
	args = sys.argv
	if len(args) > 1:
		if args[1].isdigit():
			ID = int(args[1])
			ptReader = ProThermReader(os.path.join("..", "rawdata", "ProTherm25616.dat"), quiet = True)
			record = ptReader.readRecord(ID)
			if not record:
				colortext.error("Could not read a record with number #%s" % str(ID))
			ptReader.printRecord(ID)
		else:
			if args[1] == "diff":
				flregex = re.compile("^(ProTherm)(\d+)[.]dat$", re.IGNORECASE)
				datapath = os.path.join("..", "rawdata")
				dbs = []
				for filenm in os.listdir(datapath):
					mtchs = re.match("^(ProTherm)(\d+)[.]dat$", filenm, re.IGNORECASE)
					if mtchs:
						#sqlfilename = "%s%s.sql" % (mtchs.group(1),mtchs.group(2))
						#sqlfilepath = os.path.join(datapath, sqlfilename)
						dbs.append((int(mtchs.group(2)), ProThermReader(os.path.join(datapath, filenm), mtchs.group(2))))
				dbs = sorted(dbs)
				dbs[0][1].diff(dbs[1][1], fields_to_ignore = ["E.C.NUMBER", "ION_NAME_1", "SWISSPROT_ID", "REMARKS", "NO_MOLECULE"], start_at_id = 0, end_at_id = None, showall = False)		
			if args[1] == "test":
				ptReader = ProThermReader(os.path.join("..", "rawdata", "ProTherm25616.dat"), quiet = True)
				ptReader.test()