#!/usr/bin/python2.4
# -*- coding: utf-8 -*-

import sys
import os
import re
import string
sys.path.insert(0, "..")
if __name__ == "__main__":
	sys.path.insert(0, "../common")
from string import join
import common.colortext as colortext
import ddgdbapi
import ddgobjects
from common.rosettahelper import kJtokcal, NUMBER_KJ_IN_KCAL, NUMBER_KELVIN_AT_ZERO_CELSIUS

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
	'1C6P' : {'MUTATED_CHAIN' : 'A'},
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
	'1M7T' : {'MUTATED_CHAIN' : 'A'},
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
	'1W4E' : {'MUTATED_CHAIN' : 'A'},
	'1W4H' : {'MUTATED_CHAIN' : 'A'},
	'1WIT' : {'MUTATED_CHAIN' : 'A'},
	'1WRP' : {'MUTATED_CHAIN' : 'R'},
	'1XAS' : {'MUTATED_CHAIN' : 'A'},
	'1YAL' : {'MUTATED_CHAIN' : 'A'},
	'1YCC' : {'MUTATED_CHAIN' : 'A'},
	'1YEA' : {'MUTATED_CHAIN' : 'A'},
	'1YMB' : {'MUTATED_CHAIN' : 'A'},
	'219L' : {'MUTATED_CHAIN' : 'A'},
	'2A36' : {'MUTATED_CHAIN' : 'A'},
	'2ABD' : {'MUTATED_CHAIN' : 'A'},
	'2ACE' : {'MUTATED_CHAIN' : 'A'},
	'2ACY' : {'MUTATED_CHAIN' : 'A'},
	'2ADA' : {'MUTATED_CHAIN' : 'A'},
	'2AIT' : {'MUTATED_CHAIN' : 'A'},
	'2AKY' : {'MUTATED_CHAIN' : 'A'},
	'2ASI' : {'MUTATED_CHAIN' : 'A'},
	'2B4Z' : {'MUTATED_CHAIN' : 'A'},
	'2BQA' : {'MUTATED_CHAIN' : 'A'},
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
	'2HID' : {'MUTATED_CHAIN' : 'A'},
	'2HMB' : {'MUTATED_CHAIN' : 'A'},
	'2HPR' : {'MUTATED_CHAIN' : 'A'},
	'2IFB' : {'MUTATED_CHAIN' : 'A'},
	'2IMM' : {'MUTATED_CHAIN' : 'A'},
	'2L3Y' : {'MUTATED_CHAIN' : 'A'},
	'2LZM' : {'MUTATED_CHAIN' : 'A'},
	'2NUL' : {'MUTATED_CHAIN' : 'A'},
	'2PDD' : {'MUTATED_CHAIN' : 'A'},
	'2PEC' : {'MUTATED_CHAIN' : 'A'},
	'2PGK' : {'MUTATED_CHAIN' : 'A'},
	'2PRD' : {'MUTATED_CHAIN' : 'A'},
	'2RN2' : {'MUTATED_CHAIN' : 'A'},
	'2TRT' : {'MUTATED_CHAIN' : 'A'},
	'2TS1' : {'MUTATED_CHAIN' : 'A'},
	'3BCI' : {'MUTATED_CHAIN' : 'A'},
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
	'1CYC' : {'MUTATED_CHAIN' : 'A'},
	'1FC1' : {'MUTATED_CHAIN' : 'A'},
	'1G6N' : {'MUTATED_CHAIN' : 'A'},
	'1HFY' : {'MUTATED_CHAIN' : 'A'},
	'1HFZ' : {'MUTATED_CHAIN' : 'A'},
	'1HTI' : {'MUTATED_CHAIN' : 'A'},
	'1IDS' : {'MUTATED_CHAIN' : 'A'},
	'1LRP' : {'MUTATED_CHAIN' : 'A'},
	'1MYL' : {'MUTATED_CHAIN' : 'A'},
	'1N0J' : {'MUTATED_CHAIN' : 'A'},
	'1RGG' : {'MUTATED_CHAIN' : 'A'},
	'1RN1' : {'MUTATED_CHAIN' : 'B'}, # Three identical chains A, B, C. This case is a special one - we choose chain B since while the PDB file contains identical chains, residue 45 is missing in chain A but required for records 10057 and 10058.
	'1SAK' : {'MUTATED_CHAIN' : 'A'},
	'1WQ5' : {'MUTATED_CHAIN' : 'A'},
	'1YNR' : {'MUTATED_CHAIN' : 'A'},
	'2AFG' : {'MUTATED_CHAIN' : 'A'},
	'2TRX' : {'MUTATED_CHAIN' : 'A'},
	'2ZTA' : {'MUTATED_CHAIN' : 'A'},
}

mutationsAllowedToBeStoredDespiteMissingCoordinates = set([
	965, 966, 2288, # 1WQ5 
	1577, 1578, 1579, 1580, 1603, 1604, 1605, 1606, 2012, 2013, 2054, 2055, 2066, 2073, 2074, 3166, 15947, 15948, # 1STN
	11869, # 1YCC
])


# These 726 records (333 of which are DDG records) use 2LZM as the PDB ID for T4 lysozyme. However, the mutations are performed against a cysteine-free pseudo-wildtype
# C54T-C97A (PDB ID: 1L63) described in "Control of enzyme activity by an engineered disulfide bond", Matsumura & Matthews, 1989, Science, 243, 792–794. PMID: 2916125
# http://www.jstor.org/sici?sici=0036-8075(19892)243:4892%3C792:%3E2.0.CO;2-#&origin=sfx%3Asfx
PseudoT4LysozymeCases = []
#PMID: 1854726. "Mutants T115E and Q123E and their controls were created on the wild-type template, and the rest of the mutants were created on the pseudo-wild-type template, which has two auxiliary mutations C54T and C97A"
#The measurements for T115E, T115E/K83M [and Q123E?] and the quoted values for wild-type lysozyme
#The values for N144E, N144E/K147M, K60H/L13D, K60H, K83H/A112D, K83H, S90H/Q122D, S90H and WT* were also obtained as a set under identical conditions
#
#Not all publications are included here. Those with DDG values are, those without DDG are not necessarily.
PseudoT4LysozymeCases.extend(range(368, 377 + 1))
PseudoT4LysozymeCases.extend(range(2211, 2222 + 1))
PseudoT4LysozymeCases.extend(range(13313, 13322 + 1))
PseudoT4LysozymeCases.extend(range(13973, 13984 + 1))
#PMID: 1553543 - "The mutant L133A was created by using the gene for wild-type lysozyme as a template. All of the other mutants were constructed with the gene for a pseudo-wild-type lysozyme, Cys54 -> Thr/Cys97 -* Ala (C54T/C97A, or WT*)"
PseudoT4LysozymeCases.extend(range(1021, 1024 + 1))
PseudoT4LysozymeCases.append(1026) 	#1025 is L133A
PseudoT4LysozymeCases.extend(range(13457, 13460 + 1))
PseudoT4LysozymeCases.append(13462) 	#13461 is L133A
#PMID: 1570293. "All experiments were carried out on a pseudo-wild-type lysozyme in which Cys-54 and Cys-97 in the native molecule were replaced with threonine and alanine, respectively"
PseudoT4LysozymeCases.extend(range(1057, 1064 + 1))
PseudoT4LysozymeCases.extend(range(1066, 1072 + 1))
PseudoT4LysozymeCases.extend(range(13489, 13503 + 1))
#PMID: 7507755. "All mutants are in the WT* background except where "(WT)" is explicitly indicated in the first column"
#[K16E, A41V, N116D, K147E are in the WT background]
PseudoT4LysozymeCases.extend(range(1074, 1101 + 1))
PseudoT4LysozymeCases.extend(range(13505, 13517 + 1))
#PMID: 8901549. "The 10 single-site mutants as well as various multiple-methionine mutants were constructed in cysteine-free pseudo wildtype lysozyme, hereafter identified as WT*"
PseudoT4LysozymeCases.extend(range(1102, 1112 + 1))
PseudoT4LysozymeCases.extend(range(13518, 13527 + 1))
#PMID:8218201 - "Pseudo-wild-type lysozyme, in which Cys residues at positions 54 and 97 have been replaced by Thr and Ala, respectively, was used as the reference protein and is referred to hereafter as WT*"
PseudoT4LysozymeCases.extend(range(1149, 1161 + 1))
PseudoT4LysozymeCases.extend(range(13558, 13569 + 1))
#PMID: 1911773 - "N116D, R119M, and N116D/R119M were produced on a wild-type template, while S38N, D92N, T109N, and T109D [and N144H] were produced in the "cysteine-free wild-type" (Le., C54T/C97A or WT*) background"
PseudoT4LysozymeCases.extend(range(2198, 2209 + 1))
PseudoT4LysozymeCases.extend(range(13963, 13972 + 1))
#PMID:11316887 - "All variants were constructed in the cysteine-free pseudo-wild-type lysozyme (WT*)" 
PseudoT4LysozymeCases.extend(range(11404, 11433 + 1))
#PMID:8676387 - "All variants were constructed in the cysteine-free wild-type lysozyme WT*"
PseudoT4LysozymeCases.extend(range(1113, 1120 + 1))
PseudoT4LysozymeCases.extend(range(13528, 13534 + 1))
#PMID:7869383 - "All mutant constructs made use of a pseudo wild-type T4 lysozyme gene (WT*) where Cys residues at positions 54 and 97 have been mutated to threonine and alanine, respectively"
PseudoT4LysozymeCases.extend(range(1121, 1136 + 1))
PseudoT4LysozymeCases.extend(range(13535, 13549 + 1))
#PMID:8401213 - "S117F was isolated from the background of the cysteine-free pseudo-wild-type T4 lysozyme C54T/C97A"
PseudoT4LysozymeCases.extend(range(1145, 1148 + 1))
PseudoT4LysozymeCases.extend(range(13556, 13557 + 1))
#PMID:1733941 - "The cysteine-free “pseudo wild-type” lysozyme (WT*),’ in which the 2 cysteine residues present in wild-type had been replaced in order to facilitate thermodynamic measurements (Wetzel et al., 1988; Matsumura and Matthews, 1989; Pjura  et al., 1990) was used as the reference protein."
PseudoT4LysozymeCases.extend(range(1171, 1183 + 1))
PseudoT4LysozymeCases.extend(range(13578, 13585 + 1))
#PMID:1569571 - "Mutagenesis was carried out by the uridine-labeled template method of Kunkel (1985) in the background of the pseudo wild-type cysteine-free C54T/C97A mutant of T4 lysozyme."
PseudoT4LysozymeCases.extend(range(1184, 1209 + 1))
PseudoT4LysozymeCases.extend(range(13586, 13609 + 1))
#PMID:1747370 - "Both M102K and L133D were constructed using the WT* template [C54T/C97A]" 
PseudoT4LysozymeCases.extend(range(1245, 1256 + 1))
PseudoT4LysozymeCases.extend(range(13629, 13634 + 1))
#PMID:7831309 - "the D20 substitutions, E11N, E11H, S117A/N1321, and S117A/N132M, which were made in the pseudo wild-type (WT*) background (6), and were compared to this protein"
PseudoT4LysozymeCases.extend(range(1458, 1466 + 1))
PseudoT4LysozymeCases.extend(range(13781, 13788 + 1))
#PMID:8289284 - "The site 44 mutants were constructed using a modified form of phage T4 lysozyme (referred to as WT*) in which Cys54 and Cys97 in the wild-type protein were replaced with Thr and Ala, respectively... at site 131, the additional replacements at this position were made in the wild-type rather than the WT* background"
PseudoT4LysozymeCases.extend(range(1467, 1486 + 1))
PseudoT4LysozymeCases.extend(range(13789, 13807 + 1))
#PMID:8433369 - "All mutations were made in a cysteine-free background [C54T/C97A]"
PseudoT4LysozymeCases.extend(range(1501, 1524 + 1))
PseudoT4LysozymeCases.extend(range(13817, 13838 + 1))
#PMID:1567817 - "All variants of T4 lysozyme in this study also contained the substitutions Cys54Thr and Cys97Ala... The change in melting temperature, DTm, of each mutant protein was determined relative to the Tm of the pseudo-wild-type reference protein, WT* ...DDG, was calculated from DDG =DTm.DS"
PseudoT4LysozymeCases.extend(range(1525, 1538 + 1))
PseudoT4LysozymeCases.extend(range(13839, 13850 + 1))
#PMID:8114100 - "All mutations were performed in a cysteine-free lysozyme background [C54T/C97A]"
PseudoT4LysozymeCases.extend(range(2580, 2584 + 1))
PseudoT4LysozymeCases.extend(range(14039, 14043 + 1))
#PMID:9514271 - "All mutants were constructed in the cysteine-free pseudo-wild-type lysozyme (WT*)". For I3A - "Mutant made in the wild-type protein"
PseudoT4LysozymeCases.extend(range(3290, 3312 + 1))
PseudoT4LysozymeCases.extend(range(14162, 14183 + 1))
#PMID:9541409 - "The method of Kunkel et al. (1987) was used to prepare mutants in cysteine-free (C54T, C97A) pseudo-wild-type (WT*) lysozyme... DDG is the change in free energy of unfolding of the mutant relative to WT*"
PseudoT4LysozymeCases.extend(range(3313, 3317 + 1))
PseudoT4LysozymeCases.extend(range(14184, 14188 + 1))
#PMID:2337607 - "Site-directed mutagenesis (Zoeller & Smith, 1983) was performed essentially as described by Kunkel (1985) on M13 single-strand DNA containing a derivative of the T4 lysozyme gene in which codons for cysteine residues 54 and 97 [note: a previous typo wrote 96] had been changed to encode threonine and alanine, respectively"
PseudoT4LysozymeCases.extend(range(3708, 3710 + 1))
PseudoT4LysozymeCases.extend(range(14245, 14246 + 1))
#PMID:7918421 - "A cysteine-free form of T4L, in which Cys 54 and Cys 97 have been replaced by Thr and Ala, respectively, was used for these experiments as the cysteines have been shown to complicate thermal denaturation experiments and yet to be unnecessary structurally and catalytically... Mutant proteins were constructed by Eckstein mutagenesis using an M13mp18 derivative containing a 650 base pair (bp) BamHI/HindIII fragment that encodes a cysteine-free T4L behind a twin tac promoter from plasmid pHSe54,97.TA". todo: ProTherm does not regard these DDG values as from the pseudo-wild-type. My changes should be checked.
PseudoT4LysozymeCases.extend(range(4795, 4797 + 1))
PseudoT4LysozymeCases.extend(range(14287, 14288 + 1))
#PMID:10623513 - No DDG values but "The mutants discussed here were constructed in the cysteine-free pseudo-wild-type lysozyme (WT*) (Matsumura & Matthews, 1989)"
PseudoT4LysozymeCases.extend(range(6333, 6363 + 1))
#PMID:10512706 - "The mutant lysozyme genes were constructed from the gene for pseudo-wild-type T4 lysozyme, C54T/C97A or WT*, in which the two cysteine residues were replaced with threonine and alanine. This cysteine-free lysozyme has structure and functional characteristics similar to wild-type, but displays better reversibility in thermal denaturation experiments"
PseudoT4LysozymeCases.extend(range(6578, 6589 + 1))
PseudoT4LysozymeCases.extend(range(14327, 14331 + 1))
PseudoT4LysozymeCases.extend(range(14333, 14337 + 1))
#PMID:8889173 - No DDG values but "All of the mutant proteins also had both of the two cysteines of the wild-type lysozyme mutated (C54T:C97A)"
PseudoT4LysozymeCases.extend(range(7431, 7446 + 1))
#PMID:12963380 - "The two redesigns of the C-terminal core of bacteriophage T4 lysozyme, Core-7 and Core-10, were made by iterative two-stage PCR37 using the gene for the cysteine-free (C54T/C97A) pseudo-wild-type (WTp) T4 lysozyme as the template... The double mutant V149/T152V was made in the WTp background... DDG is the change in the free energy of unfolding relative to WT*"
PseudoT4LysozymeCases.extend(range(16604, 16626 + 1))
#PMID:15340171 - "The mutants were constructed in a cysteine-free version of the T4 lysozyme gene designated as WT*"
PseudoT4LysozymeCases.extend(range(17222, 17234 + 1))
#PMID:12487988 - "Mutants were generated via the method of Kunkel et al. either in wildtype lysozyme (WT) or in the cysteine-free pseudo-wildtype protein (WT*) which includes the mutations C54T and C97A"
#Mutants with WT* as the reference protein:
#A40–49 					- N40A/K43A/S44A/E45A/L46A/D47A/K48A
#A40–49(K43L46)			- N40A/S44A/E45A/D47A/K48A
#A40–49(K43L46)/A127–132	- N40A/S44A/E45A/D47A/K48A/D127A/E128A/V131A/N132A  (17593, 17616)
#A34–49(L39K43L46) 		- T34A/K35A/S36A/P37A/S38A/N40A/S44A/E45A/D47A/K48A  	(17595)   
#A34–49(D38L39K43L46) 	- T34A/K35A/S36A/P37A/S38D/N40A/S44A/E45A/D47A/K48A  	(17597)     S38D
#A34–49(D35L39K43L46) 	- T34A/K35D/S36A/P37A/S38A/N40A/S44A/E45A/D47A/K48A  	(17598)     K35D
#A34–49(D36L39K43L46) 	- T34A/K35A/S36D/P38A/S38A/N40A/S44A/E45A/D47A/K48A* 	(17599)     S36D/P37A *ProTherm fixes P38A to the correct P37A
#A34–49(D37L39K43L46) 	- T34A/K35A/S36A/P37D/S38A/N40A/S44A/E45A/D47A/K48A 	(17600)     P37D
#A34–49(L39D40K43L46) 	- T34A/K35A/S36A/P37A/S38A/N40D/S44A/E45A/D47A/K48A  	(17601)     N40D
PseudoT4LysozymeCases.append(17593)
PseudoT4LysozymeCases.append(17616)
PseudoT4LysozymeCases.append(17595)
PseudoT4LysozymeCases.extend(range(17597, 17601 + 1))
#PMID:17400925 - "T4L* (* refers to a cysteine-free variant)... DDGeq is the difference in stability of the mutants and T4L* except where noted."
PseudoT4LysozymeCases.extend(range(22385, 22397 + 1))
#PMID:19384988 - "All of the mutant proteins were constructed in the WT background that includes cysteines at sites 54 and 97.
#todo: ProTherm does not take note of the use of WT* so this should be double-checked.
PseudoT4LysozymeCases.extend(range(25007, 25100 + 1))

# These 102 records (15 of which are DDG records) use 1LZ1 as the PDB ID for human lysozyme. However, the mutations are performed against a variant, 3SS, (PDB ID: 2BQA)
# which lacks one disulfide bond due to the double mutation C77A-C95A described in "Role of disulfide bonds in folding and secretion of human lysozyme in Saccharomyces cerevisiae", Taniyama, Y., Yamamoto, Y., Nakano, M., Kikuchi, M. & Ikehara, M., 1988, Biochem. Biophys. Res. Commun. 152, 962-967. PMID: 3288200
# http://www.sciencedirect.com/science/article/pii/S0006291X88803772
# "Substitution of Ala for Cys77 and Cys95 gave eight-fold greater secretion of a molecule with almost the same specific activity as that of the native enzyme". PMID: 3288200
# "the 3SS protein lacking a disulfide bond between Cys77 and Cys95 is destabilized by enthalpic factors". PMID: 9677301
# The human lysozyme wildtype is sometimes denoted by 4SS.
PseudoHumanLysozymeCases = []
#PMID: 9677301. "Five Ile to Val and nine Val to Ala mutants (3SS mutants) from 3SS (C77A/C95A) human lysozyme were constructed"
PseudoHumanLysozymeCases.extend(range(3384, 3398 + 1)) # 3SS thermodynamic parameters from Table 2
PseudoHumanLysozymeCases.extend(range(4247, 4315 + 1)) # 3SS thermodynamic parameters from Table 1
PseudoHumanLysozymeCases.extend(range(14189, 14202 + 1)) # 3SS DDG values from Table 3
#PMID: 10556244 mentions 3SS DDG values but the values in Table IV appear to be from 4SS mutants. Records 7046-7085 are from Table II, 7086-7102 are from Table III, 7103-7117 are from ?), 14376-14390 are from Table IV
#PMID: 12646687
PseudoHumanLysozymeCases.extend(range(16094, 16095 + 1)) # 3SS thermodynamic parameters from Table II
PseudoHumanLysozymeCases.extend(range(16096, 16097 + 1)) # 3SS thermodynamic parameters from Table III

# In these cases, the data in ProTherm is incorrect according to the publication
OverriddenEntries = []

# Missing PDB lengths.
for i in range(19104, 19151 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 71, 'PDB' : '1UZC'}))
for i in range(14973, 14999 + 1) + range(23725, 23758 + 1) + range(23764, 23791 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 238, 'PDB' : '1CHK'}))
for i in range(16592, 16600 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 435, 'PDB' : '1KFW'})) 
for i in range(16662, 16667 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 558, 'PDB' : '1W99'})) 
for i in range(16905, 16908 + 1) + range(22884, 22894 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 178,	'PDB' : '1BNL'})) 
# ** These cases have missing length for the wild-type structure, Onconase, Rana pipiens (P22069). It has been solved by X-Ray as 1ONC **
# PMID: 10913282.
for i in range(8561, 8600 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 104,	'PDB' : '1ONC'}))
# PMID: 16533040.
for i in range(19905, 19913 + 1):
	OverriddenEntries.append((i, {'LENGTH' : 104,	'PDB' : '1ONC'})) 
OverriddenEntries.append((19914, {'LENGTH' : 104,	'ASA' : None, 'PDB' : '1ONC'})) # Bad ASA record 

# ** These cases have ambiguous chain entries **
# 1WQ5. Two identical chains A, B. Only B has the information for residue 57.
# In these cases, the structural information needed for the mutations (residues A57, A62) is missing in the PDB file
# Some of the information is present in chain B (A57's corresponding residue) but I decided against remapping the residues. 
# PMIDs:2008436, 9020793
for i in [942, 963, 964, 2287, 13451]:
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'B', 'MUTATION' : 'P 1057 A', 'PDB' : '1WQ5'}))

# 1TUP. Three identical chains A, B, and C (and two DNA chains E and F) but '-' specified in ProTherm.
for i in range(3047, 3051 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1TUP'})) 

# 1G6N. Two identical chains A and B but '-' specified in ProTherm.
for i in range(3469, 3470 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1G6N'})) 

# 1AAR. Two identical chains A and B but '-' specified in ProTherm.
for i in range(5979, 5987 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'})) 
OverriddenEntries.append((2418, {'MUTATED_CHAIN' : 'A', 'PDB' : '1AAR'})) 
	
# 1FC1. Two identical chains A and B but '-' specified in ProTherm.
for i in range(3629, 3644 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1FC1'})) 
	
# 1LRP. Three identical chains A, B, and C but '-' specified in ProTherm. Todo: more to specify here.
for i in range(209, 212 + 1) + range(214, 215 + 1) + range(791, 793 + 1) + range(795, 801 + 1) + range(2170, 2190 + 1) + range(3604, 3611 + 1) + range(13412, 13420 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'})) 
for i in [13421, 13985, 13986, 14253, 14254, 14255]:
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1LRP'})) 

# 2AFG. Four identical chains A, B, C, and D but '-' specified in ProTherm
for i in range(8302, 8306 + 1) + range(14474, 14476 + 1) + range(24298, 24301 + 1): 
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '2AFG'})) 

# 1N0J. Two identical chains A and B but '-' specified in ProTherm
for i in range(14153, 14153 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A', 'PDB' : '1N0J'})) 
	
# 1OTR. Two distinct chains. There's no mutation here so there's no harm in specifying B as the chain.
for i in range(17394, 17396 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'B', 'PDB' : '1OTR'}))
	 
# PMID:15449934. Two distinct chains. All mutations in this publication are on the I chain (Eglin c).
for i in range(18378, 18385 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'I', 'PDB' : '1ACB'}))
	
# ** These cases have bad PDB IDs for the wild-type structure **
# PMID:12487987. This publication concerns MYL Arc repressor, not 1ARR. Only positions 31, 36, and 40 differ. 
for i in range(17567, 17583 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1MYL',	'PDB' : '1ARR'})) 

# These cases seem to have bad PDB IDs
# PMID: 12911302. Records 16576-16578.
# UniProt AC P83876 relates to 1PQN, 1QGV, and 1SYX. 
# 1PQN is the reduced form, 127-residue hDim1_{1-128}, deposited from this publication.
# 1QGV is the full-length, 142-residue hDim1 protein.
# R86A/K88A seems to be a mutation on the reduced form, 1PQN, rather than the full form, 1QGV
# Also, given "we report the solution structure for the reduced dominant negative form of Dim1 and compare it to the crystal structure of the oxidized full length Dim1 protein"
for i in range(16576, 16578 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1PQN', 'LENGTH' : 127, 'PDB' : '1QGV'})) 

# ** These cases have missing PDB IDs for the wild-type structure, IL6_MOUSE (P08505). It has been solved by NMR as 2L3Y where the residue IDs are off by 5 from the paper **
# PMID:9166791.
OverriddenEntries.extend([
	(2395  , {'PDB_wild'	: '2L3Y',	'PDB' : ''}),
	(2396  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'H 36 A', 'PDB' : ''}),
	(2397  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'W 39 A', 'PDB' : ''}),
	(2398  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'H 36 A, W 39 A', 'PDB' : ''}),
	(2399  , {'PDB_wild'	: '2L3Y',	'PDB' : ''}),
	(2400  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'W 39 A', 'PDB' : ''}),
	(2401  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'H 36 A, W 39 A', 'PDB' : ''}),
	(2402  , {'PDB_wild'	: '2L3Y',	'PDB' : ''}),
	(2403  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'H 36 A', 'PDB' : ''}),
	(2404  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'W 39 A', 'PDB' : ''}),
	(2405  , {'PDB_wild'	: '2L3Y',	'MUTATION' : 'H 36 A, W 39 A', 'PDB' : ''}),
])

# ** These cases have missing PDB IDs and length for the wild-type structure, Sso7d (synthetic). A PDB file 1BNZ has the same sequence and has been solved by X-Ray. **
# ** The Guerois set uses 1BF4. **
# PMID:11124040.
for i in range(10298, 10322 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1BNZ', 'LENGTH' : 64, 'MUTATED_CHAIN' : 'A', 'PDB' : ''}))

# PMID: 12144791. Records 15461 - 15500
# These cases have missing lengths, PDB IDs, and chains for the wild-type structures, HPr from Escherichia coli and Bacillus subtilis.
# "[in computer simulations] we used the pdb files 2HPR for bsHPr [Bacillus subtilis] and either 1POH or 1OPD for ecHPr [Escherichia coli]"
# 2HPR has G49 and length 87, 1POH and 1OPD have K49 and length 85. 1POH and 1OPD are not homologous - 1POH is the wildtype and has Q3 and S46 whereas 1OPD is a mutant with Q3E and S46D
# 2HPR is not the wildtype - it contains two mutations, M51V and S83C. 2HID with length 87 seems better, only having M51V. One chain of 3OQN has the wildtype sequence but is a six-chain (four unique chains) structure.
# M51V in 2HID "does not affect the function of HPr in vivo" [PMID:9336834]
# However, I chose to omit the records with bsHPr since 2HID is not the wildtype.
for i in range(15461, 15479 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1POH', 'LENGTH' : 85, 'PDB' : ''})) 
#for i in range(15480, 15486 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '2HID', 'LENGTH' : 87, 'PDB' : ''})) 
for i in range(15487, 15496 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1POH', 'LENGTH' : 85, 'PDB' : ''})) 
#for i in range(15497, 15500 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '2HID', 'LENGTH' : 87, 'PDB' : ''})) 
	
# PMID: 14756573. Records 16836-16851
# As with PMID:12144791 above, I chose to omit the records with bsHPr since I could not find a solved wildtype structure.
#for i in range(16836, 16851 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '2HID', 'LENGTH' : 87, 'PDB' : ''})) 
	
# PMID:12069590. Records 15409-15451
# These cases have missing lengths, PDB IDs, and chain for the wild-type structure, Rd-apocytochrome b562 (synthetic) ("Rd"=redesigned).
# This protein seems to match PDB file 1YYJ, the NMR solution structure of a redesigned apocytochrome b562:Rd-apocyt b562, which shares two authors with the publication.
# However, the mutations in Table 1 the reference suggest positions 16 and 41 are alanine whereas the PDB lists valine and glutamine respectively.
# However, the *text* of the reference states that valine is at position 16 and talks about the quintuple mutant (M7W/K98I/N99R/H102N/R106G) of Apocytochrome b562.
# Apocytochrome b562 has been solved by NMR as 1APC and 1YYJ has the same sequence as the quintuple mutant. All mutations in the publication match the wildtype residues bar A16 and A41 as mentioned above.
# However, 1APC has has an R at position 98, not a K as in the wildtype of the quintuple mutant.
# Apocytochrome b562 been solved by NMR as 1APC and 1YYJ has the same sequence as the quintuple mutant. All mutations in the publication match the wildtype residues bar A16 and A41 as mentioned above.
# I am using 1YYJ for all records except mutations from A16 and A41.
#for i in range(15409, 15413 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '1YYJ', 'LENGTH' : 106, 'MUTATED_CHAIN' : 'A', 'PDB' : ''})) 
	# A16 # 15414  : {'PDB_wild'	: '1YYJ', 'LENGTH' : 106, 'MUTATED_CHAIN' : 'A', 'PDB' : ''},
#for i in range(15415, 15428 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '1YYJ', 'LENGTH' : 106, 'MUTATED_CHAIN' : 'A', 'PDB' : ''})) 
	# A41 # 15429  : {'PDB_wild'	: '1YYJ', 'LENGTH' : 106, 'MUTATED_CHAIN' : 'A', 'PDB' : ''},
#for i in range(15430, 15451 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '1YYJ', 'LENGTH' : 106, 'MUTATED_CHAIN' : 'A', 'PDB' : ''})) 
	
# PMID: 15533036. Records 18311-18324. See above, PMID:12069590, reference 11 in the paper.
# "Protein expression and purification were carried out as described previously (11)" so 1YYJ may be correct again.
# 4GD7 is the quintuple mutant (W7D/L10G/L14G/V16G/I17G) described in PMID:12369818 which seems to be solved in 1YZC.
# The quintuple mutant (W7D/L10G/L14G/V16G/I17G) 1YZA of 1YYJ (itself a mutant of 1APC) may be the correct PDB ID here;
# the PDB entry for 1YZA links to a publication discussing 4GD7. 
#for i in range(18311, 18324 + 1):
#	OverriddenEntries.append((i, {'PDB_wild'	: '1YZA', 'LENGTH' : 106, 'MUTATED_CHAIN' : 'A', 'PDB' : ''})) 

# PMID:12135359. Records 15516-15518
# Thioredoxin (Human-Escherichia coli chimera), This was solved using NMR by the authors as 1M7T although the PDB ID was missing ('XXXX') in the publication.
for i in range(15516, 15518 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1M7T', 'LENGTH' : 107, 'PDB' : ''})) 

# PMID: 15769475.
# These cases have missing lengths, PDB IDs, and chain for the wild-type structure, E3BD*.
# E3BD* is a pseudo-wildtype of E3BD with the F166W mutation. The text gives the PDB ID of E3BD* as 1W4E.
# The PDB sequence matches the paper (and is from the same publication). Position 107 is valine in the PDB file and alanine in the publication sequence but in the main text is referred to as valine.
for i in range(18430, 18492 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1W4E', 'LENGTH' : 47, 'PDB' : ''})) 

# PMID: 15709759
# These cases have missing lengths, PDB IDs, and chain for the wild-type structure, APRin.
# The paper (2005) gives 1JIW (2001, X-ray) as the PDB ID for the APR-APRin complex. In 2008, the structure for APRin on its own was solved by NMR and published as 2RN4.
# I specify 1JIW here as it is solved by X-ray but we should mark 2RN4 as a "homolog" in the database.
for i in range(18674, 18686 + 1):
	OverriddenEntries.append((i, {'PDB_wild'	: '1JIW', 'LENGTH' : 106, 'PDB' : ''})) 

# ** These cases have missing length and PDB IDs for the wild-type structures, Bacillus stearothermophilus (BstHPr) and Bacillus subtilis (BsHPr). PDB IDs are given for the former in the reference. 
# PMID:15713472.
OverriddenEntries.append((18731  , {'PDB_wild'	: '1Y4Y', 'LENGTH' : 88, 'PDB' : ''}))
OverriddenEntries.append((18732  , {'PDB_wild'	: '1Y4Y', 'LENGTH' : 88, 'PDB_mutant' : '1Y51', 'PDB' : ''}))
#18733 - missing wildtype for Bacillus subtilis (BsHPr). 3OQN may be appropriate. See above.
OverriddenEntries.append((18734  , {'PDB_wild'	: '1Y4Y', 'LENGTH' : 88, 'PDB' : ''}))
OverriddenEntries.append((18735  , {'PDB_wild'	: '1Y4Y', 'LENGTH' : 88, 'PDB_mutant' : '1Y51', 'PDB' : ''}))
OverriddenEntries.append((18736  , {'PDB_wild'	: '1Y4Y', 'LENGTH' : 88, 'PDB_mutant' : '1Y51', 'PDB' : ''}))
#18737 - mBacillus subtilis (BsHPr)

# ** These cases have bad mutations **

#PMID:10623553
OverriddenEntries.extend([
	(6367  , {'MUTATION' : 'Y 68 F', 'PDB' : '1TTG'}), # Removing bad PDB residue ID corrections (meant for 1TEN)
])

#PMID:8504078
OverriddenEntries.extend([
	(2554  , {'MUTATION' : 'A 18 G', 'PDB' : '2WSY'}),
])

#PMID:9485396
OverriddenEntries.extend([
	(11869 , {'MUTATION' : 'P 76 L, K 72 M', 'PDB' : '1YCC'}), # Table 3 describes 3R22 as a double mutant (also noted in the footnote of Table 4)
])

# PMID:10956017
OverriddenEntries.extend([
	(14434 , {'MUTATION' : 'I 30 V', 'SEC.STR.' : 'Coil', 'PDB' : '1OTR'}), # Typo. I 30 V, I 36 L is the next mutation in the table. 
	(14438 , {'MUTATION' : 'I 30 F', 'SEC.STR.' : 'Coil', 'PDB' : '1OTR'}), # Typo. I 30 F, I 36 L is the next mutation in the table. 
])
	
# PMID:16503630
OverriddenEntries.extend([
	(19894 , {'MUTATION' : 'Q 19 E, Q 23 K, K 32 E, E 39 K, Q 60 K, S 65 K, E 69 K (PDB: Q28A E, Q 32A K, K 39A E, E 50A K, Q 71A K, S 76A K, E 80A K)', 'PDB' : '1AYE'}), # Typo: K 41A E given instead of K 39A E    
	(19897 , {'MUTATION' : 'Q 7 K, L 19 K, D 49 K, T 89 K (PDB: Q 808 K, L 820 K, D 850 K, T 890 K)', 'PDB' : '1TEN'}), # Missing PDB mapping.
	(19898 , {'MUTATION' : 'Q 7 K, L 19 K, D 49 K, T 89 K (PDB: Q 808 K, L 820 K, D 850 K, T 890 K)', 'PDB' : '1TEN'}), # Missing PDB mapping.
])
	
# PMID:19695265
OverriddenEntries.extend([
	(24290 , {'MUTATION' : 'V 31 I', 'PDB' : '2AFG'}), # Typo. Val31 is described as mutated to Ile in both Table 2 and throughout the text.
	(24296 , {'MUTATION' : 'K 12 V, C 83 T, C 117 V', 'PDB' : '2AFG'}), # Typo. L 12 V was entered instead of K 12 V - see Table 2.
])
	
# PMID:8142362. Entries 4489-4492 use the proper PDB numbering. Entries 14261-14264 do not.
OverriddenEntries.extend([
	(14261 , {'MUTATION' : 'R 53 E',			'PDB' : '1C5G'}),
	(14262 , {'MUTATION' : 'E 373 R',			'PDB' : '1C5G'}),
	(14263 , {'MUTATION' : 'E 373 P',			'PDB' : '1C5G'}),
	(14264 , {'MUTATION' : 'R 53 E, E 373 R',	'PDB' : '1C5G'}),
])
	
# ** These cases do not have normalized experimental conditions for parsing by my script. **
OverriddenEntries.extend([
	(11864 , {'dCp' 			 : None,'PDB' : '2FHA'}),
	(11865 , {'dCp' 			 : None,'PDB' : '2FHA'}),
	(11866 , {'dCp' 			 : None,'PDB' : '2FHA'}),
	(11867 , {'dCp' 			 : None,'PDB' : '2FHA'}),
	(24388 , {'dCp' 			 : None,'PDB' : ''}),
])
	
OverriddenEntries.extend([
	(16900  , {'m'			 : '6.86 kJ/mol/M', 'PDB' : '1C9O'}),
	
	(16093  , {'T'			 : '298 K', 	'PDB' : '2AFG'}), # Adding a space so the regex passes

	(889    , {'Tm'			 : '53.4 C',	'PDB' : '1ARR'}),
	(890    , {'Tm'			 : '67.3 C',	'PDB' : '1ARR'}),
	(5303   , {'Tm'			 : '<= 10.0 C',	'PDB' : '1YCC'}),
	(23589  , {'Tm'			 : '> 80 C', 	'PDB' : ''}),
	(25269  , {'Tm'			 : '52-54 C',	'PDB' : ''}),
	
	(17877  , {'dTm'			 : None,	'PDB' : ''}),
	(23676  , {'dTm'			 : None,	'PDB' : ''}),
	
	(14592  , {'ASA'			 : '114.3', 'PDB' : '1AM7'}),
])
	
# PMID: 15515183. Bad PDB ID and DDG calculation. The DDGs in the paper are given relative to alanine at position 33, not the wildtype lysine.
# 1UBQ also seems a better PDB ID for ubiquitin as it is an X-ray solution with just the ubiquitin chain.
OverriddenEntries.extend([
	## pH 2.25
	(19422 , {'ddG'	: "%s kJ/mol" % str( -8.1  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 A', 'PDB' : '1OTR'}),
	(19423 , {'ddG'	: "%s kJ/mol" % str( -8.8  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 E', 'PDB' : '1OTR'}),
	(19424 , {'ddG'	: "%s kJ/mol" % str( -8.7  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 F', 'PDB' : '1OTR'}),
	(19425 , {'ddG'	: "%s kJ/mol" % str(-14.2  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 G', 'PDB' : '1OTR'}),
	(19426 , {'ddG'	: "%s kJ/mol" % str( -3.5  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 I', 'PDB' : '1OTR'}),
	(19427 , {'ddG'	: "%s kJ/mol" % str( -9.7  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'wild',   'PDB' : '1OTR'}),
	(19428 , {'ddG'	: "%s kJ/mol" % str( -5.7  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 L', 'PDB' : '1OTR'}),
	(19429 , {'ddG'	: "%s kJ/mol" % str( -7.7  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 M', 'PDB' : '1OTR'}),
	(19430 , {'ddG'	: "%s kJ/mol" % str(-11.7  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 N', 'PDB' : '1OTR'}),
	(19431 , {'ddG'	: "%s kJ/mol" % str( -9.5  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 Q', 'PDB' : '1OTR'}),
	(19432 , {'ddG'	: "%s kJ/mol" % str(-10.8  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 S', 'PDB' : '1OTR'}),
	(19433 , {'ddG'	: "%s kJ/mol" % str( -9.9  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 T', 'PDB' : '1OTR'}),
	(19434 , {'ddG'	: "%s kJ/mol" % str( -6.0  -  -9.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 V', 'PDB' : '1OTR'}),
	# pH 2.5
	(19448 , {'ddG'	: "%s kJ/mol" % str( -6.1  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 A', 'PDB' : '1OTR'}),
	(19449 , {'ddG'	: "%s kJ/mol" % str( -6.9  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 E', 'PDB' : '1OTR'}),
	(19450 , {'ddG'	: "%s kJ/mol" % str( -6.6  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 F', 'PDB' : '1OTR'}),
	(19451 , {'ddG'	: "%s kJ/mol" % str(-12.8  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 G', 'PDB' : '1OTR'}),
	(19452 , {'ddG'	: "%s kJ/mol" % str( -1.5  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 I', 'PDB' : '1OTR'}),
	(19453 , {'ddG'	: "%s kJ/mol" % str( -7.7  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'wild',   'PDB' : '1OTR'}),
	(19454 , {'ddG'	: "%s kJ/mol" % str( -3.7  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 L', 'PDB' : '1OTR'}),
	(19455 , {'ddG'	: "%s kJ/mol" % str( -5.7  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 M', 'PDB' : '1OTR'}),
	(19456 , {'ddG'	: "%s kJ/mol" % str(-10.1  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 N', 'PDB' : '1OTR'}),
	(19457 , {'ddG'	: "%s kJ/mol" % str( -7.6  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 Q', 'PDB' : '1OTR'}),
	(19458 , {'ddG'	: "%s kJ/mol" % str( -9.0  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 S', 'PDB' : '1OTR'}),
	(19459 , {'ddG'	: "%s kJ/mol" % str( -8.3  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 T', 'PDB' : '1OTR'}),
	(19460 , {'ddG'	: "%s kJ/mol" % str( -4.0  -  -7.7), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 V', 'PDB' : '1OTR'}),
	# pH 2.75
	(19474 , {'ddG'	: "%s kJ/mol" % str( -3.0  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 A', 'PDB' : '1OTR'}),
	(19475 , {'ddG'	: "%s kJ/mol" % str( -3.6  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 E', 'PDB' : '1OTR'}),
	(19476 , {'ddG'	: "%s kJ/mol" % str( -3.8  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 F', 'PDB' : '1OTR'}),
	(19477 , {'ddG'	: "%s kJ/mol" % str(-10.3  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 G', 'PDB' : '1OTR'}),
	(19478 , {'ddG'	: "%s kJ/mol" % str(  1.7  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 I', 'PDB' : '1OTR'}),
	(19479 , {'ddG'	: "%s kJ/mol" % str( -4.4  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'wild',   'PDB' : '1OTR'}),
	(19480 , {'ddG'	: "%s kJ/mol" % str( -0.7  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 L', 'PDB' : '1OTR'}),
	(19481 , {'ddG'	: "%s kJ/mol" % str( -2.7  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 M', 'PDB' : '1OTR'}),
	(19482 , {'ddG'	: "%s kJ/mol" % str( -7.5  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 N', 'PDB' : '1OTR'}),
	(19483 , {'ddG'	: "%s kJ/mol" % str( -4.4  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 Q', 'PDB' : '1OTR'}),
	(19484 , {'ddG'	: "%s kJ/mol" % str( -6.4  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 S', 'PDB' : '1OTR'}),
	(19485 , {'ddG'	: "%s kJ/mol" % str( -5.2  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 T', 'PDB' : '1OTR'}),
	(19486 , {'ddG'	: "%s kJ/mol" % str( -0.6  -  -4.4), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 V', 'PDB' : '1OTR'}),
	# pH 3.0
	(19500 , {'ddG'	: "%s kJ/mol" % str(  1.0  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 A', 'PDB' : '1OTR'}),
	(19501 , {'ddG'	: "%s kJ/mol" % str(  0.3  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 E', 'PDB' : '1OTR'}),
	(19502 , {'ddG'	: "%s kJ/mol" % str(  0.7  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 F', 'PDB' : '1OTR'}),
	(19503 , {'ddG'	: "%s kJ/mol" % str( -6.6  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 G', 'PDB' : '1OTR'}),
	(19504 , {'ddG'	: "%s kJ/mol" % str(  6.0  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 I', 'PDB' : '1OTR'}),
	(19505 , {'ddG'	: "%s kJ/mol" % str( -0.1  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'wild',   'PDB' : '1OTR'}),
	(19506 , {'ddG'	: "%s kJ/mol" % str(  3.7  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 L', 'PDB' : '1OTR'}),
	(19507 , {'ddG'	: "%s kJ/mol" % str(  1.6  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 M', 'PDB' : '1OTR'}),
	(19508 , {'ddG'	: "%s kJ/mol" % str( -3.1  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 N', 'PDB' : '1OTR'}),
	(19509 , {'ddG'	: "%s kJ/mol" % str( -0.3  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 Q', 'PDB' : '1OTR'}),
	(19510 , {'ddG'	: "%s kJ/mol" % str( -2.0  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 S', 'PDB' : '1OTR'}),
	(19511 , {'ddG'	: "%s kJ/mol" % str( -0.8  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 T', 'PDB' : '1OTR'}),
	(19512 , {'ddG'	: "%s kJ/mol" % str(  3.6  -  -0.1), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 V', 'PDB' : '1OTR'}),
	# pH 3.25
	(19526 , {'ddG'	: "%s kJ/mol" % str(  5.7  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 A', 'PDB' : '1OTR'}),
	(19527 , {'ddG'	: "%s kJ/mol" % str(  4.1  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 E', 'PDB' : '1OTR'}),
	(19528 , {'ddG'	: "%s kJ/mol" % str(  4.6  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 F', 'PDB' : '1OTR'}),
	(19529 , {'ddG'	: "%s kJ/mol" % str( -2.2  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 G', 'PDB' : '1OTR'}),
	(19530 , {'ddG'	: "%s kJ/mol" % str( 10.4  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 I', 'PDB' : '1OTR'}),
	(19531 , {'ddG'	: "%s kJ/mol" % str(  3.9  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'wild',   'PDB' : '1OTR'}),
	(19532 , {'ddG'	: "%s kJ/mol" % str(  7.9  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 L', 'PDB' : '1OTR'}),
	(19533 , {'ddG'	: "%s kJ/mol" % str(  6.1  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 M', 'PDB' : '1OTR'}),
	(19534 , {'ddG'	: "%s kJ/mol" % str(  1.4  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 N', 'PDB' : '1OTR'}),
	(19535 , {'ddG'	: "%s kJ/mol" % str(  4.1  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 Q', 'PDB' : '1OTR'}),
	(19536 , {'ddG'	: "%s kJ/mol" % str(  2.4  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 S', 'PDB' : '1OTR'}),
	(19537 , {'ddG'	: "%s kJ/mol" % str(  3.4  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 T', 'PDB' : '1OTR'}),
	(19538 , {'ddG'	: "%s kJ/mol" % str(  7.8  -  3.9), 'PDB_wild' : '1UBQ', 'MUTATED_CHAIN' : 'A', 'MUTATION' : 'K 33 V', 'PDB' : '1OTR'}),
])

# PMID: 15515183. Bad PDB ID and DDG calculation. The DDGs in the paper are given relative to alanine at position 33, not the wildtype lysine.
for i in range(19409, 19421 + 1):
	OverriddenEntries.append((i, {'PDB_wild' : '1UBQ', 'PDB' : '1OTR'})) 
for i in range(19435, 19447 + 1):
	OverriddenEntries.append((i, {'PDB_wild' : '1UBQ', 'PDB' : '1OTR'})) 
for i in range(19461, 19473 + 1):
	OverriddenEntries.append((i, {'PDB_wild' : '1UBQ', 'PDB' : '1OTR'}))
for i in range(19487, 19499 + 1):
	OverriddenEntries.append((i, {'PDB_wild' : '1UBQ', 'PDB' : '1OTR'})) 
for i in range(19513, 19525 + 1):
	OverriddenEntries.append((i, {'PDB_wild' : '1UBQ', 'PDB' : '1OTR'}))

# PMID: 16922511. These cases have a missing PDB ID. 2A01 seems to be the correct PDB ID for lipid free Apolipoprotein A-I, homo sapiens.
for i in range(20165, 20181 + 1):
	if i != 20175 and i != 20181:
		OverriddenEntries.append((i, {'LENGTH' : 243,	'PDB_wild' : '2A01', 'MUTATED_CHAIN' : 'A', 'PDB' : ''})) 

# PMID: 19683006. These records are missing a chain ID. Chains A, B, C, or D are the p53 protein (tetramer complex). I chose chain A arbitrarily. 
for i in range(24275, 24282 + 1) + range(24284, 24286 + 1):
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A',	'PDB' : '2AC0'}))
OverriddenEntries.append((24283, {'MUTATED_CHAIN' : 'A', 'MUTATION' : 'M 133 L, C 141 V, Y 236 F, T 253 L', 'PDB' : '2AC0'})) # Standardizing record for parsing

# PMID: 19683006
for i in [24251, 24261, 24271]:
	OverriddenEntries.append((i, {'MUTATED_CHAIN' : 'A',	'PDB' : '2AC0'})) 

# PMID: 18077463. Records 24728-24733, 24736-24741 (DsbA, Staphylococcus aureus).
# The publication cites the PDB file for the Escherichia coli DsbA-DsbB-ubiquinone complex as 2HI7 (2006). This has an alanine at position 33 rather than the expected cysteine. ProTherm instead cites 1A23 (1998) which has DsbA on its own and C33 as expected.
# The publication cites 3BCI as the PDB file for Staphylococcus aureus DsbA and 3BD2 and 3BCK for its E96Q and T153V mutants respectively.
for i in range(24728, 24733 + 1):
	if i == 24730 or i == 24731:
		#SaDsbA T153V
		OverriddenEntries.append((i, {'PDB_wild' : '3BCI', 'LENGTH' : 186, 'PDB_mutant' : '3BCK', 'PDB' : ''}))
	elif i == 24732 or i == 24733:
		#SaDsbA E96Q
		OverriddenEntries.append((i, {'PDB_wild' : '3BCI', 'LENGTH' : 186, 'PDB_mutant' : '3BD2', 'PDB' : ''}))
	else:
		OverriddenEntries.append((i, {'PDB_wild' : '3BCI', 'LENGTH' : 186, 'PDB' : ''}))
for i in range(24736, 24741 + 1):
	if i == 24738 or i == 24739:
		#SaDsbA T153V
		OverriddenEntries.append((i, {'PDB_wild' : '3BCI', 'LENGTH' : 186, 'PDB_mutant' : '3BCK', 'PDB' : ''}))
	elif i == 24740 or i == 24741:
		#SaDsbA E96Q
		OverriddenEntries.append((i, {'PDB_wild' : '3BCI', 'LENGTH' : 186, 'PDB_mutant' : '3BD2', 'PDB' : ''}))
	else:
		OverriddenEntries.append((i, {'PDB_wild' : '3BCI', 'LENGTH' : 186, 'PDB' : ''}))

BadOrMissingMutants = [
	# ** These cases have bad PDB IDs for the mutant structure **
	# PMID: 1569557
	(90    , {'PDB_mutant'	 : '1BSD', 			'PDB' : '1BNI'}), 	# Data-entry error 
	# PMID: 8448112
	(303   , {'PDB_mutant'	 : '1SYC, 1SYD', 	'PDB' : '1STN'}),	# The two given structures are not homologous but record #2024 has the correct mutants, 1SYC and 1SYD
	# PMID: 7929430
	(13393 , {'PDB_mutant'	 : '1RBU', 			'PDB' : '2RN2'}), 	# 1RBU was given in records 15, 16, 13182, and 13183 and matches the correct sequence
	(13410 , {'PDB_mutant'	 : '1RBU', 			'PDB' : '2RN2'}), 	# 1RBU was given in records 15, 16, 13182, and 13183 and matches the correct sequence
	# PMID: 7869383
	(13535 , {'PDB_mutant'	 : None, 			'PDB' : '2LZM'}), 	# 166H is not a valid PDB ID. 166L has the T115A mutation but also has other mutations.
	# PMID: 1390669
	(1751  , {'PDB_mutant'	 : None, 			'PDB' : '4LYZ'}), 	# 1JKB is a human lysozyme mutant - 4LYZ is hen egg-white lysozyme
	(1755  , {'PDB_mutant'	 : None, 			'PDB' : '4LYZ'}), 	# 1JKB is a human lysozyme mutant - 4LYZ is hen egg-white lysozyme
	(1759  , {'PDB_mutant'	 : None, 			'PDB' : '4LYZ'}), 	# 1JKB is a human lysozyme mutant - 4LYZ is hen egg-white lysozyme
	# PMID:8784180
	(2410  , {'PDB_mutant'	 : None, 			'PDB' : '1POH'}), 	# 1POH->1OPD. 1OPD is a double mutant also containing Q3E.
	(2412  , {'PDB_mutant'	 : None, 			'PDB' : '1POH'}), 	# 1POH->1SPH. 1SPH is a bacillus subtilis HPr mutant. 1POH is an e.coli wildtype.
	(14019 , {'PDB_mutant'	 : None, 			'PDB' : '1POH'}), 	# 1POH->1OPD. 1OPD is a double mutant also containing Q3E.
	(14021 , {'PDB_mutant'	 : None, 			'PDB' : '1POH'}), 	# 1POH->1SPH. 1SPH is a bacillus subtilis HPr mutant. 1POH is an e.coli wildtype.
	
	# ** These cases have somewhat bad PDB IDs for the mutant structure **
	# PMID:7490748
	(1908  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPC. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	(1909  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPB. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	(1910  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPA. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	
	# PMID: 7937708
	# This patch seems erroneous - no mutant is specified (9822  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPC. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	# This patch seems erroneous - no mutant is specified (9823  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPB. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	# This patch seems erroneous - no mutant is specified (9824  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPA. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	
	# PMID: 7849029
	# This patch seems erroneous - no mutant is specified (14267 , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPC. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	# This patch seems erroneous - no mutant is specified(14266 , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1YPA. 2CI2 has a start tag SSVEKKPEGVNTGAGDRHN followed by an L20M mutation. The mutant also has a E78Q mutation (which may be okay - E78 in 2CI2 should be Q according to UniProt).
	(4698  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1COA. 1COA is a double mutant also containing L20M but this may be okay?
	(14280 , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1COA. 1COA is a double mutant also containing L20M but this may be okay?
	
	# PMID: 8218191
	(1845  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1COA. 1COA is a double mutant also containing L20M but this may be okay?
	(1859  , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1COA. 1COA is a double mutant also containing L20M but this may be okay?
	(13923 , {'PDB_mutant'	 : None, 			'PDB' : '2CI2'}), 	# 2CI2->1COA. 1COA is a double mutant also containing L20M but this may be okay?
	
	# Also specifying the chain ID of the mutant structure when the same chain ID from the wildtype does not exist
	# PMID: 2663837
	(2294  , {'PDB_mutant'	 : None, 	'_MUTANT_CHAIN' : 'A', 	'PDB' : '1RN1'}), 	# 1RN1->1LRA. 1LRA is a double mutant also containing Q25K.
	(13988 , {'PDB_mutant'	 : None, 	'_MUTANT_CHAIN' : 'A', 	'PDB' : '1RN1'}), 	# 1RN1->1LRA. 1LRA is a double mutant also containing Q25K.
	# PMID: 1980207
	(2570  , {'PDB_mutant'	 : None, 	'_MUTANT_CHAIN' : 'A', 	'PDB' : '1RN1'}), 	# 1RN1->1LRA. 1LRA is a double mutant also containing Q25K.
	(2572  , {'PDB_mutant'	 : None, 	'_MUTANT_CHAIN' : 'A', 	'PDB' : '1RN1'}), 	# 1RN1->1LRA. 1LRA is a double mutant also containing Q25K.
	
	# PMID: 11023787
	(9787  , {'PDB_mutant'	 : None, 			'PDB' : '1CEY'}), 	# 1CEY->1E6K. 1E6K also has A2S but this may be okay? 
	(9789  , {'PDB_mutant'	 : None, 			'PDB' : '1CEY'}), 	# 1CEY->1E6M. 1E6M also has A2S but this may be okay? 
	
	# PMID : 8535241
	(13895 , {'PDB_mutant'	 : None, 			'PDB' : '4LYZ'}), 	# 4LYZ->1LSN. 1LSN is a double mutant also containing D101N.
	
	# ** These cases have missing PDB IDs for the mutant structure **
	# PMID: 1569559
	(7385  , {'PDB_mutant'	 : '1BSA', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 69 
	(7392  , {'PDB_mutant'	 : '1BSB', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 79 
	(7393  , {'PDB_mutant'	 : '1BRI', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 80 
	(7394  , {'PDB_mutant'	 : '1BRI', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 80 
	(7397  , {'PDB_mutant'	 : '1BRJ', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 85 
	(7398  , {'PDB_mutant'	 : '1BSE', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 86 
	(7401  , {'PDB_mutant'	 : '1BAN', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 88 
	(7403  , {'PDB_mutant'	 : '1BRK', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 91 
	# PMID : 1569558
	(9982  , {'PDB_mutant'	 : '1BNS', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 58 
	(9991  , {'PDB_mutant'	 : '1BSA', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 69 
	(9994  , {'PDB_mutant'	 : '1BSB', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 79 
	(9996  , {'PDB_mutant'	 : '1BAO', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 82 
	(10003 , {'PDB_mutant'	 : '1BSC', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 84 
	(10004 , {'PDB_mutant'	 : '1BAN', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 88 
	(10010 , {'PDB_mutant'	 : '1BRH', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 49 
	# PMID: 14516751
	(16772 , {'PDB_mutant'	 : '1BAN', 			'PDB' : '1BNI'}), 	# The mutant is missing but given previously in record 88
	# PMID: 11841216
	(12817 , {'PDB_mutant'	 : '1BTI', 			'PDB' : '1BPI'}), 	# The mutant is missing but given previously in record 1728
	(12818 , {'PDB_mutant'	 : '1BPT', 			'PDB' : '1BPI'}), 	# The mutant is missing but given previously in record 1729
	(12820 , {'PDB_mutant'	 : '1FAN', 			'PDB' : '1BPI'}), 	# The mutant is missing but given previously in record 1732
	# PMID: 7540212
	(14073 , {'PDB_mutant'	 : '1BTI', 			'PDB' : '1BPI'}), 	# The mutant is missing but given previously in record 1728
	# PMID: 11254388
	(10875 , {'PDB_mutant'	 : '1HIB', 			'PDB' : '1IOB'}), 	# The mutant is missing but given previously in record 21
	# PMID: 1567879
	(2179  , {'PDB_mutant'	 : '1LLI', 			'PDB' : '1LRP'}), 	# The mutant is missing but given previously in record 2174 
	# PMID: 12646687
	(16094 , {'PDB_mutant'	 : '1IX0', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given in the publication (these mutations are of the pseudo-wildtype 2BQA) 	
	(16095 , {'PDB_mutant'	 : '1IX0', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given in the publication (these mutations are of the pseudo-wildtype 2BQA)	
	(16096 , {'PDB_mutant'	 : '1IX0', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given in the publication (these mutations are of the pseudo-wildtype 2BQA)
	(16097 , {'PDB_mutant'	 : '1IX0', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given in the publication (these mutations are of the pseudo-wildtype 2BQA)
	# PMID: 11121116
	(10838 , {'PDB_mutant'	 : '1OUA, 1LOZ', 	'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 809
	# PMID:	10556244
	(14376 , {'PDB_mutant'	 : '2HEC', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 3473 
	(14377 , {'PDB_mutant'	 : '1YAO', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 804 
	(14381 , {'PDB_mutant'	 : '1OUA, 1LOZ', 	'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 809 	
	(14382 , {'PDB_mutant'	 : '2HEE', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 3494 	
	(14383 , {'PDB_mutant'	 : '2HED', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 3478 	
	(14384 , {'PDB_mutant'	 : '1YAP', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 805 	
	
	# PMID:	9677301
	(4260  , {'PDB_mutant'	 : '1YAO', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 804 
	
	# PMID:	10556244
	(7088  , {'PDB_mutant'	 : '2HEC', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 3473 
	(7089  , {'PDB_mutant'	 : '1YAO', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 804 
	(7095  , {'PDB_mutant'	 : '2HED', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 3478 	
	(7094  , {'PDB_mutant'	 : '2HEE', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 3494 	
	(7096  , {'PDB_mutant'	 : '1YAP', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 805 	
	(7103  , {'PDB_mutant'	 : '2HEC', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 3473 
	(7104  , {'PDB_mutant'	 : '1YAO', 			'PDB' : '1LZ1'}), 	# The mutant is missing but given previously in record 804 
	(7109  , {'PDB_mutant'	 : '2HEE', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 3494 	
	(7110  , {'PDB_mutant'	 : '2HED', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 3478 	
	(7111  , {'PDB_mutant'	 : '1YAP', 			'PDB' : '1LZ1'}), 	# The mutants are missing but given previously in record 805 	
	
	# PMID: 1988046
	(4323  , {'PDB_mutant'	 : '1L19', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 13295
	(4332  , {'PDB_mutant'	 : '1L48', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1263
	(4343  , {'PDB_mutant'	 : '1L04, 1L05', 	'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1271
	(4346  , {'PDB_mutant'	 : '1L49', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1541
	
	# PMID:	1998663
	(4225  , {'PDB_mutant'	 : '1L71', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1239
	(4226  , {'PDB_mutant'	 : '1L70', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1240
	(4227  , {'PDB_mutant'	 : '1L36', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1241
	
	# PMID:	1998663
	(14258 , {'PDB_mutant'	 : '1L71', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1239
	(14259 , {'PDB_mutant'	 : '1L70', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1240
	(14260 , {'PDB_mutant'	 : '1L36', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1241
	
	# PMID:	10933506
	(8363  , {'PDB_mutant'	 : '1L44', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1030 
	(8364  , {'PDB_mutant'	 : '1L45', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1031 
	(8365  , {'PDB_mutant'	 : '1L46', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1032 
	
	# PMID:	11316887. Stored as 2LZM but actually 1L63 
	(11423 , {'PDB_mutant'	 : '237L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3303
	(11426 , {'PDB_mutant'	 : '126L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1161
	
	# PMID:	7918421
	(4796  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	
	# PMID:	8889173
	(7431  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	(7432  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	(7433  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	(7434  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	(7435  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	(7436  , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	
	# PMID:	7918421
	(14287 , {'PDB_mutant'	 : '1L68', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1059 
	
	# PMID:	12963380
	(16608 , {'PDB_mutant'	 : '1L77', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1185
	(16611 , {'PDB_mutant'	 : '235L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3302
	(16620 , {'PDB_mutant'	 : '1L77', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1185
	(16623 , {'PDB_mutant'	 : '235L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3302
	
	# PMID:	1553543
	(13460 , {'PDB_mutant'	 : '1L90', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1502
	
	# PMID:	1553543
	(1024  , {'PDB_mutant'	 : '1L90', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1502
	 
	# PMID:	10512706
	(6580  , {'PDB_mutant'	 : '1L90', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1502 
	(6586  , {'PDB_mutant'	 : '1L90', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1502 
	
	# PMID:	10512706
	(14328 , {'PDB_mutant'	 : '1L90', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1502 
	(14334 , {'PDB_mutant'	 : '1L90', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1502
	
	# PMID:	8901549
	(13521 , {'PDB_mutant'	 : '1L93', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1505
	
	# PMID:17400925
	(22386 , {'PDB_mutant'	 : '240L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3292
	(22388 , {'PDB_mutant'	 : '242L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3294
	(22390 , {'PDB_mutant'	 : '236L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3299
	(22391 , {'PDB_mutant'	 : '200L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 1023
	(22397 , {'PDB_mutant'	 : '240L', 			'PDB' : '2LZM'}), 	# The mutant is missing but given previously in record 3292
	
	# PMID: 10545167
	# Stored as 1L63
	(20184 , {'PDB_mutant'	 : '239L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3291 (wrongly as a 2LZM single mutant)
	(20186 , {'PDB_mutant'	 : '242L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3294 (wrongly as a 2LZM single mutant)
	(20189 , {'PDB_mutant'	 : '247L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3311 (wrongly as a 2LZM single mutant)
	(20190 , {'PDB_mutant'	 : '236L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3299 (wrongly as a 2LZM single mutant)
	(20192 , {'PDB_mutant'	 : '1L90', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1502 (wrongly as a 2LZM single mutant)
	(20194 , {'PDB_mutant'	 : '244L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3297 (wrongly as a 2LZM single mutant)
	(20196 , {'PDB_mutant'	 : '238L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3301 (wrongly as a 2LZM single mutant)
	(20197 , {'PDB_mutant'	 : '227L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3307 (wrongly as a 2LZM single mutant)
	(20199 , {'PDB_mutant'	 : '235L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3302 (wrongly as a 2LZM single mutant)
	(20201 , {'PDB_mutant'	 : '200L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1023 (wrongly as a 2LZM single mutant)
	(20202 , {'PDB_mutant'	 : '237L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3303 (wrongly as a 2LZM single mutant)
	(20203 , {'PDB_mutant'	 : '1L85', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1026 (wrongly as a 2LZM single mutant)
	(20216 , {'PDB_mutant'	 : '196L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1116 (wrongly as a 2LZM single mutant)
	(20237 , {'PDB_mutant'	 : '200L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1023 (wrongly as a 2LZM single mutant)
	
	# PMID:	10545167
	(20220 , {'PDB_mutant'	 : '239L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3291 (wrongly as a 2LZM single mutant)
	(20222 , {'PDB_mutant'	 : '242L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3294 (wrongly as a 2LZM single mutant)
	(20225 , {'PDB_mutant'	 : '247L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3311 (wrongly as a 2LZM single mutant)
	(20226 , {'PDB_mutant'	 : '236L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3299 (wrongly as a 2LZM single mutant)
	(20228 , {'PDB_mutant'	 : '1L90', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1502 (wrongly as a 2LZM single mutant)
	(20230 , {'PDB_mutant'	 : '244L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3297 (wrongly as a 2LZM single mutant)
	(20232 , {'PDB_mutant'	 : '238L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3301 (wrongly as a 2LZM single mutant)
	(20233 , {'PDB_mutant'	 : '227L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3307 (wrongly as a 2LZM single mutant)
	(20235 , {'PDB_mutant'	 : '235L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3302 (wrongly as a 2LZM single mutant)
	(20238 , {'PDB_mutant'	 : '237L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 3303 (wrongly as a 2LZM single mutant)
	(20239 , {'PDB_mutant'	 : '1L85', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1026 (wrongly as a 2LZM single mutant)
	(20250 , {'PDB_mutant'	 : '196L', 			'PDB' : '1L63'}), 	# The mutant is missing but given previously in record 1116 (wrongly as a 2LZM single mutant)
	
	# PMID:	12082163
	(15184 , {'PDB_mutant'	 : '1DVV', 			'PDB' : '451C'}), 	# The mutants are missing but given previously in record 9770 	
	(15193 , {'PDB_mutant'	 : '1DVV', 			'PDB' : '451C'}), 	# The mutants are missing but given previously in record 9770 	
	
	# PMID:	18957274
	(24929 , {'PDB_mutant'	 : '1DVV', 			'PDB' : '451C'}), 	# The mutants are missing but given previously in record 9770 	
	(24939 , {'PDB_mutant'	 : '1DVV', 			'PDB' : '451C'}), 	# The mutants are missing but given previously in record 9770 	

	# Bad mutants
	# PMID: 1988046. 237L and 238L are mutants of the pseudo-wildtype 1L63 and not 2LZM. The paper is mutating from 2LZM.
	(4333  , {'PDB_mutant'	 : None, 			'PDB' : '2LZM'}),
	(4340  , {'PDB_mutant'	 : None, 			'PDB' : '2LZM'}),
	# PMID:7831309 (1995). 1TLA is a mutant of the pseudo-wildtype. According to the paper, the S117F mutant was made in the wildtype background. In the reference for 1TLA (PMID:8401213, 1993), "S117F was isolated from the background of the cysteine-free pseudo-wild-type T4 lysozyme C54T/C97A"
	(13775  , {'PDB_mutant'	 : None, 			'PDB' : '2LZM'}),
	# PMID:9677301. 1OUG is a 4SS mutant. The wildtype here is the 3SS pseudo-wildtype.
	(3390   , {'PDB_mutant'	 : None, 			'PDB' : '1LZ1'}),
	(4274   , {'PDB_mutant'	 : None, 			'PDB' : '1LZ1'}),
	(4275   , {'PDB_mutant'	 : None, 			'PDB' : '1LZ1'}),
	(4276   , {'PDB_mutant'	 : None, 			'PDB' : '1LZ1'}),
	(14194  , {'PDB_mutant'	 : None, 			'PDB' : '1LZ1'}),
	# PMID:7918421. 1L68 is the mutant PDB ID for record 14287 (1L63 + S44A). This record has the mutation N68A.
	(4797   , {'PDB_mutant'	 : None, 			'PDB' : '2LZM'}),
	(14288  , {'PDB_mutant'	 : None, 			'PDB' : '2LZM'}),
	# PMID:11087397.
	# 1GB0 is the mutant for 1LZ1 + V2L as given in the publication. 1GA0 is a completely different protein.
	(9659   , {'PDB_mutant'	 : '1GB0', 			'PDB' : '1LZ1'}),
	(9660   , {'PDB_mutant'	 : '1GB0', 			'PDB' : '1LZ1'}),
	(9661   , {'PDB_mutant'	 : '1GB0', 			'PDB' : '1LZ1'}),
	(9719   , {'PDB_mutant'	 : '1GB0', 			'PDB' : '1LZ1'}),
	(9737   , {'PDB_mutant'	 : '1GB0', 			'PDB' : '1LZ1'}),
	(14499  , {'PDB_mutant'	 : '1GB0', 			'PDB' : '1LZ1'}),
	# PMID:11087397.
	# 1GB7 is the mutant for 1LZ1 + V74L as given in the publication. 1GA0 is a completely different protein.
	(9725   , {'PDB_mutant'	 : '1GB7', 			'PDB' : '1LZ1'}),
	(9742   , {'PDB_mutant'	 : '1GB7', 			'PDB' : '1LZ1'}),
	(14504  , {'PDB_mutant'	 : '1GB7', 			'PDB' : '1LZ1'}),
	# PMID:11087397.
	# 1GBX is the mutant for 1LZ1 + V110L as given in the publication. 1GA0 is a completely different protein.
	(9731   , {'PDB_mutant'	 : '1GBX', 			'PDB' : '1LZ1'}),
	(9747   , {'PDB_mutant'	 : '1GBX', 			'PDB' : '1LZ1'}),
	(14509  , {'PDB_mutant'	 : '1GBX', 			'PDB' : '1LZ1'}),
]

ddGTypos = [
	# PMID:9020793
	(970  , {'ddG_H2O'	:   "%s kcal/mol" % str(-0.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}), # Wrong sign
	
	# PMID:8378307
	(2232 , {'ddG'		:  '-3.46 kcal/mol',	'PDB' : '1VQB'}), # Incorrectly entered as -3.86
	
	# PMID:8819981
	(2508 , {'ddG'		:  '-0.06 kcal/mol',	'PDB' : '1HFY'}), # Incorrectly entered as 0.06. This is a one-off mistake for this publication so it is not included in the ddGWrongSigns dict.
	
	# PMID:8612074
	(2814 , {'ddG'		:   '0.2 kcal/mol',		'PDB' : '1BVC'}), # 2.9 - 2.7 (2.9 and 2.7 are the values of DG_NU - DG_IU for the mutant and wildtype respectively)
	
	# PMID:10079068. Wrong sign and rounding error.
	(5429 , {'ddG'		: "%s kcal/mol" % str( -0.7/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}), # This also has the wrong sign in ProTherm.
	
	# PMID:10600102. There is a mistake in the publication.
	(5982 , {'ddG'		:   '-6.2 kJ/mol',	'PDB' : '1AAR'}),
	(5983 , {'ddG'		:   '-7.0 kJ/mol',	'PDB' : '1AAR'}),
	
	# PMID:7663349. I do not know where ProTherm gets -1.4 kJ/mol (−0.3346 kcal/mol) from. It is close to the difference of the Cmid * m values i.e. 4.5*1.1 - 4.09*1.13.
	#7272 , {'ddG'		:   '%s kJ/mol' % ((1.115 * 0.41) * NUMBER_KJ_IN_KCAL),	'PDB' : '1POH'}),
	
	# PMID:10986129
	(8498  , {'ddG_H2O' 	: '-2.82 kcal/mol', 'PDB' : '1TEN'}), # Typo 
	
	# PMID:11513583. DG values entered as DDG values. Note: I am adding DDG values here. These should be double-checked.
	(11745 , {'ddG_H2O' :  '3.1 kcal/mol', 'dG_H2O' : '12.1', 'PDB' : '2TRX'}), # D26I is a stabilizing mutation
	(11746 , {'ddG_H2O' : '-3.7 kcal/mol', 'dG_H2O' : '5.3',  'PDB' : '2TRX'}), # Destabilizing mutation
	(11747 , {'ddG_H2O' : '-3.1 kcal/mol', 'dG_H2O' : '5.9',  'PDB' : '2TRX'}), # Destabilizing mutation
	(11748 , {'ddG_H2O' : '-1.4 kcal/mol', 'dG_H2O' : '7.6',  'PDB' : '2TRX'}), # Destabilizing mutation
	(11749 , {'ddG_H2O' : '-1.1 kcal/mol', 'dG_H2O' : '7.9',  'PDB' : '2TRX'}), # Destabilizing mutation
	
	# PMID:2372535
	(11893 , {'ddG'		: "%s kcal/mol" % str(( 8.6 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}), # Bad computation?
	
	# PMID:11964251 - ProTherm uses kJ/mol but does not specify them
	(13086 , {'ddG_H2O'	:   '-9.9 kJ/mol',	'PDB' : '5AZU'}),
	(13087 , {'ddG_H2O'	:   '-10.9 kJ/mol',	'PDB' : '5AZU'}),	

	# PMID:9228039
	(13204 , {'ddG'		:  '-0.46 kcal/mol', 'PDB' : '2RN2'}), # Wrong sign and rounding: "The mutant proteins A52Y, A52F, A52H, and A52K are unexpectedly unstable"
	
	# PMID: 9677301
	(14192 , {'ddG' 		: "%s kcal/mol" % str(-2.7/NUMBER_KJ_IN_KCAL),	'PDB' : '1LZ1'}), # Bad computation
	
	# PMID:12080133. DG values taken as DDG values
	(15280 , {'ddG_H2O' 		: "%s kJ/mol" % str(22.7 - 31.4),	'PDB' : '1TIT'}),
	(15281 , {'ddG_H2O' 		: "%s kJ/mol" % str(13.4 - 31.4),	'PDB' : '1TIT'}),
	
	# PMID:12215419. Also, ProTherm calls some values DDG and others DDG_H2O. Is this correct?
	(15688 , {'ddG_H2O' 	: "%s kJ/mol" % str(-7.2 + 3.5),	'PDB' : '1OTR'}), # Error in calculation.
	(15692 , {'dG' 		: "-16.4 kJ/mol", 'ddG'		: "%s kJ/mol" % str(-16.4 + 3.5),	'PDB' : '1OTR'}), # DG and DDG are swapped.
	
	# PMID:9819211
	(17873 , {'ddG_H2O'	: "-0.5 kcal/mol",	'PDB' : '1RN1'}),# Fixed copy-and-paste typo to allow parsing.
	
	# PMID:16953565
	(20134 , {'ddG'		: '-0.51 kcal/mol',			'PDB' : '1RTB'}), # minor typo (was -0.57)

	# PMID:19647749. dG_H2O value entered as ddG_H2O.
	(24386 , {'ddG_H2O' : None, 'dG_H2O' : '-17.4 kJ/mol',  'PDB' : ''}),
	
	# PMID:18077463.
	(24741, {'dG_H2O': '22.84 kJ/mol', 'PDB' : ''}), # typo as all other dG_H2O values are negated 
	
	# PMID:19565466.
	(25186 , {'ddG' : '-8.58 kJ/mol', 'PDB' : '1PIN'}), # Arithmetic error (ignored minus sign)
	(25198 , {'ddG' : '-4.00 kJ/mol', 'PDB' : '1PIN'}), # Arithmetic error
	
	# PMID:20340133
	(25658 , {'ddG_H2O'	: "-1.8 kcal/mol",	'PDB' : '1RGG'}),# Arithmetic error
	
	# PMID:20198681
	(25663 , {'ddG_H2O'	: "0.1 kcal/mol",	'PDB' : '1RGG'}),# Arithmetic error or a recalculation?
]	

ddGWrongSigns = {
	# PMID:8392867. The text explains that the stability is lowered by ~0.3 kcal/mol.
	2418 : 'ddG',
	
	# PMID:8504078. Only one case wrong out of the set.
	2566 : 'ddG',
	
	# PMID:2248951 - ProTherm appears to use the wrong sign
	3683 : 'ddG',
	3684 : 'ddG',
	3685 : 'ddG',
	
	# PMID:8652517 - ProTherm appears to use the wrong sign
	5038 : 'ddG_H2O', # N52I is a stabilizing mutation
	
	# PMID:9931001. "With the exception of L50V, all the single mutations reduce the stability of the folded state"
	5692 : 'ddG',
	5693 : 'ddG',
	5694 : 'ddG',
	5695 : 'ddG',
	5696 : 'ddG',
	5697 : 'ddG',
	5698 : 'ddG',
	5699 : 'ddG',
	5700 : 'ddG',
	5701 : 'ddG',
	5702 : 'ddG',
	5703 : 'ddG',
	5704 : 'ddG',

	#PMID:10623553. See PMID:10986129.
	6364 : 'ddG_H2O',
	6365 : 'ddG_H2O',
	6366 : 'ddG_H2O',
	6367 : 'ddG_H2O',
	6368 : 'ddG_H2O',
	
	#PMID:8917450. The Pro73->Val mutation destabilizes the folded conformation of RNase T1 by about 8.5 kJ/mol	
	8725 : 'ddG_H2O',
	8726 : 'ddG_H2O',
	8727 : 'ddG_H2O',
	8728 : 'ddG_H2O',

	# PMID:7937708 - S31A and S31G are destabilizing
	9820 : 'ddG',
	9821 : 'ddG',
	9822 : 'ddG',
	9823 : 'ddG',
	9824 : 'ddG',
	
	# PMID:7908135
	10564 : 'ddG_H2O',
	10565 : 'ddG_H2O',
	10566 : 'ddG_H2O',
	10567 : 'ddG_H2O',
	14565 : 'ddG',
	14566 : 'ddG',
	14567 : 'ddG',
	14568 : 'ddG',
	
	# PMID:8771183 - ProTherm uses the wrong sign
	11246 : 'ddG_H2O',
	11247 : 'ddG_H2O',
	11248 : 'ddG_H2O',

	# PMID:11695900 - ProTherm appears to use the wrong sign
	12235 : 'ddG_H2O', # Y14S is a destabilizing mutation
	12236 : 'ddG_H2O', # Destabilizing mutation
	12237 : 'ddG_H2O', # Destabilizing mutation

	# PMID:7918421 - ProTherm appears to use the wrong sign
	14287 : 'ddG',
	14288 : 'ddG',
	
	# PMID:12600203 - ProTherm appears to use the wrong sign			
	15807 : 'ddG_H2O',
	15808 : 'ddG_H2O',
	
	# PMID:15504413 - ProTherm appears to use the wrong sign for DDG values but the correct sign for DDG_H2O values
	17953 : 'ddG',
	17954 : 'ddG',
	17955 : 'ddG',
	17956 : 'ddG',
	17957 : 'ddG',
	17958 : 'ddG',
	17959 : 'ddG',
	17960 : 'ddG',
	17961 : 'ddG',
	17962 : 'ddG',
	17963 : 'ddG',
	17964 : 'ddG',
	17965 : 'ddG',
	17966 : 'ddG',
	17967 : 'ddG',
	17968 : 'ddG',
	17969 : 'ddG',
	17970 : 'ddG',
	17971 : 'ddG',
	17972 : 'ddG',
	17973 : 'ddG',
	17974 : 'ddG',
	17975 : 'ddG',
	17976 : 'ddG',
	17991 : 'ddG',
	17992 : 'ddG',
	17993 : 'ddG',
	17994 : 'ddG',
	17995 : 'ddG',
	17996 : 'ddG',
	17997 : 'ddG',
	17998 : 'ddG',
	17999 : 'ddG',
	18000 : 'ddG',
	18001 : 'ddG',
	18002 : 'ddG',
	18003 : 'ddG',
	18004 : 'ddG',
	18005 : 'ddG',
	18006 : 'ddG',
	18007 : 'ddG',
	18008 : 'ddG',
	
	# PMID:18189416 - ProTherm uses the wrong sign
	23154 : 'ddG_H2O',
	23155 : 'ddG_H2O',
	23156 : 'ddG_H2O',
	23157 : 'ddG_H2O',
	23158 : 'ddG_H2O',
	23159 : 'ddG_H2O',
	23160 : 'ddG_H2O',
	23161 : 'ddG_H2O',
	23162 : 'ddG_H2O',
	23163 : 'ddG_H2O',
	23164 : 'ddG_H2O',
}

RoundingErrors = [
	# PMID:1569557. Loss of precision on data entry.
	# 39 : No loss of precision.
	(40 , {'ddG'	: '-1.35 kcal/mol', 'PDB' : '1BNI'}),
	(41 , {'ddG'	: '-1.85 kcal/mol', 'PDB' : '1BNI'}),
	(42 , {'ddG'	: '-1.24 kcal/mol', 'PDB' : '1BNI'}),
	(43 , {'ddG'	: '-2.15 kcal/mol', 'PDB' : '1BNI'}),
	(44 , {'ddG'	: '-0.89 kcal/mol', 'PDB' : '1BNI'}),
	(45 , {'ddG'	: '-2.48 kcal/mol', 'PDB' : '1BNI'}),
	(46 , {'ddG'	: '-3.39 kcal/mol', 'PDB' : '1BNI'}),
	(47 , {'ddG'	: '-0.31 kcal/mol', 'PDB' : '1BNI'}),
	(48 , {'ddG'	: '-3.34 kcal/mol', 'PDB' : '1BNI'}),
	(49 , {'ddG'	: '-4.32 kcal/mol', 'PDB' : '1BNI'}),
	(50 , {'ddG'	: '-1.68 kcal/mol', 'PDB' : '1BNI'}),
	(51 , {'ddG'	:  '0.54 kcal/mol', 'PDB' : '1BNI'}),
	(52 , {'ddG'	: '-2.03 kcal/mol', 'PDB' : '1BNI'}),
	(53 , {'ddG'	: '-2.25 kcal/mol', 'PDB' : '1BNI'}),
	(54 , {'ddG'	: '-0.02 kcal/mol', 'PDB' : '1BNI'}),
	(55 , {'ddG'	: '-1.12 kcal/mol', 'PDB' : '1BNI'}),
	(56 , {'ddG'	: '-3.52 kcal/mol', 'PDB' : '1BNI'}),
	(57 , {'ddG'	: '-1.46 kcal/mol', 'PDB' : '1BNI'}),
	(58 , {'ddG'	: '-1.94 kcal/mol', 'PDB' : '1BNI'}),
	(59 , {'ddG'	: '-0.44 kcal/mol', 'PDB' : '1BNI'}),
	(60 , {'ddG'	: '-1.76 kcal/mol', 'PDB' : '1BNI'}),
	(61 , {'ddG'	: '-0.23 kcal/mol', 'PDB' : '1BNI'}),
	(62 , {'ddG'	:  '0.14 kcal/mol', 'PDB' : '1BNI'}),
	(63 , {'ddG'	: '-1.31 kcal/mol', 'PDB' : '1BNI'}),
	(64 , {'ddG'	: '-1.30 kcal/mol', 'PDB' : '1BNI'}),
	(65 , {'ddG'	: '-1.15 kcal/mol', 'PDB' : '1BNI'}),
	(66 , {'ddG'	: '-2.51 kcal/mol', 'PDB' : '1BNI'}),
	(67 , {'ddG'	: '-1.75 kcal/mol', 'PDB' : '1BNI'}),
	(68 , {'ddG'	: '-2.44 kcal/mol', 'PDB' : '1BNI'}),
	(69 , {'ddG'	: '-1.80 kcal/mol', 'PDB' : '1BNI'}),
	(70 , {'ddG'	: '-4.71 kcal/mol', 'PDB' : '1BNI'}),
	(71 , {'ddG'	: '-2.97 kcal/mol', 'PDB' : '1BNI'}),
	(72 , {'ddG'	: '-2.42 kcal/mol', 'PDB' : '1BNI'}),
	(73 , {'ddG'	: '-0.27 kcal/mol', 'PDB' : '1BNI'}),
	(74 , {'ddG'	: '-1.15 kcal/mol', 'PDB' : '1BNI'}),
	(75 , {'ddG'	: '-0.60 kcal/mol', 'PDB' : '1BNI'}),
	(76 , {'ddG'	: '-2.71 kcal/mol', 'PDB' : '1BNI'}),
	(77 , {'ddG'	:  '0.47 kcal/mol', 'PDB' : '1BNI'}),
	(78 , {'ddG'	: '-0.43 kcal/mol', 'PDB' : '1BNI'}),
	(79 , {'ddG'	: '-0.82 kcal/mol', 'PDB' : '1BNI'}),
	(80 , {'ddG'	: '-1.89 kcal/mol', 'PDB' : '1BNI'}),
	(81 , {'ddG'	: '-1.65 kcal/mol', 'PDB' : '1BNI'}),
	(82 , {'ddG'	: '-1.35 kcal/mol', 'PDB' : '1BNI'}),
	(83 , {'ddG'	: '-2.02 kcal/mol', 'PDB' : '1BNI'}),
	(84 , {'ddG'	: '-1.34 kcal/mol', 'PDB' : '1BNI'}),
	(85 , {'ddG'	: '-4.01 kcal/mol', 'PDB' : '1BNI'}),
	(86 , {'ddG'	: '-0.30 kcal/mol', 'PDB' : '1BNI'}),
	(87 , {'ddG'	: '-2.55 kcal/mol', 'PDB' : '1BNI'}),
	(88 , {'ddG'	: '-1.93 kcal/mol', 'PDB' : '1BNI'}),
	(89 , {'ddG'	: '-2.79 kcal/mol', 'PDB' : '1BNI'}),
	(90 , {'ddG'	: '-0.88 kcal/mol', 'PDB' : '1BNI'}),
	(91 , {'ddG'	: '-3.17 kcal/mol', 'PDB' : '1BNI'}),
	(92 , {'ddG'	: '-2.67 kcal/mol', 'PDB' : '1BNI'}),
	(93 , {'ddG'	:  '0.00 kcal/mol', 'PDB' : '1BNI'}),
	(94 , {'ddG'	: '-2.24 kcal/mol', 'PDB' : '1BNI'}),
	(95 , {'ddG'	: '-0.76 kcal/mol', 'PDB' : '1BNI'}),
	(96 , {'ddG'	: '-2.07 kcal/mol', 'PDB' : '1BNI'}),
	(97 , {'ddG'	: '-0.41 kcal/mol', 'PDB' : '1BNI'}),
	#2149 : No loss of precision.
	#2150 : No loss of precision.
	#2151 : No loss of precision.
	#2152 : No loss of precision.
	#2153 : No loss of precision.
	#2154 : No loss of precision.
	
	# PMID:1404369. Loss of precision on data entry.
	(171 , {'ddG' : '-0.14 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.52 * 1.91)), 'PDB'	:	'1BNI'}), #-0.1988
	(172 , {'ddG' : '-0.19 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.50 * 1.94)), 'PDB'	:	'1BNI'}), #-0.102
	(173 , {'ddG' : '-0.31 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.44 * 1.97)), 'PDB'	:	'1BNI'}), #-0.0852
	(174 , {'ddG' : '-0.35 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.41 * 2.02)), 'PDB'	:	'1BNI'}), #0.0762
	(175 , {'ddG' : '-0.41 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.38 * 1.98)), 'PDB'	:	'1BNI'}), #-0.1596
	(176 , {'ddG' : '-0.48 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.35 * 1.95)), 'PDB'	:	'1BNI'}), #-0.3495
	(177 , {'ddG' : '-0.55 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.31 * 2.08)), 'PDB'	:	'1BNI'}), #0.1328
	(178 , {'ddG' : '-0.66 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.25 * 2.00)), 'PDB'	:	'1BNI'}), #-0.332
	(179 , {'ddG' : '-0.69 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.24 * 1.98)), 'PDB'	:	'1BNI'}), #-0.4368
	(180 , {'ddG' : '-0.71 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.23 * 1.99)), 'PDB'	:	'1BNI'}), #-0.4143
	(181 , {'ddG' : '-0.78 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.19 * 1.94)), 'PDB'	:	'1BNI'}), #-0.7034
	(182 , {'ddG' : '-0.79 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.19 * 1.99)), 'PDB'	:	'1BNI'}), #-0.4939
	(183 , {'ddG' : '-0.81 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.18 * 1.95)), 'PDB'	:	'1BNI'}), #-0.681
	(184 , {'ddG' : '-0.82 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.17 * 1.95)), 'PDB'	:	'1BNI'}), #-0.7005
	(185 , {'ddG' : '-0.88 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.14 * 2.05)), 'PDB'	:	'1BNI'}), #-0.345
	(186 , {'ddG' : '-0.91 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.13 * 1.98)), 'PDB'	:	'1BNI'}), #-0.6546
	(187 , {'ddG' : '-0.98 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.09 * 1.88)), 'PDB'	:	'1BNI'}), #-1.1428
	(188 , {'ddG' : '-1.00 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (4.09 * 1.71)), 'PDB'	:	'1BNI'}), #-1.8381
	(189 , {'ddG' : '-4.08 kcal/mol', 'ddG_H2O' : "%s kcal/mol" % str(-8.832 + (2.49 * 1.89)), 'PDB'	:	'1BNI'}), #-4.1259
	
	# PMID:1870131. Loss of precision on data entry.
	(190 , {'ddG'	:  '0.39 kcal/mol', 'PDB' : '1BNI'}),
	
	# PMID:1317795. Loss of precision on conversion from kJ/mol to kcal/mol.
	(384 , {'ddG_H2O': "%s kcal/mol" % str(  1.59/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(385 , {'ddG_H2O': "%s kcal/mol" % str( -0.04/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(386 , {'ddG_H2O': "%s kcal/mol" % str(-11.00/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(387 , {'ddG_H2O': "%s kcal/mol" % str(  2.01/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(388 , {'ddG_H2O': "%s kcal/mol" % str( -2.43/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	
	# PMID:1317795. Loss of precision on conversion from kJ/mol to kcal/mol.
	(390	, {'ddG_H2O': "%s kcal/mol" % str( -9.2/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(391	, {'ddG_H2O': "%s kcal/mol" % str(  0.8/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(392	, {'ddG_H2O': "%s kcal/mol" % str(  0.4/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(393	, {'ddG_H2O': "%s kcal/mol" % str(  5.4/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(394	, {'ddG_H2O': "%s kcal/mol" % str( -4.2/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13323	, {'ddG': "%s kcal/mol" % str(-10.0/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13324	, {'ddG': "%s kcal/mol" % str( -1.9/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13325	, {'ddG': "%s kcal/mol" % str(-0.84/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13326	, {'ddG': "%s kcal/mol" % str(  2.4/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13327	, {'ddG': "%s kcal/mol" % str( -2.0/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	
	# PMID:8377205. Loss of precision on data entry.
	(510 , {'ddG'	:  '0.96 kcal/mol', 'PDB' : '1BNI'}),
	(511 , {'ddG'	:  '0.53 kcal/mol', 'PDB' : '1BNI'}),
	(512 , {'ddG'	: '-1.19 kcal/mol', 'PDB' : '1BNI'}),
	(513 , {'ddG'	:  '0.21 kcal/mol', 'PDB' : '1BNI'}),
	(514 , {'ddG'	: '-0.01 kcal/mol', 'PDB' : '1BNI'}),
	(515 , {'ddG'	: '-0.25 kcal/mol', 'PDB' : '1BNI'}),
	(516 , {'ddG'	:  '0.08 kcal/mol', 'PDB' : '1BNI'}),
	(517 , {'ddG'	: '-0.29 kcal/mol', 'PDB' : '1BNI'}),
	(518 , {'ddG'	: '-0.48 kcal/mol', 'PDB' : '1BNI'}),
	(519 , {'ddG'	:  '0.51 kcal/mol', 'PDB' : '1BNI'}),
	(520 , {'ddG'	:  '0.25 kcal/mol', 'PDB' : '1BNI'}),
	(521 , {'ddG'	:  '0.29 kcal/mol', 'PDB' : '1BNI'}),
	(522 , {'ddG'	: '-0.12 kcal/mol', 'PDB' : '1BNI'}),
	(523 , {'ddG'	: '-0.28 kcal/mol', 'PDB' : '1BNI'}),
	(524 , {'ddG'	: '-0.27 kcal/mol', 'PDB' : '1BNI'}),
	(525 , {'ddG'	: '-0.21 kcal/mol', 'PDB' : '1BNI'}),
	(526 , {'ddG'	:  '0.93 kcal/mol', 'PDB' : '1BNI'}),

	# PMID:8378307. Loss of precision on data entry.
	(544 , {'ddG'	:  '-0.68 kcal/mol', 'PDB' : '1VQB'}),
	(545 , {'ddG'	:  '-0.67 kcal/mol', 'PDB' : '1VQB'}),
	(546 , {'ddG'	:   '1.09 kcal/mol', 'PDB' : '1VQB'}),
	(547 , {'ddG'	:   '1.97 kcal/mol', 'PDB' : '1VQB'}),
	(548 , {'ddG'	:   '1.04 kcal/mol', 'PDB' : '1VQB'}),
	(549 , {'ddG'	:  '-3.49 kcal/mol', 'PDB' : '1VQB'}),
	(550 , {'ddG'	:  '-0.18 kcal/mol', 'PDB' : '1VQB'}),
	(551 , {'ddG'	:  '-2.25 kcal/mol', 'PDB' : '1VQB'}),
	(552 , {'ddG'	:  '-1.45 kcal/mol', 'PDB' : '1VQB'}),
	(553 , {'ddG'	:  '-3.21 kcal/mol', 'PDB' : '1VQB'}),
	(554 , {'ddG'	:  '-0.68 kcal/mol', 'PDB' : '1VQB'}),
	(555 , {'ddG'	:  '-2.72 kcal/mol', 'PDB' : '1VQB'}),
	(556 , {'ddG'	:  '-1.10 kcal/mol', 'PDB' : '1VQB'}),
	(557 , {'ddG'	:  '-0.62 kcal/mol', 'PDB' : '1VQB'}),
	(558 , {'ddG'	:  '-0.05 kcal/mol', 'PDB' : '1VQB'}),
	(559 , {'ddG'	:  '-5.30 kcal/mol', 'PDB' : '1VQB'}),
	(560 , {'ddG'	:  '-2.02 kcal/mol', 'PDB' : '1VQB'}),
	(561 , {'ddG'	:  '-0.67 kcal/mol', 'PDB' : '1VQB'}),
	(562 , {'ddG'	:  '-2.21 kcal/mol', 'PDB' : '1VQB'}),
	(563 , {'ddG'	:  '-2.62 kcal/mol', 'PDB' : '1VQB'}),
	(564 , {'ddG'	:   '0.50 kcal/mol', 'PDB' : '1VQB'}),
	(565 , {'ddG'	:  '-1.47 kcal/mol', 'PDB' : '1VQB'}),
	(566 , {'ddG'	:  '-4.30 kcal/mol', 'PDB' : '1VQB'}),
	(567 , {'ddG'	:   '0.76 kcal/mol', 'PDB' : '1VQB'}),
	(568 , {'ddG'	:   '1.63 kcal/mol', 'PDB' : '1VQB'}),
	(569 , {'ddG'	:   '1.23 kcal/mol', 'PDB' : '1VQB'}),
	(570 , {'ddG'	:  '-1.50 kcal/mol', 'PDB' : '1VQB'}),
	(571 , {'ddG'	:  '-0.66 kcal/mol', 'PDB' : '1VQB'}),
	(572 , {'ddG'	:   '0.47 kcal/mol', 'PDB' : '1VQB'}),
	
	# PMID:7764048. Loss of precision on conversion from kJ/mol to kcal/mol.
	(723 , {'ddG_H2O'	: "%s kcal/mol" % str((40.0 - 38.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(724 , {'ddG_H2O'	: "%s kcal/mol" % str((39.4 - 38.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(725 , {'ddG_H2O'	: "%s kcal/mol" % str((39.2 - 38.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(726 , {'ddG_H2O'	: "%s kcal/mol" % str((37.5 - 38.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(728 , {'ddG_H2O'	: "%s kcal/mol" % str((37.9 - 35.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(729 , {'ddG_H2O'	: "%s kcal/mol" % str((35.1 - 35.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(730 , {'ddG_H2O'	: "%s kcal/mol" % str((35.7 - 35.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(731 , {'ddG_H2O'	: "%s kcal/mol" % str((34.2 - 35.5)/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	
	# PMID:2000379. Removing 'precision' of -3.94 here since I do not see where the extra decimal point of precision in ProTherm comes from.
	(868 , {'ddG' : '-3.9 kcal/mol', 'PDB' : '1VQB'}), 
	
	# PMID:9020793. Loss of precision on conversion from kJ/mol to kcal/mol.
	(961  , {'ddG_H2O'	:   "%s kcal/mol" % str(-8.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(962  , {'ddG_H2O'	:   "%s kcal/mol" % str( 0.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(963  , {'ddG_H2O'	:   "%s kcal/mol" % str(-1.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(964  , {'ddG_H2O'	:   "%s kcal/mol" % str(-1.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(965  , {'ddG_H2O'	:   "%s kcal/mol" % str(-2.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(966  , {'ddG_H2O'	:   "%s kcal/mol" % str( 1.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(967  , {'ddG_H2O'	:   "%s kcal/mol" % str(-6.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(968  , {'ddG_H2O'	:   "%s kcal/mol" % str(-2.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(969  , {'ddG_H2O'	:   "%s kcal/mol" % str(-3.8/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	# 970 has the wrong sign and is fixed above
	(971  , {'ddG_H2O'	:   "%s kcal/mol" % str(-7.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(972  , {'ddG_H2O'	:   "%s kcal/mol" % str( 2.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(2286 , {'ddG_H2O'	:   "%s kcal/mol" % str(-8.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(2287 , {'ddG_H2O'	:   "%s kcal/mol" % str(-2.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(2288 , {'ddG_H2O'	:   "%s kcal/mol" % str(-1.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(2289 , {'ddG_H2O'	:   "%s kcal/mol" % str(-8.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(2290 , {'ddG_H2O'	:   "%s kcal/mol" % str(-4.4/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(2291 , {'ddG_H2O'	:   "%s kcal/mol" % str(-4.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	 
	# PMID:2217161. Loss of precision on conversion from kJ/mol to kcal/mol.
	(2513  , {'ddG'	:   "%s kcal/mol" % str((2.8 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2514  , {'ddG'	:   "%s kcal/mol" % str((4.3 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2515  , {'ddG'	:   "%s kcal/mol" % str((1.6 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2516  , {'ddG'	:   "%s kcal/mol" % str((6.2 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2517  , {'ddG'	:   "%s kcal/mol" % str((2.3 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2518  , {'ddG'	:   "%s kcal/mol" % str((3.4 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2519  , {'ddG'	:   "%s kcal/mol" % str((5.1 - 1.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	
	
	# PMID:3409879. Loss of precision on conversion from kJ/mol to kcal/mol.
	(2745  , {'ddG_H2O'	:   "%s kcal/mol" % str((15.0 - 20.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2746  , {'ddG_H2O'	:   "%s kcal/mol" % str((15.0 - 20.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	(2747  , {'ddG_H2O'	:   "%s kcal/mol" % str((22.0 - 20.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1IGV'}),
	
	# PMID:1765074. Loss of precision on conversion from kJ/mol to kcal/mol.
	(2749  , {'ddG_H2O'	: "%s kcal/mol" % str((12.5 - 16.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '3PGK'}),
	
	# PMID:9020874. Loss of precision on conversion from kJ/mol to kcal/mol.
	(2776  , {'ddG_H2O'	: "%s kcal/mol" % str((13.4 - 23.8)/NUMBER_KJ_IN_KCAL), 'PDB' : '1AXB'}),
	(2777  , {'ddG_H2O'	: "%s kcal/mol" % str((18.0 - 21.7)/NUMBER_KJ_IN_KCAL), 'PDB' : '1AXB'}),
	
	# PMID:9215576. Loss of precision on conversion from kJ/mol to kcal/mol.
	(3112 , {'ddG_H2O'	: "%s kcal/mol" % str( 4.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1CYO'}),
	(3113 , {'ddG_H2O'	: "%s kcal/mol" % str(-7.8/NUMBER_KJ_IN_KCAL), 'PDB' : '1CYO'}),
	(14154 , {'ddG'		: "%s kcal/mol" % str( 3.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1CYO'}),
	(14155 , {'ddG'		: "%s kcal/mol" % str(-7.8/NUMBER_KJ_IN_KCAL), 'PDB' : '1CYO'}),
	(14156 , {'ddG'		: "%s kcal/mol" % str(-11.8/NUMBER_KJ_IN_KCAL), 'PDB' : '1CYO'}),

	# PMID:9588945. Loss of precision on conversion from kJ/mol to kcal/mol.
	(3469 , {'ddG'		: "%s kcal/mol" % str(-1.1/NUMBER_KJ_IN_KCAL), 'PDB' : '1G6N'}),
	(3470 , {'ddG'		: "%s kcal/mol" % str( 0.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1G6N'}),
	
	# PMID:9533624. Loss of precision on conversion from kJ/mol to kcal/mol.
	(3519 , {'ddG_H2O'	: "%s kcal/mol" % str(-3.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(3520 , {'ddG_H2O'	: "%s kcal/mol" % str(-6.4/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(3521 , {'ddG_H2O'	: "%s kcal/mol" % str(-9.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(3522 , {'ddG_H2O'	: "%s kcal/mol" % str( 0.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(14241 , {'ddG'		: "%s kcal/mol" % str(-2.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(14242 , {'ddG'		: "%s kcal/mol" % str(-5.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(14243 , {'ddG'		: "%s kcal/mol" % str(-7.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	(14244 , {'ddG'		: "%s kcal/mol" % str( 1.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1CSP'}),
	
	# PMID:10079068. Loss of precision on conversion from kJ/mol to kcal/mol.
	(5424 , {'ddG'		: "%s kcal/mol" % str(  1.4/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5425 , {'ddG'		: "%s kcal/mol" % str(  7.2/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5426 , {'ddG'		: "%s kcal/mol" % str(  8.0/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5427 , {'ddG'		: "%s kcal/mol" % str(  2.1/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5428 , {'ddG'		: "%s kcal/mol" % str( 19.3/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	#5429 has the wrong sign as well
	(5430 , {'ddG'		: "%s kcal/mol" % str( 10.3/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5431 , {'ddG'		: "%s kcal/mol" % str(  0.6/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5432 , {'ddG'		: "%s kcal/mol" % str(  6.0/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5433 , {'ddG'		: "%s kcal/mol" % str(  1.1/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	(5434 , {'ddG'		: "%s kcal/mol" % str(  8.9/-NUMBER_KJ_IN_KCAL), 'PDB' : '1AG2'}),
	
	# PMID:10388847. Loss of precision on conversion from kJ/mol to kcal/mol.
	(5541 , {'ddG_H2O'		: "%s kcal/mol" % str((22.9 - 27.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1P2P'}),
	(5542 , {'ddG_H2O'		: "%s kcal/mol" % str((17.7 - 27.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1P2P'}),
	(5543 , {'ddG_H2O'		: "%s kcal/mol" % str((13.8 - 27.2)/NUMBER_KJ_IN_KCAL), 'PDB' : '1P2P'}),
	
	# PMID:7549876. Loss of precision on conversion from kcal/mol to kJ/mol.
	(7253 , {'ddG_H2O'	: "-4.3 kcal/mol", 'PDB' : '3MBP'}),
	(7254 , {'ddG_H2O'	: "-1.9 kcal/mol", 'PDB' : '3MBP'}),
	
	# PMID:8043610. Loss of precision on conversion from kcal/mol to kJ/mol.
	(7257 , {'ddG_H2O'	: "-6.9 kcal/mol",  'PDB' : '1B0O'}),
	
	# PMID:8795042. Loss of precision on conversion on conversion from kJ/mol to kcal/mol.
	(11772 , {'ddG'	: "%s kcal/mol" % str(-5.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11773 , {'ddG'	: "%s kcal/mol" % str(-6.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11774 , {'ddG'	: "%s kcal/mol" % str(-2.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11775 , {'ddG'	: "%s kcal/mol" % str(-5.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11776 , {'ddG'	: "%s kcal/mol" % str(-9.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11777 , {'ddG'	: "%s kcal/mol" % str(-5.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11778 , {'ddG'	: "%s kcal/mol" % str(-5.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11779 , {'ddG'	: "%s kcal/mol" % str(-7.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11780 , {'ddG'	: "%s kcal/mol" % str(-4.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	(11781 , {'ddG'	: "%s kcal/mol" % str(-4.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1WQ5'}),
	
	# PMID:8539253. Loss of precision on conversion on conversion from kJ/mol to kcal/mol.
	(11860 , {'ddG'	: "%s kcal/mol" % str((61.1 - 71.7)/NUMBER_KJ_IN_KCAL), 'PDB' : '1ROP'}),
	(11861 , {'ddG'	: "%s kcal/mol" % str((46.1 - 71.7)/NUMBER_KJ_IN_KCAL), 'PDB' : '1ROP'}),
	
	# PMID:9336842. Loss of precision on conversion from kJ/mol to kcal/mol.
	(11863 , {'ddG'	: "%s kcal/mol" % str( 2.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1BTA'}),
	
	# PMID:2372535. Loss of precision on conversion on conversion from kJ/mol to kcal/mol.
	(11887 , {'ddG'	: "%s kcal/mol" % str((-1.5 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}),
	(11888 , {'ddG'	: "%s kcal/mol" % str(( 0.9 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}),
	(11889 , {'ddG'	: "%s kcal/mol" % str(( 6.0 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}),
	(11890 , {'ddG'	: "%s kcal/mol" % str((-0.47 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}),
	(11891 , {'ddG'	: "%s kcal/mol" % str(( 6.0 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}),
	(11892 , {'ddG'	: "%s kcal/mol" % str(( 5.1 - 6.0)/NUMBER_KJ_IN_KCAL), 'PDB' : '1STN'}),
	#11893 Seems to have been entered incorrectly so is corrected above
	
	# PMID:9228039. Loss of precision on data entry.
	(13199 , {'ddG'	:   '1.88 kcal/mol', 'PDB' : '2RN2'}),
	(13200 , {'ddG'	:   '1.67 kcal/mol', 'PDB' : '2RN2'}),
	(13201 , {'ddG'	:   '1.31 kcal/mol', 'PDB' : '2RN2'}),
	(13202 , {'ddG'	:   '0.76 kcal/mol', 'PDB' : '2RN2'}),
	(13203 , {'ddG'	:   '0.49 kcal/mol', 'PDB' : '2RN2'}),
	#13204 also has the wrong sign
	(13205 , {'ddG'	:  '-0.82 kcal/mol', 'PDB' : '2RN2'}),
	(13206 , {'ddG'	:  '-1.19 kcal/mol', 'PDB' : '2RN2'}),
	(13207 , {'ddG'	:  '-1.52 kcal/mol', 'PDB' : '2RN2'}),
	(13208 , {'ddG'	:  '-1.64 kcal/mol', 'PDB' : '2RN2'}),
	(13209 , {'ddG'	:  '-1.76 kcal/mol', 'PDB' : '2RN2'}),
	(13210 , {'ddG'	:  '-1.79 kcal/mol', 'PDB' : '2RN2'}),
	(13211 , {'ddG'	:  '-1.85 kcal/mol', 'PDB' : '2RN2'}),
	(13212 , {'ddG'	:  '-2.31 kcal/mol', 'PDB' : '2RN2'}),
	(13213 , {'ddG'	:  '-2.71 kcal/mol', 'PDB' : '2RN2'}),
	(13214 , {'ddG'	:  '-3.59 kcal/mol', 'PDB' : '2RN2'}),
	(13215 , {'ddG'	:  '-5.93 kcal/mol', 'PDB' : '2RN2'}),
	
	# PMID:8955106. Loss of precision on data entry.
	(13216 , {'ddG'	:  '1.96 kcal/mol', 'PDB' : '2RN2'}),
	(13217 , {'ddG'	:  '4.05 kcal/mol', 'PDB' : '2RN2'}),
	(13218 , {'ddG'	:  '0.96 kcal/mol', 'PDB' : '2RN2'}),
	(13219 , {'ddG'	:  '2.67 kcal/mol', 'PDB' : '2RN2'}),
	(13220 , {'ddG'	:  '2.37 kcal/mol', 'PDB' : '2RN2'}),
	(13221 , {'ddG'	:  '0.28 kcal/mol', 'PDB' : '2RN2'}),
	(13222 , {'ddG'	: '-0.28 kcal/mol', 'PDB' : '2RN2'}),
	(13223 , {'ddG'	: '-0.22 kcal/mol', 'PDB' : '2RN2'}),
	(13224 , {'ddG'	:  '1.57 kcal/mol', 'PDB' : '2RN2'}),
	(13225 , {'ddG'	:  '1.08 kcal/mol', 'PDB' : '2RN2'}),
	(13226 , {'ddG'	:  '0.11 kcal/mol', 'PDB' : '2RN2'}),
	(13227 , {'ddG'	:  '1.84 kcal/mol', 'PDB' : '2RN2'}),
	(13228 , {'ddG'	:  '1.96 kcal/mol', 'PDB' : '2RN2'}),
	(13229 , {'ddG'	: '-0.74 kcal/mol', 'PDB' : '2RN2'}),
	(13230 , {'ddG'	:  '2.42 kcal/mol', 'PDB' : '2RN2'}),
	(13231 , {'ddG'	:  '1.12 kcal/mol', 'PDB' : '2RN2'}),
	(13232 , {'ddG'	:  '0.70 kcal/mol', 'PDB' : '2RN2'}),
	(13233 , {'ddG'	:  '0.44 kcal/mol', 'PDB' : '2RN2'}),
	(13234 , {'ddG'	: '-0.06 kcal/mol', 'PDB' : '2RN2'}),
	(13235 , {'ddG'	: '-0.06 kcal/mol', 'PDB' : '2RN2'}),
	(13236 , {'ddG'	:  '0.23 kcal/mol', 'PDB' : '2RN2'}),
	(13237 , {'ddG'	: '-0.35 kcal/mol', 'PDB' : '2RN2'}),
	(13238 , {'ddG'	: '-0.09 kcal/mol', 'PDB' : '2RN2'}),
	(13239 , {'ddG'	:  '0.53 kcal/mol', 'PDB' : '2RN2'}),
	(13240 , {'ddG'	: '-0.09 kcal/mol', 'PDB' : '2RN2'}),
	(13241 , {'ddG'	:  '0.85 kcal/mol', 'PDB' : '2RN2'}),
	
	# PMID:8251481. Loss of precision on data entry.
	(13271 , {'ddG'	:  '-0.55 kcal/mol', 'PDB' : '1BVC'}),
	(13272 , {'ddG'	:  '-0.56 kcal/mol', 'PDB' : '1BVC'}),
	(13273 , {'ddG'	:   '0.04 kcal/mol', 'PDB' : '1BVC'}),
	(13274 , {'ddG'	:  '-1.33 kcal/mol', 'PDB' : '1BVC'}),
	(13275 , {'ddG'	:  '-0.64 kcal/mol', 'PDB' : '1BVC'}),
	(13276 , {'ddG'	:  '-1.14 kcal/mol', 'PDB' : '1BVC'}),
	(13277 , {'ddG'	:  '-1.84 kcal/mol', 'PDB' : '1BVC'}),
	(13278 , {'ddG'	:   '0.63 kcal/mol', 'PDB' : '1BVC'}),
	(13279 , {'ddG'	:   '0.93 kcal/mol', 'PDB' : '1BVC'}),
	(13280 , {'ddG'	:  '-0.12 kcal/mol', 'PDB' : '1BVC'}),
	(13281 , {'ddG'	:  '-1.10 kcal/mol', 'PDB' : '1BVC'}),
	(13282 , {'ddG'	:  '-1.12 kcal/mol', 'PDB' : '1BVC'}),
	(13283 , {'ddG'	:   '0.12 kcal/mol', 'PDB' : '1BVC'}),
	(13284 , {'ddG'	:  '-1.72 kcal/mol', 'PDB' : '1BVC'}),
	(13285 , {'ddG'	:  '-2.37 kcal/mol', 'PDB' : '1BVC'}),
	(13286 , {'ddG'	:  '-0.02 kcal/mol', 'PDB' : '1BVC'}),
	(13287 , {'ddG'	:   '0.00 kcal/mol', 'PDB' : '1BVC'}),
	(13288 , {'ddG'	:  '-0.10 kcal/mol', 'PDB' : '1BVC'}),
	(13289 , {'ddG'	:  '-1.18 kcal/mol', 'PDB' : '1BVC'}),
	(13290 , {'ddG'	:  '-1.54 kcal/mol', 'PDB' : '1BVC'}),
	(13291 , {'ddG'	:  '-0.79 kcal/mol', 'PDB' : '1BVC'}),
	(13292 , {'ddG'	:  '-2.25 kcal/mol', 'PDB' : '1BVC'}),
	(13293 , {'ddG'	:  '-0.80 kcal/mol', 'PDB' : '1BVC'}),
	
	# PMID:1854726. Loss of precision on data entry of some records.
	(13308 , {'ddG'	:   '0.04 kcal/mol', 'PDB' : '2LZM'}),
	
	# PMID:8358293. Loss of precision on data entry of some records.
	(13328 , {'ddG'	: '1.05 kcal/mol', 'PDB' : '1BVC'}),
	(13329 , {'ddG'	: '0.75 kcal/mol', 'PDB' : '1BVC'}),
	(13330 , {'ddG'	: '0.59 kcal/mol', 'PDB' : '1BVC'}),
	(13331 , {'ddG'	: '0.16 kcal/mol', 'PDB' : '1BVC'}),
	(13332 , {'ddG'	: '-0.67 kcal/mol', 'PDB' : '1BVC'}),
	(13333 , {'ddG'	: '-0.26 kcal/mol', 'PDB' : '1BVC'}),
	(13334 , {'ddG'	: '-0.26 kcal/mol', 'PDB' : '1BVC'}),
	(13335 , {'ddG'	: '-0.44 kcal/mol', 'PDB' : '1BVC'}),
	(13336 , {'ddG'	: '-0.41 kcal/mol', 'PDB' : '1BVC'}),
	(13337 , {'ddG'	: '-1.12 kcal/mol', 'PDB' : '1BVC'}),
	(13338 , {'ddG'	: '-1.78 kcal/mol', 'PDB' : '1BVC'}),
	(13339 , {'ddG'	: '-1.45 kcal/mol', 'PDB' : '1BVC'}),
	(13340 , {'ddG'	: '-1.41 kcal/mol', 'PDB' : '1BVC'}),
	(13341 , {'ddG'	: '-1.60 kcal/mol', 'PDB' : '1BVC'}),
	(13342 , {'ddG'	: '-1.92 kcal/mol', 'PDB' : '1BVC'}),
	
	# PMID:8125123. Loss of precision on conversion from kJ/mol to kcal/mol.
	(13351 , {'ddG'	: "%s kcal/mol" % str(3.7/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13352 , {'ddG'	: "%s kcal/mol" % str(8.1/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13353 , {'ddG'	: "%s kcal/mol" % str(3.6/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13354 , {'ddG'	: "%s kcal/mol" % str(5.5/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13355 , {'ddG'	: "%s kcal/mol" % str(4.5/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13356 , {'ddG'	: "%s kcal/mol" % str(4.5/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13357 , {'ddG'	: "%s kcal/mol" % str(4.7/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13358 , {'ddG'	: "%s kcal/mol" % str(5.3/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13359 , {'ddG'	: "%s kcal/mol" % str(6.3/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13360 , {'ddG'	: "%s kcal/mol" % str(6.3/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13361 , {'ddG'	: "%s kcal/mol" % str(-0.4/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13362 , {'ddG'	: "%s kcal/mol" % str(2.2/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13363 , {'ddG'	: "%s kcal/mol" % str(3.0/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13364 , {'ddG'	: "%s kcal/mol" % str(2.0/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13365 , {'ddG'	: "%s kcal/mol" % str(1.1/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13366 , {'ddG'	: "%s kcal/mol" % str(0.5/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13367 , {'ddG'	: "%s kcal/mol" % str(1.3/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13368 , {'ddG'	: "%s kcal/mol" % str(2.9/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13369 , {'ddG'	: "%s kcal/mol" % str(4.8/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	(13370 , {'ddG'	: "%s kcal/mol" % str(3.7/NUMBER_KJ_IN_KCAL), 'PDB' : '2RN2'}),
	
	# PMID:7473760. Loss of precision on conversion from kJ/mol to kcal/mol.
	(13422 , {'ddG'	: "%s kcal/mol" % str(-1.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13423 , {'ddG'	: "%s kcal/mol" % str(-5.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13424 , {'ddG'	: "%s kcal/mol" % str(-4.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13425 , {'ddG'	: "%s kcal/mol" % str(-2.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13426 , {'ddG'	: "%s kcal/mol" % str(-3.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	
	# PMID:9010773. Loss of precision on conversion from kJ/mol to kcal/mol.
	(13427 , {'ddG'	: "%s kcal/mol" % str(-15.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	
	# PMID:9020766. Loss of precision on conversion from kJ/mol to kcal/mol.
	(13428 , {'ddG'	: "%s kcal/mol" % str(-6.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13429 , {'ddG'	: "%s kcal/mol" % str(-1.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13430 , {'ddG'	: "%s kcal/mol" % str(-3.1/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13431 , {'ddG'	: "%s kcal/mol" % str(-4.1/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13432 , {'ddG'	: "%s kcal/mol" % str(-1.1/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13433 , {'ddG'	: "%s kcal/mol" % str( 2.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13434 , {'ddG'	: "%s kcal/mol" % str(-6.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13435 , {'ddG'	: "%s kcal/mol" % str(-5.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(13436 , {'ddG'	: "%s kcal/mol" % str(-3.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	
	# PMID:8433369. Table 6 has more precise DDG values at pH 3 than Table 2.
	(13817 , {'ddG'	: '-5.00 kcal/mol', 'PDB' : '2LZM'}),
	(13818 , {'ddG'	: '-2.27 kcal/mol', 'PDB' : '2LZM'}),
	(13819 , {'ddG'	: '-1.40 kcal/mol', 'PDB' : '2LZM'}),
	(13820 , {'ddG'	: '-0.75 kcal/mol', 'PDB' : '2LZM'}),
	(13821 , {'ddG'	: '-0.36 kcal/mol', 'PDB' : '2LZM'}),
	(13822 , {'ddG'	: '-3.53 kcal/mol', 'PDB' : '2LZM'}),
	(13823 , {'ddG'	: '-1.78 kcal/mol', 'PDB' : '2LZM'}),
	(13824 , {'ddG'	: '-0.49 kcal/mol', 'PDB' : '2LZM'}),
	(13825 , {'ddG'	:  '0.20 kcal/mol', 'PDB' : '2LZM'}),
	(13826 , {'ddG'	: '-0.81 kcal/mol', 'PDB' : '2LZM'}),
	(13827 , {'ddG'	: '-8.30 kcal/mol', 'PDB' : '2LZM'}),
	
	# PMID:9677301. Loss of precision on conversion from kJ/mol to kcal/mol.
	(14189 , {'ddG'	: "%s kcal/mol" % str(-1.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14190 , {'ddG'	: "%s kcal/mol" % str(-5.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14191 , {'ddG'	: "%s kcal/mol" % str(-3.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	#14192 was entered incorrectly so is corrected above
	(14193 , {'ddG'	: "%s kcal/mol" % str(-4.1/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14194 , {'ddG'	: "%s kcal/mol" % str(-6.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14195 , {'ddG'	: "%s kcal/mol" % str(-1.8/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14196 , {'ddG'	: "%s kcal/mol" % str(-4.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14197 , {'ddG'	: "%s kcal/mol" % str(-3.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14198 , {'ddG'	: "%s kcal/mol" % str(-2.4/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14199 , {'ddG'	: "%s kcal/mol" % str( 0.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14200 , {'ddG'	: "%s kcal/mol" % str(-7.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14201 , {'ddG'	: "%s kcal/mol" % str(-6.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14202 , {'ddG'	: "%s kcal/mol" % str(-4.7/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	
	# PMID:9649316. Loss of precision on conversion from kJ/mol to kcal/mol.
	(14203 , {'ddG'	: "%s kcal/mol" % str(-2.1/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14204 , {'ddG'	: "%s kcal/mol" % str(-0.8/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14205 , {'ddG'	: "%s kcal/mol" % str( 0.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14206 , {'ddG'	: "%s kcal/mol" % str(-4.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14207 , {'ddG'	: "%s kcal/mol" % str(-1.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14208 , {'ddG'	: "%s kcal/mol" % str(-1.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	
	# PMID:9685334. Loss of precision on conversion from kJ/mol to kcal/mol.
	(14224 , {'ddG'	: "%s kcal/mol" % str(-4.97/NUMBER_KJ_IN_KCAL), 'PDB' : '1HUE'}),
	(14225 , {'ddG'	: "%s kcal/mol" % str( 1.72/NUMBER_KJ_IN_KCAL), 'PDB' : '1HUE'}),
	(14226 , {'ddG'	: "%s kcal/mol" % str(-3.43/NUMBER_KJ_IN_KCAL), 'PDB' : '1HUE'}),
	(14227 , {'ddG'	: "%s kcal/mol" % str(-0.54/NUMBER_KJ_IN_KCAL), 'PDB' : '1HUE'}),
	(14228 , {'ddG'	: "%s kcal/mol" % str(  0.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1HUE'}),
	
	# PMID:9398521. Loss of precision on conversion from kJ/mol to kcal/mol.
	(14235 , {'ddG'	: "%s kcal/mol" % str(-10.6/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14236 , {'ddG'	: "%s kcal/mol" % str(-15.5/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14237 , {'ddG'	: "%s kcal/mol" % str( -7.2/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14238 , {'ddG'	: "%s kcal/mol" % str(-11.3/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14239 , {'ddG'	: "%s kcal/mol" % str( -3.9/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),
	(14240 , {'ddG'	: "%s kcal/mol" % str(-16.0/NUMBER_KJ_IN_KCAL), 'PDB' : '1LZ1'}),	
]

# Check for collisions within patch arrays
for patch_array in [OverriddenEntries, BadOrMissingMutants, ddGTypos, RoundingErrors]:
	colortext.write("Checking patch array: ", "silver")
	okay = True
	duplicatecount = {}
	for tpl in patch_array:
		duplicatecount[tpl[0]] = duplicatecount.get(tpl[0], 0) + 1 
	for k, v in duplicatecount.iteritems():
		if v > 1:
			okay = False
			colortext.error("Duplicate patch found for record #d." % k)
	if okay:
		colortext.message("okay.")

# Merge the patch datasets into one dict and check that no collisions occur
PDBTagSet = set(["PDB"])
MergedPatchSet = {}
for patch_array in [OverriddenEntries, BadOrMissingMutants, ddGTypos, RoundingErrors]:
	for tpl in patch_array:
		ID = tpl[0]
		patchrecord = tpl[1]
		if MergedPatchSet.get(ID):
			assert(patchrecord["PDB"] == MergedPatchSet[ID]["PDB"])
			assert(set(patchrecord.keys()).intersection(set(MergedPatchSet[ID].keys())) == PDBTagSet)
			# Merge
			for k, v in patchrecord.iteritems():
				MergedPatchSet[ID][k] = v
		else:
			MergedPatchSet[ID] = {}
			for k, v in patchrecord.iteritems():
				MergedPatchSet[ID][k] = v

# Make sure that the patch dict does not overwrite the ddGWrongSigns dict
for ID, patchrecord in MergedPatchSet.iteritems():
	if ID in ddGWrongSigns:
		assert(ddGWrongSigns[ID] not in patchrecord.keys())

# PMID: 7507755. 1L63 is missing N163 and L164 so we use 219L instead.
PseudoLysozyme163Cases = [13516, 13517]
for ID in PseudoLysozyme163Cases:
	PseudoT4LysozymeCases.remove(ID)
	assert(not(MergedPatchSet.get(ID)))
	MergedPatchSet[ID] = {'PDB_wild'	 : '219L', 			'PDB' : '2LZM'}

for ID in PseudoT4LysozymeCases:
	if MergedPatchSet.get(ID):
		assert(MergedPatchSet[ID]['PDB'] == '2LZM')
		assert(not(MergedPatchSet[ID].get('PDB_wild')))
		MergedPatchSet[ID]['PDB_wild'] = '1L63'
	else:
		MergedPatchSet[ID] = {'PDB_wild'	 : '1L63', 			'PDB' : '2LZM'}

for ID in PseudoHumanLysozymeCases:
	if MergedPatchSet.get(ID):
		assert(MergedPatchSet[ID]['PDB'] == '1LZ1')
		assert(not(MergedPatchSet[ID].get('PDB_wild')))
		MergedPatchSet[ID]['PDB_wild'] = '2BQA'
	else:
		MergedPatchSet[ID] = {'PDB_wild'	 : '2BQA', 			'PDB' : '1LZ1'}

# In the records below where the wildtype for Bacillus subtilis (BsHPr) was entered incorrectly, 3OQN may be an appropriate PDB file as one chain has the correct wildtype sequence.
recordsWithUnresolvedMissingData = []
# PMID:8448200. Records 203-206. http://www.sciencedirect.com/science/article/pii/016748389390133C#
# The protein is human cyanomet myoglobin.
# The UniProt AC P02144 in ProTherm corresponds to PDB file 3RGK but that file is a K45R mutant of human cyanomet myoglobin. 
recordsWithUnresolvedMissingData.extend(range(203, 206 + 1))
# PMID:2765493. Records 3875-3883.
# Again, the protein is human cyanomet myoglobin.
# The UniProt AC P02144 in ProTherm corresponds to PDB file 3RGK but that file is a K45R mutant of human cyanomet myoglobin. 
recordsWithUnresolvedMissingData.extend(range(3875, 3883 + 1))
# PMID:1581299. Records 4215, 4216.
# This structure with UniProt AC P00912 does not appear to have been solved yet.
recordsWithUnresolvedMissingData.append(4215)
recordsWithUnresolvedMissingData.append(4216)
# PMID:8755725. Records 4670, 4671, 4672.
# This structure with UniProt AC P67881 does not appear to have been solved yet.
recordsWithUnresolvedMissingData.extend(range(4670, 4672 + 1))
# PMID:7727438. Records 5090, 5091.
# This structure with UniProt AC P19756 does not appear to have been solved yet.
recordsWithUnresolvedMissingData.append(5090)
recordsWithUnresolvedMissingData.append(5091)
# PMID:8621402. Records 5241-5244.
# HIV-1 protease. This structure has no UniProt AC given. Has it been solved?
recordsWithUnresolvedMissingData.extend(range(5241, 5244 + 1))
# PMID:10884358. Records 8380-8383.
# Acylphosphatase (Synthetic). This structure has no UniProt AC given. Has it been solved?
recordsWithUnresolvedMissingData.extend(range(8380, 8383 + 1))
# PMID:10873472. Records 8384-8401.
# Coiled-coil protein (Synthetic). De novo design which does not appear to have been solved.
recordsWithUnresolvedMissingData.extend(range(8384, 8401 + 1))
# PMID:8580838. Records 8655, 14481
# I did not add these records as the best PDB match for bsHPr is 2HID which contains the mutation M51V. 2HPR has M51V and S83C.
recordsWithUnresolvedMissingData.append(8655)
recordsWithUnresolvedMissingData.append(14481)
# PMID:8161701. Record 8659.
# The mutation is of a hybrid mutant, nuclease conA S28G, of wild-type Staphylococcal nuclease A. No PDB ID given.
recordsWithUnresolvedMissingData.append(8659)
# PMID:9265621. Records 8819-8824.
# No PDB/Uniprot ID given. The PMD says that the source is Bacillus amyloliquefaciens and the sequence is:
# AGKSNGEKKYIVGFKQTMSTMSAAKKKDVISEKGGKVQKQFKYVDAASATLNEKAVKELKKDPSVAYVEEDHVAHAY
# 1SPB seems the closest match and has two authors in common.
recordsWithUnresolvedMissingData.extend(range(8819, 8824 + 1))
# PMID:8580841. Records 9138-9142
# Thioredoxin (Rhodobacter sphaeroides), UniProt AC P08058. UniProt has no PDB mapping. Has it been solved?
recordsWithUnresolvedMissingData.extend(range(9138, 9142 + 1))
# PMID:8539251. Records 12366-12398
# TEM beta-lactamase from an unknown source. This protein seems to be described in previous paper PMID:7990143. Is there a corresponding PDB file?
recordsWithUnresolvedMissingData.extend(range(12366, 12398 + 1))
# PMID:10421435. Records 12944,12945.
# There is no mutation here. The UniProt AC Q96HL2 had no corresponding PDB ID.
recordsWithUnresolvedMissingData.extend(range(12944, 12945 + 1))
# PMID:9685334. Records 14224-14233.
# The UniProt AC P08821 had no corresponding PDB ID.
recordsWithUnresolvedMissingData.extend(range(14224, 14233 + 1))
# PMID: 12144791. Records 15480-15486, 15497-15500
# I did not add these records as the best PDB match for bsHPr is 2HID which contains the mutation M51V. 2HPR has M51V and S83C.
recordsWithUnresolvedMissingData.extend(range(15480, 15486 + 1))
recordsWithUnresolvedMissingData.extend(range(15497, 15500 + 1))
# PMID: 10079065. Records 16294-16298.
# The structure does not appear to have been solved ("the combination of N40(fragment 1-40) and C16 (fragment 41-56)")
recordsWithUnresolvedMissingData.extend(range(16294, 16298 + 1))
# PMID: 10079065. Records 16471-16474.
# The structure does not appear to have been solved.
recordsWithUnresolvedMissingData.extend(range(16471, 16474 + 1))
# PMID: 12974622. Record 16701.
# The closest PDB match is 2KYC which is a point mutation, C72S, of chicken parvalbumin 3 (CPV3):
# "Because the wild-type CPV3 sequence includes a solvent-exposed cysteine at position 72, rather than the consensus serine, the protein readily forms disulfide-linked dimers and trimers in the absence of reductant. We have chosen to work with the C72S variant to avoid the experimental complications attendant to this behavior. This sequence substitution has no discernible impact on divalent ion affinity."
# However, this mutation here *is* the C72S mutation. Also, is the published DDG between rat and chicken proteins?
recordsWithUnresolvedMissingData.append(16701)
# PMID: 14978309. Records 16817-16832.
# These peptide structures do not seem to have been solved.
recordsWithUnresolvedMissingData.extend(range(16817, 16832 + 1))
# PMID: 14756573. Records 16836-16851
# As with PMID:12144791 above, I chose to omit the records with bsHPr since I could not find a solved wildtype structure.
recordsWithUnresolvedMissingData.extend(range(16836, 16851 + 1))
# PMID: 15023058. Records 16941-16944
# These artificial Antennafinger (Ant-F) proteins do not seem to have been solved.
recordsWithUnresolvedMissingData.extend(range(16941, 16944 + 1))
# PMID: 15504416. Records 17916-17925
# These records involve Staphylococcal nuclease from various sources. A quick PDB search does not show solved structures from these sources.
recordsWithUnresolvedMissingData.extend(range(17916, 17925 + 1))
# PMID:12069590. Records 15409-15451
# I am not completely sure that 1YYJ is the correct PDB ID due to the difference at position 98 and the described mutation.  
recordsWithUnresolvedMissingData.extend(range(15409, 15451 + 1))
# PMID: 15533036. Records 18311-18324. See reference 11, PMID:12069590
# I am not completely sure that 1YZA is the correct PDB ID due to the uncertainty with 1YYJ for records 15409-15451. 
recordsWithUnresolvedMissingData.extend(range(18311, 18324 + 1))
# PMID:16566582. Records 19960-19965
# I did not add these records as the best PDB match for bsHPr is 2HID which contains the mutation M51V. 2HPR has M51V and S83C.
recordsWithUnresolvedMissingData.extend(range(19960, 19965 + 1))
# PMID: 16705642. Records 20013-20048.
# These records involve a synthetic Eglin C variant, Eglin C-F10W. A quick PDB search does not show solved structures from this source.
recordsWithUnresolvedMissingData.extend(range(20013, 20048 + 1))
# PMID: 16762367. Records 20104-20111.
# These records involve Family 10 Xylanase, Cellvibrio mixtus (CmXyn10B).
# The paper gives the PDB ID as 1UR1. "The crystal structure of the double mutant A334V/G348D was solved to 2.4A resolution. Crystals were obtained only in the presence of AX_2 and no crystal of the quadruple mutant A26C/A334V/G348D/L380C was generated."
# 2CNC appears to be one of the mutants mentioned in the paper but neither PDB ID matches the wild type w.r.t. the residue positions.
recordsWithUnresolvedMissingData.extend(range(20104, 20111 + 1))
# PMID: 16501226. Records 20479-20567.
# These records involve the human FynSH3 domain. A quick PDB search does not show solved structures matching the wildtype amino acids I28 and V55 in the publication. 
recordsWithUnresolvedMissingData.extend(range(20479, 20567 + 1))
# PMID: 17188709. Records 20665-21892.
# Missing the wildtype for the synthetic CspB protein (CspB-TB).
# These mutations are from wildtype CspB-Bs (1CSP). The other mutations in this range (21040-21388) are of the 6H-CspB-Bs* His-tagged variant (6H-WT*) or of the variant with the His-tag removed (WT*) so I did not include them as the DG values differ and we have no associated solved structures.
CspBBsRecords = set([
	21042, 21043, 21046, 21050, 21074, 21084, 21087, 21092,
	21099, 21100, 21103, 21107, 21131, 21141, 21144, 21149,
	21159, 21160, 21163, 21168, 21192, 21202, 21205, 21210,
	21217, 21218, 21221, 21226, 21250, 21260, 21263, 21268,
	21275, 21276, 21279, 21284, 21308, 21318, 21321, 21326,
	21333, 21334, 21337, 21342, 21366, 21376, 21379, 21384,])
recordsWithUnresolvedMissingData.extend([i for i in range(20665, 21892 + 1) if i not in CspBBsRecords]) 
# PMID: 16799151. Records 22700-22713.
# These records involve HFV protease, human. A quick PDB search does not show solved structures. 7HVP is a HIV-1 complex.
recordsWithUnresolvedMissingData.extend(range(22700, 22713 + 1))
# PMID: 19683006. Some records in the range 24245-24286.
# These records involve the P53 DNA binding domain from various sources.
# The human source has a PDB ID (2AC0) but at some some of the other sources (Drosophila melanogaster, Caenorhabditis elega, Xenopus laevis, Reflectopallium marm, Mus domestica, Bovine, Mouse, Goat) do not.
recordsWithUnresolvedMissingData.extend([24265, 24272, 24273, 24274])
recordsWithUnresolvedMissingData.extend(range(24267, 24270 + 1))
# PMID: 18586378. Records 24714-24717. 
# These records involve the human GTPase effector domain. A quick PDB search does not show solved structures from this source. 
recordsWithUnresolvedMissingData.extend(range(24714, 24717 + 1))
# PMID: 18840434. Records 25005-25006. 
# These records involve a fusion of the amyloidogenic Abeta42 species to GFP. A quick PDB search does not show solved structures from this source. 
recordsWithUnresolvedMissingData.extend(range(25005, 25006 + 1))
recordsWithUnresolvedMissingData = set(recordsWithUnresolvedMissingData)

# Records with the wrong PMID
badPublicationReferences = {}
badPublicationReferences[15502] = 12079391
for recordID in range(13376, 13381 + 1):
	badPublicationReferences[recordID] = 8390295
for recordID in range(15714, 15733 + 1):
	badPublicationReferences[recordID] = 12473461
for recordID in range(12725, 12730 + 1):
	badPublicationReferences[recordID] = 11714927
for recordID in range(14301, 14305 + 1):
	badPublicationReferences[recordID] = 11714927
#for recordID in range(14571, 14574 + 1):
#	badPublicationReferences[recordID] = 11714927

# Records with the wrong secondary structure
badSecondaryStructure = dict.fromkeys(
	[2747, 4611, 11869, 12310, 12701, 12702, 12979, 12980, 12982, 12983, 16895, 19886, 19887, 19888, 19889, 19893, 22231, 24335]
	+ range(15529, 15534 + 1)
	+ range(24921, 24929 + 1)
	+ range(24931, 24939 + 1)
	+ range(24962, 24964 + 1) + [24966]
	+ range(24968, 24983 + 1)
	+ range(24985, 25000 + 1)
	,True)

# Records with the wrong ASA
badASA = dict.fromkeys(
	[11869, 14408, 14409, 14413, 14414, 14434, 14438, 14450, 16895, 19886, 19887, 19888, 19889, 19893, 22231]
	+ range(15529, 15534 + 1)
	, True)


# Set of fixed records
fixedRecordIDs = set()
for dset in [MergedPatchSet, ddGWrongSigns, badPublicationReferences, badSecondaryStructure, badASA]:
	fixedRecordIDs = fixedRecordIDs.union(dset.keys())

# In these cases, the protein is elongated at the 67th position. This is more than a mutation so I ignore it. 	
skipTheseCases = [12156, 12159, 12161, 12188, 12191, 12193, 12218, 12220, 14468]

# In these cases of PMID:16503630, the PDB ID in the publication is wrong - this is not human ubiquitin.
skipTheseCases.extend([19883, 19884, 19885, 19886, 19887, 19888, 19889])

# In this case, the mutation is an insertion mutation. This is more than a mutation so I ignore it. 	
skipTheseCases.extend([17241])

# todo: These cases need to be checked
skipTheseCases.extend([17756, 18137, 4596])

# In this case, the wild type is Ile but while the paper writes 'V22A', ProTherm writes 'I 22 A (PDB: I 23 A; PIR: I 1785 A)'. Do they mutate to Val and then to Ala in the paper? 
skipTheseCases.extend([18114])

# PMID: 16922511. Deletion mutants
skipTheseCases.extend([20175, 20181])

# PMID:19647749. dG_H2O value entered as ddG_H2O.
skipTheseCases.append(24386)

#These cases fail parsing the mutation line - need to write a new regex
#skipTheseCases.extend([19893,19894,19895])

# In this case, I couldn't map the paper's structure well enough onto 1YCC. The mutations mentioned e.g. A57 N->I do not correspond with the PDB structure (attributed to a later paper).
skipTheseCases.append(11817)

# In this case, the structural information needed for the mutation (residues A62) is missing in the PDB file
# A57 is also missing in chain A but is present in chain B as A1057.
skipTheseCases.extend([13452])

# The DDG value here is a limit rather than a value (the paper states that the mutant (G107V) "was too unstable to yield an accurate determination of DDG_H20".
skipTheseCases.append(2043)

skipTheseCases = set(skipTheseCases)

# Make sure that we do not repeat records to be skipped in the recordsWithUnresolvedMissingData set
assert(len(set(skipTheseCases).intersection(recordsWithUnresolvedMissingData)) == 0)

# Mutations involving cysteine which have a different format involving bridges (S-H or S-S)
CysteineMutationCases = [13663, 13664, 13677, 13678]

# Mutations with different parsing requirements and their regexes
multimapCases1 = [16597, 16598, 16599, 16600, 19897, 19898, 19893, 19894, 22383, 22384]
mmapCases1regex = re.compile("PDB:(.+[,]{0,1})+\)")

multimapCases2Ranges = ((17681, 17687), (17692, 17698), (18104, 18105), (18108, 18136), (18138, 18175))
multimapCases2 = []
for m2r in multimapCases2Ranges:
	multimapCases2.extend(range(m2r[0], m2r[1] + 1))
#multimapCases2 = range(17681, 17687 + 1) + range(17692, 17698 + 1) + range(18104, 18105 + 1) +
#				  range(18108, 18136 + 1) + range(18138, 18175 + 1)
mmapCases2regex = re.compile("^.*PDB:(.*);PIR.*$")	

		
multimapCases3Ranges = ((3352, 3383), (14215, 14223), (8470, 8504), (12235, 12237),
					(12308, 12310), (15402, 15408), (16251, 16253), (16255, 16257), 
					(16259, 16261), (16263, 16265))
multimapCases3 = []
for m3r in multimapCases3Ranges:
	multimapCases3.extend(range(m3r[0], m3r[1] + 1))
multimapCases3.extend([6366, 6368, 10384, 16991, 17678, 17679, 17680, 17689, 17690, 17691])
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
	"J AM CHEM SOC 115, 8523-8526 (1993) "		: ("PMID", 99999999999), # No PMID for this article. "Phospholipase A2 Engineering.  10.  The Aspartate...Histidine Catalytic Diad Also Plays an Important Structural Role."  Y. Li and M.-D. Tsai, J. Am. Chem. Soc. 115, 8523-8526 (1993).
	"J MOL BIOL 351, 402-416 (2005)"			: ("PMID", 16002092),
	"PROTEIN SCI 16, 227-238 (2007)"			: ("PMID", 17189482),
} 


# todo: Unused data
PMIDReferencesWhichICouldNotAccess = ['PMID:14529489']
DuplicatedRecords = [
	(1163, 13570),
	(8911, 14482),
	(12193, 14468),
	(13410, 13182), # Maybe more from this publication
	(13393, 13183), # Maybe more from this publication
	(8363, 13477), # Publication states that the duplication but the T values do not match
	(8364, 13478), # Publication states that the duplication but the T values do not match
	(8365, 13479), # Publication states that the duplication but the T values do not match
]
# These publications have DDG values but none are stored in ProTherm.
MissingPublicationsPMIDsWithDDGValuesForReview = [10350481, 10504240, 10623513]

def getDDGUnitsUsedInDB(ddGDB):
	results = ddGDB.locked_execute('SELECT ID, DGUnitUsedInProTherm FROM Source')
	unitsUsed = {}
	for r in results:
		unitsUsed[r["ID"]] = r["DGUnitUsedInProTherm"]
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
		25844 : "2012-07-31",
	}
	
	def __init__(self, infilepath, ddGDB = None, quiet = False, skipIndexStore = False):
		
		if not ddGDB:
			ddGDB = ddgdbapi.ddGDatabase()
		self.ddGDB = ddGDB
		
		mtchs = re.match(".*(ProTherm)(\d+)[.]dat$", infilepath, re.IGNORECASE)
		if mtchs:
			lastrecord = int(mtchs.group(2))
		if lastrecord in [23581, 25616, 25844]:
			# These fields of ProTherm records cannot be empty for our purposes
			self.requiredFields = ["NO.", "PDB_wild", "LENGTH", "MUTATION", "MUTATED_CHAIN"]
			self.quiet = quiet
			# For bad data
			self.iCodeRecords = iCodeRecords
			self.singleChainPDBs = singleChainPDBs
			self.identicalChainPDBs = identicalChainPDBs
			self.MergedPatchSet = MergedPatchSet
			self.ddGWrongSigns = ddGWrongSigns
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
			self.mutationsAllowedToBeStoredDespiteMissingCoordinates = mutationsAllowedToBeStoredDespiteMissingCoordinates
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
		self.singleErrors["either_ddG"] = 0
		
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
	
	def createGenerator(self, require_one_of = [], specific_cases = None, args = None, start_index = None, end_index = None):
		'''This function creates a generator for iterating over the set of records.
		It was motivated by me duplicating similar code in various scripts which use the reader.
		 
		specific_cases is a list of specific cases to test.
		args is expected to be a list/tuple of up to two values: the start index and the end index.
		args is ignore if specific_cases is given.'''
		
		if not self.test():
			raise Exception("The ProTherm reader self-test failed.")
		
		cases = specific_cases
		if not cases:
			if args:
				if args[0].isdigit():
					if len(args) > 1 and args[1].isdigit():
						cases = [c for c in self.list_of_available_keys if int(args[0]) <= c <= int(args[1])]
						if not self.quiet:
							colortext.message("[Starting from record %s]" % args[0])
					else:
						cases = [c for c in self.list_of_available_keys if int(args[0]) <= c]
						if not self.quiet:
							colortext.message("[Starting from record %s]" % args[0])
			else:
				cases = self.list_of_available_keys
		else:
			assert(isinstance(cases, list))
			cases = [c for c in cases if c in self.list_of_available_keys]
			
		if not self.quiet:
			colortext.message("\nParsing ProTherm")
			colortext.printf("|" + ("*" * (int(len(self.list_of_available_keys)/1000) - 2)) + "|")
		thousands_read = 0
		for ID in cases:
			#Progress meter
			if not self.quiet:
				if ID/1000 > thousands_read:
					colortext.write("." * ((ID/1000) - thousands_read), "green")
					colortext.flush()
					thousands_read += (ID/1000) - thousands_read
			
			if start_index != None and ID < start_index:
				continue
			if end_index != None and ID > end_index:
				continue
			
			# Skip bad cases in ProTherm
			if ID in self.skipTheseCases:
				continue
			
			record = self.readRecord(ID)		
			
			fixed_record = self.readRecord(ID)
			self.fixRecord(ID, fixed_record, quiet = True)
			
			# Skip records where none of the require_one_of keys has a corresponding value 
			found_one_required_key = False or require_one_of == []
			for recordkey in require_one_of:
				if record[recordkey] != None:
					found_one_required_key = True
					break
			if not found_one_required_key:
				for recordkey in require_one_of:
					if fixed_record[recordkey] != None:
						found_one_required_key = True
						break
			if not found_one_required_key:
				continue
			
			yield (ID, record)
			
		if not self.quiet:
			colortext.write("." * ((len(self.list_of_available_keys)/1000) - thousands_read), "green")
			colortext.flush()
			print("")

	def testRounding(self):
		'''This function checks that the rounding table does not introduce large changes in DDG values due to input error.'''

		colortext.write("Testing rounding table for errors: ", color = "silver")
		allowed_drift = 0.0501 # The maximum amount in kcal/mol that the fixed DDG value is allowed to deviate from the original DDG value
		founderror = False
		for tpl in RoundingErrors:
			ID = tpl[0]
			patch_dict = tpl[1]
			assert(patch_dict.get("ddG") != None or patch_dict.get("ddG_H2O") != None)
			
			record = self.readRecord(ID)
			if patch_dict.get("ddG") != None:
				original_ddG = self.getDDGInKcal(ID, record, useRosettaConvention = False)
				record["ddG"] = patch_dict["ddG"]
				patched_ddG = self.getDDGInKcal(ID, record, useRosettaConvention = False)
				if abs(original_ddG - patched_ddG) > allowed_drift:
					colortext.write("\nRecord %d: The original rounded DDG value %s kcal/mol differs from the fixed value %s kcal/mol by %s kcal/mol (DDG in ProTherm convention)" % (ID, original_ddG, patched_ddG, abs(original_ddG - patched_ddG)), color="red")
					founderror = True
			if patch_dict.get("ddG_H2O") != None:
				original_ddG = self.getDDGH2OInKcal(ID, record, useRosettaConvention = False)
				record["ddG_H2O"] = patch_dict["ddG_H2O"]
				patched_ddG = self.getDDGH2OInKcal(ID, record, useRosettaConvention = False)
				if abs(original_ddG - patched_ddG) > allowed_drift:
					colortext.error("\nRecord %d: The original rounded DDG_H2O value %s kcal/mol differs from the fixed value %s kcal/mol by %s kcal/mol (DDG in ProTherm convention)" % (ID, original_ddG, patched_ddG, abs(original_ddG - patched_ddG)), color="red")
					founderror = True
		if not founderror:
			colortext.message("passed")
		else:
			print("")
			
	def test(self):
		self.testRounding()
		MutationO = ddgobjects.Mutation
		success = True
		expected_results = {
		# PLAIN
			# L 121 A, A 129 M, F 153 L
			1163  : [MutationO('L',  '121', 'A'), MutationO('A',  '129', 'M'), MutationO('F',  '153', 'L')],
			# S 31 G, E 33 A, E 34 A
			1909  : [MutationO('S',   '31', 'G'), MutationO('E',   '33', 'A'), MutationO('E',   '34', 'A')],
			#Q 15 I, T 16 R, K 19 R, G 65 S, K 66 A, K 108 R
			2227  : [MutationO('Q',   '15', 'I'), MutationO('T',   '16', 'R'), MutationO('K',   '19', 'R'), MutationO('G',   '65', 'S'), MutationO('K',   '66', 'A'), MutationO('K',  '108', 'R')],
			# I 6 R, T 53 E, T 44 A
			12848 : [MutationO('I',    '6', 'R'), MutationO('T',   '53', 'E'), MutationO('T',   '44', 'A')],
			# E 128 A, V 131 A, N 132 A
			17608 : [MutationO('E',  '128', 'A'), MutationO('V',  '131', 'A'), MutationO('N',  '132', 'A')],
		# PLAIN BUT BADLY FORMATTED
			11146 : [MutationO('G',   '23', 'A'), MutationO('G',   '25', 'A')],
			21156 : [MutationO('E',    '3', 'R'), MutationO('F',   '15', 'A'), MutationO('D',   '25', 'K')],
		# MUTATION (PDB: MUTATION; PIR MUTATION), repeat
			#S 15 R (PDB: S 14 R; PIR: S 262 R), H 19 E (PDB: H 18 E; PIR: H 266 E), N 22 R ( PDB: N 21 R; PIR: N 269 R)
			14215 : [MutationO('S',   '14', 'R'), MutationO('H',   '18', 'E'), MutationO('N',   '21', 'R')],
		# MUTATION_LIST (PDB: MUTATION_LIST; PIR MUTATION_LIST)
			# N 51 H, D 55 H (PDB: N 47 H, D 51 H; PIR: N 135 H, D 139 H)
			17681 : [MutationO('N',   '47', 'H'), MutationO('D',   '51', 'H')],
		# MUTATION_LIST (PDB: MUTATION_LIST)
			# A 2 K, L 33 I (PDB: A 1 K, L 32 I)
			22383 : [MutationO('A',    '1', 'K'), MutationO('I',   '23', 'V')],
			16597 : [MutationO('G',   '92', 'P')],
			19893 : [MutationO('Q',  '11A', 'E'), MutationO('H',  '53A', 'E'), MutationO('S',  '76A', 'K'), MutationO('M',  '78A', 'K'), MutationO('D',  '81A', 'K')],
		}
		
		errors = []
		colortext.write("Testing %d mutation formats: " % len(expected_results), "silver")	
		for ID, expected in sorted(expected_results.iteritems()):
			numerrormsgs = len(errors)
			mutations = self.getMutations(ID)
			if expected and (expected != mutations):
				errors.append("Record %d: Expected %s, got %s." % (ID, expected, mutations))
			if (not badSecondaryStructure.get(ID)):
				mutationsWithMissingPosition = [mutation for mutation in mutations if mutation.SecondaryStructurePosition == None]
				if mutationsWithMissingPosition:
					errors.append("%s: \n" % record["MUTATION"])
					errors.append("Record %d: Read %d mutations but %d have no secondary structure positions." % (ID, len(mutations), len(mutationsWithMissingPosition)))
					errors.append(str(mutations))
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
			5980  : 1.1 / NUMBER_KJ_IN_KCAL, # 1.1 kJ/mol
			12144 : -3.9 / NUMBER_KJ_IN_KCAL, #-3.9 kJ/mol
			13535 : -0.14, #-0.14
			14451 : -20.7 / NUMBER_KJ_IN_KCAL, # -20.7 kJ/mol
			14751 : 3.35, # 3.35
			17849 : 12.7, # 12.7 kcal/mol
			17858 : 0.06, # 0.06 kcal/mol
			20129 : 0.08, # 0.08 
			21100 : 9.5 / NUMBER_KJ_IN_KCAL, # 9.5 kJ/mol,
			22261 : -0.4 / NUMBER_KJ_IN_KCAL, # -0.4 kJ/mol
			23706 : -2.8, # -2.8
			25089 : -3.3, # -3.3
		}
		errors = []
		colortext.write("Testing %d ddG conversions: " % len(expected_results), "silver")	
		for ID, expectedDDG in sorted(expected_results.iteritems()):
			numerrormsgs = len(errors)
			record = self._getRecord(ID, None)
			referenceID = self.getReference(ID, record)
			#record["dbReferencePK"] = "PMID:%s" % referenceID
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
			try:
				mts = self.getMutations(ID)
			except:
				pass
			if mts and len(mts) == 1:
				#if str(mts[0]["ResidueID"]) == "32" and  m["WildTypeAA"] == 'A':
				#	print(ID, mts, record["ddG"], record["PDB_wild"])
				if str(mts[0].ResidueID) == "45" and mts[0].WildTypeAA == 'A' and mts[0].MutantAA == 'G':
					print(ID, mts, record["ddG"], record["PDB_wild"])
			
			if (pdbID == None) or (pdbID != None and record["PDB_wild"] == pdbID):
				#if mts:
				#	print(ID, mts)
				if not(hasddG) or (hasddG == True and record["ddG"] != None):
					if (numMutations == None) or (mts and numMutations == len(mts)):
						foundResidueID = False
						if residueIDs:
							for m in mts:
								if m.ResidueID in residueIDs and m.WildTypeAA == 'A':
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
			if self.MergedPatchSet.get(ID) and self.MergedPatchSet[ID].get("Tm"):
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
				# hack: We usually expect a float value so we hack the odd cases
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
				# hack: for special case
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
	
	def fixRecord(self, ID, record = None, quiet = False):
		'''Override bad data in ProTherm.'''
		record = self._getRecord(ID, record)
		
		# Patch for typos
		if 24390 <= ID <= 24420 and record["ddG_H2O"] and record["ddG_H2O"].find("kal/mol") != -1:
			record["ddG_H2O"] = record["ddG_H2O"].replace("kal/mol", "kcal/mol")
		
		passed = True
		
		singleChainPDBs = self.singleChainPDBs
		identicalChainPDBs = self.identicalChainPDBs
		
		# Apply the override patches
		MergedPatchSet = self.MergedPatchSet
		if MergedPatchSet.get(ID):
			if record.get('PDB_wild'):
				if record["PDB_wild"] != MergedPatchSet[ID]["PDB"]:
					raise colortext.Exception("Error in MergedPatchSet table: Record %d. Read '%s' for PDB_wild, expected '%s'." % (ID, record["PDB_wild"], MergedPatchSet[ID]["PDB"]))
			for k, v in MergedPatchSet[ID].iteritems():
				if k != "PDB":
					record[k] = v
	
		# Fix the chain when it is possible to be ambiguous up to homology
		if record["PDB_wild"]:
			pdbID = record["PDB_wild"].upper()
			if not(MergedPatchSet.get(ID)) or not(MergedPatchSet[ID].get("MUTATED_CHAIN")):
				# Update PDB IDs so long as we don't have them in the override dict 
				if singleChainPDBs.get(pdbID):
					for k,v in singleChainPDBs[pdbID].iteritems():
						record[k] = v
				elif identicalChainPDBs.get(pdbID):
					for k,v in identicalChainPDBs[pdbID].iteritems():
						record[k] = v
		
		# Reverse the signs of DDG/DDG_H2O values where wrong
		ddGWrongSigns = self.ddGWrongSigns
		if ddGWrongSigns.get(ID):
			tokens = record[ddGWrongSigns[ID]].split(" ")
			ddGvalue = tokens[0]
			newddGvalue = str(-float(ddGvalue)) # This should always pass
			#print(ddGvalue, newddGvalue)
			record[ddGWrongSigns[ID]] =  "%s %s" % (newddGvalue, join(tokens[1:], " "))
		
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
				if not quiet:
					colortext.error("No measure mapping for %d" % ID)
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
					if not quiet:
						colortext.error("No method mapping for %d" % ID)
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
		if not(record["ddG"] or record["ddG_H2O"]):
			missingFields.append("either_ddG")
		
				
		if not quiet:
			if not passed:
				if ID not in recordsWithUnresolvedMissingData:
					colortext.error("#%d: Could not fix record; Fields %s are missing" % (ID, str(missingFields)))
				
		return passed
		
	def getMutations(self, ID, record = None):
		'''Returns a list of ddgobjects.Mutation objects.'''
		record = self._getRecord(ID, record)
			
		mutations = []
			
		mutationline = record["MUTATION"].split()
		#print(len(mutationline))
		cline = join(mutationline, "")
		if ID in self.CysteineMutationCases: # Hack for input which includes information on the bonds generated from mutation to cysteine
			if not mutationline[1].isdigit():
				# todo: This will fail with insertion codes. I just want to see when they are added to make sure parsing works.
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
						if not residueID[-1].isalpha():
						#if (not residueID[:-1].isdigit()) or (not residueID[-1].isalpha()):
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
						# todo: This will fail with insertion codes. I just want to see when they are added to make sure parsing works.
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
						# todo: This will fail with insertion codes. I just want to see when they are added to make sure parsing works.
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
				return None
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
					raise Exception("An exception occurred parsing the mutation %s in record %d." % (cline, ID))
				 
		else:
			raise Exception("We need to add a case to handle this mutation string: '%s'." % cline)
		
		mutobjects = []
		for mutation in mutations:
			mutobjects.append(ddgobjects.Mutation(mutation["WildTypeAA"], mutation["ResidueID"], mutation["MutantAA"]))
			
		mutation_locations = []
		if (not badSecondaryStructure.get(ID)) and record["SEC.STR."]:
			mutation_locations = record["SEC.STR."].split(",")
			numlocations = len(mutation_locations)
			if not numlocations == len(mutobjects):
				raise Exception("The mutations '%s' do not have corresponding locations '%s' in record %d." % (str(mutobjects), mutation_locations, ID))
			for i in range(numlocations):
				n_location = self.secondary_structure_values.get(mutation_locations[i].strip())
				if not n_location:
					raise Exception("Bad secondary structure information '%s' for the mutation in record %d." % (record["SEC.STR."], ID))
				mutobjects[i].SecondaryStructurePosition = n_location 
		
		ASAs = []
		if ID == 1927:
			# Hack for this case
			ASAs = [None, 205.4]
		elif (not badASA.get(ID)) and record["ASA"]:
			ASAs = record["ASA"].split(",")
			numASAs = len(ASAs)
			if not numASAs == len(mutobjects):
				raise Exception("The mutations '%s' do not have corresponding ASAs '%s' in record %d." % (str(mutobjects), ASAs, ID))
			for i in range(numASAs):
				try:
					n_ASA = float(ASAs[i].strip())
				except :
					raise Exception("Bad ASA information '%s' for the mutation in record %d." % (record["ASA"], ID))
				mutobjects[i].AccessibleSurfaceArea = n_ASA 
		
		return mutobjects
	
	def getDDGH2OInKcal(self, ID, record = None, useRosettaConvention = False, getDDGH2OInstead = False):
		return self.getDDGInKcal(ID, record, useRosettaConvention, getDDGH2OInstead = True)
	
	def getDDGInKcal(self, ID, record = None, useRosettaConvention = False, getDDGH2OInstead = False):
		record = self._getRecord(ID, record)
		
		referenceID = self.getReference(ID, record)
		if not referenceID:
			raise Exception("Error processing reference: ID %d, %s" %  (ID, record["REFERENCE"]))
		dbReferencePK = "PMID:%s" % referenceID
		
		ddGKey = "ddG"
		if getDDGH2OInstead:
			ddGKey = "ddG_H2O"
		
		ddGline = record[ddGKey]
		#if getDDGH2OInstead and not ddGline:
			# todo: I should handle this better
		#	return None
		if ddGline.find("kJ/mol") == -1 and ddGline.find("kcal/mol") == -1:
			try:
				x = float(ddGline)
			except:
				colortext.error("Error processing %s: ID %d, %s" % (ddGKey, ID, record["ddG"]))
			if self.ddGUnitsUsed.get(dbReferencePK):
				unitsUsed = self.ddGUnitsUsed[dbReferencePK]
				assert(unitsUsed[-1] != '?') # todo: assertion after removing legacy code - this assertion can be removed on the next commit 
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
				raise Exception("Excepting parsing %s." % ddGKey)
			else:
				colortext.printf("No %s unit specified for new record: ID %d, %s, publication ID='%s'" % (ddGKey, ID, ddGline, dbReferencePK), "cyan")
				raise Exception("Excepting parsing %s." % ddGKey)
			mutations = self.getMutations(ID, record)
			mutations = join([join(map(str, m),"") for m in mutations], ",")
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
					refPMID = refPMID[0:-1] # todo: Hack for record 25600
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
			if not record.get(h):
				self.missingExpData[h] = self.missingExpData.get(h, 0) + 1
				record[h] = None
			fielddata = record[h]
			if record[h]:
				self.maxDBfieldlengths[h] = max(self.maxDBfieldlengths.get(h, 0), len(fielddata))
			if fielddata == "Unknown":
				fielddata = None
			if fielddata and h == "METHOD":
				assert(record['METHODS'])
				for i in range(len(record['METHODS'])):
					ExperimentalConditions["%s%d" % (exp2DBfield[h], i + 1)] = record['METHODS'][i]
				for i in range(len(record['METHODS']), 2):
					ExperimentalConditions["%s%d" % (exp2DBfield[h], i + 1)] = None
			elif fielddata and h == "MEASURE":
				assert(record['MEASURES'])
				for i in range(len(record['MEASURES'])):
					ExperimentalConditions["%s%d" % (exp2DBfield[h], i + 1)] = record['MEASURES'][i]
				for i in range(len(record['MEASURES']), 3):
					ExperimentalConditions["%s%d" % (exp2DBfield[h], i + 1)] = None
			elif fielddata and h in self.numericexpfields:
				try:
					fielddata = float(fielddata)
				except:
					if h == "T":
						# Values are stored in degrees Celsius
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
				colortext.error("Record %s not found" % str(recordnumber))
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
		if fields_to_ignore:
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
						if myrecord["ddG"] or theirrecord["ddG"] or myrecord["ddG_H2O"] or theirrecord["ddG_H2O"]:
							records_with_ddG_values[id] = True
						if not fieldname in myrecord.keys():
							colortext.error("Record %s is missing field %s in %s." % (str(id), fieldname, self.infilepath))
						elif not fieldname in theirrecord.keys():
							colortext.error("Record %s is missing field %s in %s." % (str(id), fieldname, secondDB.infilepath))
						elif myrecord[fieldname] != theirrecord[fieldname]:
							print(id, fieldname, myrecord[fieldname], theirrecord[fieldname])
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
					
		#existingIDs = getIDsInDB(self.ddGDB, source = "ProTherm-2008-09-08-23581")
		for i in added_info:
			if showall:
				colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
				#if i[0] not in existingIDs:
				#	colortext.printf("\tIgnored as record %d not in our database." % i[0], "cyan")
				if not records_with_ddG_values.get(i[0]):
					colortext.printf("\tIgnored as record %d has no ddG value." % i[0], "cyan")
				if not fields_of_interest.get(i[1]):
					colortext.printf("\tIgnored as field '%s' is not of interest." % i[1], "cyan")
			else:
				#if i[0] in existingIDs and
				if records_with_ddG_values.get(i[0]) and fields_of_interest.get(i[1]):
					colortext.message("Record %10.d: %s (%s). Value '%s' added." % (i[0], i[1], field_descriptions[i[1]], i[2]))			
		for i in deleted_info:
			if showall:
				colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
				#if i[0] not in existingIDs:
				#	colortext.printf("\tIgnored as record %d not in our database." % i[0], "cyan")
				if not records_with_ddG_values.get(i[0]):
					colortext.printf("\tIgnored as record %d has no ddG value." % i[0], "cyan")
				if not fields_of_interest.get(i[1]):
					colortext.printf("\tIgnored as field '%s' is not of interest." % i[1], "cyan")
			else:
				#if i[0] in existingIDs and
				if records_with_ddG_values.get(i[0]) and fields_of_interest.get(i[1]):
					colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
		for i in changed_info:
			if showall:
				colortext.error("Record %10.d: %s (%s). Value '%s' deleted." % (i[0], i[1], field_descriptions[i[1]], i[2]))
				#if i[0] not in existingIDs:
				#	colortext.printf("\tIgnored as record %d not in our database." % i[0], "cyan")
				if not records_with_ddG_values.get(i[0]):
					colortext.printf("\tIgnored as record %d has no ddG value." % i[0], "cyan")
				if not fields_of_interest.get(i[1]):
					colortext.printf("\tIgnored as field '%s' is not of interest." % i[1], "cyan")
			else:
				#if i[0] in existingIDs and 
				if records_with_ddG_values.get(i[0]) and fields_of_interest.get(i[1]):
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
			ptReader = ProThermReader(os.path.join("..", "rawdata", "ProTherm", "ProTherm25844.dat"), quiet = True)
			record = ptReader.readRecord(ID)
			if not record:
				colortext.error("Could not read a record with number #%s" % str(ID))
			ptReader.printRecord(ID)
		else:
			if args[1] == "diff":
				flregex = re.compile("^(ProTherm)(\d+)[.]dat$", re.IGNORECASE)
				datapath = os.path.join("..", "rawdata", "ProTherm")
				dbs = []
				for filenm in sorted(os.listdir(datapath)):
					mtchs = re.match("^(ProTherm)(\d+)[.]dat$", filenm, re.IGNORECASE)
					if mtchs:
						#sqlfilename = "%s%s.sql" % (mtchs.group(1),mtchs.group(2))
						#sqlfilepath = os.path.join(datapath, sqlfilename)
						dbs.append((int(mtchs.group(2)), ProThermReader(os.path.join(datapath, filenm))))
				dbs = sorted(dbs)
				#dbs[-2][1].diff(dbs[-1][1], fields_to_ignore = ["E.C.NUMBER", "ION_NAME_1", "SWISSPROT_ID", "REMARKS", "NO_MOLECULE"], start_at_id = 0, end_at_id = None, showall = False)		
				dbs[-2][1].diff(dbs[-1][1], fields_to_ignore = [], start_at_id = 0, end_at_id = None, showall = False)		
			if args[1] == "test":
				ptReader = ProThermReader(os.path.join("..", "rawdata", "ProTherm25844.dat"), quiet = True)
				ptReader.test()