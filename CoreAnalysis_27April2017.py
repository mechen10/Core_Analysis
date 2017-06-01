#!/bin/bash/python

# Built under QIIME 1.9.0 and python 2.7.3

# Script for determining core microbiomes in different treatment groups

# So if you have multiple treatment groups with envrinmental samples, there are a few 'groups' of OTUs that will exist:
	# 1. OTUs that are everywhere-- including in the outgroup/environmental samples
	# 2. OTUs that are only on treatment groups (eg. surface-specific microbes across all latex shapes)
	# 3. OTUs that are only on specific treatment groups (eg. specific to finely branched/bladed, etc)
# The metadata should have an outgroup column and a treatment type column. They can be called anything you want.
	# HOWEVER, the treatment type column should have 'x' for outgroup members.
# Thus, you want to parse out all of these options.

import argparse
import os
import sys
import copy


#################################

parser = argparse.ArgumentParser(
	description="Bins and classifies OTUs according to salinity specialization")
parser.add_argument(
	'-i',
	'--OTU_table',
	help = "OTU table from QIIME containing all samples",
	required = True,)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'Mapping file for OTU table',
	required = True)
parser.add_argument(
	'-t',
	'--treatment_groups',
	help = 'Treatment groups; must match headers in metadata file. Can have multiple, but must be comma-separated',
	type = str,
	required = True)
parser.add_argument(
	'--outgroup',
	help = 'Name of header indicated treatment groups vs outgroup-- not required. If not included, default is to not calculate outgroup analysis',
	required = False,
	default = 'False')
parser.add_argument(
	'--outgroup_remove',
	help = 'Name of type to remove in header-- not required. If not included, default is to not calculate outgroup analysis',
	required = False,
	default = 'False')
parser.add_argument(
	'-T',
	'--threshold',
	help = 'Percent threshold for OTU to be observed within group to be considered "core" [Default: 0.9]',
	required = False,
	default = .90)
	
	
args = parser.parse_args()

OTUFP = args.OTU_table
metadataFP = args.metadata
treatment_groups = args.treatment_groups.split(",")
outgroup = args.outgroup
outgroupRemove = args.outgroup_remove
threshold = args.threshold


#################################
# FUNCTIONS

# LOAD DATA

def loadOTUTable(OTUFP): # Loads and converts to relative abundance
	# Make OTU Table
	os.system('biom convert -i ' + OTUFP + ' --to-tsv --header-key taxonomy -o OTU_Table_text.txt')
	# Load OTU Table
	biomOpen = open("OTU_Table_text.txt", 'r') # This is in Unix format
	biomTemp = []
	for i in biomOpen:
		tempLine = i.strip()
		tempLineSplit = tempLine.split("\t")
		biomTemp.append(tempLineSplit)
	biomOpen.close()
	# Get rid of last column
	OTUTable = {} # Now make dictionary
	taxaIDs = {} # Make taxonomy reference 
	for lineN in range(len(biomTemp)):
		if lineN == 0: # This is first line; skip this
			pass
		elif lineN == 1: # This is the header line
			headers = biomTemp[lineN][1:len(biomTemp[lineN])-1]
		else:
			OTUTable[str(biomTemp[lineN][0])] = {}
			taxaIDs[str(biomTemp[lineN][0])] = biomTemp[lineN][len(biomTemp[lineN])-1]
			for abund in range(len(biomTemp[lineN][1:])-1):
				OTUTable[str(biomTemp[lineN][0])][headers[abund]] = biomTemp[lineN][1:][abund]
	# Get total counts for each site
	totalCounts = {}
	for h in range(len(headers)):
		totalCounts[headers[h]] = 0
	for OTU in OTUTable:
		for h in range(len(headers)):
			totalCounts[headers[h]] += float(OTUTable[OTU][str(headers[h])])
	# Convert to relative abundance
	for OTU in OTUTable.keys():
		tempOTUlist = OTUTable[OTU].copy()
		for sites in OTUTable[OTU].keys():
			tempOTUlist[sites] = float(OTUTable[OTU][sites])/float(totalCounts[sites])
		OTUTable[OTU] = tempOTUlist
	os.system("rm OTU_Table_text.txt")
	return OTUTable,taxaIDs # Output is a 2-layer dictionary; first is OTU IDs and second is samples. Also, one-layer dict with taxa IDs
	
def loadMetadata(metadataFP):
	metadataOpen = open(metadataFP, 'U') # U is for 'Universal read'-- automatically turns into Unix LF
	metadataTemp = []
	for i in metadataOpen:
		lineTemp = i.strip()
		lineTempSplit = lineTemp.split("\t")
		metadataTemp.append(lineTempSplit)
	metadataOpen.close()
	metadata = {}
	metadataSites = []
	for lineN in range(len(metadataTemp)):
		if lineN == 0:
			headerList = metadataTemp[lineN]
			for headerName in metadataTemp[lineN]:
				metadata[headerName] = {}
		else:
			for i in range(1,len(metadataTemp[lineN])):
				metadataSites.append(metadataTemp[lineN][0])
				sortHeader = headerList[i]
				metadata[sortHeader][metadataTemp[lineN][0]] = metadataTemp[lineN][i]
	return metadata # output is 2-layer dictionary: first is Metadata and second is samples
	
def getCore(OTUTable, metadata, Treatment, outgroupRemove):
	global threshold
	# Get Treatment
	outgroupMETA = metadata[Treatment]
	# Get treatment types
	allGroups = []
	for i in outgroupMETA.keys():
		allGroups.append(outgroupMETA[i])
	treatmentGroups = set(allGroups)
	treatmentGroups.remove(outgroupRemove)
	OVERALLCORE = {}
	# Loop through groups
	for treatment in treatmentGroups:
		samplesToDiscard = []
		for i in outgroupMETA.keys(): # Get Treatment members and save as sampleToKeep
			if outgroupMETA[i] != treatment:
				samplesToDiscard.append(i)
		# Copy and filter OTUTable
		outgroupOTUTable = copy.deepcopy(OTUTable)
		for sample in samplesToDiscard:
			for OTU in outgroupOTUTable.keys():
				del outgroupOTUTable[OTU][sample]
		# Then, calculate OVERALL core
		totalSamples = len(outgroupOTUTable[outgroupOTUTable.keys()[0]])
		OVERALLSUBCORE = []
		for OTU in outgroupOTUTable.keys():
			appearances = 0
			for sample in outgroupOTUTable[OTU]:
				if outgroupOTUTable[OTU][sample] != 0:
					appearances += 1
			if appearances/float(totalSamples) > float(threshold):
				OVERALLSUBCORE.append(OTU)
		OVERALLCORE[treatment] = OVERALLSUBCORE
	return OVERALLCORE

def commonCore(dict): # takes dictionary output from getCore and finds common elements between all elements
	# Get overall cores
	treatmentcores = []
	for i in dict:
		for j in dict[i]:
			for k in dict[i][j]:
				treatmentcores.append(k)
	treatmentcores = set(treatmentcores)
	treatmentcoresFINAL = []
	for OTU in treatmentcores:
		commonOTUTF = []
		for i in dict:
			for j in dict[i]:
				commonOTUTF.append(OTU in dict[i][j])
		if all(commonOTUTF):
			treatmentcoresFINAL.append(OTU)
	return treatmentcoresFINAL

def getTaxaIDs(OTUID):
	global taxaIDs
	taxonomyList = {}
	for i in OTUID:
		taxonomyList[i] = taxaIDs[i]
	return taxonomyList
	
	
	
# def getCore(OTUTable, metadata, Treatment, outgroupRemove):
# 	global threshold
# 	# Get Treatment
# 	outgroupMETA = metadata[Treatment]
# 	# Get treatment types
# 	treatmentGroups = set(outgroupMETA)
# 	treatmentGroups.remove(outgroupRemove)
# 	
# 	# Loop through groups
# 	samplesToDiscard = []
# 	for i in outgroupMETA.keys(): # Get Treatment members and save as sampleToKeep
# 		if outgroupMETA[i] == outgroupRemove:
# 			samplesToDiscard.append(i)
# 	# Copy and filter OTUTable
# 	outgroupOTUTable = OTUTable.copy()
# 	for sample in samplesToDiscard:
# 		for OTU in outgroupOTUTable.keys():
# 			del outgroupOTUTable[OTU][sample]
# 	# Then, calculate OVERALL core
# 	totalSamples = len(outgroupOTUTable[outgroupOTUTable.keys()[0]])
# 	
# 
# 	OVERALLCORE = []
# 	for OTU in outgroupOTUTable.keys():
# 		appearances = 0
# 		for sample in outgroupOTUTable[OTU]:
# 			if outgroupOTUTable[OTU][sample] != 0:
# 				appearances += 1
# 		print appearances/float(totalSamples)
# 		if appearances/float(totalSamples) >= threshold:
# 			OVERALLCORE.append(OTU)
# 	return OVERALLCORE
# 
# def getTaxaIDs(OTUID):
# 	global taxaIDs
# 	taxonomyList = {}
# 	for i in OTUID:
# 		taxonomyList[i] = taxaIDs[i]
# 	return taxonomyList

#################################

os.system("mkdir COREOTUs")

OTUTable,taxaIDs = loadOTUTable(OTUFP)
print "DONE LOADING OTU TABLE"
metadata = loadMetadata(metadataFP)
print "DONE LOADING METADATA"

# First, calculate OVERALL core
OVERALLCORE = getCore(OTUTable, metadata, outgroup, outgroupRemove)
overallCore = OVERALLCORE[OVERALLCORE.keys()[0]]

# Then, get secondary core
cores = {}
for i in treatment_groups:
	cores[i] = getCore(OTUTable, metadata, i, 'X')

treatmentcores = commonCore(cores)


# # Remove overallcore from treatmentcore
# uniqueTreatmentCores = [x for x in treatmentcores if not x in overallCore]

# # Remove overallCore from individual cores
# uniqueCores = {}
# for i in cores:
# 	uniqueCores[i] = {}
# 	for j in cores[i]:
# 		uniqueCores[i][j] = [x for x in cores[i][j] if not x in overallCore]
# 	

# Don't want to get rid of un-unique cores; keep these
uniqueTreatmentCores = copy.deepcopy(treatmentcores)
uniqueCores = copy.deepcopy(cores)
	
#################################
# Convert to taxonomy and print out as list
overallCoreWrite = ''
for i in overallCore:
	overallCoreWrite += taxaIDs[i] + "\n"
open("./COREOTUs/OverallCore.txt", 'w').write(overallCoreWrite)

TreatmentCoresWrite = ''
for i in uniqueTreatmentCores:
	TreatmentCoresWrite += taxaIDs[i] + "\n"
open("./COREOTUs/TreatmentCores.txt", 'w').write(TreatmentCoresWrite)

for i in uniqueCores:
	for j in uniqueCores[i]:
		tempToWrite = ''
		for k in uniqueCores[i][j]:
			tempToWrite += taxaIDs[k] + "\n"
		open(str("./COREOTUs/" +i + "_" + j + ".txt"), 'w').write(tempToWrite)
		
	
#################################
		
toWrite = "OverallCore"
for i in overallCore:
	toWrite += "\t" + i
toWrite += "\n"

toWrite = "AllTreatmentsCore"
for i in uniqueTreatmentCores:
	toWrite += "\t" + i
toWrite += "\n"

for i in uniqueCores:
	for j in uniqueCores[i]:
		toWrite += i + ":" + j
		for k in uniqueCores[i][j]:
			toWrite += "\t" + k
		toWrite += "\n"

coreOTUs = open("./COREOTUs/CORE_OTUS.txt", 'w')	
coreOTUs.write(toWrite)
coreOTUs.close()

#################################
# Print relative abundance OTU table

first = True
toWrite = "#OTU ID"
for row in OTUTable:
	if first:
		for sample in OTUTable[row].keys():
			toWrite += "\t" + sample
		first = False
		toWrite += "\n"
	toWrite += row
	for abund in OTUTable[row]:
		toWrite += "\t" + str(OTUTable[row][abund])
	toWrite += "\n"

open("./COREOTUs/OTU_Table_relabund.txt", 'w').write(toWrite)

#################################
# Print taxaIDLegend

toWrite = "OTUID\ttaxonomy\n"
for i in taxaIDs.keys():
	toWrite += i + "\t" + str(taxaIDs[i]) + "\n"
open("./COREOTUs/taxaIDLegend.txt", 'w').write(toWrite)

