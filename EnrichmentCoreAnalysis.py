#!/macqiime/bin/python

# V2017-05-29
# This takes the output from UniqueOTUsAcrossGroups or Batch Core Analysis and compares to a control to see how much enrichment there is

import argparse
import os
import re
import math as m
import numpy as np


#################################

parser = argparse.ArgumentParser(
	description="Takes output from UniqueOTUsAcrossGroups and compares to a control to find enrichment of OTUs")
parser.add_argument(
	'-i',
	'--inputOTUs',
	help = 'Comma separated list of uniqueOTU or coreOTU files-- output from UniqueOTUsAcrossGroups.py or BatchCoreAnalysis.py',
	type = str,
	required = True)
parser.add_argument(
	'-t',
	'--OTUTable',
	help = 'OTUTable, biom format',
	type = str,
	required = True)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'Metadata',
	type = str,
	required = True)
parser.add_argument(
	'-c',
	'--control',
	help = 'Group:Treatment for which to be treated as control comparison',
	type = str,
	required = True)
parser.add_argument(
	'-r',
	'--replicates',
	help = 'Metadata column that has replicates to compare. If not included, will do overall calculation [Default: NONE]',
	type = str,
	required = False,
	default = 'NONE')
parser.add_argument(
	'-o',
	'--outputFolder',
	help = 'Output folder [Default: Enrichment_of_unique_OTUs]',
	required = False,
	default = 'Enrichment_of_unique_OTUs')
	
args = parser.parse_args()

inputFP = args.inputOTUs
inputList = inputFP.split(",")
OTUFP = args.OTUTable
metadataFP = args.metadata
control = args.control
controlSplit = control.split(":")
rep = args.replicates

outputFolder = args.outputFolder

os.system("mkdir " + outputFolder)

#################################

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

#################################


OTUTable,taxaID = loadOTUTable(OTUFP)
metadata = loadMetadata(metadataFP)


#################################
# Calculate Enrichment

if rep != 'NONE':

	# Make uniqueOTUs file
	uniqueOTUs = {} # Final dictionary
	sampleNames = {} # For holding sample names
	for file in inputList: # makes 3-layer dict; first layer is treatment and second is OTUs, third is treatment (value is abundance)
		nameTemp = file.replace(".txt","")
		name = re.sub("^.*/","", nameTemp)
		uniqueOTUs[name] = {}
		sampleNames[name] = {}
		tempFile = open(file, 'U')
		first = True
		for i in tempFile:
			newLine = i.strip()
			newLine = newLine.split("\t")
			if first:
				sampleNames[name] = newLine[0:]
				first = False
			else:
				uniqueOTUs[name][newLine[0]] = {}
				for j in range(len(sampleNames[name])):
					replicateN = [metadata[rep][s] for s in metadata[rep] if sampleNames[name][j] == s]
					uniqueOTUs[name][newLine[0]][replicateN[0]] = float(newLine[j+1])
		file.close()
	# Make list of control sample names
	allConSamples = {}
	for sample in metadata[controlSplit[0]]: # Get all control sample names
		if metadata[controlSplit[0]][sample] == controlSplit[1]:
			allConSamples[metadata[rep][sample]] = sample

	# Make file with both abundance values
	OTUCompare = {} # dictionary with only rel.abundance values
	for treatment in uniqueOTUs:
		otuList = uniqueOTUs[treatment].keys()
		OTUCompare[treatment] = {}
		for OTU in otuList:
			OTUCompare[treatment][OTU] = {}
			abundCon = {}
			for nCon in allConSamples:
				abundCon[nCon] = float(OTUTable[OTU][allConSamples[nCon]])
			OTUCompare[treatment][OTU]['Control'] = abundCon
			OTUCompare[treatment][OTU]['Treatment'] = uniqueOTUs[treatment][OTU] 
			# output is dictionary (1) treatment (2) OTU (3) control or treatment (4) abundance

	# Make sure all have same number of replicates
	for treatment in OTUCompare:
		for OTU in OTUCompare[treatment]:
			for replicate in OTUCompare[treatment][OTU]['Treatment']:
				if replicate in OTUCompare[treatment][OTU]['Control']:
					pass
				else:
					del OTUCompare[treatment][OTU]['Treatment'][replicate]


	# Compare means of control and treatment for each OTU
	for treatment in OTUCompare:
		for OTU in OTUCompare[treatment]:
			OTUCompare[treatment][OTU]['FoldChangeList'] = []
			OTUCompare[treatment][OTU]['FoldChangeAve'] = []
			OTUCompare[treatment][OTU]['sd'] = []
			for replicate in OTUCompare[treatment][OTU]['Treatment']:
				treatVal = OTUCompare[treatment][OTU]['Treatment'][replicate]
				conVal = OTUCompare[treatment][OTU]['Control'][replicate]
				if conVal == 0 and treatVal > 0:
					OTUCompare[treatment][OTU]['FoldChangeList'].append(float("inf"))
				elif treatVal == 0:
					OTUCompare[treatment][OTU]['FoldChangeList'].append(float("-inf"))
				else:
					tempVal = m.log(treatVal/conVal,2)
					if tempVal < 1e-10:
						tempVal = 0
					OTUCompare[treatment][OTU]['FoldChangeList'].append(tempVal)
			OTUCompare[treatment][OTU]['FoldChangeAve'] = np.mean(OTUCompare[treatment][OTU]['FoldChangeList'])
			OTUCompare[treatment][OTU]['sd'] = np.std(OTUCompare[treatment][OTU]['FoldChangeList'])

else:

	# Make uniqueOTUs file
	uniqueOTUs = {} # Final dictionary
	sampleNames = {} # For holding sample names
	for file in inputList: # makes 3-layer dict; first layer is treatment and second is OTUs, third is treatment (value is abundance)
		nameTemp = file.replace(".txt","")
		name = re.sub("^.*/","", nameTemp)
		uniqueOTUs[name] = {}
		sampleNames[name] = {}
		tempFile = open(file, 'U')
		first = True
		for i in tempFile:
			newLine = i.strip()
			newLine = newLine.split("\t")
			if first:
				sampleNames[name] = newLine[0:]
				first = False
			else:
				uniqueOTUs[name][newLine[0]] = []
				for j in newLine[1:]:
					uniqueOTUs[name][newLine[0]].append(float(j))

	# Make list of control sample names
	allConSamples = []
	for sample in metadata[controlSplit[0]]: # Get all control sample names
		if metadata[controlSplit[0]][sample] == controlSplit[1]:
			allConSamples.append(sample)

	# Make file with both abundance values
	OTUCompare = {} # dictionary with only rel.abundance values
	controlDone = False
	for treatment in uniqueOTUs:
		otuList = uniqueOTUs[treatment].keys()
		OTUCompare[treatment] = {}
		for OTU in otuList:
			OTUCompare[treatment][OTU] = {}
			abundCon = []
			for nCon in allConSamples:
				abundCon.append(float(OTUTable[OTU][nCon]))
			OTUCompare[treatment][OTU]['Control'] = abundCon
			OTUCompare[treatment][OTU]['Treatment'] = uniqueOTUs[treatment][OTU] 
			# output is dictionary (1) treatment (2) OTU (3) control or treatment (4) abundance

	# Compare means of control and treatment for each OTU
	for treatment in OTUCompare:
		for OTU in OTUCompare[treatment]:
			treatVal = np.mean(OTUCompare[treatment][OTU]['Treatment'])
			conVal = np.mean(OTUCompare[treatment][OTU]['Control'])
			if treatVal == 0:
				OTUCompare[treatment][OTU]['FoldChangeAve'] = float("-inf")
			elif conVal == 0:
				OTUCompare[treatment][OTU]['FoldChangeAve'] = float("inf")
			else:
				OTUCompare[treatment][OTU]['FoldChangeAve'] = m.log(treatVal/conVal,2)
#################################
# Add taxa IDs

for treatment in OTUCompare:
	for OTU in OTUCompare[treatment]:
		OTUCompare[treatment][OTU]['taxonomy'] = taxaID[OTU]
				

#################################
# Print results

if rep != 'NONE':
	orderheadings = ['Control','Treatment','FoldChangeList','FoldChangeAve','sd','taxonomy']
	for treatment in OTUCompare:
		toWrite = ''
		for title in orderheadings:
			toWrite += "\t" + title
		toWrite += "\n"
		for OTU in OTUCompare[treatment]:
			toWrite += OTU
			for title in orderheadings:
				if isinstance(OTUCompare[treatment][OTU][title],dict):
					toWrite += "\t" + str(OTUCompare[treatment][OTU][title].values())
				else:
					toWrite += "\t" + str(OTUCompare[treatment][OTU][title])
			toWrite += "\n"
		open(outputFolder + "/" + treatment + ".txt", 'w').write(toWrite)
		
else:
	orderheadings = ['Control','Treatment','FoldChangeAve','taxonomy']
	for treatment in OTUCompare:
		toWrite = ''
		for title in orderheadings:
			toWrite += "\t" + title
		toWrite += "\n"
		for OTU in OTUCompare[treatment]:
			toWrite += OTU
			for title in orderheadings:
				if isinstance(OTUCompare[treatment][OTU][title],dict):
					toWrite += "\t" + str(OTUCompare[treatment][OTU][title].values())
				else:
					toWrite += "\t" + str(OTUCompare[treatment][OTU][title])
			toWrite += "\n"
		open(outputFolder + "/" + treatment + ".txt", 'w').write(toWrite)

