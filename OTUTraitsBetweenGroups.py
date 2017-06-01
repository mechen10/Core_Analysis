#!/macqiime/bin/python

# V2017-05-31

# This script takes output from Enrichment Core Analysis and compares how many are shared between treatments
# Rows: all OTUs
# Columns: each treatment + taxonomy at end
# You pick which trait you want to include in table; originally used for fold-change

import argparse
import re
import os

#################################

parser = argparse.ArgumentParser(
	description="Takes multiple trait tables with OTUs as first column, taxonomy as last column, and traits in middle columns and compares one trait across all tables.")
parser.add_argument(
	'-i',
	'--inputOTUs',
	help = 'Comma separated list of trait table (eg output from EnrichmentCoreAnalysis.py)',
	type = str,
	required = True)
parser.add_argument(
	'-c',
	'--column',
	help = 'Comma separated list of metadata columns that has traits to compare. If not included, will do all columns. [Default: NONE]',
	type = str,
	required = False,
	default = 'NONE')
parser.add_argument(
	'-o',
	'--outputFolder',
	help = 'Output folder [Default: Comparison_across_treatments]',
	required = False,
	default = 'Comparison_across_treatments')
	
args = parser.parse_args()

inputFP = args.inputOTUs
inputList = inputFP.split(",")

columnraw = args.column
columnList = columnraw.split(",")

outputFolder = args.outputFolder
os.system('mkdir ' + outputFolder)

#################################
# #TESTING
# 
# inputList = ['/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/WATERCORE_Enrichment_of_unique_OTUs/NereotestExNWater_test.txt','/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/WATERCORE_Enrichment_of_unique_OTUs/NereotestMastWater.txt','/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/WATERCORE_Enrichment_of_unique_OTUs/NereotestNereoMastWater.txt','/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/WATERCORE_Enrichment_of_unique_OTUs/NereotestNereoWater.txt']
# columnList = ['FoldChangeAve']
#################################


def loadTable(tableFP): # Loads and converts to relative traitance
	tableOpen = open(tableFP, 'U') # This is in Unix format
	tableTemp = []
	for i in tableOpen:
		tempLine = i.strip()
		tempLineSplit = tempLine.split("\t")
		tableTemp.append(tempLineSplit)
	tableOpen.close()
	# Get rid of last column
	OTUTable = {} # Now make dictionary
	taxaIDs = {} # Make taxonomy reference 
	headers = []
	for lineN in range(len(tableTemp)):
		if lineN == 0: # This is the header line
			headers = tableTemp[lineN]
		else:
			OTUTable[str(tableTemp[lineN][0])] = {}
			taxaIDs[str(tableTemp[lineN][0])] = tableTemp[lineN][len(tableTemp[lineN])-1]
			for trait in range(len(tableTemp[lineN][1:])-1):
				OTUTable[str(tableTemp[lineN][0])][headers[trait]] = tableTemp[lineN][1:][trait]
	return OTUTable,taxaIDs # Output is a 2-layer dictionary; first is OTU IDs and second is samples. Also, one-layer dict with taxa IDs

def printTableFromDictionary(dictionary, sampleNames, output):
	toPrint = ''
	for i in sampleNames:
		toPrint += '\t' + i
	toPrint += '\n'
	for row in dictionary:
		toPrint += str(row)
		for column in sampleNames:
			toPrint += "\t" + str(dictionary[row][column])
		toPrint += "\n"
	open(str(output)+".txt", 'w').write(toPrint)
	print "DONE"
	
#################################

allTables = {}
allTaxaID = {}
for file in inputList: # makes 3-layer dict; first layer is treatment and second is OTUs, third is treatment (value is abundance)
	nameTemp = file.replace(".txt","")
	name = re.sub("^.*/","", nameTemp)
	allTables[name], allTaxaID[name] = loadTable(file)
	
# Get all OTUs
allOTUs = []
for treatment in allTaxaID:
	allOTUs += allTaxaID[treatment].keys()
allOTUs = list(set(allOTUs))
	
# Now, iterate through Traits, OTUs and through samples. Then add taxonomy on at end
traitsTable = {}
for trait in columnList:
	traitsTable[trait] = {}
	for OTU in allOTUs:
		traitsTable[trait][OTU] = {}
		for treatment in allTables:
			if OTU in allTables[treatment].keys():
				print traitsTable
				print allTables[treatment][OTU]
				traitsTable[trait][OTU][treatment] = allTables[treatment][OTU][trait]
				traitsTable[trait][OTU]['taxonomy'] = allTaxaID[treatment][OTU]
			else:
				traitsTable[trait][OTU][treatment] = 'NA'
			
# Now, print the results.
# Get sample name order
allheaders = []
for trait in [traitsTable.keys()[0]]:
	for OTU in [traitsTable[trait].keys()[0]]:
		allheaders = traitsTable[trait][OTU].keys()
allheaders.remove('taxonomy')
allheaders.append('taxonomy')

for trait in traitsTable:
	printTableFromDictionary(traitsTable[trait], allheaders, outputFolder + "/" + trait)


