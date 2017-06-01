#!/macqiime/bin/python

import argparse
import re
import copy
import os

# V2017-05-29
# This script takes output from BatchCoreAnalysis and see which OTUs are unique (on top of being core)


#################################

parser = argparse.ArgumentParser(
	description="Takes output from BatchCoreAnalysis and sees which core OTUs are unique to each treatment")
parser.add_argument(
	'-c',
	'--cores',
	help = 'Comma separated list of core files-- output from BatchCoreAnalysis',
	type = str,
	required = True)
parser.add_argument(
	'-o',
	'--outputFolder',
	help = 'Output folder [Default: Unique_OTUs_Across_Groups]',
	required = False,
	default = 'Unique_OTUs_Across_Groups')
	
args = parser.parse_args()

coresFP = args.cores
coresList = coresFP.split(",")

outputFolder = args.outputFolder

os.system("mkdir " + outputFolder)

#################################

def printTableFromDictionary(dictionary, sampleNames, output):
	toPrint = ''
	for i in sampleNames:
		toPrint += '\t' + i
	toPrint += '\n'
	for row in dictionary:
		toPrint += str(row)
		for element in dictionary[row]:
			toPrint += "\t" + str(element)
		toPrint += "\n"
	open(str(output)+".txt", 'w').write(toPrint)
	print "DONE"
	
#################################

cores = {}
sampleNames = {}
for file in coresList: # makes 2-layer dict; first layer is cores and second is OTUs, third is abundance
	nameTemp = file.replace(".txt","")
	name = re.sub("^.*/","", nameTemp)
	cores[name] = {}
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
			cores[name][newLine[0]] = newLine[1:]
	
uniqueOTUs = {}
for sample in cores:
	uniqueOTUs[sample] = []
	outGroupSamples = cores.keys()
	outGroupSamples.remove(sample) # Get rid of outgroup
	# Now get all OTUs unique to ingroup and outgroup
	ingroup = cores[sample].keys()
	outgroup = []
	for s in outGroupSamples:
		for OTU in cores[s]:
			outgroup.append(OTU)
	outgroup = list(set(outgroup))
	uniqueOTUs[sample] = list(set(ingroup) - set(outgroup))
	
coresFilt = copy.deepcopy(cores)
for sample in cores:
	for OTU in cores[sample]:
		if OTU not in uniqueOTUs[sample]:
			del coresFilt[sample][OTU]
		
for i in coresFilt:
	printTableFromDictionary(coresFilt[i], sampleNames[i], outputFolder + '/' + i)