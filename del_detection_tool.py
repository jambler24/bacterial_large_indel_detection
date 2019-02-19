# Importing functions
from subprocess import call, Popen, PIPE

import argparse

import numpy

import sys

import shutil

import csv

import copy

import os

import itertools

import networkx as nx

from itertools import islice

import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import pylab as pl
import math

from scipy.optimize import curve_fit


def input_parser(file_path, parse_as='default'):
	if file_path[-3:] == ".fa" or file_path[-6:] == ".fasta":
		input_file = open(file_path, "r")
		output_list = []
		# set variables
		sequence_details = ""
		sequence = ""

		for line in input_file:
			if line[0] == ">":
				if len(sequence_details) > 0:
					sequence_details = sequence_details.strip()
					sequence = sequence.strip()
					sequence = sequence.replace("\n","")
					gene_ID_dict = {"gene_details" : sequence_details[1:], "DNA_seq" : sequence}
					output_list.append(gene_ID_dict)
					sequence_details = ""
					sequence = ""
				sequence_details = line
			else:
				sequence += line
		
		sequence_details = sequence_details.strip()
		
		sequence = sequence.strip()
		sequence = sequence.replace("\n","")
		sequence = sequence.replace(" ","")
		
		gene_ID_dict = {"gene_details" : sequence_details[1:], "DNA_seq" : sequence}

		output_list.append(gene_ID_dict)

		return output_list

	if file_path[-4:] == ".bim":
		input_file = open(file_path, "r")
		return input_file
		
	if file_path[-4:] == ".csv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=',')
		data_matrix = list(data_table)
		result=numpy.array(data_matrix)
		return result

	if file_path[-4:] == ".ssv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=';')
		data_matrix = list(data_table)
		result=numpy.array(data_matrix)
		return result
		
	if file_path[-4:] == ".txt":
		list_of_dicts  = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts
		
	if file_path[-4:] == ".vcf":
		list_of_dicts  = []
		# Deal with random info at start
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if line[0:2] != "##" and line[0] == "#":
				vcf_headder_line = line.split('\t')
				vcf_headder_line[0] = vcf_headder_line[0][1:]
				vcf_headder_line[-1] = vcf_headder_line[-1].strip()

			if not line.startswith('#'):
				entries = line.split('\t')
				entry_dict = {vcf_headder_line[0]:entries[0], vcf_headder_line[1]:entries[1], vcf_headder_line[2]:entries[2], vcf_headder_line[3]:entries[3], vcf_headder_line[4]:entries[4], vcf_headder_line[5]:entries[5], vcf_headder_line[6]:entries[6], vcf_headder_line[7]:entries[7], vcf_headder_line[8]:entries[8], vcf_headder_line[9]:entries[9].strip(), 'ORIGIN':entry_label} 

				list_of_dicts.append(entry_dict)
		return list_of_dicts
	
	if file_path[-5:] == ".diff":
		list_of_dicts  = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts

	if file_path[-5:] == "kbone":
		backb_listOlists = []
		in_file = open(file_path, 'r')
		for line in in_file:
			curr_list = []
			line = line.strip()
			curr_list = line.split('\t')
			backb_listOlists.append(curr_list)
		return backb_listOlists

	if file_path[-4:] == ".gtf":
		list_of_lists = []
		in_file = open(file_path, 'r')
		for line in in_file:
			entries = line.split('\t')
			entries[8] = entries[8].strip('\n')
			#entries_extra_info = entries[8].split(';')

			if '; ' in entries[8]:
				entries_extra_info = entries[8].split('; ')
			else:
				entries_extra_info = entries[8].split(';')

			extra_info_dict = {}

			for info_byte in entries_extra_info:

				if info_byte != '#' and len(info_byte) > 1:

					info_byte = info_byte.split(' ')

					extra_info_dict[info_byte[0]] = info_byte[1]

					entries[8] = extra_info_dict

			list_of_lists.append(entries)


		return list_of_lists

	if file_path[-5:] == ".gff3" or file_path[-4:] == ".gff":
		list_of_lists  = []
		in_file = open(file_path, 'r')
		entry_label = file_path
		
		if parse_as == 'default':
			for line in in_file:
				if not line.startswith('#'):
					entries = line.split('\t')
					if len(entries) > 5:
						entries[8] = entries[8].strip('\n')

						entries_extra_info = entries[8].split(';')

						extra_info_dict = {}

						for info_byte in entries_extra_info:
							info_byte = info_byte.split('=')
							extra_info_dict[info_byte[0]] = info_byte[1]

						entries[8] = extra_info_dict
				
						list_of_lists.append(entries)
		if parse_as == 'gtf':
			# Reformatting as gtf
			for line in in_file:
				if not line.startswith('#'):
					entries = line.split('\t')
					if len(entries) > 5:
						entries[8] = entries[8].strip('\n')

						entries_extra_info = entries[8].split(';')

						extra_info_dict = {}
						gtf_info_dict = {}

						for info_byte in entries_extra_info:
							info_byte = info_byte.split('=')
							extra_info_dict[info_byte[0]] = info_byte[1]

						if 'locus_tag' in extra_info_dict.keys():

							gtf_info_dict['gene_id'] = extra_info_dict['locus_tag']

							entries[8] = gtf_info_dict
					
							list_of_lists.append(entries)


		return list_of_lists

	if file_path[-4:] == ".gff OLD":
		list_of_dicts  = []
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if not line.startswith('#'):
				entries = line.split('\t')
				entries[8] = entries[8].strip('\n')
				entries_extra_info = entries[8].split(';')

				if entries[2] == 'gene':
					
					NOTE = ''
					for add_info in entries_extra_info:
						if 'locus_tag' in add_info:
							LOCUS = add_info[10:]
						if 'Name'in add_info:
							SYMBOL = add_info[5:]
						if 'Note'in add_info:
							NOTE = add_info[5:]



						#row['LOCUS'] row['SYMBOL'] row['SYNOYM'] row['LENGTH']  row['START'] row['STOP']  row['STRAND']  row['NAME']  row['CHROMOSOME']  row['GENOME ONTOLOGY']  row['ENZYME CODE']  row['KEGG']  row['PATHWAY']  row['REACTION'] row['COG'] row['PFAM']  row['OPERON'] 

					entry_dict = {'CHROM':entries[0], 'LOCUS':LOCUS, 'START':entries[3], 'STOP':entries[4], 'STRAND':entries[6], 'SYMBOL':SYMBOL, 'INFO':entries_extra_info, 'NOTE':NOTE}
					list_of_dicts.append(entry_dict)
		return list_of_dicts


def getPerBaseCoverage(bamFilePath):

	bwa_build_call = [bedtoolsPath, 'genomecov', '-ibam', bamFilePath, '-d']

	outfileCatch = open(temp_folder + 'outStream.txt', 'w')

	return call(bwa_build_call, stdout=outfileCatch)


def getNegStrandPerBaseCoverage(bamFilePath):

	#posStrand_call = [bedtoolsPath, 'genomecov', '-ibam', bamFilePath, '-d', '-strand', '+']

	#posStrandCatch = open('outStreamPosStrand.txt', 'w')

	#return call(posStrand_call, stdout=posStrandCatch)

	negStrand_call = [bedtoolsPath, 'genomecov', '-ibam', bamFilePath, '-d', '-strand', '-']

	negStrandCatch = open(temp_folder + 'outStreamNegStrand.txt', 'w')

	return call(negStrand_call, stdout=negStrandCatch)


def getPosStrandPerBaseCoverage(bamFilePath):

	#posStrand_call = [bedtoolsPath, 'genomecov', '-ibam', bamFilePath, '-d', '-strand', '+']

	#posStrandCatch = open('outStreamPosStrand.txt', 'w')

	#return call(posStrand_call, stdout=posStrandCatch)

	posStrand_call = [bedtoolsPath, 'genomecov', '-ibam', bamFilePath, '-d', '-strand', '-']

	posStrandCatch = open(temp_folder + 'outStreamPosStrand.txt', 'w')

	return call(posStrand_call, stdout=posStrandCatch)


def formatCoverageFiles():

	impGenome = input_parser(refGenome)

	genomeLength = len(impGenome[0]["DNA_seq"])

	print genomeLength

	posStrand = open(temp_folder + "outStreamPosStrand.txt", 'r')
	negStrand = open(temp_folder + "outStreamNegStrand.txt", 'r')
	bothStrand = open(temp_folder + "outStream.txt", 'r')

	posIndex = range(1, genomeLength + 1)
	#posIndex = range(1, 11)

	outDF = pd.DataFrame(np.nan, index=posIndex, columns=['both', 'pos', 'neg'])

	for line in posStrand:
		line = line.replace('\n', '')
		line = line.split('\t')
		outDF['neg'][int(line[1])] = float(line[2])

	for line in negStrand:
		line = line.replace('\n', '')
		line = line.split('\t')
		outDF['pos'][int(line[1])] = float(line[2])

	for line in bothStrand:
		line = line.replace('\n', '')
		line = line.split('\t')
		outDF['both'][int(line[1])] = float(line[2])




	return outDF


def getBaseCovStats(bamFilePath):

	getPerBaseCoverage(bamFilePath)
	getNegStrandPerBaseCoverage(bamFilePath)
	getPosStrandPerBaseCoverage(bamFilePath)

	perBaseMappingDF = formatCoverageFiles()

	# Cleaning up

	print "Add temp file removal here"

	return perBaseMappingDF


def getBaseCovStatsStranded(bothbamFilePath, forwardbamFilePath, reversebamFilePath):

	getPerBaseCoverage(bothbamFilePath)
	getNegStrandPerBaseCoverage(forwardbamFilePath)
	getPosStrandPerBaseCoverage(reversebamFilePath)

	perBaseMappingDF = formatCoverageFiles()

	# Cleaning up

	print "Add temp file removal here"

	return perBaseMappingDF


def plotSubRegion(startPos, stopPos, dataframe, fileName='aNewFig', windowSize=0):

	subDF = dataframe[startPos - windowSize : stopPos + windowSize]

	ax = subDF.plot()

	fig = ax.get_figure()

	fig.savefig(temp_folder + fileName + '.pdf')


def convertPosListToRanges(inList):
	for a, b in itertools.groupby(enumerate(inList), lambda (x, y): y - x):
		b = list(b)
		yield b[0][1], b[-1][1]


def findMappingPeaksOld(dataframe, threshold, windowSize=5, stepSize=3, col='both'):

	posCount = 0
	rowCount = 1

	returnPosList = []

	numberOfWindows = dataframe.shape[0] / stepSize

	#print numberOfWindows


	while posCount < numberOfWindows:
		windowList = []
		windowCount = 0
		while windowCount < windowSize:
			if math.isnan(dataframe[col][[rowCount + windowCount]]) == False:
				print '---'
				print rowCount + windowCount
				print dataframe.iloc[rowCount]

				windowList.append(float(dataframe[col][[rowCount + windowCount]]))

			windowCount += 1

		print windowList
		exit()

		if len(windowList) > 0:
			windAve = sum(windowList) / float(len(windowList))
		else:
			windAve = 0

		if windAve > threshold:
			#print 'Window frame'
			#print windAve
			#print rowCount

			windowPosList = range(rowCount, rowCount + stepSize)
			returnPosList = returnPosList + windowPosList


		rowCount = rowCount + stepSize
		#print rowCount
		posCount += 1


	uniqueReturnPosList = list(set(returnPosList))

	return uniqueReturnPosList
	#dataframe[]


def findMappingPeaks(dataframe, threshold, windowSize=5, stepSize=3, col='both'):

	posCount = 0
	rowCount = 1

	returnPosList = []

	numberOfWindows = dataframe.shape[0] / stepSize

	#print numberOfWindows
	for index, row in dataframe.iterrows():
		if float(row[col]) > float(threshold):
			returnPosList.append(index)


	return returnPosList
	#dataframe[]


def findMappingTroughs(dataframe, threshold, windowSize=5, stepSize=3, col='both'):

	posCount = 0
	rowCount = 1

	returnPosList = []

	numberOfWindows = dataframe.shape[0] / stepSize

	#print numberOfWindows
	for index, row in dataframe.iterrows():
		if float(row[col]) < float(threshold):
			returnPosList.append(index)


	return returnPosList
	#dataframe[]


def removeBothStrandRegions(annoFile, dataframe):


	anno_obj = input_parser(annoFile)

	removePosList = []
	
	for line in anno_obj:
		#print line
		if int(line[4]) != 4418548:

			featureRange = range(int(line[3]), int(line[4]))
			removePosList = removePosList + featureRange 

	uniqueRemovePosList = list(set(removePosList))

	dataframe = dataframe.drop(uniqueRemovePosList)

	return dataframe


def calcMappingdepthDist(dataframe, figName, col='both'):

	depthList = dataframe[col].tolist()
	print "Mean"
	print np.mean(depthList)
	print "Median"
	print np.median(depthList)
	print "Skew"
	print stats.skew(depthList)

	print 'Sorting'
	print len(depthList)
	srtdepthList = sorted(depthList)

	print 'Sorted'
	print len(srtdepthList)

	fit = stats.norm.pdf(srtdepthList, np.mean(srtdepthList), np.std(srtdepthList))

	pl.plot(srtdepthList,fit,'-o')

	pl.hist(srtdepthList,normed=True)

	print 'Plot'
	pl.savefig(temp_folder + figName)

	return np.mean(depthList)


def removeRegion(startPos, stopPos, dataframe):

	posList = range(startPos, stopPos + 1)

	indexList = dataframe.index

	for aPos in posList:
		if aPos in indexList:
			dataframe = dataframe.drop([aPos])

	return dataframe


def outParseToGFF(windowPassList, gffName, minDpth, chromName, minWindow=50):

	outFile = open(temp_folder + gffName + '_potential_deletions.gff3', 'w')

	annoSource = 'delPred'

	orientation = '.'

	outFile.write('##gff-version 3\n')

	sortedWindowPassList = sorted(windowPassList)

	windowFrames = list(convertPosListToRanges(sortedWindowPassList))

	count = 1
	for line in windowFrames:
		if line[1] - line[0] > minWindow:

			startPos = line[0]
			stopPos = line[1]

			additionalInfo = 'MinDepth=' + str(minDpth) + ";ID=del_region_" + str(count) + ";Name=del_region_" + str(count) + ";gene_name=del_region_" + str(count)


			outLineString = '\t'.join([chromName, annoSource, 'del_region', str(startPos), str(stopPos), '.', orientation, '.', additionalInfo, '\n'])

			outFile.write(outLineString)

			count += 1


def outParseDeletionsToGFF(dataframe, windowPassList, gffName, minDpth, chromName, minWindow=50):

	outFile = open(temp_folder + gffName + '_potential_deletions.gff3', 'w')

	annoSource = 'delPred'

	orientation = '.'

	sigmoid_curve_rate_threshold = 2

	outFile.write('##gff-version 3\n')

	sortedWindowPassList = sorted(windowPassList)

	windowFrames = list(convertPosListToRanges(sortedWindowPassList))


	# Scan for edges

	delWindow = 50

	count = 1
	for line in windowFrames:

		if line[1] - line[0] > minWindow:

			#print "\nDel region: " + str(count)
			startPos = line[0]
			stopPos = line[1]

			# Extract data for left region

			ydata_left = np.array(list(dataframe[startPos - delWindow:startPos + delWindow]['both']))
			xdata_left = np.array(np.arange(0,len(ydata_left),1))

			#print ydata_left
			#print xdata_left

			if np.mean(ydata_left) > 10:

				is_left_edge = check_for_edge_window(xdata_left, ydata_left, sigmoid_curve_rate_threshold, 'left', plot_curve=True, del_count=count)

			else:
				is_left_edge = False

			#print 'Result left:'
			#print is_left_edge

			# Extract data for right region

			ydata_right = np.array(list(dataframe[stopPos - delWindow:stopPos + delWindow]['both']))
			xdata_right = np.array(np.arange(0,len(ydata_right),1))

			#print ydata_right
			#print xdata_right

			if np.mean(ydata_right) > 10:

				is_right_edge = check_for_edge_window(xdata_right, ydata_right, sigmoid_curve_rate_threshold, 'right', plot_curve=True, del_count=count + 100)

			else:
				is_right_edge = False

			#print 'Result right:'
			#print is_right_edge

			if is_left_edge is True or is_right_edge is True:

				if is_left_edge is True:
					edgeLR = 'LeftDeletion'
				if is_right_edge is True:
					edgeLR = 'RightDeletion'
				if is_left_edge is True and is_right_edge is True:
					edgeLR = 'LeftRightDeletion'


				additionalInfo = 'MinDepth=' + str(minDpth) + ";ID=del_region_" + str(count) + ";Name=del_region_" + str(count) + ";gene_name=del_region_" + str(count) + ';edges=' + edgeLR


				outLineString = '\t'.join([chromName, annoSource, 'del_region', str(startPos), str(stopPos), '.', orientation, '.', additionalInfo, '\n'])

				outFile.write(outLineString)

			count += 1


def findsRNAregions(bamFilePath, annoFile, covSD, annoOutName, windowSize=50, stepSize=5):

	print "Making converage DF -", annoOutName
	RESdataframe = getBaseCovStats(bamFilePath)

	RESdataframe.to_pickle(temp_folder + annoOutName + '_allStrandIGR.pkl')

	print "Removing coding regions -", annoOutName
	IRGdataFrame = removeBothStrandRegions(annoFile, RESdataframe)

	print IRGdataFrame

	print "Getting Mapping stats"
	meanDepth = calcMappingdepthDist(IRGdataFrame, annoOutName, col='both')

	covThreshold = meanDepth * covSD

	print 'Coverage Threshold'
	print covThreshold

	print "Finding peaks", annoOutName
	peakList = findMappingPeaks(IRGdataFrame, covThreshold, windowSize=windowSize, stepSize=stepSize, col='both')

	chrName = getChromName(annoFile)

	print "Formatting for GFF3 -", annoOutName
	outParseToGFF(peakList, annoOutName, covThreshold, chrName)


def getDominantStrand(dataframe, startPos, stopPos):

	#print dataframe[int(startPos):int(stopPos)].sum()


	bothValues = dataframe[int(startPos):int(stopPos)].sum()['both']
	posValues = dataframe[int(startPos):int(stopPos)].sum()['pos']
	negValues = dataframe[int(startPos):int(stopPos)].sum()['neg']


	if posValues > negValues:
		strandSymbol = '+'
	else:
		strandSymbol = '-'

	return strandSymbol


def getChromName(annotationFile):

	temp_gff = input_parser(annotationFile)

	return temp_gff[0][0] 


def combineAllGff(GFFdir):

	posList = []



	for file in os.listdir(GFFdir):
		if file.endswith(".gff3"):

			openFile = open(GFFdir + file, 'r')

			for line in openFile:
				if line[0] != '#':
					lineList = line.split('\t')
					posList = posList + range(int(lineList[3]), int(lineList[4]) + 1)


	uniqueReturnPosList = list(set(posList))

	outParseToGFF(uniqueReturnPosList, 'Combined_sRNA_potentials', '2SD', 'NZ_CP012090.1')


def seqValStats(geneID, annoFile, genomeFile, e_val=0.05):

	resDict = {'ID':geneID}

	cwd = '/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/mappingTest/'

	anno_obj = input_parser(annoFile)
	genome_obj = input_parser(genomeFile)

	geneEnt = ''

	for geneEnt in anno_obj:
		if geneEnt[8]['Name'] == geneID:
			geneInfo = geneEnt

	geneSeq = genome_obj[0]['DNA_seq'][int(geneInfo[3]) - 1:int(geneInfo[4])]


	# Get flanking genes



	# Get uniqueness
	e_val = str(e_val)

	database_name = cwd +  'tempBlastDB'

	call(["makeblastdb", "-in", genomeFile, "-out", database_name, "-dbtype", "nucl"])

	temp_fasta = open(cwd + 'blasta.fa','w')
	temp_fasta.write(">temp_fasta_file\n")
	temp_fasta.write(geneSeq + "\n")
	temp_fasta.close()

	full_directory = cwd + "blasta.fa"

	blast_result_raw = Popen(["blastn", "-db", database_name, "-query", full_directory, "-evalue", e_val, "-outfmt", "7"], stdout=PIPE)
	
	blast_result = blast_result_raw.stdout.read()

	blast_list = blast_result.split('\n')

	passCount = 0

	for line in blast_list:
		if len(line) > 0:
			if line[0] != '#':
				listLine = line.split('\t')
				passCount += 1


	resDict['sequence'] = geneSeq
	resDict['Uniquness'] = passCount

	return resDict


# A set of functions fer detecting deletions

def sigmoid(x, x0, k, a, c):
	y = a / (1 + np.exp(-k * (x - x0))) + c
	return y


def check_for_edge_window(x_np_array, y_np_array, grad_threshold, orientation, plot_curve=False, del_count=0, window_size=6):

	is_edge = False

	if orientation == 'left':

		y_np_array = y_np_array[::-1]


	# sliding window funct

	position = 0

	while position < len(x_np_array) - window_size*2:

		window_list_1 = y_np_array[position:position + window_size]
		window_list_2 = y_np_array[position + window_size:position + window_size*2]

		ave_1 = reduce(lambda x, y: x + y, window_list_1) / len(window_list_1)
		ave_2 = reduce(lambda x, y: x + y, window_list_2) / len(window_list_2)

		gradient = (ave_2 - ave_1) / window_size

		if gradient > grad_threshold:
			is_edge = True

		position += 1

	if plot_curve is True:

		pl.plot(x_np_array, y_np_array, 'o', label='data')
		pl.ylim(0, 120.05)
		pl.legend(loc='upper left')
		pl.grid(True)
		if is_edge is False:
			pl.title('Negative')
		else:
			pl.title('Positive')

		fig_filename = "edge_" + str(del_count) + '.png'
		pl.savefig(fig_filename)
		pl.clf()

	return is_edge


def check_for_edge(x_np_array, y_np_array, rate_threshold, orientation, plot_curve=False, del_count=0):


	#xdata = np.array([1.0,2.0,3.0,4.0, 5.0, 6.0, 7.0, 8.0 ,9.0, 10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0])

	#pos_edge = np.array([20,21,24,21,23,25,21,2,1,0,2,1,0,2,1,0,0,1])
	#neg_edge = np.array([22.0, 21.0, 22.0, 18.0, 15.0, 10.0, 9.0, 5.0, 4.0, 4.0, 3.0, 2.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0])[::-1]


	#xdata = np.array([0.0,   1.0,  3.0,  4.3,  7.0,   8.0,   8.5, 10.0,  12.0, 14.0])
	#ydata = np.array([20,21,24,21,23,25,21,2,1,0,2,1,0,2,1,0,0,1])[::-1]
	#ydata = np.array([22.0, 21.0, 22.0, 18.0, 15.0, 10.0, 9.0, 5.0, 4.0, 4.0, 3.0, 2.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0])[::-1]

	if orientation == 'left':

		print 'array flipped'
		y_np_array = y_np_array[::-1]
		print 'flipped array'
		print y_np_array
		print '\n'
		print x_np_array


	popt = 0
	pcov = 0

	popt, pcov = curve_fit(sigmoid, x_np_array, y_np_array)

	if debug == True:
		print "Fit:"

		# This is the midpoint of the curve
		print "x0 =", popt[0]

		# This variable is the rate of change, what we are wanting.
		print "k  =", popt[1]

		#This is the top right y height of the curve
		print "a  =", popt[2]

		#This is the bottom left y height of the curve
		print "c  =", popt[3]
		print "Asymptotes are", popt[3], "and", popt[3] + popt[2]


	if plot_curve == True:

		x = np.linspace(-1, 120, 50)
		y = sigmoid(x, *popt)


		pl.plot(x_np_array, y_np_array, 'o', label='data')
		pl.plot(x,y, label='fit')
		pl.ylim(0, 100.05)
		pl.legend(loc='upper left')
		pl.grid(True)
		if abs(popt[1]) > rate_threshold:
			pl.title('Negative')
		else:
			pl.title('Positive')

		fig_filename = "edge_" + str(del_count) + '.png'
		pl.savefig(fig_filename)
		pl.clf()


	if abs(popt[1]) > rate_threshold:
		return True
	else:
		return False


def findDeletions(bamFilePath, annoFile, minCov, annoOutName, windowSize=50, stepSize=5):
	# minCov is minimum coverage 
	print "Making converage DF -", annoOutName
	RESdataframe = getBaseCovStats(bamFilePath)

	RESdataframe.to_pickle(temp_folder + annoOutName + '_allStrandIGR.pkl')

	print "Getting Mapping stats"
	meanDepth = calcMappingdepthDist(RESdataframe, annoOutName, col='both')

	covThreshold = minCov

	print 'Coverage Threshold'
	print covThreshold

	print "Finding deletions", annoOutName
	peakList = findMappingTroughs(RESdataframe, covThreshold, windowSize=windowSize, stepSize=stepSize, col='both')

	chrName = getChromName(annoFile)
	print "Get chromosome name"
	print chrName

	print "Formatting no-coverage regions for GFF3 -", annoOutName
	outParseToGFF(peakList, annoOutName + '_potential', covThreshold, chrName)

	print "Formatting deletions for GFF3"
	outParseDeletionsToGFF(RESdataframe, peakList, annoOutName, covThreshold, chrName)






'''
resDFIGR = pd.read_pickle('/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/Sample_507stationary1_allStrandIGR.pkl')

print resDFIGR

filtDF = removeBothStrandRegions(refAnno, resDFIGR)

print filtDF

filtDF.to_pickle('xTESTallStrandFiltered.pkl')


filtDF = pd.read_pickle('xTESTallStrandFiltered.pkl')
	
chrName = getChromName(refAnno)

print chrName



peakPosList = findMappingPeaks(filtDF, 300, windowSize=50, stepSize=5, col='both')


outParseToGFF(peakPosList, 'xTESTnew', 300, chrName, minWindow=50)
	
'''

#combineAllGff("/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/")

# ------------------------------------------------------------------------------------------
'''

filtDF = pd.read_pickle('/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/Sample_5527earlylog2_allStrandIGR.pkl')

plotSubRegion(674274, 674337, filtDF, fileName='sRNA_34', windowSize=100)

quit()


print "COUNT STUFF"

count_list = []
var_list = []

cuffdiff_file = "/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/mappingTest/cuffdiffRun/isoforms.count_tracking"

cuffDF = pd.read_table(cuffdiff_file, index_col=0)

for headItem in cuffDF.columns.values.tolist():
	if "count" in headItem:
		count_list.append(headItem)


for headItem in cuffDF.columns.values.tolist():
	if "_var" in headItem:
		var_list.append(headItem)


cuffDF['count_avg'] = cuffDF[count_list].mean(axis=1)
cuffDF['var_avg'] = cuffDF[var_list].mean(axis=1)

print cuffDF['count_avg']
cuffDFListcount = cuffDF['count_avg'].tolist()
from math import log
cuffDFListcount = [log(y,10) for y in cuffDFListcount]
cuffDFListcount = sorted(cuffDFListcount)

print cuffDF['var_avg'] 
cuffDFListvar = cuffDF['var_avg'].tolist()


xLab = cuffDF.index.values.tolist()

xlabPos = range(0,len(xLab))

result = cuffDF.sort('count_avg', ascending=False)
print '------------------'
print result['count_avg']

passCount = 0
for countVal in result['count_avg']:
	if countVal > 200:
		passCount +=1

sRNAanno = '/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/mappingTest/Combined_sRNA_potentials_potential_sRNA.gff3'

ordered_sRNA_list = result.index.values.tolist()

result['uniqueness'] = ''
result['sequence'] = ''

outsRNA_seq = open('/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/mappingTest/sRNAseq.fa', 'w')

for sRNA_name in ordered_sRNA_list:

	sRNA_seq_dict = seqValStats(sRNA_name, sRNAanno, refGenome)

	print sRNA_seq_dict

	fasta_head = ">" + sRNA_name + '\n'

	outsRNA_seq.write(fasta_head)

	outsRNA_seq.write(sRNA_seq_dict['sequence'] + '\n')

	result['uniqueness'][sRNA_name] = sRNA_seq_dict['Uniquness']
	result['sequence'][sRNA_name] = sRNA_seq_dict['sequence']



print result

result.to_csv('/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/mappingTest/sRNAsumary_apple.csv', sep=';')
'''

# ------------------------------------------------------------------------------------------

samtoolsPath = 'samtools'

bedtoolsPath = '/Users/panix/Dropbox/Programs/tools/bedtools2/bin/bedtools'

inBamDir = "/Users/panix/W-148_mapping/"

refGenome = "/Users/panix/Desktop/W-148-StudyVersion/W_148_NCBI.fa"
refAnno = "/Users/panix/Desktop/W-148-StudyVersion/W_148_NCBI.gff3"


#plt.bar(xlabPos, cuffDFListcount, color="w")
#plt.savefig('/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/mappingTest/noLogcountsRNAsorted.svg', format='svg', dpi=1200)

temp_folder = '/Users/panix/W-148_mapping/temp'
#bothbamFilePath = '/Volumes/External/CIDRI_Data/TASH_TERRY_DATA/del_detection/inBamFiles/sample_2_condition_genomeic_replicate_1_filtered_sorted.bam'
#forwardbamFilePath = '/Volumes/External/CIDRI_Data/TASH_TERRY_DATA/del_detection/inBamFiles/sample_2_condition_genomeic_replicate_1_filtered_sorted_fwd.bam'
#reversebamFilePath = '/Volumes/External/CIDRI_Data/TASH_TERRY_DATA/del_detection/inBamFiles/sample_2_condition_genomeic_replicate_1_filtered_sorted_rev.bam'

#candidateGFF = '/Volumes/HDD/Genomes/M_tuberculosis/sRNA/sRNAdiscovery/Combined_sRNA_potentials_potential_sRNA.gff3'

#resDF = getBaseCovStatsStranded(bothbamFilePath, forwardbamFilePath, reversebamFilePath)

#resDF.to_pickle(temp_folder + 'allStrand.pkl')

'''

RESdataframe = pd.read_pickle('/Volumes/External/CIDRI_Data/TASH_TERRY_DATA/del_detection/temp_filessample_2_condition_genomeic_replicate_1_filtered_sorted_allStrandIGR.pkl')

#peakList = findMappingTroughs(RESdataframe, 2, windowSize=50, stepSize=5, col='both')
#print peakList

#outParseToGFF(peakList, 'outAnnoTest_1', 2, 'TestChr')


xdata = np.array([1.0,2.0,3.0,4.0, 5.0, 6.0, 7.0, 8.0 ,9.0, 10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0])

pos_edge = np.array([20,21,24,21,23,25,21,2,1,0,2,1,0,2,1,0,0,1])
neg_edge = np.array([22.0, 21.0, 22.0, 18.0, 15.0, 10.0, 9.0, 5.0, 4.0, 4.0, 3.0, 2.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0])


debug = True

#buggy_data = np.array([   0. ,   0.   , 0.  ,  0.  ,  0.   , 0.  ,  0. ,   0.  ,  0. , 0. ,   0.  ,  0. ,   0.  ,  0.  ,  0. ,   0.  ,  0.  ,  0.  ,  0. ,   0. ,   0.  , 0.   , 0.  ,  0. ,   0. ,   0. ,   0.  ,  0.  ,  0.  ,  0.  ,  0. ,   0.  ,  0. ,   1.  ,  1. ,  64. ,  67. ,  71. ,  71.  , 71.  , 71. ,  71. ,  71. ,  72. ,  73.  , 77. ,  77. ,  78.  , 79. ,  80.  , 80.  , 82. ,  84. ,  85. ,  85.  , 85.  , 87.  , 91. ,  93. ,  96.  , 96. ,  96. ,  99. ,  99. , 101. , 102. , 103. , 105. , 102. , 102. , 101. , 103. , 102. , 100. , 100. ,  98. ,  96. ,  95. ,  95. ,  94. ,  94. ,  95. ,  96. ,  97. ,  96.])
#xdata = np.array(np.arange(0,len(buggy_data),1))

#print check_for_edge_window(xdata, buggy_data, 4, 'right', plot_curve=True)

# outParseDeletionsToGFF(dataframe, windowPassList, gffName, minDpth, chromName, minWindow=50)

quit()

'''

sampleDict = {}

for file in os.listdir(inBamDir):
	if file.endswith(".bam"):
		sampleDict[file.strip('.bam')] = os.path.join(inBamDir, file)


for exp in sampleDict.keys():
	print exp
	print sampleDict[exp]

	findDeletions(sampleDict[exp], refAnno, 2, exp, windowSize=50, stepSize=5)


'''
filtDF = pd.read_pickle(temp_folder + 'allStrand.pkl')

good_sRNA = ['sRNA_187','sRNA_47','sRNA_42','sRNA_211','sRNA_95','sRNA_170','sRNA_31','sRNA_214','sRNA_141','sRNA_213','sRNA_99','sRNA_144','sRNA_223','sRNA_135','sRNA_8','sRNA_3','sRNA_72','sRNA_127','sRNA_4','sRNA_62','sRNA_23','sRNA_174','sRNA_56','sRNA_164','sRNA_30','sRNA_158','sRNA_40','sRNA_78','sRNA_64','sRNA_22','sRNA_39','sRNA_175','sRNA_190','sRNA_151','sRNA_138','sRNA_186','sRNA_119','sRNA_71','sRNA_172','sRNA_45','sRNA_202','sRNA_92','sRNA_32','sRNA_98','sRNA_85','sRNA_5','sRNA_27','sRNA_117','sRNA_77','sRNA_163','sRNA_126','sRNA_112','sRNA_131','sRNA_153','sRNA_109','sRNA_220']

# Import old gtf

annoGTF = input_parser(candidateGFF, parse_as='default')

print len(annoGTF)

out_gtf = open(temp_folder + 'filterpass_sRNA_candidates_Stranded.gtf', 'w')

for line in annoGTF:

	if line[8]['Name'] in good_sRNA:

		strand = getDominantStrand(filtDF, line[3], line[4])

		line[6] = strand

		annoLine = line[0] + '\t' + line[1] + '\t' + 'exon' + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + 'gene_id ' + line[8]['Name'] + '; ' + 'transcript_id ' + line[8]['Name'] + ';' +  '\n'

		out_gtf.write(annoLine)


		print line


#plotSubRegion(4107239, 4107538, filtDF, fileName='sRNA_211', windowSize=100)


quit()
'''

'''
sampleDict = {}

for file in os.listdir(inBamDir):
	if file.endswith(".bam"):
		sampleDict[file.strip('.bam')] = os.path.join(inBamDir, file)


for exp in sampleDict.keys():
	print exp
	print sampleDict[exp]

	findsRNAregions(sampleDict[exp], refAnno, 2, exp, windowSize=50, stepSize=5)


'''