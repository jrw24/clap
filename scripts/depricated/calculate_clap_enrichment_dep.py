__author__ = 'Jamie Wangen'

"""
This script will calculate enrichments from a CLAP/CLIP dataset 
	based on enrichment realative to an input dataset

The input arguments are:
	Input compressed numpy arrays with clap reads mapped to single nucleotide genomic coordinates
	Capture compressed numpy arrays with clap reads mapped to single nucleotide genomic coordinates
	Bedfile with genome annotation
	desired binSize to use for calculating enrichment windows


Strategy for enrichment calculation:
	1) iterate through normalized np array input, 1 gene at a time
	2) count reads in 100bp bins (or binSize bp) within each reads
	3) calculate the median value of counts per bin for that gene
	4) discard genes with no counts / median count of 0
		** still not totally sure what the best way to manage zero bins will be
		** currently with numpy, N/0 == Inf and 0/0 == NaN
	5) determine the max() of read counts or medians for input bins
	6) save these to ordered dictionaries
	7) iterate through normalized np array capture, 1 gene at a time
	8) count reads in same size bins 
	9) calculate fold enrichment for each bin of capture to input
	10) calculate statistical test that will determine significance of enrichment difference


To Do List (230411):
	- addition of pseudocounts for low expressed input genes
		- how should we handle Inf?
	- make qc plots of reads per window across the genome
	- incorporate strandedness
		- this will require minor modifactions to npz array generation
		- have the code, just need to uncomment and change to np.divide()


"""



import sys
import os
import multiprocessing
import pandas as pd
import numpy as np 
from collections import OrderedDict
from datetime import datetime
from pathos.pools import ProcessPool
from scipy.stats import hypergeom 
import argparse
import importlib
import subprocess

parser= argparse.ArgumentParser()
parser.add_argument('--samp', help= 'name of sample being processed')
parser.add_argument('--libset', help= 'path to libsettings file')
parser.add_argument('--binSize', help = 'desired size for calculating bins')
parser.add_argument('--threadNumb', help= 'number of threads to use')
args = parser.parse_args()

### load libset into namespace 
libset = importlib.machinery.SourceFileLoader('libset', args.libset).load_module()
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)

def load_bedfile(bedInput):
		"""
		take an input bedfile, and load this into a pandas data.frame
		"""
		# bedInput = '/Users/jamiewangen/Desktop/HPC-Mount-guttman/genomes/mouse/GRCm38_p6/appris_pricipal_gtf/gencode.vM25.appris_principal.unique.bed'
		b = pd.read_csv(bedInput, sep="\t", header=None)

		### check column length here:

		b.columns = [   
			'chrom', 
			'chromStart', 
			'chromEnd', 
			'name', 
			'score', 
			'strand', 
			'thickStart', 
			'thickEnd', 
			'itemRgb', 
			'blockCount', 
			'blockSizes', 
			'blockStarts', 
			'gene_id', 
			'gene', 
			'gene_type'
		]
		return b 

class clapEnrichment(object):
	"""
	Create a class to manage bed enrichment calculation

	_____________
	Output enriched table:
		Write sufficiently enriched peaks to an output data.frame

	_____________
	Output chromsome table:
		Want to write out a tsv table for all bins across the chromosome

	"""


	def __init__(self, inputNpzArray, captureNpzArray, bed, binSize, outRootDir, 
			samp, inputTotalReads, captureTotalReads, **kwargs):
		self.inputNpzArray = inputNpzArray
		self.captureNpzArray = captureNpzArray
		self.bed = bed
		self.binSize = int(binSize)
		self.outRootDir = outRootDir
		self.samp = samp 
		self.inputTotalReads = inputTotalReads
		self.captureTotalReads = captureTotalReads
		self.__dict__.update(kwargs)
	
	def make_outDir(self):
		if not os.path.exists(self.outRootDir):	os.makedirs(self.outRootDir)


	def start_end_arrays(self, start, end, binSize):
		"""
		given starting position, ending position, and binSize
			determine start positions and end positions for a given gene
		"""

		### Check to see if gene is divisible by binSize
		###	If not, then need to set the final window to overlap the penultimate window

		if (end-start) % binSize == 0:
			### start array
			startArray = np.arange(start,end,binSize)
			### end array
			endArray = np.arange(start+binSize,end+binSize,binSize)
		else:
			### start array		   
			startArrayPartial = np.arange(start,end-binSize,binSize)
			startArrayLast = np.arange(end-binSize,end,binSize)
			startArray = np.concatenate([startArrayPartial, startArrayLast])
			### end array
			endArrayPartial = np.arange(start+binSize,end,binSize)
			endArrayLast = np.arange(end,end+binSize,binSize)
			endArray = np.concatenate([endArrayPartial, endArrayLast])

		return (startArray, endArray)

	def start_end_arrays_chrom(self, start, end, binSize):
		"""
		given starting position, ending position, and binSize
			get the bins across an entire chromosome
			do not include the final bin 
		"""

		### start array
		startArray = np.arange(start,end,binSize)
		### end array
		endArray = np.arange(start+binSize,end+binSize,binSize)

		return (startArray, endArray)

	def bin_array_counts(self, npArray):
		"""
		Take an input 1d numpy array, and slice this into 2d arrays of binSize
			Then calculate the sum of entries within each sliced array

			If the array is not evenly divisible by binsize, 
				calculate an additional array from the -binSize to the end of the array

			example:
				binsize = 3
				in = [1,2,3,4,5,6,7,8]

				bins: 	
					[1,2,3]
					[4,5,6]
					[6,7,8] <-- overlap at value 6, window here will be calculated twice

		Currently using overlap of final bin to include the end of the gene
			possible alternative to search downstream for next gene, and include 'unannotated' region

		Other possibility is to just calculate these windows for the entire chromosome, 
			agnostic of predetermined gene positions
		
		"""
		
		### check to see if final window need to be caluculated (array length / binsize remainder == 0)
		if npArray.size % self.binSize == 0:

			binnedArray = npArray.reshape(-1, self.binSize)				
			binCounts = binnedArray.sum(axis=1)
		
		else:

			binnedArray = npArray[:(npArray.size // self.binSize) * self.binSize].reshape(-1, self.binSize)
			binPartialCounts = binnedArray.sum(axis=1)
			finalBin = npArray[-self.binSize:] ## overlap penultimate bin to include terminus of gene
			finalBinCounts = finalBin.sum()

			binCounts = np.concatenate([binPartialCounts, [finalBinCounts]])
			
		return binCounts
			
	def bin_array_counts_chrom(self, npArray):
		"""
		Take an input 1d numpy array, and slice this into 2d arrays of binSize
			Then calculate the sum of entries within each sliced array

		This is for getting counts across an entire chrom
			incomplete bins at the end will be excluded
		
		"""
		binnedArray = npArray.reshape(-1, self.binSize)				
		binCounts = binnedArray.sum(axis=1)
			
		return binCounts
			

	def median_floor_array(self, npArray):
		"""
		look up the median value across an array
			and set values to the median if they are less than the median
		"""

		medianVal = np.median(npArray)
		print(medianVal, '<--- median Val')
		npArray[npArray < medianVal] = medianVal
		
		return npArray

	def calc_gene_enrichment_single_chrom(self, chrom):
		"""
		for a given chromosome specified by chrom
			take the input and capture arrays, and calculate enrichment scores
			only for each gene in the bedfile on that chromosome

		"""
		### set error messages to ingore warning for division by zero
		### these will be set to 'Inf' and filtered later if desired
		np.seterr(divide='ignore') 
		
		### slice the bed data.frame to include only the chromosome of interest
		print('slicing bedChrom')
		bedChrom = self.bed.loc[self.bed['chrom']==chrom]

		print('loading numpy arrays')
		### load in the input npzArray into memory for the chromosome of interest
		inputArray = self.inputNpzArray[chrom]
		captureArray = self.captureNpzArray[chrom]

		### create an list to store output data.frames for each gene in the chromosome
		outDfListChrom = []

		### define the output column names for each data.frame
		outDfColNames = [
				'chrom',
				'chromStartWindow',
				'chromEndWindow',
				'gene',
				'inputReadCount',
				'inputFloorCount',
				'captureReadCount',
				'enrichValue',
				]

		for row in bedChrom.index:

			gene = bedChrom.loc[row,'gene']
			start = int(bedChrom.loc[row, 'chromStart'])
			end = int(bedChrom.loc[row, 'chromEnd'])

			### slice the arrays to get only the region containing the gene
			iarray = inputArray[start:end]
			carray = captureArray[start:end]

			### get the start positions and end positions over the gene
			startArray, endArray = self.start_end_arrays(start, end, self.binSize)

			### get the counts from the input file
			inputCounts = self.bin_array_counts(iarray)

			### set the minimum value over the inputCounts to the median value of input counts
			floorInputCounts = self.median_floor_array(inputCounts)

			### calculate counts from the caputre file
			captureCounts = self.bin_array_counts(carray)

			### calculate enrichments of the caputure sample relative to the floored input sample
			enrichArray = captureCounts/floorInputCounts

			### build arrays with the chrom and gene information of the same length
			chrArray = np.repeat(chrom, startArray.size)
			geneArray = np.repeat(gene, startArray.size)

			### store all of the arrays in a single list
			outArrays = [
				chrArray, 
				startArray, 
				endArray, 
				geneArray, 
				inputCounts,
				floorInputCounts,
				captureCounts,
				enrichArray 
				] 

			### zip these into a dictionary and convert to a dataframe
			geneDf = pd.DataFrame.from_dict(dict(zip(outDfColNames, outArrays)))

			### append this data.frame to the list of data.frames for that chromosome
			outDfListChrom.append(geneDf)

		### concatnate the list of data.frames into a single dataframe per chromosome
		dfOutChrom = pd.concat(outDfListChrom)

		return dfOutChrom

	def calc_chrom_enrichment_single_chrom(self, chrom):
		"""
		for a given chromosome specified by chrom
			take the input and capture arrays, and calculate enrichment scores
			across the entire chromosome

		"""
		### set error messages to ingore warning for division by zero
		### these will be set to 'Inf' and filtered later if desired
		np.seterr(divide='ignore') 

		print('loading numpy arrays')
		### load in the input npzArray into memory for the chromosome of interest
		inputArray = self.inputNpzArray[chrom]
		captureArray = self.captureNpzArray[chrom]

		### create an list to store output data.frames for each gene in the chromosome
		outDfListChrom = []

		### define the output column names for each data.frame
		outDfColNames = [
				'chrom',
				'chromStartWindow',
				'chromEndWindow',
				'inputReadCount',
				'inputFloorCount',
				'captureReadCount',
				'enrichValue',
				]

		start = 0
		end = int(inputArray.size)

		### get the start positions and end positions over the gene
		startArray, endArray = self.start_end_arrays(start, end, self.binSize)

		### get the counts from the input file
		inputCounts = self.bin_array_counts(inputArray)

		### set the minimum value over the inputCounts to the median value of input counts
		floorInputCounts = self.median_floor_array(inputCounts)

		### calculate counts from the capture file
		captureCounts = self.bin_array_counts(captureArray)

		### calculate enrichments of the caputure sample relative to the floored input sample
		enrichArray = captureCounts/floorInputCounts

		### build arrays with the chrom and gene information of the same length
		chrArray = np.repeat(chrom, startArray.size)

		### store all of the arrays in a single list
		outArray = [
			chrArray, 
			startArray, 
			endArray,  
			inputCounts,
			floorInputCounts,
			captureCounts,
			enrichArray 
			] 

		### zip these into a dictionary and convert to a dataframe
		dfOutChrom = pd.DataFrame.from_dict(dict(zip(outDfColNames, outArray))) 
		return dfOutChrom


	def gene_enrichments_iterate_through_chroms(self):
		
		"""
		for iterating through chromosomes one at time to build the gene-level output dataframe
		"""

		### define output data.frame file
		dfOutFile = '%s/%s_clapGeneEnrichments_gene.tsv.gz' % (self.outRootDir, self.samp)

		### define a list of valid chromosomes, in order they appear in the bedfile
		### must start with 'chr'
		chromList = []
		for entry in self.bed['chrom']:
			if entry not in chromList and entry[0:3]=='chr':
				chromList.append(entry)
		print(f'valid chromosomes --> {chromList}')

		### define OrderedDict that will store the output dataframes

		dfOutDict = OrderedDict()

		for chrom in chromList:
			print(f'starting analysis for {chrom}')
			dfOutChrom = self.calc_gene_enrichment_single_chrom(chrom)
			dfOutDict[chrom] = dfOutChrom

		### after collecting dataframes for each chromosome, build into master dataframe
		dfOut = pd.concat(list(dfOutDict.values()))

		### reset the index so all values are unique
		dfOut.reset_index(drop=True, inplace=True)

		print('writing data.frame %s' % dfOutFile)

		dfOut.to_csv(dfOutFile, sep='\t', compression='gzip')
		return dfOut 

	def chrom_enrichments_iterate_through_chroms(self):
		
		"""
		for iterating through chromosomes one at time to build the chrom-level output dataframe
		"""

		### define output data.frame file
		dfOutFile = '%s/%s_clapGeneEnrichments_chrom.tsv.gz' % (self.outRootDir, self.samp)

		### define a list of valid chromosomes, in order they appear in the bedfile
		### must start with 'chr'
		chromList = []
		for entry in self.bed['chrom']:
			if entry not in chromList and entry[0:3]=='chr':
				chromList.append(entry)
		print(f'valid chromosomes --> {chromList}')

		### define OrderedDict that will store the output dataframes

		dfOutDict = OrderedDict()

		for chrom in chromList:
			print(f'starting analysis for {chrom}')
			dfOutChrom = self.calc_chrom_enrichment_single_chrom(chrom)
			dfOutDict[chrom] = dfOutChrom

		### after collecting dataframes for each chromosome, build into master dataframe
		dfOut = pd.concat(list(dfOutDict.values()))

		### reset the index so all values are unique
		dfOut.reset_index(drop=True, inplace=True)

		print('writing data.frame %s' % dfOutFile)

		dfOut.to_csv(dfOutFile, sep='\t', compression='gzip')
		return dfOut

	def filter_enriched_peaks_chrom(self, df, name):
		"""
		take an input data.frame for the whole geneome
			and extract the most enriched peaks only
		
		perform hypergeometric test to obtain 
			P value for statistical significance

		look up overlapping genes from the bedFile and 
			add these to the table

		write out the most enriched genes

		"""


		print('--- log transforming enrchments ---',datetime.now())
		### calculate log10 enrichment values for each window
		df['log10enrich']=np.log(df['enrichValue'])

		### replace all 'Inf' values with 'NaN'
		### and then drop all 'NaN' values in the table
		df.replace([np.inf, -np.inf], np.nan, inplace=True)
		df.dropna(axis=0, inplace=True)

		### --- perform hypergeometric test ---
		"""
		discrete probability function that requires integers and not floats!
			Can not draw 6.666 jellybeans out a jar...
			Therefore, we must convert back to integer numbers 
				of reads from normalized reads

		R:
		read.stats$phyper.bis <- 
			phyper(
				chip.counts = capReads,
				m.per.window = windowReads,
				n.per.window = outsideReads,
				k = captot
				)

		Python:
		hp = hypergeom.sf(
			k=capReads, ## number of successes
			M=libraryTot, ## population size (also called N)
			n=captot, ## successes in populatin (also called K)
			N=windowReads ## sample size (also called n)
			)

		"""

		print('--- calculating hypergeometric tests ---', datetime.now())
		### retrieve total mapped reads for input (in millions of reads)
		with open(self.inputTotalReads,'r') as f:
			inputTotFloat = float(f.read())
		inputTot = int(np.round(inputTotFloat * 1E6)) ## convert back to total reads

		with open(self.captureTotalReads,'r') as f:
			capTotFloat = float(f.read())
		capTot = int(np.round(capTotFloat * 1E6)) ## np.round rounds correctly here, int() does not always

		libraryTot = inputTot + capTot 

		for row in df.index:
			### convert to raw read values
			inputReads = int(np.round(df.loc[row,'inputReadCount']*inputTotFloat))
			capReads = int(np.round(df.loc[row,'captureReadCount']*capTotFloat))
			### sum total reads in window (total number of jellybeans)
			windowReads = inputReads + capReads
			### calc remaining reads not in window (jellybeans still in the jar)
			outsideReads = libraryTot - windowReads

			### calculate a p-value using the hypergeometric distribution
			### sf ('survival function') is 1-cdf (opposite of the cumulative distribution function)
			hp = hypergeom.sf(
				k=capReads, ## number of successes
				M=libraryTot, ## population size - total reads
				n=capTot, ## total number of reads in capture
				N=windowReads ## total number of reads in window
				)

			### finally, add the p-value to the data.frame
			df.at[row,'hypergeom_p'] = hp 

		### --- Bed file gene overlaps --- ###
		print('--- looking up gene overlaps ---', datetime.now())
		for row in df.index:
			
			chrom = df.loc[row,'chrom']
			s = int(df.loc[row,'chromStartWindow'])
			e = int(df.loc[row,'chromEndWindow'])

			bedchr = self.bed.loc[self.bed['chrom']==chrom]
			bedChrom = bedchr[(bedchr['chromStart'] <= s) & (bedchr['chromEnd'] >= e )]

			# bedChrom = self.bed.loc[(self.bed['chrom']==chrom) & (self.bed['chromStart'] <= s) & (self.bed['chromEnd'] >= e )]
	
			geneList = []
			if len(bedChrom) >= 1:
				for g in bedChrom.index:
					bedgene = bedChrom.loc[g,'gene']
					if bedgene not in geneList:
						geneList.append(bedgene)
			
				joinGene = ("_".join(geneList))			    
				df.at[row,'gene'] = joinGene
			else:
				df.at[row,'gene'] = 'none'

		### sort dataframe by the most enriched genes (logscale) and write to output table
		print('--- writing output enrichment table ---', datetime.now())
		df.sort_values('log10enrich', ascending=False, inplace=True)
		df.reset_index(drop=True,inplace=True)

		dfOutFile = '%s/%s_%s_topEnriched_chrom.tsv.gz' % (self.outRootDir, self.samp, name)
		df.to_csv(dfOutFile, sep='\t', compression='gzip')




def main():
	
	print('starting enrichment analysis', datetime.now())

	clapEnrichDir = f'{rootDir}/clapEnrichment'

	### load the bedfile as a pandas data.frame
	bed = load_bedfile(bedFile)

	print(bed.head())

	### make generator objects for loading numpy arrays
	inputSamp = input_sample[0]
	inputNpPath = f'{rootDir}/npArrayDir/{inputSamp}/{inputSamp}_secondReadGenomeRPM.npz'
	captureNpPath = f'{rootDir}/npArrayDir/{args.samp}/{args.samp}_secondReadGenomeRPM.npz'

	inputNpzArray = np.load(inputNpPath)
	captureNpzArray = np.load(captureNpPath)

	inputTotalReads= f'{rootDir}/npArrayDir/{inputSamp}/{inputSamp}_npReadsPerMillion.txt'
	captureTotalReads= f'{rootDir}/npArrayDir/{args.samp}/{args.samp}_npReadsPerMillion.txt'

	clap = clapEnrichment(
		inputNpzArray=inputNpzArray,
		captureNpzArray=captureNpzArray,
		bed=bed,
		outRootDir=clapEnrichDir,
		binSize=args.binSize,
		samp=args.samp,
		inputTotalReads=inputTotalReads,
		captureTotalReads=captureTotalReads
		)


	df = clap.chrom_enrichments_iterate_through_chroms()
	clap.filter_enriched_peaks_chrom(df, 'chrom')

	dfgene = clap.gene_enrichments_iterate_through_chroms()
	clap.filter_enriched_peaks_chrom(dfgene, 'gene')

	print('done', datetime.now())


if __name__ == '__main__':
	main()








