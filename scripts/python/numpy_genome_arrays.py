"""
Input: 
	- Bam file of aligned reads

Outputs:
	- numpy arrays for each chromosome, stored in a unified folder
	- bigwig file of aligned reads

Workflow for calculating genome coverage:

	1) -- Iterate through a single chromosome in the genome,
		Retrieve read positions for 2nd reads --
		** Use multiprocessing here for parallelization **
			Add 2nd read stop positions to NP array
				all, pos, and negative strands
			Add score of 1 for each position in in read for coverage array
				all only for now
		Write these *4* arrays to compressed np.array files
		Write the total valid reads counted for each chromsome

	2) -- combine single numpy arrays into a compressed dictionary archive --
		Load each chromosome np.array into memory and add to a dictionary
		Write this to an output file
		delete the individual chromosome arrays

	3) -- Perform read normalization --
		take the sum of reads in the np.array
		Write the total_validReads.txt file

		Load each np.array and divide each position by the normalizer value
		write this compressed dictionary to an output file

	4) -- Bedgraph output file --
		Write each numpy array to a bedgraph file 

	5) -- Bigwig output file --
		Convert the bedgraph file into a bigwig file 
		Delete the bedgraph files to convserve storage space

"""

import argparse
import importlib
import subprocess
import sys
import os
import twobitreader
import glob
import numpy as np 
from pathos.pools import ProcessPool
import pysam
import csv
from datetime import datetime

parser= argparse.ArgumentParser()
parser.add_argument('--samp', help= 'name of sample being processed')
parser.add_argument('--libset', help= 'path to libsettings file')
parser.add_argument('--threadNumb', help= 'number of threads to use')
args = parser.parse_args()

### load libset into namespace 
libset = importlib.machinery.SourceFileLoader('libset', args.libset).load_module()
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)


class npArrayBuilder(object):

	def __init__(self, inBam, outDir, genome, samp, **kwargs):
		self.inBam = inBam
		self.outDir = outDir
		self.genome = genome
		self.samp = samp 
		self.__dict__.update(kwargs)


	def build_coverage_array_singleChrom(self, chrom):
		"""

		"""

		bamfile = pysam.AlignmentFile(self.inBam, "rb")

		# wrong_strand_reads = 0
		valid_reads = 0
		valid_read_pos = 0
		valid_read_neg = 0
		# reads_out_of_bounds = 0
		not_proper_pair = 0
		qc_fail = 0

		print("building arrays for %s" % (chrom)) 
		# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))
		### define np arrays that are the length of each chromosome
		### create seperate arrays for 2nd reads, and for coverage
		
		## -- both strands -- ##
		chromSecReadArray = np.zeros(len(self.genome[chrom]))
		chromCovArray = np.zeros(len(self.genome[chrom]))
		# chromCovReadLenNormArray = np.zeros(len(genome[chrom]))

		## postivie strand -- ##
		chromSecReadArrayPos = np.zeros(len(self.genome[chrom]))
		# chromCovArrayPos = np.zeros(len(genome[chrom]))
		# chromCovReadLenNormArrayPos = np.zeros(len(genome[chrom]))

		## negative strand -- ##
		chromSecReadArrayNeg = np.zeros(len(self.genome[chrom]))
		# chromCovArrayNeg = np.zeros(len(genome[chrom]))
		# chromCovReadLenNormArrayNeg = np.zeros(len(genome[chrom]))


		### iterate over each read in the bamfile that is found on that chromosome
		for read in bamfile.fetch(chrom): 

			if not read.is_proper_pair:
				not_proper_pair+=1
				continue
			if read.is_qcfail:
				qc_fail+=1
				continue

			if read.is_read2: #correct strand for clap data is read2 alignments
				if read.is_reverse:
					read_strand = "-"
				else:
					read_strand = "+"
				
				# if not read_strand == geneFeatures['strand']:
				# 	wrong_strand_reads +=1
				# 	continue
				
				if read_strand == "+":

					readGenomicCoord = read.reference_start ## this is 0 based, add 1 for gff coords

					### add to second read array
					chromSecReadArray[readGenomicCoord] +=1 ## add count of 1 for position of second read
					chromSecReadArrayPos[readGenomicCoord] +=1

					#### build read position array relative to genomic chromosome position
					read_array = []
					counter = 0
					read_len = len(read.seq)
					for tup in read.cigartuples:
						if tup[0] == 0:
							for i in range(tup[1]):
								# if readGenomicCoord+counter < len(geneArrayCov): ## prevent reads that run off end of transcript
								# chromCovArray[readGenomicCoord+counter]+=1
								read_array.append(readGenomicCoord+counter) ### add this position to the read array
								counter +=1
						elif tup[0] == 1: ## insertion to reference
							read_len +=1
						elif tup[0] == 2: ## deletion to reference
							read_len -=1
						elif tup[0] == 3: ## spliced junction
							counter+=tup[1] ## add spliced distance
						elif tup[0] ==4: ## softclip
							pass ## softclips should not be included in read_genomic_coord start position
						else:
							pass ## ignore other cigar entries
					### add each position in the read to the chromosome array
					for pos in read_array:
						chromCovArray[pos] +=1
						# chromCovArrayPos[pos] +=1
					# rl_norm = 1.0/len(read_array)

					### add each position in read to chromosome array, divide by read_length
					# 	### each read has a sum of 1.0 now
					# for pos in read_array:
					# 	chromCovReadLenNormArray[pos] += rl_norm
					# 	chromCovReadLenNormArrayPos[pos] += rl_norm

					valid_reads += 1
					valid_read_pos += 1

				if read_strand == "-":
					readGenomicCoord = read.reference_end-1 ## this is 1 based for actual position, subtract 1 for zero based
					
					### add to second read array
					chromSecReadArray[readGenomicCoord] +=1 ## add count of 1 for second reads
					chromSecReadArrayNeg[readGenomicCoord] +=1

					### build read position array relative to genomic chromosome position
					read_array = []
					counter = 0
					read_len = len(read.seq)
					for tup in read.cigartuples:
						# print(tup[0])
						if tup[0] == 0:
							for i in range(tup[1]):
								read_array.append(readGenomicCoord+counter)
								# if tr_coord+counter >= 0: ## prevent reads that run off start of transcript
								# 	# print(i, tr_coord, counter)
								# 	geneArrayCov[tr_coord+counter]+=1
								counter -=1
					#             print(i)
						elif tup[0] == 1: ## insertion to reference
							read_len +=1
						elif tup[0] == 2: ## deletion to reference
							read_len -=1
						elif tup[0] == 3: ## spliced junction
							counter-=tup[1] ## add spliced distance
						elif tup[0] ==4: ## softclip
							pass ## softclips should not be included in read_genomic_coord start position
						else:
							pass ## ignore other bam entries
					### add each position in the read to the chromosome array
					for pos in read_array:
						chromCovArray[pos] +=1
						# chromCovArrayNeg[pos] +=1
					# rl_norm = 1.0/len(read_array)

					### add each position in read to chromosome array, divide by read_length
						### each read has a sum of 1.0 now
					# for pos in read_array:
					# 	chromCovReadLenNormArray[pos] += rl_norm
					# 	chromCovReadLenNormArrayNeg[pos] += rl_norm

					valid_reads += 1
					valid_read_neg += 1

		### convert read length normalized array to reads per million (including only number of valid reads)
		# print("valid reads == %s" % (valid_reads))
		# readsPerMillion = float(valid_reads)/1E6 ## calculate reads per million normalization value

		# chromRpmRlNormArray = [x/readsPerMillion for x in chromCovReadLenNormArray]

		print("finished chrom array for %s" % chrom)


		## write total valid reads to an output file
		valid_reads_file = "%s/%s_%s_validReads.txt" % (self.outDir, chrom, self.samp)
		valid_reads_file_pos = "%s/%s_%s_validReadsPos.txt" % (self.outDir, chrom, self.samp)
		valid_reads_file_neg = "%s/%s_%s_validReadsNeg.txt" % (self.outDir, chrom, self.samp)
		
		# with open(valid_reads_file, 'w') as f:
		# 	f.write("%s" % valid_reads)

		# with open(valid_reads_file_pos, 'w') as f:
		# 	f.write("%s" % valid_read_pos)

		# with open(valid_reads_file_neg, 'w') as f:
		# 	f.write("%s" % valid_read_neg)

		### write single arrays to numpy array files
		secReadOut = "%s/%s_%s_secondReadGenome.npy" % (self.outDir, self.samp, chrom)
		# secReadDictRpmOut = "%s/%s_secondReadGenomeDictRpm.npz" % (outDir, samp)
		covOut = "%s/%s_%s_covGenome.npy" % (self.outDir, self.samp, chrom)
		# covDictRlNormOut = "%s/%s_covGenomeDictRlNormOut.npz" % (outDir, samp)
		# covDictRpmRlNormOut = "%s/%s_covDictRpmRlNormOut.npz" % (outDir, samp)
		secReadPosOut = "%s/%s_%s_secondReadGenomePos.npy" % (self.outDir, self.samp, chrom)
		secReadNegOut = "%s/%s_%s_secondReadGenomeNeg.npy" % (self.outDir, self.samp, chrom)

		print("outfile paths: \n %s \n %s \n %s \n %s" % (secReadOut, covOut, secReadPosOut, secReadNegOut))

		np.save(secReadOut, chromSecReadArray)
		np.save(covOut, chromCovArray)
		np.save(secReadPosOut, chromSecReadArrayPos)
		np.save(secReadNegOut, chromSecReadArrayNeg)



def concat_np_array_single(outDir, samp, chromList, arraySuffix):
	"""
	After assigning reads to numpy arrays, 
		load these into memory, 
		add them to a dictionary,
		write the dictionary as a compressed npz archive,
		delete the uncompressed numpy array binary file

	Naming structure:
		samp_chrom_secondReadGenome.npy
		samp_chrom_covGenome.npy
		samp_chrom_secReadGenomePos.npy
		samp_chrom_secReadGenomeNeg.npy
	"""
	npArrayDict = dict()
	npArrayDictOutfile = '%s/%s_%s.npz' % (outDir, samp, arraySuffix)
	for chrom in chromList:
		chromArrayFile = '%s/%s_%s_%s.npy' % (outDir, samp, chrom, arraySuffix)
		chromArray = np.load(chromArrayFile)
		npArrayDict[chrom] = chromArray
		subprocess.Popen('rm %s' % (chromArrayFile), shell=True).wait()
	np.savez_compressed(npArrayDictOutfile, **npArrayDict)

def normalize_np_array(outDir, samp, chromList, arraySuffix):
	"""
	calculate the total reads present in a chromosome array
		divide each position in the array by reads-per-million normalizer

	write the normalized RPM array

	"""

	npzFile = '%s/%s_%s.npz' % (outDir, samp, arraySuffix) ## reads to be normalized
	npzFileAll = '%s/%s_secondReadGenome.npz' % (outDir, samp) ## use all reads for norm
	rpmFile = '%s/%s_npReadsPerMillion.txt' % (outDir, samp) ## write read normalizer to txt file
	npzRpmFile = '%s/%s_%sRPM.npz' % (outDir, samp, arraySuffix)
	npArrayDictAll = np.load(npzFileAll)

	totalReads = 0
	for chrom in chromList:
		chromReads = npArrayDictAll[chrom].sum()
		totalReads += chromReads
	readsPerMillion = totalReads/1E6
	print('reads per million = %s' % (readsPerMillion))

	## write readsPerMillion to a text file
	with open(rpmFile, 'w') as f:
		f.write('%s' % (readsPerMillion))
	f.close()

	npArrayDictRpm = dict()
	npArrayDict = np.load(npzFile)

	for chrom in chromList:
		npArrayDictRpm[chrom] = np.divide(npArrayDict[chrom], readsPerMillion)
	np.savez_compressed(npzRpmFile, **npArrayDictRpm)

def bedGraph_writer_genome(npArray, outBedGraph):
	"""
	npArray --> numpy array with reads assigned to genomic positions
	outBedGraph --> output file to write the bedgraph
	"""

	print("--> writing bedGraph to: ", outBedGraph)
	chromArrayDict = np.load(npArray) ### load the numpy array dictionary

	with open(outBedGraph, 'w') as bedGraph:
		bgwriter = csv.writer(bedGraph, delimiter='\t')
		
		for chrom in chromArrayDict:
			print("writing chrom: %s" % (chrom))
			print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))

			chromArray = chromArrayDict[chrom]
			curPos = -1
			counter = 0
			val = chromArray[0] ## set to the starting value
			startPos = float(curPos+1) ## 

			for i in chromArray:
				curPos +=1 
				
				if i != val:
					endPos = float(curPos)
					row = [chrom, int(startPos), int(endPos), val]
					bgwriter.writerow(row)
					startPos = endPos
					val = i
				counter += 1
			if curPos != chromArray[len(chromArray)-1]:
				finalPos = float(curPos+1)
				finalVal = chromArray[-1]

				finalRow = [chrom, int(startPos), int(finalPos), finalVal]
				bgwriter.writerow(finalRow)
	print("--> finished writing bedGraph")

def npz_array_to_bedgraph_and_bw_mp(outDir, chrSizes):
	"""
	Use wangenUtils to write each npz archive to a bedgraph and bigwig file
	"""

	# npArrayFile = "%s/"

	npz_list = glob.glob('%s/*.npz' % (outDir))
	print(npz_list)

	outbg_list = []
	for npz in npz_list:
		outbg = npz[:-4]+'.bedgraph'
		outbg_list.append(outbg)

	bgtuple = list(zip(npz_list,outbg_list))

	# bedgraph_mp_handler(bgtuple)

	# for npz in npz_list:
	# 	outbg = npz[:-4]+'.bedgraph'
	# 	jrw.bedGraph_writer_genome(npArray=npz , outBedGraph=outbg)

	print('starting mp handler for bedgraphs')
	pool = ProcessPool(nodes=len(bgtuple))
	results = pool.map(bedgraph_mp_handler, bgtuple)

	return results

def bedgraph_to_bigwig_conversion(npArray, inBedGraph):
	outBigWig = inBedGraph[:-9]+".bw"
	# print("bigwig location: ", outBigWig )

	### use kentUtils here for now
	bedGraph_to_bigWig = "bedGraphToBigWig %s %s %s" % (
		inBedGraph, chrSizes, outBigWig)
	print(bedGraph_to_bigWig)
	subprocess.Popen(bedGraph_to_bigWig, shell=True).wait()
	remove_bedgraph = 'rm %s' % (inBedGraph)
	print(remove_bedgraph)
	subprocess.Popen(remove_bedgraph, shell=True).wait()

def bedgraph_mp_handler(bgtuple):

	bedGraph_writer_genome(npArray=bgtuple[0], outBedGraph=bgtuple[1])
	bedgraph_to_bigwig_conversion(npArray=bgtuple[0], inBedGraph=bgtuple[1])

def build_coverage_arrays_mp_hanlder(chromList, inBam, outDir, genome, samp):

	npArray = npArrayBuilder(
		inBam = inBam,
		outDir = outDir,
		genome = genome,
		samp= samp
		)


	print('starting mp handler')
	pool = ProcessPool(nodes=len(chromList))
	results = pool.map(npArray.build_coverage_array_singleChrom, chromList)

	return results

def validChroms(genomeFile, genomeFormat='twobit', chromPrefix='chr', chromPrefixChars=3):
	"""
	given a genome in the 2bit format,
		look up valid chromsome names given 'chr' prefix

	return a list of valid chromsome names
	"""
	if genomeFormat == 'twobit':
		genome = twobitreader.TwoBitFile(genomeFile) # create handler to open the 2bit file

		chromList = []

		for chrom in genome:
			if chrom[0:chromPrefixChars] == chromPrefix:
				chromList.append(chrom)

	return chromList


if __name__ == '__main__':
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))

	genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file
	chromList = validChroms(twobitfile)
	inBam = f'{rootDir}/{alignDir}/{genomeName}/star/dedupe/{args.samp}_dedupe.sorted.bam'
	outDir = f'{rootDir}/npArrayDir/{args.samp}'
	if not os.path.exists(outDir):	os.makedirs(outDir)


	build_coverage_arrays_mp_hanlder(chromList, inBam, outDir, genome, args.samp)

	print(chromList)

	# build_coverage_array_singleChrom('chr1', args.inBam, args.outDir, genome, args.samp)

	print('!!! DONE building single arrays!!!')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))
	print(' Concatinating arrays ')

	"""
	secondReadGenome.npy
	covGenome.npy
	secReadGenomePos.npy
	secReadGenomeNeg.npy
	"""

	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'secReadGenome')
	concat_np_array_single(outDir, args.samp, chromList, arraySuffix='secondReadGenome')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'covGenome')
	concat_np_array_single(outDir, args.samp, chromList, arraySuffix='covGenome')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'pos')
	concat_np_array_single(outDir, args.samp, chromList, arraySuffix='secondReadGenomePos')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'neg')
	concat_np_array_single(outDir, args.samp, chromList, arraySuffix='secondReadGenomeNeg')
	print('done concating arrays')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'))

	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'secReadGenome-norm')
	normalize_np_array(outDir, args.samp, chromList, arraySuffix='secondReadGenome')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'secReadGenome-norm-pos')
	normalize_np_array(outDir, args.samp, chromList, arraySuffix='secondReadGenomePos')
	# print(datetime.today().strftime('%Y-%m-%d-%H:%M:%S'), 'secReadGenome-norm-neg')
	normalize_np_array(outDir, args.samp, chromList, arraySuffix='secondReadGenomeNeg')

	npz_array_to_bedgraph_and_bw_mp(outDir, chrSizes)

