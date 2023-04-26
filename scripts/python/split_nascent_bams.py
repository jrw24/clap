"""
given a deduplicated bamfile,
	split the reads into nascent and matrue
		"nascent" reads overlap an intron in either read1 or read2
		"mature" reads fully map to exons only

"""

import argparse
import importlib
import subprocess
import sys
import os
import ast
import pysam
import twobitreader
from Bio import Seq
import pandas as pd
pd.set_option('display.max_columns', None)
import numpy as np
from collections import OrderedDict
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

### define inputs and outputs
bam = "%s/%s/%s/star/dedupe/%s_dedupe.sorted.bam" % (rootDir, alignDir, genomeName, args.samp)
outDir = '%s/%s/%s/star/nascent' % (rootDir, alignDir, genomeName)
outBamPrefix = '%s/%s.split' % (outDir, args.samp)


### load in genome files:
def load_genomes(twobitfile, bedInput):
	genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file
	b = pd.read_csv(bedInput, sep="\t", header=None)
	b.columns = [   'chrom', 
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

	return genome, b

### define exon arrays:
def build_exon_array(chooseGene, bedDF, genome):
	"""
	output dictionaries with geneFeatures, exon positions, and intron positions 
	"""
	g = chooseGene
	df = bedDF.loc[bedDF['gene']==g] ## slice the bed data.frame to just the chosen gene
	print(df)
	if len(df) > 1: ## if multiple entries for gene
		df = df.iloc[[0]] ## just take the first entry

	### define gene specific variables
	s = int(df['chromStart'])
	e = int(df['chromEnd'])
	chrom = df['chrom'].to_list()[0]
	blockStarts = ast.literal_eval(df['blockStarts'].to_list()[0])
	blockSizes = ast.literal_eval(df['blockSizes'].to_list()[0])
	strand = df['strand'].to_list()[0]

	### convert blockStarts and blockSizes to tuples if not already tuples
	### handles the case where blockStarts = 0

	if not type(blockStarts) == tuple:
		blockStarts = tuple([blockStarts])
	if not type(blockSizes) == tuple:
		blockSizes = tuple([blockSizes])


	### check that each exon has an entry for sizes
	# assert len(blockStarts) == len(blockSizes), "check exon lenghts are identical"

	### define genomic values for chosen gene
	genomic_length = e-s
	gseq = genome[chrom][s:e] ## retrieve the genomic sequence information
	gseq_rc = str(Seq.Seq(gseq).reverse_complement())

	### define features of the gene and store in a dict
	geneFeatures = dict()
	geneFeatures['start'] = s
	geneFeatures['end'] = e 
	geneFeatures['chrom'] = chrom 
	geneFeatures['strand'] = strand
	geneFeatures['genomicLength'] = genomic_length
	geneFeatures['genomicSeq'] = gseq 
	geneFeatures['genomicSeqRC'] = gseq_rc
	geneFeatures['gene'] = chooseGene

	### transcript relative coordinates

	exonRegions = []
	intronRegions = []

	### for positive strand

	## Define Exons
	for i in range(len(blockStarts)):

		singleExon = [int(blockStarts[i]), int(blockStarts[i]+blockSizes[i])]
		exonRegions.append(singleExon)

	exonCount = len(exonRegions)

	## Define Introns
	for i in range(len(exonRegions)-1):

		intStart = exonRegions[i][1]
		intEnd = exonRegions[i+1][0]
		intronRegions.append([intStart, intEnd])

	## Build regions on the negative strand
	exonRegionsNeg = []
	intronRegionsNeg = []

	for i in exonRegions[::-1]:
		# print([abs(i[1]-len(gseq)), abs(i[0]-len(gseq))])
		exonRegionsNeg.append([abs(i[1]-len(gseq)), abs(i[0]-len(gseq))])
	for i in intronRegions[::-1]:
		# print([abs(i[1]-len(gseq)), abs(i[0]-len(gseq))])
		intronRegionsNeg.append([abs(i[1]-len(gseq)), abs(i[0]-len(gseq))])

	return exonRegions, intronRegions, exonRegionsNeg, intronRegionsNeg, geneFeatures



def define_nascent_mature_reads(genome, bedDF, bam):

	"""
	iterate through every gene in the annotation file and extract all reads that are found within the gene
		filter reads to ensure they are on the correct strand, have a valid mate, ect.

	for each mate pair, compute the coverage values along the pre-mRNA transcript
		store this numpy coverage array in a dictionary
		each gene will have dicttionary, with an entry for each valid mate-pair 

	then compute the coverage of each mate pair within the introns of the gene
		if this is greater than zero, add the read name to the nascent(intron containing) list
		if reads only span exons, add this to the mature (exon only) list

	return the list of these genes so that they can be written to new output files
	"""
	
	### open bamfile for reading
	bamfile = pysam.AlignmentFile(bam, "rb")

	### create output dictionaries
	nascentRNA = OrderedDict()
	matureRNA = OrderedDict()

	wrong_strand_reads = 0
	valid_reads = 0
	reads_out_of_bounds = 0

	for gene in bedDF['gene'].tolist():

		exonRegions, intronRegions, exonRegionsNeg, intronRegionsNeg, geneFeatures = build_exon_array(gene, bedDF, genome)

		geneReads = OrderedDict()

		for read in bamfile.fetch(geneFeatures['chrom'], geneFeatures['start'], geneFeatures['end']): 
			
			### basic quality control filtering
			if not read.is_proper_pair:
				continue
			if read.is_qcfail:
				continue

			### check that reads are on the appropriate strand
			### 2nd reads should be same as genes, 1st reads should be opposite for CLAP
			if read.is_read2: #correct strand for clap data is read2 alignments
				if read.is_reverse:
					read_strand = "-"
				else:
					read_strand = "+"
				if not read_strand == geneFeatures['strand']:
					wrong_strand_reads +=1
					continue
			else:
				if read.is_reverse:
					read_strand = "-"
				else:
					read_strand = "+"
				if read_strand == geneFeatures['strand']:
					wrong_strand_reads +=1
					continue
			
			### add the read name do read dictionary, and instantiate a numpy vector the length of the pre-mRNA
			if read.query_name not in geneReads:
				geneReads[read.query_name] = np.zeros(geneFeatures['genomicLength'])
			
			### identify the starting position of the read on the 5' end
				### this is valid irrespective of the strand of the read
			readGenomicCoord = read.reference_start ## this is 0 based, add 1 for gff coords
			
			### discard read if it falls outside the bounds of the annotated gene
			if (readGenomicCoord < geneFeatures['start']) or (readGenomicCoord >= geneFeatures['end']):
				reads_out_of_bounds +=1
				continue

			### calculate local (transcript relative) positions for each read 
			### use the length of the read along with the cigar strings to add densities of 1.0 for each position along transcript
			tr_coord = readGenomicCoord-geneFeatures['start']
			counter = 0
			read_len = len(read.seq)
			for tup in read.cigartuples:
				if tup[0] == 0: ## match for the reference
					for i in range(tup[1]):
						if tr_coord+counter < len(geneReads[read.query_name]): ## prevent reads that run off end of transcript
							geneReads[read.query_name][tr_coord+counter]+=1
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

			#### coverage above
			valid_reads += 1

		### assign each read pair to either the nascentRNA or the matureRNA dictionary
		for r in geneReads:
			intronReadCoverage = 0.0

			for intron in intronRegions:
				intrCov = geneReads[r][intron[0]:intron[1]] ## take slice for each intron
				intronReadCoverage += intrCov.sum() ## add the sum of he values to the intron coverage

			if intronReadCoverage > 0.0:
				nascentRNA[r] = intronReadCoverage
			elif intronReadCoverage == 0.0:
				matureRNA[r] = intronReadCoverage
			else:
				print("coverage not assigned for %s" % (r))
	# print("nascent: ", nascentRNA)
	# print("mature: ", matureRNA)
	return nascentRNA, matureRNA

def write_output_bams(bam, nascentRNA, matureRNA, outBamPrefix):

	### open bamfile for reading
	bamfile = pysam.AlignmentFile(bam, "rb")
	nascentBam = pysam.AlignmentFile("%s.nascent.bam" % outBamPrefix, "wb", template=bamfile)
	matureBam = pysam.AlignmentFile("%s.mature.bam" % outBamPrefix, "wb", template=bamfile)

	nascent_reads = 0
	mature_reads = 0
	unassigned_reads = 0

	for read in bamfile.fetch():
		if read.query_name in nascentRNA:
			nascentBam.write(read)
			nascent_reads += 1
		elif read.query_name in matureRNA:
			matureBam.write(read)
			mature_reads += 1
		else:
			unassigned_reads += 1

	bamfile.close()
	nascentBam.close()
	matureBam.close()

	sort_cmnd_nascent = "samtools sort -@ 4 -o %s.nascent.sorted.bam %s.nascent.bam" % (outBamPrefix, outBamPrefix)
	index_cmnd_nascent = "samtools index %s.nascent.sorted.bam" % (outBamPrefix)

	subprocess.Popen(sort_cmnd_nascent, shell=True).wait()
	subprocess.Popen(index_cmnd_nascent, shell=True).wait()

	sort_cmnd_mature = "samtools sort -@ 4 -o %s.mature.sorted.bam %s.mature.bam" % (outBamPrefix, outBamPrefix)
	index_cmnd_mature = "samtools index %s.mature.sorted.bam" % (outBamPrefix)

	subprocess.Popen(sort_cmnd_mature, shell=True).wait()
	subprocess.Popen(index_cmnd_mature, shell=True).wait()

	### remove unsorted bams:

	rm_nascent = 'rm %s.nascent.bam' % (outBamPrefix)
	rm_mature = 'rm %s.mature.bam' % (outBamPrefix)

	subprocess.Popen(rm_nascent, shell=True).wait()
	subprocess.Popen(rm_mature, shell=True).wait()

	print("nascent reads = %s" % nascent_reads)
	print("mature reads = %s" % mature_reads)
	print("unassigned_reads = %s" % unassigned_reads)



def main():
	sys.stdout.write(str(datetime.now())+'\n')
	genome, bedDF = load_genomes(twobitfile, bedFile)
	nascentRNA, matureRNA = define_nascent_mature_reads(genome, bedDF, bam)
	write_output_bams(bam, nascentRNA, matureRNA, outBamPrefix)
	sys.stdout.write(str(datetime.now())+'\n')



if __name__ == '__main__':
	main()

