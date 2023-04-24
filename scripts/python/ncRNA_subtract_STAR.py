"""
Trim read adapter using skewer

"""

import argparse
import importlib
import subprocess

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

def ncRNA_subtract_STAR(samp, threadNumb):
	"""
	Deplete noncoding RNA using STAR prior to mapping reads to genome
	using --alignIntronMax 1 to prevent splicing
	set to 2 mismatches using --outFilterMismatchNmax 
	"""

	ncRNA_subSTAR_outpath = '%s/%s' % (rootDir, readProcDir)
	ncRNA_subSTAR_outfile = '%s/%s' % (ncRNA_subSTAR_outpath, samp)

	fq1_in = "%s/%s/%s-trimmed-pair1.fastq.gz" % (rootDir, readProcDir, samp)
	fq2_in = "%s/%s/%s-trimmed-pair2.fastq.gz" % (rootDir, readProcDir, samp)
	
	ncRNA_sub_star = 'STAR \
					  --runThreadN %s \
					  --genomeDir %s \
					  --readFilesIn %s %s \
					  --readFilesCommand gunzip -c \
					  --outReadsUnmapped Fastx \
					  --limitBAMsortRAM 20000000000 \
					  --outSAMtype BAM SortedByCoordinate \
					  --outFilterMismatchNmax 2 \
					  --outFilterMultimapNmax 100 \
					  --alignIntronMax 1 \
					  --outWigType wiggle read1_5p \
					  --outFileNamePrefix %s' % (
		threadNumb, ncRNAstarGenome, fq1_in, fq2_in, ncRNA_subSTAR_outfile)

	print(ncRNA_sub_star)
	subprocess.Popen(ncRNA_sub_star, shell=True).wait()
	
	index_command = 'samtools index %sAligned.sortedByCoord.out.bam' % (ncRNA_subSTAR_outfile)
	subprocess.Popen(index_command, shell=True).wait()

	compression_command = 'pigz -p %s -c %sUnmapped.out.mate1 > %s_no_ncRNA_1.fastq.gz' % (threadNumb, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
	compression_command2 = 'pigz -p %s -c %sUnmapped.out.mate2 > %s_no_ncRNA_2.fastq.gz' % (threadNumb, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
	subprocess.Popen(compression_command, shell=True).wait()
	subprocess.Popen('rm %sUnmapped.out.mate1' % (ncRNA_subSTAR_outfile), shell=True).wait()
	subprocess.Popen(compression_command2, shell=True).wait()
	subprocess.Popen('rm %sUnmapped.out.mate2' % (ncRNA_subSTAR_outfile), shell=True).wait() 

def main():
	ncRNA_subtract_STAR(args.samp, args.threadNumb)

if __name__ == '__main__':
	main()


