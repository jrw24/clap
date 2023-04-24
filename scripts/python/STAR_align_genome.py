"""
Trim read adapter using skewer

"""

import argparse
import importlib
import subprocess
import os

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

def STAR_align_genome(samp, threadNumb):
	"""
	Align reads to the genome using STAR
	Make sure to allow for soft-clipping to permit alignment of reads with A or T added by superscript III
	45 seems to be the max number of supported threads for alignment
	only allowing single alignments with outFilterMultimapNmax 1

	Not writing wig files here for simplicity - this comes later after numpy array building
	"""
	if pairedEndReads == True:
		star_infile_1 = '%s/%s/%s_no_ncRNA_1.fastq.gz' % (rootDir, readProcDir, samp)
		star_infile_2 = '%s/%s/%s_no_ncRNA_2.fastq.gz' % (rootDir, readProcDir, samp)
		star_out_path = '%s/%s/%s/star/%s' % (rootDir, alignDir, genomeName, samp) # dir for each sample
		star_out_file = '%s/%s' % (star_out_path, samp)
		if not os.path.exists(star_out_path): os.makedirs(star_out_path)
		
		star_command = 'STAR \
						--runThreadN %s \
						--genomeDir %s \
						--readFilesIn %s %s \
						--readFilesCommand gunzip -c \
						--outSAMtype BAM SortedByCoordinate \
						--alignSJDBoverhangMin 1 \
						--alignSJoverhangMin 8 \
						--outFilterType BySJout \
						--outFilterMultimapNmax 1 \
						--outFileNamePrefix %s \
						--outReadsUnmapped Fastx' % (
			threadNumb, starGenome, star_infile_1, star_infile_2, star_out_file)
		print(star_command)
		subprocess.Popen(star_command, shell=True).wait()
		
		## index bam outputs
		index_genome = 'samtools index %sAligned.sortedByCoord.out.bam' % (star_out_file)
		print(index_genome)
		subprocess.Popen(index_genome, shell=True).wait()
		print("STAR alignment has finished successfully \n")

def main():
	STAR_align_genome(args.samp, args.threadNumb)

if __name__ == '__main__':
	main()



