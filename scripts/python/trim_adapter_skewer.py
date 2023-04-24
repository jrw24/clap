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

def remove_adapter(samp, threadNumb):
	'''
	*** This function uses skewer to trim a standard illumina adapter for paired end reads ***
	-x specifies the 3' adapter on on the read
	-l specifies the minimum post-trimming read length
	-L specifies the maximum post-trimming read length
	-Q specifies the minimum quality score 
	-m pe specifies paired end trimming - uses default adapter seqs
	For TruSeq Kits:
	Illumina Read 1 Adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	Illumina Read 2 Adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
	'''
	# input fastq files:
	
	if pairedEndReads == True:
		fastq1 = '%s/%s*%s.fastq.gz' % (fastqPath, samp, fqSuffix1)
		fastq2 = '%s/%s*%s.fastq.gz' % (fastqPath, samp, fqSuffix2)

		trimming_out_path = '%s/%s' % (rootDir, readProcDir)
		# if not os.path.exists(trimming_out_path): os.makedirs(trimming_out_path)
		trimming_out_file = '%s/%s' % (trimming_out_path, samp)

		trimming_out_filegz_1 = '%s/%s-trimmed-pair1.fastq.gz' % (trimming_out_path, samp)
		trimming_out_filegz_2 = '%s/%s-trimmed-pair2.fastq.gz' % (trimming_out_path, samp)
		
		# if os.path.isfile(trimming_out_filegz_1) and os.path.isfile(trimming_out_filegz_2):
		# 	print("adapter previously removed")
		# else: 
		# 	# leaving this without a defined adapter for now...
		# 	# skewer is using -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC, -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
		read_trimming = 'skewer -m pe -Q 10 -l 15 --quiet -o %s --threads %s %s %s' % (
						trimming_out_file, threadNumb, fastq1, fastq2)
		print(read_trimming)
		subprocess.Popen(read_trimming, shell=True).wait()
		
		compression_command_1 = 'pigz -p %s %s-trimmed-pair1.fastq' % (threadNumb, trimming_out_file)
		compression_command_2 = 'pigz -p %s %s-trimmed-pair2.fastq' % (threadNumb, trimming_out_file)
		print(compression_command_1)
		subprocess.Popen(compression_command_1, shell=True).wait()
		print(compression_command_2)
		subprocess.Popen(compression_command_2, shell=True).wait()
		print("adapter trimming finished")

def main():
	remove_adapter(args.samp, args.threadNumb)

if __name__ == '__main__':
	main()


