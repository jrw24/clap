"""
Deduplicate bam files using picard MarkDuplicates
	
	Alternate options should be used ideally if reads have UMIs

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


def deduplicate_bams(samp, threadNumb):
	"""
	use picard's MarkDuplicates function to remove duplicate reads from aligned bam files
	"""

	picardPath = "/groups/guttman/software/picard.2.18.7/picard.jar" # path to Picard

	bamOutDir = "%s/%s/%s/star/dedupe" % (rootDir, alignDir, genomeName)

	bamIn = "%s/%s/%s/star/%s/%sAligned.sortedByCoord.out.bam" % (
		rootDir, alignDir, genomeName, samp, samp)

	bamOut = "%s/%s_dedupe.bam" % (bamOutDir, samp)
	metricFile = "%s/%s_dedupeMetrics.txt" % (bamOutDir, samp)

	dedupe_cmnd = "java -jar %s MarkDuplicates REMOVE_DUPLICATES=true I=%s O=%s M=%s" % (
		picardPath, bamIn, bamOut, metricFile)

	print(dedupe_cmnd)
	subprocess.Popen(dedupe_cmnd, shell=True).wait()

	bamSorted = "%s/%s_dedupe.sorted.bam" % (bamOutDir, samp)
	sort_cmnd = "samtools sort -@ %s -o %s %s" % (threadNumb, bamSorted, bamOut)
	index_cmnd = "samtools index %s" % (bamSorted)
	rm_cmnd = "rm %s" % (bamOut) ## clean up

	print(sort_cmnd)
	subprocess.Popen(sort_cmnd, shell=True).wait()
	print(index_cmnd)
	subprocess.Popen(index_cmnd, shell=True).wait()
	print(rm_cmnd)
	subprocess.Popen(rm_cmnd, shell=True).wait()

def main():
	deduplicate_bams(args.samp, args.threadNumb)

if __name__ == '__main__':
	main()
