"""
Run Jar file for clap peak calling


"""


"""
JAR arguments;
	args[0] = sample bam
	args[1] = input bam
	args[2] = bin size
	args[3] = genes (bed)
	args[4] = save
	args[5] = reverse strand (optional, default=false)
"""

import argparse
import importlib
import subprocess

parser= argparse.ArgumentParser()
parser.add_argument('--samp', help= 'name of sample being processed')
parser.add_argument('--libset', help= 'path to libsettings file')
parser.add_argument('--scriptDir', help = 'directory with workflow script')
parser.add_argument('--threadNumb', help= 'number of threads to use')
args = parser.parse_args()

### load libset into namespace 
libset = importlib.machinery.SourceFileLoader('libset', args.libset).load_module()
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)

def run_jar(jarFile, capBam, inputBam, binSize, bedFile, jarOut, reverseStand):

	jar_cmnd = f'java -jar {jarFile} {capBam} {inputBam} {binSize} {bedFile} {jarOut}'
	print(jar_cmnd)
	subprocess.Popen(jar_cmnd, shell=True).wait()

def main():
	print('running JAR analysis')

	jarFile = f'{args.scriptDir}/scripts/java/CLAPAnalysis_V7.jar'
	inputSamp = input_sample[0]

	capBam = f'{rootDir}/{alignDir}/{genomeName}/star/dedupe/{args.samp}_dedupe.sorted.bam'
	inputBam = f'{rootDir}/{alignDir}/{genomeName}/star/dedupe/{inputSamp}_dedupe.sorted.bam' 
	clapEnrichDir = f'{rootDir}/clapEnrichment'
	jarOut = f'{clapEnrichDir}/{args.samp}_peaks_jar.tsv'
	reverseStand = False


	run_jar(jarFile, capBam, inputBam, enrichBinSize, bedFile, jarOut, reverseStand)


if __name__ == '__main__':
	main()
