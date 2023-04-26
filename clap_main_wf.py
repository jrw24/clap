"""
Main workflow for running analysis on a CLAP dataset
"""

### import requried packages here ###
import os
import sys
import subprocess 
import importlib
import yaml
from simple_slurm import Slurm 

scriptDir = os.path.dirname(os.path.realpath(__file__)) ## directory containing this script
sys.path.append("%s/libset" % scriptDir) ## add libsettings dir to path
from clap_libset import * ## import libsettings from py file into namespace
libset = f'{scriptDir}/libset/clap_libset.py'

print(rootDir)
print(samplelist)

"""
Data Processing Steps:

To Do List:
	- temp files to scratch directory
	- standardize standard out and error messages
	- get conda environment working
	- save scripts to analysi directory
		- done, but should build in a check to see if scripts already exist
	- peak calls may not work if input sample is not done... 
		may need to restrategize this function 

I) Preprocessing
	- install required python packages using conda and .yml file
		Python Packages: (as of April 2023)
			- Python version: 3.8.5

			- numpy 1.21.2
			- pandas 1.3.3
			- pysam 0.17.0
			- twobitreader 3.1.7
			- pathos 0.2.7
			- scipy 
			- simple_slurm 
			...

		Command line tools:
			- pigz 
			- skewer 
			- STAR 
			- samtools 
			- picard '/groups/guttman/software/picard.2.18.7/picard.jar'
			- kentUtils bedGraphToBigWig 
			...

	- generate noncoding RNA genome and STAR indexes
	- generate appris-principal BED file from gencode GTF file
	- 


II) Read processing
	1) Trim adapter using skewer
	2) Subtract noncoding RNA using STAR
	3) Align reads to genome using STAR
	4) Deduplicate reads using Picard
	5) Build genome Numpy Arrays and bigwigs
	6) 

"""


def make_dir(outDir):
	if not os.path.exists(outDir):	os.makedirs(outDir)

def make_out_dirs(rootDir, readProcDir, alignDir, genomeName):
	"""
	make output directories if they do not exist
	"""
	scratchDir = '/central/scratch/%s/%s/%s' % (os.environ['USER'], rootDir, experiment) 
	readProcOutDir = f'{rootDir}/{readProcDir}' 
	alignOutDir = f'{rootDir}/{alignDir}/{genomeName}'
	starOutDir = f'{alignOutDir}/star'
	dedupeDir = f'{starOutDir}/dedupe'
	nascentDir = f'{starOutDir}/nascent'
	jobDir = f'{rootDir}/jobs'
	logDir = f'{rootDir}/logs'
	npArrayDir = f'{rootDir}/npArrayDir'
	clapEnrichDir = f'{rootDir}/clapEnrichment'
	outScriptDir = f'{rootDir}/workflow'

	make_dir(scratchDir)
	make_dir(readProcOutDir)
	make_dir(alignOutDir)
	make_dir(starOutDir)
	make_dir(dedupeDir)
	make_dir(nascentDir)
	make_dir(jobDir)
	make_dir(logDir)
	make_dir(npArrayDir)
	make_dir(clapEnrichDir)
	make_dir(outScriptDir)

def copy_scripts(rootDir, scriptDir):
	"""
	create a copy of all scripts using during the analysis for reference
	
	scriptDir is where this script is run 
	rootDir is where all output files will be written

	should build in a check to see if scripts exist probably
	"""

	cp_cmnd = f'cp -rf {scriptDir} {rootDir}/workflow'
	subprocess.Popen(cp_cmnd, shell=True).wait()


def set_slurm_settings(samp, logDir, slurm_settings_yaml):
	"""
	given a sample and slurm settings in a yaml file
		submit a job for the sample
	"""

	slurm_settings = yaml.safe_load(open(slurm_settings_yaml))
	slurm_settings['output'] = slurm_settings['output'].replace('{logDir}', logDir)
	slurm_settings['output'] = slurm_settings['output'].replace('{samp}', samp)
	slurm_settings['error'] = slurm_settings['error'].replace('{logDir}', logDir)
	slurm_settings['error'] = slurm_settings['error'].replace('{samp}', samp)

	slurm_cmnd = Slurm()
	for key, val in slurm_settings.items():
		slurm_cmnd.add_arguments(key, val)
	return slurm_cmnd

def write_slurm_job(samp, jobfile, slurm_settings, input_command):
	"""
	once slurm settings have been assigned, write the jobfile to output file
	"""
	jobDir = f'{rootDir}/jobs'
	jobOutFile = f'{jobDir}/{samp}_{jobfile}.job.sh'

	with open(jobOutFile, 'w') as f:
		f.write(slurm_settings.arguments())
		f.write(input_command+'\n')

def ntasks_lookup(slurm_set):
	"""
	Given a simple_slurm 
	"""

	slurmArgs = slurm_set.arguments() ## string of arguments
	slurmArgsList = slurmArgs.split('\n')
	for i in slurmArgsList:
		if i[0:16] == '#SBATCH --ntasks':
			taskn = i[-1]
	return taskn

def slurm_single_sample(samp):
	"""
	take a single sample, passed as an object to this function

	1) set python script settings
	2) define slurm commands
	3) precheck for output files 
	4) set slurm dependencies
	5) submit sbatch jobs

	"""
	print(samp)

	commands = dict()
	commands['trim_adapter_cmnd'] = f'python \
		{scriptDir}/scripts/python/trim_adapter_skewer.py --samp {samp} --libset {libset} --threadNumb ' 
	commands['ncRNA_subtract_cmnd'] = f'python \
		{scriptDir}/scripts/python/ncRNA_subtract_STAR.py --samp {samp} --libset {libset} --threadNumb ' 
	commands['star_align_cmnd'] = f'python \
		{scriptDir}/scripts/python/STAR_align_genome.py --samp {samp} --libset {libset} --threadNumb ' 
	commands['deduplicate_bams_cmnd'] = f'python \
		{scriptDir}/scripts/python/deduplicate_bams.py --samp {samp} --libset {libset} --threadNumb ' 
	commands['split_nascent_cmnd'] = f'python \
		{scriptDir}/scripts/python/split_nascent_bams.py --samp {samp} --libset {libset} --threadNumb '
	commands['numpy_genome_arrays'] = f'python \
		{scriptDir}/scripts/python/numpy_genome_arrays.py --samp {samp} --libset {libset} --threadNumb '
	# commands['calculate_clap_enrichment'] = f'python \
	# 	{scriptDir}/scripts/python/calculate_clap_enrichment.py --samp {samp} --libset {libset} \
	# 	--binSize {enrichBinSize} --threadNumb '
	# commands['run_clap_jar'] = f'python \
	# 	{scriptDir}/scripts/python/run_clap_jar.py --samp {samp} --libset {libset} \
	# 	--scriptDir {scriptDir} --threadNumb '

	### define settings for slurm commands ###
	slurmSetDir = f'{scriptDir}/scripts/slurm/yaml' ## dir with yaml files
	jobDir = f'{rootDir}/jobs' ## output of jobfiles, and error logs
	logDir = f'{rootDir}/logs'

	### read in slurm settings from yaml files ###
	trim_adapter_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_trim_adapter_skewer.yaml')
	trim_task_numb = ntasks_lookup(trim_adapter_slurm)

	ncRNA_subtract_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_ncRNA_subtract_STAR.yaml')
	ncRNA_task_numb = ntasks_lookup(ncRNA_subtract_slurm)

	star_align_genome_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_STAR_align_genome.yaml')
	star_task_numb = ntasks_lookup(star_align_genome_slurm)

	dedeuplicate_bams_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_deduplicate_bams.yaml')
	dedupe_tasks_numb = ntasks_lookup(dedeuplicate_bams_slurm)

	split_nascent_bams_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_split_nascent.yaml')
	nascnet_tasks_numb = ntasks_lookup(split_nascent_bams_slurm)

	numpy_arrays_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_numpy_genome_arrays.yaml')
	npArray_tasks_numb = ntasks_lookup(numpy_arrays_slurm)

	numpy_nascent_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_numpy_genome_arrays_nascent.yaml')
	numpy_mature_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_numpy_genome_arrays_mature.yaml')

	# calc_clap_enrich_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_calc_clap_enrichments.yaml')
	# calcEnrich_tasks_numb = ntasks_lookup(calc_clap_enrich_slurm)

	# run_clap_jar_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_clapJar.yaml')
	# clap_jar_tasks_numb = ntasks_lookup(run_clap_jar_slurm)

	### modify commands to add thread number ###
	commands['trim_adapter_cmnd'] = commands['trim_adapter_cmnd']+trim_task_numb
	commands['ncRNA_subtract_cmnd'] = commands['ncRNA_subtract_cmnd']+ncRNA_task_numb
	commands['star_align_cmnd'] = commands['star_align_cmnd']+star_task_numb
	commands['deduplicate_bams_cmnd'] = commands['deduplicate_bams_cmnd']+dedupe_tasks_numb
	commands['split_nascent_cmnd'] = commands['split_nascent_cmnd']+nascnet_tasks_numb
	commands['numpy_genome_arrays'] = commands['numpy_genome_arrays']+npArray_tasks_numb
	commands['numpy_genome_arrays_nascent'] = commands['numpy_genome_arrays']+' --split nascent'
	commands['numpy_genome_arrays_mature'] = commands['numpy_genome_arrays']+' --split mature'


	# commands['calculate_clap_enrichment'] = commands['calculate_clap_enrichment']+calcEnrich_tasks_numb
	# commands['run_clap_jar'] = commands['run_clap_jar']+clap_jar_tasks_numb

	### submit sbatch jobs to scheduler ###
	jobID_trim_adapt = trim_adapter_slurm.sbatch(commands['trim_adapter_cmnd']) ## submit the job!
	write_slurm_job(samp, 'trim_adapter', trim_adapter_slurm, commands['trim_adapter_cmnd'])

	ncRNA_subtract_slurm.set_dependency('afterok:%s' % (jobID_trim_adapt))
	jobID_ncRNA_subtract = ncRNA_subtract_slurm.sbatch(commands['ncRNA_subtract_cmnd'])
	write_slurm_job(samp, 'ncRNA_subtract', ncRNA_subtract_slurm, commands['ncRNA_subtract_cmnd'])
	
	star_align_genome_slurm.set_dependency('afterok:%s' % (jobID_ncRNA_subtract))
	jobID_star_align_genome = star_align_genome_slurm.sbatch(commands['star_align_cmnd'])
	write_slurm_job(samp, 'star_align_genome', star_align_genome_slurm, commands['star_align_cmnd'])

	dedeuplicate_bams_slurm.set_dependency('afterok:%s' % (jobID_star_align_genome))
	jobID_deduplicate_bams = dedeuplicate_bams_slurm.sbatch(commands['deduplicate_bams_cmnd'])
	write_slurm_job(samp, 'dedeuplicate_bams', dedeuplicate_bams_slurm, commands['deduplicate_bams_cmnd'])

	split_nascent_bams_slurm.set_dependency('afterok:%s' % (jobID_deduplicate_bams))
	jobID_split_nascent = split_nascent_bams_slurm.sbatch(commands['split_nascent_cmnd'])
	write_slurm_job(samp, 'split_nascent_bams', split_nascent_bams_slurm, commands['split_nascent_cmnd'])

	numpy_arrays_slurm.set_dependency('afterok:%s' % (jobID_deduplicate_bams))
	jobID_numpy_arrays = numpy_arrays_slurm.sbatch(commands['numpy_genome_arrays'])
	write_slurm_job(samp, 'numpy_genome_arrays', numpy_arrays_slurm, commands['numpy_genome_arrays'])

	## nascent
	numpy_nascent_slurm.set_dependency('afterok:%s' % (jobID_split_nascent))
	jobID_numpy_nascent = numpy_nascent_slurm.sbatch(commands['numpy_genome_arrays_nascent'])
	write_slurm_job(samp, 'numpy_nascent_arrays', numpy_nascent_slurm, commands['numpy_genome_arrays_nascent'])

	## mature
	numpy_mature_slurm.set_dependency('afterok:%s' % (jobID_numpy_nascent)) ## run last
	jobID_numpy_mature = numpy_mature_slurm.sbatch(commands['numpy_genome_arrays_mature'])
	write_slurm_job(samp, 'numpy_mature_arrays', numpy_mature_slurm, commands['numpy_genome_arrays_mature'])

	return jobID_numpy_mature


	#### Requires input as well #### 


	# if samp not in input_sample:
	# 	# calc_clap_enrich_slurm.set_dependency('afterok:%s' % (jobID_numpy_arrays))
	# 	jobID_calc_enrich = calc_clap_enrich_slurm.sbatch(commands['calculate_clap_enrichment'])
	# 	write_slurm_job(samp, 'calculate_clap_enrichment', calc_clap_enrich_slurm, commands['calculate_clap_enrichment'])
	# 	return jobID_calc_enrich

	# 	# run_clap_jar_slurm.set_dependency('afterok:%s' % (jobID_numpy_arrays))
	# 	jobID_clap_jar = run_clap_jar_slurm.sbatch(commands['run_clap_jar'])
	# 	write_slurm_job(samp, 'run_clap_jar', run_clap_jar_slurm, commands['run_clap_jar'])

def slurm_capture_only(samp, jobList):

	"""
	analysis for capture compared to input_sample
		requires both capture and input to have finished processing before running

	set jobList = 0 to ignore afterok settings

	jobList contains jobID's of final processing jobs,
	use these for afterok
	"""
	slurmSetDir = f'{scriptDir}/scripts/slurm/yaml' ## dir with yaml files
	jobDir = f'{rootDir}/jobs' ## output of jobfiles, and error logs
	logDir = f'{rootDir}/logs'

	if not jobList == 0:
		jl = [str(x) for x in jobList]
		jobOkList = ':'.join(jl)

	### set commands to run
	commands = dict()
	commands['calculate_clap_enrichment'] = f'python \
		{scriptDir}/scripts/python/calculate_clap_enrichment.py --samp {samp} --libset {libset} \
		--binSize {enrichBinSize} --threadNumb '
	commands['run_clap_jar'] = f'python \
		{scriptDir}/scripts/python/run_clap_jar.py --samp {samp} --libset {libset} \
		--scriptDir {scriptDir} --threadNumb '

	### define slurm settings for each command
	calc_clap_enrich_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_calc_clap_enrichments.yaml')
	calcEnrich_tasks_numb = ntasks_lookup(calc_clap_enrich_slurm)

	run_clap_jar_slurm = set_slurm_settings(samp, logDir, f'{slurmSetDir}/slurm_clapJar.yaml')
	clap_jar_tasks_numb = ntasks_lookup(run_clap_jar_slurm)

	### add ntasks for commands
	commands['calculate_clap_enrichment'] = commands['calculate_clap_enrichment']+calcEnrich_tasks_numb
	commands['run_clap_jar'] = commands['run_clap_jar']+clap_jar_tasks_numb

	### submit jobs using sbatch
	if not jobList == 0:
		calc_clap_enrich_slurm.set_dependency('afterok:%s' % (jobOkList))
	jobID_calc_enrich = calc_clap_enrich_slurm.sbatch(commands['calculate_clap_enrichment'])
	write_slurm_job(samp, 'calculate_clap_enrichment', calc_clap_enrich_slurm, commands['calculate_clap_enrichment'])

	if not jobList == 0:
		run_clap_jar_slurm.set_dependency('afterok:%s' % (jobOkList))
	jobID_clap_jar = run_clap_jar_slurm.sbatch(commands['run_clap_jar'])
	write_slurm_job(samp, 'run_clap_jar', run_clap_jar_slurm, commands['run_clap_jar'])


def main():

	make_out_dirs(rootDir, readProcDir, alignDir, genomeName)
	copy_scripts(rootDir, scriptDir)

	### processing of all samples

	jobList = []
	for samp in samplelist:
		print(samp, 'starting')
		finalJob = slurm_single_sample(samp)
		jobList.append(finalJob)

	print(jobList)

	for samp in capture_samples:
		slurm_capture_only(samp, jobList=jobList)

	### processing of capture samples only

if __name__ == '__main__':
	main()

