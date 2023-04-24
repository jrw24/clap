"""
Inputs for RNAseq sample processing
"""
# import os

## set paths to input fastq files and output data files
# libsetDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
rootDir = "/groups/guttman/jamie/data/sharp/clap/plathClap"
fastqPath =  "%s/fastq" % (rootDir)
fqSuffix1 = "1"
fqSuffix2 = "2"
#####

### Library input Settings for experiment and samples:
experiment = '220417'

### Annotation Files: 
genomeName = 'mm10_default_ncRNAsub' 	# add this for storing alignment files to unique directories
# starGenome = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/star_mm10_default'
starGenome = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/star_mm10_default'
twoBitGenome = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/mm10.2bit'
ncRNAstarGenome = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/star_mm10_ncRNA'
salmonGenome = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/salmon/salmon_index'
trXgene = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/gencode.vM25.trXgene.csv'
trFile = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/gencode.vM25.transcripts.fa.gz'
GTFfile = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/gencode.vM25.annotation.gtf' 
twobitfile = "/groups/guttman/jamie/genomes/mouse/GRCm38_p6/GRCm38.primary_assembly.genome.2bit"
genome_fasta = "/groups/guttman/jamie/genomes/mouse/GRCm38_p6/GRCm38.primary_assembly.genome.fa"
chrSizes = "/groups/guttman/jamie/genomes/mouse/GRCm38_p6/GRCm38.primary_assembly.genome.sizes"
bedFile = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/appris_pricipal_gtf/gencode.vM25.appris_principal_expanded.bed'


### Alignment variables
pairedEndReads = True
multiMap = False
readProcDir = 'read_proc_ncSub'
alignDir = 'alignments_ncSub'
enrichBinSize = '100' ## bin size to use for calculating genome enrichments
"""
samples MUST have 3 "_"'s in the name to be compatible with workflow'
schema == number_name_condition_replicate   ex: 1_wt_ut_1
"""
capture_samples = [ 
'1_sharp_flSHARP_A',
'2_sharp_flSHARP_B',
'3_sharp_dRRM_B',
'4_sharp_dIDR_A',
'5_sharp_dIDR_B'
]
input_sample = ['input_merged']
samplelist = capture_samples + input_sample



