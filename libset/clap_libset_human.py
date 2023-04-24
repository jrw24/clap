"""
Inputs for RNAseq sample processing
"""
# import os

## set paths to input fastq files and output data files
# libsetDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
rootDir = "/groups/guttman/jamie/data/sharp/clap/hek1_bams/clap_main"
fastqPath =  "/groups/guttman/jamie/data/sharp/clap/hek1_bams"
fqSuffix1 = "1"
fqSuffix2 = "2"
#####

### Library input Settings for experiment and samples:
experiment = '230413'

### Annotation Files: 
genomeName = 'hg38_default_ncRNAsub' 	# add this for storing alignment files to unique directories
# starGenome = '/groups/guttman/jamie/genomes/mouse/GRCm38_p6/star_mm10_default'
starGenome = '/groups/guttman/jamie/genomes/human/GRCh38/star_gtf_gencodeV30annotation'
twoBitGenome = '/groups/guttman/jamie/genomes/human/GRCh38/hg38.2bit'
ncRNAstarGenome = '/groups/guttman/jamie/genomes/human/GRCh38/star_hg38_ncRNA'
salmonGenome = '/groups/guttman/jamie/genomes/human/GRCh38/salmon/salmon_index'
trXgene = '/groups/guttman/jamie/genomes/human/GRCh38/gencode.v30.trXgene.csv'
trFile = '/groups/guttman/jamie/genomes/human/GRCh38/gencode.v30.transcripts.fa.gz'
GTFfile = '/groups/guttman/jamie/genomes/human/GRCh38/gencode.v30.annotation.gtf' 
twobitfile = '/groups/guttman/jamie/genomes/human/GRCh38/hg38.2bit'
genome_fasta = '/groups/guttman/jamie/genomes/human/GRCh38/hg38.fa'
chrSizes = "/groups/guttman/jamie/genomes/human/GRCh38/hg38.sizes"
bedFile = '/groups/guttman/jamie/genomes/human/GRCh38/appris_pricipal_gtf/gencode.v30.appris_principal.bed'

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
"SHARP_HEK293T",
]
input_sample = ['Input_HEK293T']
samplelist = capture_samples + input_sample




