"""
testing for clap_enrichments

"""

## import requried packages here ###
import os
import sys
import subprocess 
import importlib

# run_cmnd

scriptDir = '/home/jwangen/scripts/clap/clap_main_jrw'
sys.path.append(f'{scriptDir}/libset')
from clap_libset import * ## import libsettings from py file into namespace
libset = f'{scriptDir}/libset/clap_libset.py'

print(rootDir)
print(samplelist)
print(libset)

samp = samplelist[-1]
threadNumb = '1'

run_enrich = f'python {scriptDir}/scripts/python/calculate_clap_enrichment.py \
	--samp {samp} --libset {libset} --binSize {enrichBinSize} --threadNumb {threadNumb}'

print(run_enrich)

subprocess.Popen(run_enrich, shell=True).wait()