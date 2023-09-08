#!/usr/bin/env python3

# Script to automate the download of refseq using  sample-taxon tuples ( for ANI calculation)
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import yaml
import itertools
import subprocess
import argparse
import os, gzip,shutil

"""
Business logic
This program downloads resef genomes per sample using tuple sampleId_taxonId as input
The input is a text file issued from Kraken species ID
Usage:
 python get_sampleTaxonId.py --sample_taxon 

"""

parser = argparse.ArgumentParser(prog = 'refseqDwld_sampleTaxonId.py [Options]')
#parser.add_argument('--sample_file', type = argparse.FileType('r'), help = 'provide the sample file or its path', required=True)
parser.add_argument('--sample_taxon', type = argparse.FileType('r'), help =' provide the species name', required=True)

args = parser.parse_args()

samplesTaxId = args.sample_taxon

# Opening file tuple sample-taxonId

samptax_list = samplesTaxId.readlines() # list of tuples
samptax_list =[line.rstrip('\n') for line in samptax_list]
for sgt in samptax_list:
    tp= eval(sgt) # treating list each item as tuple
    #Sample directories respective to each taxon genomes
    subprocess.run('mkdir -p '+ tp[0], shell = True, check = True)
    subprocess.run('ncbi-genome-download bacteria --species-taxids ' + tp[1] + ' -F fasta -o ' + tp[0] + ' --parallel 16 ', shell = True, check = True)
    subprocess.run('mv ' + tp[0] + '/refseq/bacteria/*/*.fna.gz ' + tp[0] + '/', shell = True, check = True)                                                                                       
    subprocess.run("find "+tp[0]+ "/ -iname '*fna*' -exec rename .fna.gz .fasta.gz '{}' \;", shell = True, check = True)
    subprocess.run("gzip -d " + tp[0]+ "/*.fasta.gz",shell = True, check = True)  
    #subprocess.run('rm -rf ' + tp[0] + '/refseq/', shell = True, check = True)
    #cd tp[0] && gzip -d $(find ./ -type f -name '*.gz') a way to unzip
    shutil.rmtree(tp[0] + '/refseq', ignore_errors=True) # remove subdirectories refseq
    # Decompressing or unzip gz files
    for item in os.listdir(tp[0]): # loop through items in dir
      if item.endswith(".gz"): # check for ".gz" extension
          gz_name = os.path.abspath(item) # get full path of files
          file_name = (os.path.basename(gz_name)).rsplit('.',1)[0] #get file name for file within
          with gzip.open(gz_name,"rb") as f_in, open(file_name,"wb") as f_out:
              shutil.copyfileobj(f_in, f_out)
          os.remove(gz_name) # delete zipped file                                                                                           
#gzip -d *.fasta.gz  
