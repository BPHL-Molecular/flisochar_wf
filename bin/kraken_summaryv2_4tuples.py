#!/usr/bin/env python3
# 
# Script to parse kraken output
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import yaml
import itertools
import subprocess
import argparse
import os
import glob
import subprocess # copy file
import shutil
parser = argparse.ArgumentParser(prog = 'kraken_summaryv2_3tuples.py [Options]')
#parser.add_argument('--kreport_file', type = argparse.FileType('r'), help = 'provide the sample file or its path', required=True)
parser.add_argument('--reppath', type=str,help= 'paste path to sample.report file',required=True)
parser.add_argument('--out', type=str,help= 'paste txt file to write ID report file',required=True)
#parser.add_argument('--tuplefile', type=str,help= 'paste file to write the tuple sample-specName file',required=True )
parser.add_argument('--tupleGfile', type=str,help= 'paste file to write the tuple sample-Genus file',required=True )
#parser.add_argument('--threads', default = 8, dest = 'threads', help =' specify the number of threads, (default:%(default)s)')
#threads = str(args.threads)
args = parser.parse_args()
kreport_path = args.reppath
# May use either one based on accuracy of specied ID
#specl=args.tuplefile # sample_species
specl_g=args.tupleGfile # sample-genus
# Species list

#speclist=open(specl, 'w')
samp_genlist=open(specl_g, 'w')
# adding report to a file
out_rep=args.out

#Report file
report = open(out_rep, 'w') # user chose a name for the file
header = ['SampleID','SpeciesID_kraken', 'kraken_percent']
report.write('\t'.join(map(str,header)) + '\n')

#Parse Kraken2 output
lfiles=os.listdir(kreport_path)

# list of report files 
lrepfile=[krep for krep in lfiles if krep.endswith('.report')]

# Opening each report files
for srep in lrepfile:
    with open(kreport_path +'/'+srep,'r') as fhrep: lines=fhrep.readlines()

    samp_report=srep.split(".")
    
    sampID=samp_report[0] # sample ID from sampleID.report
    # fetching through each line
    for l in lines:
        l_parse= l.lstrip().rstrip().split("\t")
        percent = l_parse[0]
        tax_level = l_parse[3]
        tax = l_parse[5].lstrip()
        gen=tax.split(' ')[0] 
        if tax_level == 'S':
            break
    results = [sampID,tax, percent]
    report.write('\t'.join(map(str,results)) + '\n')
    #samp_spec=[(sampID,tax)]
    #speclist.write('\t'.join(map(str,samp_spec)) + '\n')
    samp_gen = [(sampID,gen)]
    samp_genlist.write('\t'.join(map(str,samp_gen)) + '\n') 
report.close()
samp_genlist.close()


