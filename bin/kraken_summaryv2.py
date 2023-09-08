#!/usr/bin/env python3

# Script to parse kraken output
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
#import yaml
import itertools
import subprocess
import argparse
import os
import glob
import subprocess # copy file
import shutil
parser = argparse.ArgumentParser(prog = 'kraken_summary.py [Options]')
#parser.add_argument('--kreport_file', type = argparse.FileType('r'), help = 'provide the sample file or its path', required=True)
parser.add_argument('--reppath', type=str,help= 'paste path to sample.report file')
parser.add_argument('--out', type=str,help= 'paste file to write report file')
args = parser.parse_args()
kreport_path = args.reppath

# adding report file to a directory
subprocess.run('mkdir -p krakenTempDir', shell = True, check = True)
out_rep=args.out
subprocess.run('cp '+ kreport_path+ ' krakenTempDir', shell = True, check = True)
#Report file

report = open(out_rep +'.txt', 'w') # user chose a name for the file
header = ['SampleID','SpeciesID_kraken', 'kraken_percent']
report.write('\t'.join(map(str,header)) + '\n')
#Parse Kraken2 output
lfiles=os.listdir('krakenTempDir')

# list of report files 
lrepfile=[krep for krep in lfiles if krep.endswith('.report')]

# Opening each report files
for srep in lrepfile:
    with open('krakenTempDir/'+srep,'r') as fhrep: lines=fhrep.readlines()

    samp_report=srep.split(".")
    sampID=samp_report[0] # sample ID from sampleID.report
    # fetching through each line
    for l in lines:
        l_parse= l.lstrip().rstrip().split("\t")
        percent = l_parse[0]
        tax_level = l_parse[3]
        tax = l_parse[5].lstrip()
        if tax_level == 'S':
            break
    results = [sampID,tax, percent]
    report.write('\t'.join(map(str,results)) + '\n')


