#!/usr/bin/env python3

# Script to run the prokaryotic genome annotation pipeline pgap (for NCBI pgap annotation tool)
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import os
import itertools
import subprocess
import argparse
import shutil
import collections

"""
Business Logic:
This script runs pgap to confirm the taxomic classification to the samples using ANI

To run the script, use the following command:
python pgap_tax_v1.py --ppgroup_dir 
"""


parser = argparse.ArgumentParser(prog = 'pgap_run_v2.py [Options]')
#parser.add_argument('--asb_dir', type=str, help = 'provide the assembly files or its path', required=True)
parser.add_argument('--ppgroup_dir', type=str, help =' provide the directory path that contains assembly and metada', required=True)
#parser.add_argument('--threads', default = 8, dest = 'threads', help =' specify the number of threads, (default:%(default)s)')
args = parser.parse_args()

#asbl_dir = args.asb_dir
#threads = str(args.threads)
pg_dir = args.ppgroup_dir
# os.path.join(path, "User/Desktop", "file.txt")    
#aslist=os.listdir(asbl_dir)
pgdirlist=os.listdir(pg_dir)
pgdirlist.sort()
#suf =".fa"
#path = os.getcwd()
#metdirlist.sort()
#aslist.sort()
#print(aslist, len(aslist))
#print(metdirlist, len(metdirlist))
#ll =[l.split(".")[0] for l in aslist]
#print(pgdirlist)
for sample_id in pgdirlist:
    subprocess.run('${PGAP_INPUT_DIR}./pgap.py -n --no-internet --ignore-all-errors --container-path ${PGAP_INPUT_DIR}pgap_2022-04-14.build6021.sif --taxcheck-only -o ' + sample_id+'_taxCheck '+ pg_dir +'/'+ sample_id+'/input_'+ sample_id + '.yaml --docker singularity', shell = True, check = True)


