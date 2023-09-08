#!/usr/bin/env python3

# Script to run the prokaryotic genome annotation pipeline pgap (for NCBI pgap annotation tool)
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import os
import itertools
import subprocess
import argparse
import shutil
import collections
parser = argparse.ArgumentParser(prog = 'pgap_run_v3.py [Options]')
#parser.add_argument('--asb_dir', type=str, help = 'provide the assembly files or its path', required=True)
parser.add_argument('--ppgroup_dir', type=str, help =' provide the assembly and metada directory path', required=True)
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
#${PGAP_INPUT_DIR}./pgap.py -n --ignore-all-errors -D singularity -o resul pgap_dirmetdir_test/JBI20000547/input_JBI20000547.yaml
for sample_id in pgdirlist:
    #subprocess.run('./${PGAP_INPUT_DIR}pgap.py -n --no-internet --ignore-all-errors -o ' + sample_id+'_results '+ pg_dir +'/'+ sample_id+'/input_'+ sample_id+'.yaml', shell = True, check = True)
    subprocess.run('${PGAP_INPUT_DIR}./pgap.py -n --ignore-all-errors -D singularity -o ' + sample_id+'_results '+ pg_dir +'/'+ sample_id+'/input_'+ sample_id + '.yaml', shell = True, check = True)
"""
subprocess.run('${PGAP_INPUT_DIR}./pgap.py -n --no-internet --ignore-all-errors --container-path ${PGAP_INPUT_DIR}pgap_2022-04-14.build6021.sif -o ' + sample_id+'_results '+ p\
g_dir +'/'+ sample_id+'/input_'+ sample_id + '.yaml --docker singularity --cpus '+threads, shell = True, check = True)
for asf in aslist:
    #print(asf.split(".")[0])
    source = os.path.join(path, asbl_dir, asf)
    target = os.path.join(path, meta_dir, asf.split(".")[0])
    if asf.split(".")[0] in metdirlist:
        #print(asf, met)
        #source = asbl_dir+'/'+asf
        #source = os.path.join(path, asbl_dir, asf)
        #target = meta_dir+'/'+asf.split(".")[0]+'/
       #target = os.path.join(path, meta_dir, asf.split(".")[0])
        #subprocess.run("copy", asf, met, shell=True, check = True)
        #shutil.copyfile(source, target)
        subprocess.run('cp '+ source + ' ' + target , shell = True, check = True)
        print(os.listdir(target))

      
"""
