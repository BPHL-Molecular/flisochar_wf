#!/usr/bin/env python3

# Script to run pyani Avere Nucleotide Identity (ANI)
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import os
import itertools
import subprocess
import argparse
import shutil
import collections
parser = argparse.ArgumentParser(prog = 'ani_Msrun_v02.py [Options]')
#parser.add_argument('--asb_dir', type=str, help = 'provide the assembly files or its path', required=True)
parser.add_argument('--asbrefseq_dir', type=str, help =' provide the assembly and resef directory path', required=True)
#parser.add_argument('--aniOutDir', type=str, help =' output directory path', required=True)
#parser.add_argument('--threads', default = 8, dest = 'threads', help =' specify the number of threads, (default:%(default)s)')
args = parser.parse_args()

#asbl_dir = args.asb_dir
#threads = str(args.threads)
asr_dir = args.asbrefseq_dir
#out_dir = args.aniOutDir
# os.path.join(path, "User/Desktop", "file.txt")    
#aslist=os.listdir(asbl_dir)
pgdirlist=os.listdir(asr_dir)
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
    #subprocess.run('mkdir ' +sample_id + '_ani', shell = True, check = True)
    subprocess.run('average_nucleotide_identity.py -i ' + asr_dir+'/'+sample_id + ' -o ' + sample_id +'_ani -m ANIm -g --gformat png,pdf', shell = True, check = True)

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
