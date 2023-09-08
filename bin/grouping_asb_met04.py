#!/usr/bin/env python3

# Script to automate the copy of assembly files for each sample genome (for NCBI pgap annotation tool)
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import os
import itertools
import subprocess
import argparse
import shutil
import collections
parser = argparse.ArgumentParser(prog = 'grouping_asb_met01.py [Options]')
parser.add_argument('--asb_dir', type=str, help = 'provide the assembly files or its path', required=True)
parser.add_argument('--met_dir', type=str, help =' provide the met dir path', required=True)

args = parser.parse_args()
asbl_dir = args.asb_dir
meta_dir = args.met_dir
# os.path.join(path, "User/Desktop", "file.txt")    
aslist=os.listdir(asbl_dir)
metdirlist=os.listdir(meta_dir)
suf =".fa"
path = os.getcwd()
metdirlist.sort()
aslist.sort()
#print(aslist, len(aslist))
#print(metdirlist, len(metdirlist))
ll =[l.split(".")[0] for l in aslist]
#print(ll)
for (ab, met) in zip(aslist, metdirlist):
        #srs = ab
        #tgt=  met
        if ab.split(".")[0] ==met: 
            #print(srs)
            subprocess.run('cp ' + asbl_dir +'/' +ab + ' ' + meta_dir+'/' + met +'/', shell = True, check = True)
        else:
            subprocess.run('mv ' + meta_dir+'/' + met +'/'+ab,  shell = True, check = True)            
"""
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
