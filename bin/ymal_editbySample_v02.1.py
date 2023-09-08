#!/usr/bin/env python3

# Script to automate the creation of input and submol yaml files for each sample genome ( for NCBI pgap annotation tool)
# Tassy J-S Bazile, Bioinformatician APHL-Fellow (2021-2023)
import yaml
import itertools
import subprocess
import argparse
parser = argparse.ArgumentParser(prog = 'ymal_editbySample_v01.py [Options]')
#parser.add_argument('--sample_file', type = argparse.FileType('r'), help = 'provide the sample file or its path', required=True)
parser.add_argument('--sample_species', type = argparse.FileType('r'), help =' provide the species name', required=True)
parser.add_argument('--input_yaml', dest="input_yml", help =' provide the input template *.yaml', required=True, default ='input_template.yaml')
parser.add_argument('--input_submol', dest="input_sub", help =' provide the input submol template *.yaml', required=True)
parser.add_argument('--input_ext', dest="input_exten", type=str, help =' provide the extension in fa or fasta', required=True)

args = parser.parse_args()
#sample_f = args.sample_file
#print(sample_f.readlines())
samplesGenus = args.sample_species
input_template = args.input_yml
submol_template = args.input_sub
input_ex=args.input_exten
 
# Getting the sample genomes in a list
#with open(sample_f) as sampl:sampl_list = sampl.readlines()
#sampl_list = sample_f.readlines()
#sampl_list = [line.strip() for line in sampl_list]
#sampl_list =[line.rstrip('\n') for line in sampl_list]
#sampl_list.sort()
#print(sampl_list) checked

#Getting the samples identified genus
#with open("samplesGenus") as identgen:genus_list = identgen.readlines()
sgenus_list = samplesGenus.readlines() # list of tuples
sgenus_list =[line.rstrip('\n') for line in sgenus_list]
for sgt in sgenus_list:
    tp= eval(sgt) # treating list each item as tuple
    #print(tp[1])

    subprocess.run('mkdir -p '+ tp[0], shell = True, check = True)
    # Editing the input yaml file
    with open(input_template) as inp:
        y = yaml.safe_load(inp)
        y['fasta']['location'] = tp[0] + '.' + input_ex
        y['submol']['location'] = 'submol_'+tp[0]+'.yaml'
    
    # Saving the new input as a file
    with open( tp[0]+"/input_" + tp[0] + ".yaml", "w") as samplinp_ymal:
        yaml.dump(y,samplinp_ymal, default_flow_style=False, sort_keys=False)
        #print(yaml.dump(y, default_flow_style=False, sort_keys=False)) checked
   
    # Editing submol yaml file
    with open(submol_template) as sb:
        z = yaml.safe_load(sb)
        z['organism']['genus_species'] = tp[1] 
    # Saving new submol as a file    
    with open(tp[0]+"/submol_"+tp[0]+".yaml", "w") as samplesub_ymal:
        yaml.dump(z,samplesub_ymal, default_flow_style=False, sort_keys=False)

