#!/usr/bin/env python
# Parsing output for mash on multiple samples
import os
import re
import pandas as pd
import numpy as np
#import xlsxwriter
import itertools
from functools import reduce
import argparse
"""
This module parses Mash results to an excel sheet, text file and  generates giving genus, species, distance and accession number
"""
parser = argparse.ArgumentParser(prog = 'parsing_outMashv2MS3tryA3.py [Options]')
parser.add_argument('--mreppath', type=str,help= 'paste path to the directory of mash top 10 report file', required=True)
#parser.add_argument('--samppath', type=str,help= 'paste path to sample file', required=True) 
args = parser.parse_args()
mreport_path = args.mreppath
#sample_path = args.samppath

def mash_report(path_to_file,s):
    with open(path_to_file, 'r') as mash: top_hit = mash.readline()
    # distance           
    distance = top_hit.split()[2]
    accession = top_hit.split('-')[5]
    gen_space = top_hit.split()[0] # split on space
    gen_str = gen_space.split('-')[7] # split on dashes
    gen_str01 =gen_str.split('_') # split on underscore
    gen_str02 =[it for it in gen_str01 if it.strip()]# remove space in the list
    genus = gen_str02[0]
    species_ = gen_str02[1]
    species = species_.split('.')[0]# split on period
  
    #Representing report in a dataframe wherein each sample is a row(record) 
    #First, assign the results to a list for each sample
    mash_result=[]
    mash_result.append(s)
    mash_result.append(genus)
    mash_result.append(species)
    mash_result.append(distance)
    mash_result.append(accession)
    mash_result[1:3] = [' '.join(mash_result[1:3])] # rapproaching genus and species
    # List of all lists (these are the rows for the dataframe)
    return mash_result

# Appending result to list for datafradef result_list(mash_result):
data_samples =[]
def result_list(mash_result):  
    data_samples.append(mash_result)
    return data_samples # list has duplicates
# Unique list
new_sampdata =[]
def uniq_sample(data_samples):
    for x in data_samples:
        if x not in new_sampdata:
            new_sampdata.append(x)
    return new_sampdata
# Data frame
def result_to_df(new_sampdata, saved_path):
    #df=pd.DataFrame(data=new_sampdata, columns=['Sample_ID', 'Genus', 'Species', 'Distance', 'Accession_Number'])
    df=pd.DataFrame(data=new_sampdata, columns=['Sample_ID', 'SpeciesID_Mash', 'Mash_Distance', 'Accession_Number'])
    df.index = np.arange(1, len(df)+1)# having index start from 1 to number of row in df
    return df
# Converting data frame into csv amd text files
def df_to_csv(df):
    df.to_csv('mash_MSreport_tab.csv')
    df.to_csv('mash_MSreport_tab.txt', index=None, sep='\t')
    return
       
# Testing the module
def main():
    saved_path = mreport_path # path to mash output
    report_list =os.listdir(saved_path) # list of file in dir   
    report_list.sort()
    lsamp = [ mrep.split('_')[0] for mrep in report_list] # list of samples
    for s in lsamp:
            rep_file = s +  '_distances_top10.tab'
            saved_path = mreport_path 
            path_to_file = os.path.join(saved_path, rep_file)
            pre_result= mash_report(path_to_file,s) # Function 1 opens mash report, produces len(sample_list) lists
            list_of_results= result_list(pre_result) # function 2 making a list of  all lists
    # replaces lsit_list.appenf()
    datafr1 = result_to_df(list_of_results, saved_path) # creates a df 
    df_to_csv(datafr1)
    #print(datafr)
    print(datafr1)
    
if __name__ == "__main__":
    main()

