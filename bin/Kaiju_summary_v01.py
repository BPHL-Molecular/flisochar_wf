#!/usr/bin/env python
# Kaiju_summary_v01.py

#sys.path.append('/usr/sammy/') to get appen the path of a module
import os
import pandas as pd
import numpy as np
import re
import argparse

"""
Business Logic:
This script parses through the output of Kaiju report for species identication into csv and txt files

Usage:
python Kaiju_summary_v01.py --kjreppath

"""

#import xlsxwriter
parser = argparse.ArgumentParser(prog = 'Kaiju_summary_v01.py [Options]')
parser.add_argument('--kjreppath', type=str,help= 'paste path to mash top 10 report file', required=True)

args = parser.parse_args()
kjreport_path = args.kjreppath

spec_con=[] #  to concatenate df into list 
CFI_list=[]
rel_abundList=[]
species_list=[]
# Parsing Kraken2 output
def kaiju_report(path_to_file,s):
    #Opening Kaijou report
    with open(path_to_file, 'r') as kaiReport: lines = kaiReport.readlines() # open each related file in midas output
    firstline = lines[1] # skip the header or lines[0]
    lin_parse = firstline.strip().split('\t')# split tab ->list of columns produced
    rel_abund = lin_parse[1]
    #species_part= lin_parse[0].split('_')[0:2] #['Bacillus', 'niacini', '62369'] then  1st two ['Bacillus', 'niacini'] ie genus and species
    species= lin_parse[4] # genus and species
    # Assigning items to list
    rel_abundList.append(rel_abund) 
    species_list.append(species)
    CFI_list.append(s)
    return  CFI_list, species_list, rel_abundList
    
#Data frame for the report using the 4 lists(cfi, prcent, txlevel, scname) representing columns
def data_Frame(CFI_list, species_list, rel_abundList):
    #Merging the for three lists (column variables)
    list_tuples = list(zip(CFI_list, species_list, rel_abundList))
    #pandasDataFrame the zipped lists for converting them
    df=pd.DataFrame(list_tuples, columns=["Sample_ID", "SpeciesID_Kaiju", "Kaiju_Percent"])
    return df # len(samples_list) dataframes,in this case 4

def df_to_csv(df):
    df.to_csv('kaiju_MSreport_tab.csv')
    df.to_csv('kaiju_MSreport_tab.txt', index=None, sep='\t')
    return

# Report file JBI20000547_kaiju_summary_esp.tsv
# Creating a main function, testing the module
def main():
    saved_path = kjreport_path
    rep_list = os.listdir(kjreport_path)
    sampList = [repka.split('_')[0] for repka in rep_list]
    sampList.sort()  
    cfi = []
    sp_co =[]
    for s in sampList:
        rep_file = s +'_kaiju_summary_esp.tsv'
        path_to_file = os.path.join(saved_path, rep_file) 
        CFI_list, species_list, rel_abundList= kaiju_report(path_to_file,s) # function 1 opnes report
        dfr = data_Frame(CFI_list, species_list, rel_abundList) #Function 2 creates df 
        dfr.index = np.arange(1, len(dfr.index)+1)
    df_to_csv(dfr) # Function 3 creates text files
    print(dfr)
# Managing how the program will run (whther a module or standalone)
if __name__ == "__main__":
    main()



