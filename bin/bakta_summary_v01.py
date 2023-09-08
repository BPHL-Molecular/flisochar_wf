#!/usr/bin/env python3
#bakta_summary_v01.py
# Author" Tassy J. Bazile Bioinformatics Fellow 2021-2023
#sys.path.append('/usr/sammy/') to get append the path of a module
import os
import pandas as pd
import numpy as np
import re
import xlsxwriter
import csv
from Bio import SeqIO
from collections import Counter
import argparse

"""
Business Logic:
This script parses through the output of bakta annotation report for species identication into csv and txt files

Usage:
python bakta_summary_v01.py --kjreppath

"""

#import xlsxwriter
parser = argparse.ArgumentParser(prog = 'bakta_summary_v01.py [Options]')
parser.add_argument('--bkreppath', type=str,help= 'paste path to mash top 10 report file', required=True)

args = parser.parse_args()
bkreport_path = args.bkreppath

def pgap_report(s,path_to_file):
      
    with open(path_to_file, "r") as annot: ancont=annot.readlines()
    rep_list= [ x.replace("\n", "") for x in ancont]
    rep_list = [ x.strip() for x in rep_list] # remove whitespace
    rep_list_sl = rep_list[9:24] # list of items like k:v from raw report
    rep_list_sl = [i for i in rep_list_sl if i] # remove empty items
    dict_rep_ = {i.split(":")[0]:i.split(":")[1] for i in rep_list_sl} # converting list into a dictionary
    df_features = pd.DataFrame(list(dict_rep_.items())) # Making a  Dataframe
    df_featuresT =  df_features.transpose() # to having features as colunms
    header_row = df_featuresT.iloc[0] # converting row into column header
    df_features2 = pd.DataFrame(df_featuresT.values[1:], columns =  header_row)  
    df_features2.insert(0,"SampleID", s, True)  # adding "SampleID" column at the position zero
    # renaming columns
    df_features2.rename(columns={'ncRNA regions': 'ncRNA_regions', 'CRISPR arrays': 'CRISPR_arrays' , 'signal peptides': 'signal_peptides'}, inplace=True)
    return df_features2
    
# Converting report into txt and excel files
def report_csv(df_features2):
    df_features2.to_csv('bakta_annotation_report.txt', sep = ' ', index=True, header=True, quoting = csv.QUOTE_NONE, escapechar = ' ')
    with pd.ExcelWriter('bakta_annotation_report.xlsx', engine = 'xlsxwriter') as writer: # writer is to add more than one sheet in the workbook
        df_features2.to_excel(writer, sheet_name='Annotation_Report', startrow=1)
    return

# Main function to test out the module
def main():
    df_list=[]
    sampList = os.listdir(bkreport_path) # list of dir for each sample
    sampleList = [rep.split("_")[0] for rep in sampList] 
    sampleList.sort()
    for s in sampleList:
        saved_path = bkreport_path + '/' + s +'_bakta'
        rep_file = s +'.txt'
        path_to_file = os.path.join(saved_path, rep_file)
        feature_df = pgap_report(s, path_to_file)
        df_list.append(feature_df)
        feature_dfn = pd.concat(df_list)
        feature_dfn.index = np.arange(1, len(feature_dfn) + 1)
        report_csv(feature_dfn)  
    #print(feature_dfn)
    

if __name__ == "__main__":
    main()

