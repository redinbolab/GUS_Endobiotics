import sys
import os
from os.path import exists, join
import shutil
import pandas as pd
import numpy as np
from scipy import stats
import scipy
import math
from itertools import chain, combinations

############ Setup ############
# - conda/mamba create -n fitFinder python numpy pandas scipy
# - conda/mamba create -n fitFinder python numpy pandas scipy "blas=*=*accelerate*" # apple silicon

############ Usage ############
# - conda/mamba activate fitFinder
# - python omics_ttest_v2.py {mgx/mpx} {out_filename.xlsx}
# - python omics_ttest_v2.py 

############ Reformat Abundance Data for Regression Analysis + Function Definitions ############

# Definitions
# Calculate t test between 2 sample groups for a given GUS_group
def ttest(intensity_df, sample_grouping_df, sample_grouping_col, intensity_col, headerSet):
    intensity_df[sample_grouping_col] = intensity_df[sample_grouping_col].astype('string')
    # map vars
    group_0 = intensity_df[intensity_df[sample_grouping_col]=='0'][intensity_col]
    group_1 = intensity_df[intensity_df[sample_grouping_col]=='1'][intensity_col]
    ttest_result = scipy.stats.ttest_ind(group_0,group_1,equal_var=False) # perform t test
    t_test_tmp_summary_df = pd.DataFrame({'feature_group':[headerSet],'sample_grouping_col':[sample_grouping_col],'ttest_pvalue':[ttest_result.pvalue]}) # gather results
    return(t_test_tmp_summary_df)

def mod_log2(a):
    a = float(a)
    if a == 0:
        return(0)
    else:
        return(math.log2(a))

def pivot_df(in_df,index_col,in_col,value_col):

    def get_colname(x, in_col):
        if x == in_col:
            return(x)
        outname = '_'.join([x,in_col])
        return(outname)
    
    in_df = in_df[[index_col,in_col,value_col]].drop_duplicates()
    pivot_df = in_df.pivot(index=index_col,columns=[in_col],values=value_col).reset_index()
    new_cols = [get_colname(x, index_col) for x in pivot_df.columns]
    pivot_df.columns = new_cols
    cols_to_drop = [in_col,value_col]
    in_df = in_df.drop(columns=cols_to_drop)
    merged_df = in_df.merge(pivot_df,how='inner').drop_duplicates()

    return(merged_df)

# Data modification
import re
def remove_trailing_number_after_underscore(s):
    s = str(s)
    # Pattern to match a string ending with an underscore followed by a number
    pattern = r'_(\d+)$'
    # Search for the pattern in the string
    match = re.search(pattern, s)
    # If the pattern is found, remove the trailing underscore and the number
    if match:
        # The start position of the matched pattern
        start_pos = match.start()
        # Return the string up to the start position of the matched pattern
        return s[:start_pos]
    # If no match is found, return the original string
    else:
        return s

############ Perform T-tests ############
if __name__ == '__main__':
    
    omic_mode = sys.argv[1]
    ofnm = sys.argv[2]

    if omic_mode not in ['mgx','mpx']:
        print('please list mgx or mpx for the omic mode')
        sys.exit()
    if '.xlsx' not in ofnm:
        print('please list a .xlsx as the outfile name')
        sys.exit()

    ### import data
    # rate classifications
    sample_grouping_in = 'rates_z_score_binned.csv'
    sample_grouping_df = pd.read_csv(sample_grouping_in)
    sample_grouping_df['sample_ID'] = sample_grouping_df['sample_ID'].map(lambda a: str(a).replace('D','')) 

    log2 = False # set initial log2 status

    # set input colnames
    sample_ID_col = 'sample_ID'
    intensity_col = 'value'
    struct_class_col = 'Struct_Class'

    groupColList = ['Phyla','Class','Order','Family','Genus','Species']

    if omic_mode == 'mpx':
        groupColList = ['merged_'+k for k in groupColList]
        header_col = 'Uniprot_ID'
        intensity_col = 'intensity'
        intensity_csv = 'metaproteomic_abundance.csv'

        intensity_df = pd.read_csv(intensity_csv)
        intensity_df['merged_Species'] = intensity_df['merged_Species'].map(lambda a: remove_trailing_number_after_underscore(a))
        intensity_df['merged_Class'] = intensity_df['merged_Class'].map(lambda a: a.replace('Clostridia_A','Clostridia'))
        log2 = True

    if omic_mode == 'mgx':
        intensity_col = 'value'
        header_col = 'variable'
        intensity_csv = 'metagenomic_abundance.csv' # mgx
        intensity_df = pd.read_csv(intensity_csv)

    # initialize vals, lists
    ttest_out_df_list = []
    min_samples = 7 # set min samples needed for a fit
    headerList = list(pd.unique(intensity_df[header_col]))

    # trim df to only relevant columns
    cols = [sample_ID_col,intensity_col,header_col]+groupColList
    intensity_df = intensity_df[cols]
    
    taxLevel_DF_list = []

    for groupCol in groupColList:
        print(groupCol)

    for groupCol in groupColList:
        outFileCt = 0

        ttest_summary_list = []

        cols = [sample_ID_col,intensity_col,header_col]+[groupCol]

        intensity_df_groupCol = intensity_df[cols]

        groupCol_level_list = list(pd.unique(intensity_df_groupCol[groupCol])) # get unique entries in current grouping column
        # print(groupCol_level_list)

        groupsChecked = []

        for groupCol_level in groupCol_level_list:

            # # shape df by matching to current seq's group class
            group_working_df = intensity_df_groupCol[intensity_df_groupCol[groupCol] == groupCol_level]
            # print('inGroup',group_working_df) 

            # sum abundances by matching sample_ID to prepare for regression
            if groupCol != 'Species':
                group_working_df = group_working_df.groupby(['sample_ID',groupCol],as_index=False)[intensity_col].sum()
            else:
                pass

            # log transform if vals have not been normalized
            if log2 == True:
                group_working_df[intensity_col] = group_working_df[intensity_col].map(lambda a: math.log2(a))

            # filter 0s
            group_working_df = group_working_df[group_working_df[intensity_col]>0]

            # if filtered group is in >= min_samples, run regression
            new_n = len(group_working_df[sample_ID_col])
            new_nStr = str(new_n)

            if (new_n >= min_samples):

                # mutate data
                group_working_df[sample_ID_col] = group_working_df[sample_ID_col].map(lambda a: str(a))
                sample_grouping_df[sample_ID_col] = sample_grouping_df[sample_ID_col].map(lambda a: str(a))
                # print(group_working_df,'\n',sample_grouping_df)
                
                group_working_df = pd.merge(group_working_df, sample_grouping_df, on = "sample_ID", how='right').fillna(0) # merge, fill na

                # mutate col to match sample_ID as needed
                group_working_df[sample_ID_col] = group_working_df[sample_ID_col].map(lambda a: str(a))

                # define sample grouping colnames
                sample_grouping_df_cols = [k for k in sample_grouping_df.columns if k!='sample_ID']

                # run t test on each prot grouping
                for sample_grouping_col in sample_grouping_df_cols:
                    to_summary_df = ttest(group_working_df, sample_grouping_df, sample_grouping_col, intensity_col, groupCol_level)
                    # print(to_summary_df)
                    ttest_summary_list.append(to_summary_df)

                # concatenate ttest results
                if len(ttest_summary_list) > 0:
                    concat_out_df = pd.concat(ttest_summary_list).reset_index(drop=True).sort_values('ttest_pvalue',ascending=True)
                    # print(concat_out_df)
                    ttest_out_df_list.append([concat_out_df,groupCol])
                    
    # export results
    ttest_excel_out_name = ofnm
    with pd.ExcelWriter(ttest_excel_out_name) as writer:
        sheetOrder = groupColList
        for name in sheetOrder:
            for x in range(len(ttest_out_df_list)):
                DF = ttest_out_df_list[x][0]
                taxLevel = ttest_out_df_list[x][1]
                if taxLevel == name:
                    DF.to_excel(writer, sheet_name=taxLevel,index=False)