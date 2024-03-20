import sys
import os
from os.path import exists, join
import shutil
import pandas as pd
import numpy as np
from scipy import stats
import math
from itertools import chain, combinations
import time

############ Setup ############
# - conda/mamba create -n fitFinder python numpy pandas scipy
# - conda/mamba create -n fitFinder python numpy pandas scipy "blas=*=*accelerate*" # apple silicon

############ Usage ############
# - conda/mamba activate fitFinder
# - python fitFinder_metaproteomics_v4.2.py {RMSD_threshold, float or NA}
# - python fitFinder_metaproteomics_v4.2.py 2
# - python fitFinder_metaproteomics_v4.2.py NA

############ Reformat Abundance Data for Regression Analysis + Function Definitions ############

##### define fxns
def mod_log2(a):
    a = float(a)
    if a == 0:
        return(0)
    else:
        return(math.log2(a))

# initialize powerset function
def powerset(iterable, ):
    # "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)" # explanation of function's general behavior
    s = list(iterable)
    # return chain.from_iterable(combinations(s, r) for r in range(len(s)+1)) # default behavior - consider groups of max size equal to the size of the whole set
    return chain.from_iterable(combinations(s, r) for r in range(setSize+1)) # consider groups of max size equal to a manual limit (N unique proteins)

# get number of sim cols above thresh for each header
def rmsd_cols_above_thresh(header,df,header_col,simThresh,sim_outCol_list):
    header_df = df[df[header_col]==header]
    aboveThresh = 0
    for col in header_df:
        if col in sim_outCol_list:
            col_sim = header_df[col].iloc[0]
            if col_sim >= simThresh:
                aboveThresh+=1
    return(aboveThresh)

def filter_df_by_sim(tmp_df, headerStr, sim_outCol_list):
    tmp_header_list = list(pd.unique(tmp_df[header_col]))

    rmsd_cols_to_keep = [x for x in tmp_df.columns if 'RMSD_' in x]
    rmsd_cols_to_keep = [x for x in rmsd_cols_to_keep if any(header in x for header in tmp_header_list)]
    all_cols_to_keep = cols+rmsd_cols_to_keep
    tmp_df = tmp_df[all_cols_to_keep]          
                    
    for col in tmp_df.columns:
        if col in sim_outCol_list:  
            tmp_df = tmp_df.copy() # avoid warning msgs              
            tmp_df['above_thresh'] = tmp_df[header_col].map(lambda a: rmsd_cols_above_thresh(a,tmp_df,header_col,simThresh,sim_outCol_list))
            max_above_thresh = tmp_df['above_thresh'].max()

            if max_above_thresh > 0:
                tmp_df = tmp_df[tmp_df['above_thresh']<max_above_thresh]                                                   
                if headerStr not in list(pd.unique(tmp_df[header_col])):
                    return()

            # get RMSD cols, subset df
            tmp_header_list = list(pd.unique(tmp_df[header_col]))
            rmsd_cols_to_keep = [x for x in tmp_df.columns if 'RMSD_' in x]
            rmsd_cols_to_keep = [x for x in rmsd_cols_to_keep if any(header in x for header in tmp_header_list)]
            all_cols_to_keep = cols+rmsd_cols_to_keep
            tmp_df = tmp_df[all_cols_to_keep]

    return(tmp_df)

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

############ Perform Regressions ############
if __name__ == '__main__':
    # start timer
    start = time.time()

    # set similarity thresh
    simThresh = sys.argv[1]
    
    # Import concatenated rates
    rate_csv = 'GUS_rates.csv'
    rate_df = pd.read_csv(rate_csv)
    
    
    intensity_csv = 'GUS_metaproteomics.csv' # import intensity sheet - metaproteomics
    intensity_df = pd.read_csv(intensity_csv)

    out_suffix = ''

    # set input colnames
    sample_ID_col = 'sample_ID'
    intensity_col = 'intensity'
    struct_class_col = 'Struct_Class'
    header_col = 'Uniprot_ID'

    ##### initialize vals, lists
    outFileCt = 0
    headerIter = 0
    
    rateList = []
    uniqueSets_goodN = []
    uniqueSets_full = []
    summary_df_list = []

    setSize = 20 # set max number of unique sequences to check in each grouping

    # stats thresholds - change as desired
    rsq_thresh = .7
    p_val_thresh = .05
    slope_thresh = 0
    min_samples = 5 # set min samples needed for a fit
    #####

    for col in rate_df.columns:
        if col != 'sample_ID':
            strings_to_check = ['activity','Activity']
            if any(s in col for s in strings_to_check):
                rate_df[col+'_log2'] = rate_df[col].map(lambda a: mod_log2(a))
                rate_df.pop(col)

    intensity_df['Uniprot_ID'] = intensity_df['Uniprot_ID'].map(lambda a: a.replace('S-','S'))

    # filter intensity_df by only those samples which are in rate_df['sample_ID']
    rates_sampleList = rate_df['sample_ID'].to_list()
    intensity_df = intensity_df[intensity_df['sample_ID'].isin(rates_sampleList)]

    # trim df to only relevant columns, make sure it's in expected order
    cols = [sample_ID_col,intensity_col,header_col,struct_class_col]
    intensity_df = intensity_df[cols]

    dropList = ['sample_ID']
    rate_df_trunc = rate_df.drop(columns=dropList)

    if simThresh != 'NA':
        try:
            simThresh = float(simThresh) # get similarity threshold
        except ValueError as e:
            # Handle ValueError (e.g., input cannot be converted to float)
            print("Error: Invalid input. Please enter a valid float RMSD similarity threshold between 0 and inf.")
            print("Error below: \n",e)
            sys.exit(1)  # Exit with a non-zero status code to indicate error
        except Exception as e:
            # Handle any other exceptions not caught above
            print("An error occurred:", str(e))
            sys.exit(1)  # Exit with a non-zero status code to indicate error
            
        rmsd_csv = 'MPX_GUS_RMSD.csv'
        rmsd_df = pd.read_csv(rmsd_csv)

        rmsd_df = rmsd_df.set_index('Header')

        newHeaderList = [str(x) for x in intensity_df[header_col].to_list()]

        rmsd_df = rmsd_df.drop(columns=[col for col in rmsd_df if col not in newHeaderList])

        # get list of headers from rmsd_df
        headerList = rmsd_df.index.to_list()
        headerStr = headerList[headerIter]

        dropRowIndexList = [row for row in rmsd_df.index if row not in newHeaderList]
        rmsd_df = rmsd_df.drop(index=dropRowIndexList)

        path = '_'.join(['./'+intensity_csv.replace('.csv',''),rate_csv.replace('.csv',''),'RMSD',str(simThresh),intensity_col,str(rsq_thresh)])+out_suffix
        outfnm = '_'.join([intensity_csv.replace('.csv',''),rate_csv.replace('.csv',''),'RMSD',str(simThresh),intensity_col,str(rsq_thresh)])+out_suffix+'_summary.csv'
        summary_outfnm = join(path,outfnm)
        outFilePath = join('./',path)

        # delete output path if it exists, re-create
        path_exists = exists(path)
        if path_exists == True:
            shutil.rmtree(path)
        os.mkdir(path)

        for headerIter in range(len(headerList)):
            headerStr = headerList[headerIter]
            dfHeaderList = list(intensity_df[header_col])
            header_index = dfHeaderList.index(headerStr)

            # retrieve loop class of current header
            headerLoop = intensity_df.iloc[header_index, 3]
            # shape df by matching to current seq's class
            class_working_df = intensity_df[intensity_df[struct_class_col] == headerLoop]

            sim_outCol_list = []

            # get similarity matrix
            for x in class_working_df[header_col]:
                x_sim_list = []
                outCol = 'RMSD_'+x
                x_index = headerList.index(x)
                for y in class_working_df[header_col]:
                    y_index = headerList.index(y)
                    seqSim = rmsd_df.iloc[x_index, y_index]
                    x_sim_list.append(seqSim)
                class_working_df = class_working_df.copy()
                class_working_df[outCol] = x_sim_list
                sim_outCol_list.append(outCol)

            # get initial n, pre-filtering
            n_int = len(pd.unique(class_working_df['sample_ID']))
            n = str(n_int)

            # filter to only seqs which meet thresh to current header
            header_sim_col = 'RMSD_'+headerStr
            class_working_df = class_working_df[class_working_df[header_sim_col]<=simThresh]

            if n_int >=4:
                uniqueHeaderList = class_working_df[header_col].unique()

                # for each set in powerset, run a filtering process
                for group in powerset(uniqueHeaderList):
                    groupList = list(group) # convert iterator to list
                    groupSet = set(group) # convert tuple to set
                    break_var = False
                    # check break_var and check that group has not been checked
                    if (break_var == False) and (groupSet not in uniqueSets_full):
                        uniqueSets_full.append(groupSet)
                        # set limit of subset size to consider
                        if len(group)<setSize+1:
                            currentList = group
                            # eliminate empty set from consideration, check for header
                            if groupSet and (headerStr in groupSet):
                                # shape df by current iter list in power set, removing all others
                                new_class_working_df = class_working_df[class_working_df[header_col].apply(lambda x: any([k in x for k in currentList]))].reset_index(drop=True)
                                test = new_class_working_df.copy()
                                filtered_df = filter_df_by_sim(test, headerStr, sim_outCol_list)
                                
                                if len(filtered_df) > 0:
                                    uniqueEntries = filtered_df[header_col].unique() # get number of unique entries in group
                                    uniqueEntries_set = set(uniqueEntries)
                                    uniqueEntriesStr = str(len(uniqueEntries))
                                    
                                    class_list = list(filtered_df[struct_class_col].unique()) # get list of classes in grouping (if not grouping by class)

                                    output_df = filtered_df.copy()
                                    output_df = output_df.reset_index(drop=True).sort_values(['sample_ID'],ascending=True)

                                    # sum abundances by matching sample_ID to prepare for regression
                                    filtered_df = filtered_df.groupby([sample_ID_col], as_index=False)[intensity_col].sum()
                                    new_n = len(filtered_df[sample_ID_col])
                                    new_nStr = str(new_n)
                                    fitFound = False

                                    # if filtered group in powerset has min_samples or more unique sample points, run regression
                                    if (new_n >= min_samples) and (uniqueEntries_set not in uniqueSets_goodN):
                                        uniqueSets_goodN.append(uniqueEntries_set)
                                        duplicate_entries = [k for k in uniqueSets_goodN if len([x for x in uniqueSets_goodN if x == k])>=2]
                                        regression_df = filtered_df
                                        regression_df[intensity_col] = regression_df[intensity_col].map(lambda a: math.log2(a)) # log transform if desired (mpx data only, not mgx)
                                        new_n = len(regression_df[sample_ID_col])
                                        new_nStr = str(new_n)

                                        # if filtered group in powerset has at least 1 element present in >= min_samples, run regression
                                        if (new_n >= min_samples):
                                            merged_df = pd.merge(regression_df, rate_df, on = "sample_ID", how='inner') # merge regression and rate df
                                            substrateList = list(rate_df_trunc.columns) # get substrate list
                                            # initialize vals and lists
                                            fitsFound = 0
                                            rsq_list = []
                                            pvalue_list = []
                                            slope_list = []
                                            spearmanr_list = []
                                            class_list = []
                                            hit_cols = []
                                            # calculate regression for each column in rate df
                                            for column in rate_df_trunc:
                                                # select cols for regression
                                                y = merged_df[column]
                                                X = merged_df[intensity_col]
                                                # get results of linear regression
                                                slope, intercept, r_value, p_value, std_err = stats.linregress(X,y)
                                                spearmanr = stats.spearmanr(X,y).statistic
                                                rsq_value = r_value**2
                                                # export values to lists
                                                rsq_list.append(rsq_value)
                                                pvalue_list.append(p_value)
                                                slope_list.append(slope)
                                                spearmanr_list.append(spearmanr)
                                                class_list.append(headerLoop)
                                                # identify hit columns and update counts
                                                if (p_value <= p_val_thresh) and (rsq_value >= rsq_thresh) and (slope>slope_thresh):
                                                    fitFound = True
                                                    fitsFound += 1
                                                    hit_cols.append(column)

                                        if fitFound == True:
                                            # update outfile count
                                            outFileCt += 1
                                            outFileCt_str = str(outFileCt)
                                            # export details of current grouping
                                            outList = [uniqueEntries, headerStr, headerLoop]
                                            outList_df = pd.DataFrame({ 'Set info': outList})
                                            # make df from this grouping's stats
                                            statsDict = {'substrate':substrateList, 'rsq_val': rsq_list, 'p_val': pvalue_list, 'slope': slope_list,'spearman_val': spearmanr_list}
                                            stats_df =  pd.DataFrame(statsDict)
                                            # filter stats_df to only positive results with a minimum rsq thresh (could change to spearman)
                                            paredStatsDF = stats_df[(stats_df.slope >= 0) & (stats_df.rsq_val >= .7)].reset_index(drop=True)
                                            # id best fit by p value
                                            bestFitSubstrate = str(paredStatsDF.loc[paredStatsDF['p_val'].idxmin(),'substrate'])

                                            for hit in hit_cols:
                                                # get hit stats
                                                hit_stats = paredStatsDF[paredStatsDF['substrate']==hit].reset_index(drop=True)
                                                hit_rsq_value = hit_stats['rsq_val'][0]
                                                hit_p_value = hit_stats['p_val'][0]
                                                hit_slope = hit_stats['slope'][0]       
                                                hit_spearmanr = hit_stats['spearman_val'][0]
                                                hit_spearmanr_str = str(hit_spearmanr)     
                                                hit_p_valueStr = str(round(hit_p_value,4))
                                                hit_rsq_valueStr = str(round(hit_rsq_value,4))
                                                # pass results to summary df output
                                                to_summary_df = pd.DataFrame({'result_num':[outFileCt_str],'unique_sample_count':[new_nStr],'unique_seq_count':[len(uniqueEntries)],'unique_sequences':[';'.join(uniqueEntries)],'Struct_Class':[headerLoop],'substrate':[hit], 'spearman_val':[hit_spearmanr], 'rsq_val': [hit_rsq_valueStr], 'p_val': [hit_p_valueStr], 'slope': [hit_slope]})
                                                summary_df_list.append(to_summary_df)
                                            
                                            # get stats for best fit by p val
                                            good_rsq_value = paredStatsDF.loc[paredStatsDF['p_val'].idxmin(),'rsq_val']
                                            good_p_value = paredStatsDF.loc[paredStatsDF['p_val'].idxmin(),'p_val']
                                            good_slope = paredStatsDF.loc[paredStatsDF['p_val'].idxmin(),'slope']        
                                            good_spearmanr = paredStatsDF.loc[paredStatsDF['p_val'].idxmin(),'spearman_val']
                                            good_spearmanr_str = str(good_spearmanr)                    
                                            bestFitSubstrate = str(paredStatsDF.loc[paredStatsDF['p_val'].idxmin(),'substrate'])
                                            good_p_valueStr = str(round(good_p_value,4))
                                            good_rsq_valueStr = str(round(good_rsq_value,4))
                                            # cat, export df for given run
                                            concatList = [merged_df, outList_df, stats_df, output_df]
                                            exportDF = pd.concat(concatList, axis=1)
                                            run_outfnm = '_'.join([outFileCt_str,uniqueEntriesStr,new_nStr,str(fitsFound),bestFitSubstrate,headerLoop,good_spearmanr_str,good_p_valueStr])
                                            exportDF.to_csv(join(path,run_outfnm+".csv"), sep=',', index=False)
        # export summary
        if len(summary_df_list) > 0:
            exportDF = pd.concat(summary_df_list, axis=0).reset_index(drop=True).sort_values('p_val',ascending=True)
            # pivot wider
            export_df_pivot = pivot_df(exportDF,'unique_sequences','substrate','p_val')
            export_df_pivot.to_csv(summary_outfnm,index=False)

        else:
            pass

    ##### by struct class, no similarity
    elif simThresh == 'NA':
        headerList = intensity_df[header_col].unique()

        path = '_'.join(['./'+intensity_csv.replace('.csv',''),rate_csv.replace('.csv',''),'all-class',intensity_col])+out_suffix
        outfnm = '_'.join([intensity_csv.replace('.csv',''),rate_csv.replace('.csv',''),'all-class',intensity_col])+out_suffix+'_byClass_summary.csv'
        summary_outfnm = join(path,outfnm)
        outFilePath = join('./',path)

        # delete output path if it exists, re-create
        path_exists = exists(path)
        if path_exists == True:
            shutil.rmtree(path)
        os.mkdir(path)
        
        loopsChecked = []
        for headerIter in range(len(headerList)):
            uniqueSets_full = []
            headerStr = headerList[headerIter]
            dfHeaderList = list(intensity_df[header_col])
            header_index = dfHeaderList.index(headerStr)
            # get class of current header
            headerLoop = intensity_df.iloc[header_index, 3]
            if headerLoop not in loopsChecked:
                loopsChecked.append(headerLoop)

                class_working_df = intensity_df

                # shape df by matching to current seq's loop class
                class_working_df = intensity_df[intensity_df[struct_class_col] == headerLoop]
                group = list(pd.unique(class_working_df[header_col]))

                # sum abundances by matching sample_ID to prepare for regression
                class_working_df = class_working_df.groupby(['sample_ID'], as_index=False)[intensity_col].sum()

                regression_df = class_working_df

                # log transform if vals have not been normalized
                regression_df[intensity_col] = regression_df[intensity_col].map(lambda a: math.log2(a))

                merged_df = pd.merge(regression_df, rate_df, on = "sample_ID", how='inner')
                fitFound = False
                fitsFound = 0
                rsq_list = []
                pvalue_list = []
                slope_list = []
                spearmanr_list = []
                class_list = []

                substrateList = list(rate_df_trunc.columns)

                hit_cols = []

                for response_col in rate_df_trunc:
                    y = merged_df[response_col]
                    X = merged_df[intensity_col]
                    slope, intercept, r_value, p_value, std_err = stats.linregress(X,y)
                    spearmanr = stats.spearmanr(X,y).statistic
                    rsq_value = r_value**2
                    rsq_list.append(rsq_value)
                    pvalue_list.append(p_value)
                    slope_list.append(slope)
                    spearmanr_list.append(spearmanr)
                    class_list.append(headerLoop)

                    p_valueStr = str(round(p_value,4))
                    rsq_valueStr = str(round(rsq_value,4))
                    if (p_value <= p_val_thresh) and (rsq_value >= rsq_thresh) and (slope>slope_thresh):
                        fitFound = True
                        fitsFound += 1
                        hit_cols.append(response_col)
            
                outList = [headerLoop]
                outList_df = pd.DataFrame({ 'Set info': outList})
                statsDict = {'substrate':substrateList, 'rsq_val': rsq_list, 'p_val': pvalue_list, 'slope': slope_list,'spearman_val': spearmanr_list}
                stats_df =  pd.DataFrame(statsDict)

                good_rsq_value = stats_df.loc[stats_df['p_val'].idxmin(),'rsq_val']
                good_rsq_valueStr = str(round(good_rsq_value,4))
                good_p_value = stats_df.loc[stats_df['p_val'].idxmin(),'p_val']
                good_p_valueStr = str(round(good_p_value,4))
                good_slope = stats_df.loc[stats_df['p_val'].idxmin(),'slope']  
                good_spearmanr = stats_df.loc[stats_df['p_val'].idxmin(),'spearman_val']
                good_spearmanr_str = str(good_spearmanr)    

                bestFitSubstrate = str(stats_df.loc[stats_df['p_val'].idxmin(),'substrate'])
                concatList = [merged_df, outList_df, stats_df]
                exportDF = pd.concat(concatList, axis=1)
                
                outFileCt += 1
                outFileCt_str = str(outFileCt)

                # export details of current grouping
                outList = [group, headerStr, headerLoop]

                outList_df = pd.DataFrame({ 'Set info': outList})

                statsDict = {struct_class_col:class_list,'substrate':substrateList, 'rsq_val': rsq_list, 'p_val': pvalue_list, 'spearman_val': spearmanr_list,'slope': slope_list}
                stats_df =  pd.DataFrame(statsDict)
                summary_df_list.append(stats_df)

                run_outfnm = '_'.join([outFileCt_str, str(fitsFound), bestFitSubstrate,headerLoop,good_spearmanr_str, good_p_valueStr])
                exportDF.to_csv(join(path,run_outfnm+".csv"), sep=',', index=False)

            # export summary
            if len(summary_df_list) > 0:
                exportDF = pd.concat(summary_df_list, axis=0).sort_values('p_val',ascending=True)
                # pivot wider
                export_df_pivot = pivot_df(exportDF,'Struct_Class','substrate','p_val')
                export_df_pivot.to_csv(summary_outfnm,index=False)
            else:
                pass

        ### final regressions plotting total GUS abundance against rates 
        allGUS_DF = intensity_df
        uniqueEntriesStr = str(len(allGUS_DF[header_col].unique()))
        uniqueEntries = allGUS_DF[header_col].unique()

        # output df used before grouping by sample_ID
        outputDF = allGUS_DF

        # sum abundances by matching sample_ID to prepare for regression
        allGUS_DF = allGUS_DF.groupby(['sample_ID'], as_index=False)[intensity_col].sum()

        new_n = len(allGUS_DF['sample_ID'])
        new_nStr = str(new_n)
        regression_df = allGUS_DF
        
        regression_df[intensity_col] = regression_df[intensity_col].map(lambda a: math.log2(a))

        merged_df = pd.merge(regression_df, rate_df, on = "sample_ID", how='inner')
        fitFound = False
        fitsFound = 0
        rsq_list = []
        pvalue_list = []
        slope_list = []
        spearmanr_list = []

        substrateList = list(rate_df_trunc.columns)

        hit_cols = []

        for response_col in rate_df_trunc:
            y = merged_df[response_col]
            X = merged_df[intensity_col]
            slope, intercept, r_value, p_value, std_err = stats.linregress(X,y)
            spearmanr = stats.spearmanr(X,y).statistic
            rsq_value = r_value**2
            rsq_list.append(rsq_value)
            pvalue_list.append(p_value)
            slope_list.append(slope)
            spearmanr_list.append(spearmanr)

            p_valueStr = str(round(p_value,4))
            rsq_valueStr = str(round(rsq_value,4))
            if (p_value <= p_val_thresh) and (rsq_value >= rsq_thresh) and (slope>slope_thresh):
                fitFound = True
                fitsFound += 1
                hit_cols.append(response_col)
    
        outList_df = pd.DataFrame({ 'Set info': ['all_GUS']})
        statsDict = {'substrate':substrateList, 'rsq_val': rsq_list, 'p_val': pvalue_list,'spearman_val': spearmanr_list, 'slope': slope_list}
        stats_df =  pd.DataFrame(statsDict)
        
        good_rsq_value = stats_df.loc[stats_df['p_val'].idxmin(),'rsq_val']
        good_rsq_valueStr = str(round(good_rsq_value,4))
        good_p_value = stats_df.loc[stats_df['p_val'].idxmin(),'p_val']
        good_p_valueStr = str(round(good_p_value,4))
        good_slope = stats_df.loc[stats_df['p_val'].idxmin(),'slope']  
        good_spearmanr = stats_df.loc[stats_df['p_val'].idxmin(),'spearman_val']
        good_spearmanr_str = str(good_spearmanr)
        bestFitSubstrate = str(stats_df.loc[stats_df['p_val'].idxmin(),'substrate'])

        stats_df = stats_df.sort_values(['p_val'])

        stats_df.to_csv(join(path,"allGUS_stats.csv"), sep=',', index=False)
        concatList = [merged_df, outList_df, stats_df]
        exportDF = pd.concat(concatList, axis=1)
        outFileCt += 1
        outFileCt_str = str(outFileCt)
        
        to_summary_df = pd.DataFrame({struct_class_col:[headerLoop],'best_substrate':[bestFitSubstrate], 'best_spearman':[good_spearmanr], 'best_rsq': [good_rsq_valueStr], 'best_p_val': [good_p_valueStr], 'best_slope': [good_slope], 'num_fits_at_thresh':[str(fitsFound)]})

        summary_df_list.append(to_summary_df)

        exportDF.to_csv(join(path,"allGUS_details.csv"), sep=',', index=False)

    end = time.time()
    totalTime = round((end - start)/60,2)