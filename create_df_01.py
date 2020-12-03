#!/usr/bin/env python3
#coding: utf-8

""" This script creates the dataframe from the merge between manifest and igs1 files. 
Then the dataframe is trimmed down to the relevant numeric features
to create the finished dataframe, df_num. 

Author: James Sacco 
Date: 04/14/2020
Contact email: jsacco001@gmail.com
"""


# import modules
import numpy as np 
import pandas as pd
import csv

import sys

from datetime import datetime, date, time, timedelta

'''Set important variables'''
cellType = 'Bulk PBMC'
stim = 'a-CD3' # Stimulus.in.Readout

def manifest_reader(fname):
    '''Read sample manifests'''
    manifest = pd.read_csv(fname, parse_dates=['MNFD','MNF01', 'MNFTM'], infer_datetime_format=['MNFTM'])

    return manifest

def igs1_reader(fname):
    '''Read IGS1 file'''
    cd3 = pd.read_csv(fname, index_col=False, infer_datetime_format=['Sample.Date'])

    # check that only one col has your stimulus
    assert(len(cd3.columns[cd3.isin(['a-CD3']).sum() > 0]) == 1)

    # rename some cols for easier processing
    cd3.rename(columns = {'Patient.ID': 'PATNUM', 'Visit': 'VISIT', 'Mean.Spot.Count': 'counts'}, inplace =True)

    # filter for 'type of cell' column
    typecell_col = ''.join([col for col in [col for col in cd3.columns if "Cells" in col] if "Type" in col]) 

    # find and filter stimulus column. chained conditions
    stim_col = ''.join(cd3.columns[cd3.isin(['a-CD3']).sum() > 0])
    print(stim_col) 

    # shape filtered on a-CD3
    cd3[cd3[stim_col] == stim].shape

    # shape filtered on a-CD3
    cd3[cd3[stim_col] == stim].shape

    # shape filtered on desired type of cell (Bulk PBMC)
    cd3[cd3[typecell_col] == cellType].shape

    # isolate Bulk PBMC
    cd3 = cd3[cd3[typecell_col] == cellType]

    # check bulk PBMC isolated
    cd3[typecell_col].value_counts()

    # isolate a-CD3
    cd3 = cd3[cd3[stim_col] == stim]

    # drop all from cd3 except stimulus.in.readout, patnum, visit, and sample.id, type of cells, and mean spot count ('counts')
    cd_keepers = ['PATNUM', 'VISIT', 'Sample.ID', 'Sample.Date', typecell_col, stim_col, 'counts']
    print(cd_keepers)

    cd3 = cd3[cd_keepers]

    cd3.shape

    # number of unique patients
    len(cd3.PATNUM.unique())
    
    return cd3

def drop_mnf_cols(fname)
    '''Merge MNF1 and MNF2 on common columns'''
    mnf_keepers = ['PATNUM', 'MNFD', 'MNFTM', 'VISIT', 'MNFTPT', 'MNFALID', 'MNFCNTRY','MNFSITE', 'MNFSPRDT', 'MNFBIOM', 'MNFTUBEN', 'MNF01', 'MNFREFID', 'MNF06', 'MNF14', 'MNF15', 'MNF16', 'MNF17']

    return fname[mnf_keepers]

def merge_mnfs(mnf1, mnf2):
    '''Merge MNF1 and MNF2'''
    merged = mnf1_dropped.append(mnf2_dropped, ignore_index=True)
    return merged



def time_format(merged):
    '''Format time columns'''
    # Derive Turn-Around Time for Sample Manifests
    # Make a collection time column 
    # Create collection column and reformat
    merged['MNFTM']  = merged.MNFTM.dt.time

    merged['MNFTM'] = pd.to_timedelta(merged['MNFTM'].astype(str))

    merged["Collection"] = pd.to_datetime(merged.pop('MNFD')) + pd.to_timedelta(merged.pop('MNFTM').astype(str))

    # create TAT column
    merged["TAT"] = merged["MNF01"] - merged["Collection"]
    merged['TAT'] = merged['TAT']/np.timedelta64(1,'h')

    return merged

def biom_filter(dataframe):
    '''Filter MNFBIOM (Biological Matrix)'''
    merged = dataframe[merged['MNFBIOM']== 'PBMC']

    # Merge Manifest Dataframe with CD3 (assay performance dataframe
     merged['PATNUM'].isin(cd3['PATNUM']).value_counts()

    print(pd.merge(cd3, merged, how = 'left', on=['PATNUM', 'VISIT']).shape)

    print(pd.merge(cd3, merged, how = 'left', on=['PATNUM', 'VISIT']).isnull().mean() * 100)

    return merged

def joined(merged, cd3):
    '''Create merged df with cd3 and manifests'''
    drop_cols = ['Type.of.Cells','MNFTPT', 'MNFALID', 'MNFCNTRY', 'MNFSITE', 'MNFSPRDT','MNFBIOM', 'MNFTUBEN', 'MNF01', 'MNFREFID', 'MNF16', 'MNF17', 'Collection']

    combined = cd3.set_index(['PATNUM', 'VISIT']).join(merged.set_index(['PATNUM', 'VISIT']), lsuffix='_caller', rsuffix='_other')

    comb_drop = combined.drop(drop_cols, axis = 1)

    return comb_drop

def bin_counts(dataframe):
    '''Binarize Mean Spot Counts into above/below 500'''
    # No counts above 1300. Upper limit of detection.
    bins = [-1, 500, 1300]
    labels = [1,0]
    df['binned'] = pd.cut(df['counts'], bins=bins, labels=labels)
    try:
        assert min(df['binned']) > 0
    except:
        print('Negative spot counts found!')

    print("Binned class imbalance: \n", df.binned.value_counts())

    print("All patnums, visits, and samples aligned!")

    # check for null binned values
    try:
        assert len(df[df["binned"].isnull()]) == 0
    except:
        print('Null binned values detected!')


def make_numeric(df):
    '''Create numeric dataframe with the binned y category'''
    df_num = df.select_dtypes(include=[np.number])

    return df_num

def reorg_df(df):
    '''Move counts column to last index position'''
    # Rename columns 
    df.rename(columns = {'MNF06': 'viability','MNF14': 'min_temp', 'MNF15': 'max_temp'}, inplace =True)

    cols_at_end = ['counts']
    df_num = df[[c for c in df if c not in cols_at_end] + [c for c in cols_at_end if c in df]]

    return df_num

def main(manifest1, manifest2, igs1, outfile):

    mnf1 = manifest_reader(manifest1)
    mnf2 = manifest_reader(manifest2)

    cd3 = igs1_reader()

    mnf1_d = drop_mnf_cols(mnf1)
    mnf2_d = drop_mnf_cols(mnf2)

    samples_merged = merge_mnfs(mnf1_d, mnf2_d)

    samples_time = time_format(samples_merged)

    samples_biom = biom_filter(samples_time)

    dataframe = joined(samples_biom, cd3)

    counts_binned = bin_counts(dataframe)

    numeric_df = make_numeric(counts_binned)

    combinedFeats = reorg_df(numeric_df)

    combinedFeats.to_csv(outfile, sep = ',', header = True, index = False)

    print('Dataframe created.' + '\n' + 'Done!')

if __name__ = "__main__":
    # Parse command line args
    if len(sys.argv) == 5:
        manifest1 = sys.argv[1]
        manifest2 = sys.argv[2]
        igs1 = sys.arg[3]
        outfile = sys.argv[4]
        main(manifest1 = infile1, manifest2 = infile2, igs1 = igs1, outfile = outfile)
    else:
        print("Invalid number of arguments. Please provide:\n2 input sample manifests, 1 IGS1 file, and output file name")

