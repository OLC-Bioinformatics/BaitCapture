#!/usr/bin/env python

"""
Parse a sorted BAM file to generate two tables:
a table of alignment summary statistics ('_stats.tsv'), including fold-coverages
(fold_cov) and proportions of target lengths covered by mapped reads (prop_cov),
and a table of read depths ('_depth.tsv') for each bp position of each target.

Dependencies: python=3.10.8, pandas=1.5.2, samtools=1.16.1
Other package versions may work but are untested.
"""

__author__ = 'Liam Brown'
__email__ = 'liam.brown@inspection.gc.ca'

import os
import sys
import argparse
import subprocess
import pandas as pd

#-------------------------------------------------------------------------------
# parse_arguments()
#-------------------------------------------------------------------------------

def parse_arguments():
    """
    Parse command-line arguments.

    :returns args: List of parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description = """
        Parse a sorted BAM file to generate two tables:
        a table of alignment summary statistics ('_stats.tsv'), including 
        fold-coverages (fold_cov) and proportions of target lengths covered by
        mapped reads (prop_cov), and a table of read depths ('_depth.tsv') for
        each bp position of each target.  
        """)

    # Required arguments
    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-i', '--input', type = str, required = True,
        help = """
        Path to sorted BAM file to process.
        """)
    
    # Optional arguments
    optional_args = parser.add_argument_group('Optional')
    optional_args.add_argument('-o', '--output', type = str, required = False,
        default = 'sample',
        help = """
        Path and base name of files to save as tab-separated tables
        ('[output]_stats.tsv', '[output]_depth.tsv').
        Default: 'sample'
        """)

    # If no arguments provided:
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

#-------------------------------------------------------------------------------
# Other functions
#-------------------------------------------------------------------------------

def compute_alignment_stats(bamfile):
    """
    Compute alignment summary statistics and store in a Pandas dataframe.

    :param bamfile: Path to the sorted BAM file to process.
    :type bamfile: str
    :returns idxstats_df: Alignment summary statistics.
    :rtype idxstats_df: <class 'pandas.core.frame.DataFrame'>
    """
    with open(os.devnull, 'w') as dev_null:
        out = subprocess.check_output(['samtools', 'idxstats', bamfile])

    out = out.decode()
    parsed = [line.split('\t') for line in out.splitlines()]
    idxstats_df = pd.DataFrame(parsed, 
        columns = ['target', 'seqlen', 'mapped_segments', 'unmapped_segments'])
    # Drop the '*' row from the dataframe
    idxstats_df.drop(idxstats_df[idxstats_df['target'] == '*'].index, 
        inplace = True)
    idxstats_df['seqlen'] = pd.to_numeric(idxstats_df['seqlen'])    

    return idxstats_df

def compute_depth(bamfile):
    """
    Summarize read depth for each target and store in a Pandas dataframe.

    :param bamfile: Path to the sorted BAM file to process.
    :type bamfile: str
    :returns depth_df: Per base pair alignment depths.
    :rtype depth_df: <class 'pandas.core.frame.DataFrame'>
    :returns total_depth_df: Per target alignment depths.
    :rtype total_depth_df: <class 'pandas.core.frame.DataFrame'> 
    """
    with open(os.devnull, 'w') as dev_null:
        out = subprocess.check_output(['samtools', 'depth', '-a', args.input])
    out = out.decode()
    parsed = [line.split('\t') for line in out.splitlines()]
    depth_df = pd.DataFrame(parsed,
        columns = ['target', 'position', 'depth'])
    depth_df['depth'] = pd.to_numeric(depth_df['depth'])
    total_depth_df = depth_df.groupby(['target'])['depth'].sum()

    return depth_df, total_depth_df

def compute_len_cov(depth_df):
    """
    Compute the number of bp of each target that are covered by at least one
    mapped read.

    :param depth_df: Per base pair alignment depths.
    :type depth_df: <class 'pandas.core.frame.DataFrame'>
    :returns len_cov_df: Number of base pairs covered by at least one mapped
    read.
    :rtype len_cov_df: <class 'pandas.core.frame.DataFrame'> 
    """
    # If read depth at a position is greater than 0, set to 1
    depth_df.loc[depth_df.depth > 0, 'depth'] = 1
    len_cov_df = depth_df.groupby(['target'])['depth'].sum()
    len_cov_df.name = 'len_cov'

    return len_cov_df

def join_dfs(idxstats_df, total_depth_df, len_cov_df):
    """
    Join the idxstats_df, total_depth_df, and len_cov_df dataframes, and compute
    the proportions of target lengths covered by mapped reads (prop_cov), and
    the fold-coverages of mapped reads to targets (fold_cov).

    :param idxstats_df: Alignment summary statistics.
    :type idxstats_df: <class 'pandas.core.frame.DataFrame'>
    :param total_depth_df: Per target alignment depths.
    :type total_depth_df: <class 'pandas.core.frame.DataFrame'> 
    :param len_cov_df: Number of base pairs covered by at least one mapped
    read.
    :type len_cov_df: <class 'pandas.core.frame.DataFrame'>
    :returns joined_df: idxstats_df, total_depth_df, and len_cov_df dataframes
    joined together with additional columns for fold-coverage and proportion of
    target length covered.
    :rtype joined_df: <class 'pandas.core.frame.DataFrame'>
    """
    # Join dataframes
    joined_df = idxstats_df.join(total_depth_df, on = 'target').fillna(0)
    joined_df = joined_df.join(len_cov_df, on = 'target').fillna(0)
    # Compute prop_cov and fold_cov
    joined_df['prop_cov'] = joined_df['len_cov'] / joined_df['seqlen']
    joined_df['fold_cov'] = joined_df['depth'] / joined_df['seqlen']
    # Drop unnecessary columns
    joined_df.drop(['mapped_segments', 'unmapped_segments'], axis = 1,
        inplace = True)

    return joined_df

def write_dfs(joined_df, depth_df, outbasename):
    """
    Save the dataframes as tab-separated tables (.tsv).

    :param joined_df: idxstats_df, total_depth_df, and len_cov_df dataframes
    joined together with additional columns for fold-coverage and proportion of
    target length covered.
    :type joined_df: <class 'pandas.core.frame.DataFrame'>   
    :param depth_df: Per base pair alignment depths.
    :type depth_df: <class 'pandas.core.frame.DataFrame'>
    :param outbasename: Path and base name of files to save as tab-separated
    tables.
    """
    joined_df.to_csv(f"{outbasename}_stats.tsv", sep = '\t', index = False)
    depth_df.to_csv(f"{outbasename}_depth.tsv", sep = '\t', index = False)

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):

    idxstats_df = compute_alignment_stats(bamfile = args.input)
    depth_df, total_depth_df = compute_depth(bamfile = args.input)
    len_cov_df = compute_len_cov(depth_df)
    joined_df = join_dfs(idxstats_df, total_depth_df, len_cov_df)
    write_dfs(joined_df, depth_df, outbasename = args.output)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)