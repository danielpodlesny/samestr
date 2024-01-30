import pandas as pd
import numpy as np
from functools import reduce
from os.path import basename
import re

import logging

LOG = logging.getLogger(__name__)


def identify_taxonomic_profile_database(file_fn):

    profile_tool = ''
    profile_version = ''

    with open(file_fn, 'r') as f:
        lines = f.readlines()
    header = [l.split('#')[1].strip() for l in lines if l.startswith('#')]

    if header[0].startswith('mpa_'):
        profile_tool = 'MetaPhlAn'
        profile_version = header[0].strip()
    elif any(['motus version' in h for h in header]):
        profile_tool = 'mOTUs'
        profile_str = [h for h in header if 'motus version' in h][0]
        profile_str = [p.strip() for p in profile_str.split('|') if 'motus version' in p][0]
        profile_version = profile_str.replace('motus version','').strip()
    
    if profile_tool == '' or profile_version == '':
        LOG.error('Could not identify data source from taxonomic profile: %s' % file_fn)
        exit(1)
    else:
        LOG.info('Database source and version identified from taxonomic profile: %s, %s' % (profile_tool, profile_version))
        return profile_tool, profile_version

def read_metaphlan_profile(file_fn):
    """
    Read and parse a table to a dataframe, dropping the 'NCBI_tax_id' column and
    keeping only the 'clade' and 'relative_abundance' columns, and adding the
    'file_fn' column to the dataframe.
    """
    with open(file_fn, 'r') as f:
        lines = f.readlines()

    # Determine the format of the input file
    if lines[0].startswith('#SampleID\tMetaphlan_Analysis'):
        data_lines = lines[1:]
        col_names = ['clade', 'relative_abundance']
    else:
        data_lines = [line for line in lines if not line.startswith('#')]
        col_names = ['clade', 'NCBI_tax_id', 'relative_abundance']

    # Parse the data
    data = [line.strip().split('\t#')[0].split('\t')[:3] for line in data_lines] 
    df = pd.DataFrame(data, columns=col_names)
    df['file_fn'] = basename(file_fn)

    # Drop the 'NCBI_tax_id' column if present
    if 'NCBI_tax_id' in df.columns:
        df = df.drop('NCBI_tax_id', axis=1)

    # Split lineage by pipe, keeping only the last level
    df['clade'] = df['clade'].str.split('|').str[-1]

    # Filter for t__ (SGB-level) only and assign the result back to df
    df = df[df['clade'].str.startswith('t__')]

    return df

def read_motus_profile(file_fn):
    """
    Read and parse a table to a dataframe, dropping the header and
    keeping only the 'clade' and 'relative_abundance', and adding the
    'file_fn' column to the dataframe.
    """
    with open(file_fn, 'r') as f:
        lines = f.readlines()

    # Determine the format of the input file
        # #consensus_taxonomy     sample1 sample2
        # #consensus_taxonomy     sample1
        # #consensus_taxonomy     NCBI_tax_id     sample1
        # Leptospira alexanderi [ref_mOTU_v3_00001]       100053  1
    lines = [line for line in lines if not line.startswith('# ')]
    header = lines[0].split('\t')
    if header[0] != '#consensus_taxonomy' or len(header) > 2 and header[1] != 'NCBI_tax_id':
        LOG.error('Failed parsing taxonomy file: %s. Could not identify headers or merged profile.' % file_fn)
        exit(1)

    # header check
    if len(header) == 2:
        col_names = ['clade', 'relative_abundance']
    elif len(header) == 3:
        col_names = ['clade', 'NCBI_tax_id', 'relative_abundance']
    else:
        LOG.error('Failed parsing taxonomy file: %s. Too few/many columns.' % file_fn)
        exit(1)

    # Parse the data
    data_lines = lines[1:]
    data = [line.strip().split('\t')[:3] for line in data_lines]
    df = pd.DataFrame(data, columns=col_names)

    # Reduce clade name to mOTUs ID or keep the entire string if no ID is found
    def extract_or_keep_original(clade_name):
        match = re.search(r'\[([^\]]+)\]$', clade_name)
        return match.group(1) if match else clade_name
    df['clade'] = df['clade'].apply(extract_or_keep_original)

    # Convert 'relative_abundance' to numeric and filter rows where it's greater than 0
    df['relative_abundance'] = pd.to_numeric(df['relative_abundance'], errors='coerce')
    df = df[df['relative_abundance'] > 0]

    # Add file name as a column
    df['file_fn'] = basename(file_fn)

    # Drop the 'NCBI_tax_id' column if it's present
    df = df.drop('NCBI_tax_id', axis=1, errors='ignore')

    return df

def concat_vertical_taxonomic_profile(read_function, file_fns: list = None):
    """
    Read and concatenate (vertical) multiple taxonomic profile files into a single dataframe.
    """
    dfs = [read_function(file_fn) for file_fn in file_fns]
    return pd.concat(dfs, ignore_index=True)


def merge_horizontal_taxonomic_profile(read_function, file_fns: list = None):
    """
    Read and merge (horizontal) multiple taxonomic profile files into a single dataframe,
    so that the column name for Abundance is the file name.
    """
    dfs = [read_function(file_fn) for file_fn in file_fns]

    # set index
    dfs = [df.set_index('clade') for df in dfs]

    # rename columns to file name
    dfs = [df.rename(columns={'relative_abundance': df['file_fn'][0]}).drop(
        columns=['file_fn']) for df in dfs]

    # join the dataframes
    dfs = pd.concat(dfs, axis=1, join='outer').fillna(0)

    # bring the index to a 'clade' column
    dfs = dfs.reset_index().rename(columns={'index': 'clade'})

    return dfs

def get_clade_profile(ifn, db_source=['MetaPhlAn', 'mOTUs']):
    """
    Read taxonomic profile from `ifn` list,
    return all clades and relative abundances as pandas df
    """
    if db_source == 'MetaPhlAn':
        fun = read_metaphlan_profile
    elif db_source == 'mOTUs':
        fun = read_motus_profile

    return merge_horizontal_taxonomic_profile(read_function=fun, file_fns=ifn)

def get_clade_profile_dict(ifn, db_source=['MetaPhlAn', 'mOTUs']):
    """
    Read taxonomic profile from `ifn` list,
    return all clades and relative abundances as dict
    """
    profile = get_clade_profile(ifn=ifn, db_source=db_source)

    return profile.set_index('clade')[profile.columns[1]].to_dict()


def binarize_taxonomic_profile(df):
    """
    Take a taxonomic profile dataframe and return a binary dataframe
    """
    # convert to binary > numeric
    return df.applymap(lambda x: int(bool(x)))


def get_taxon_cooccurrence(wide_taxonomic_profile):
    """
    Take a horizontally merged, filtered, and taxonomy annotated profile matrix
    and return a dataframe with the co-occurrence of taxa 
    for each sample pair and taxonomic level
    this does not return self-co-occurrences and only one side of the matrix
    """
    # convert to binary > numeric, index on first column
    df = binarize_taxonomic_profile(wide_taxonomic_profile)

    # calculate co-occurrence
    cooc = np.dot(df.T, df)

    # convert to long format
    cooc_long = pd.DataFrame(cooc, columns=df.columns,
                             index=df.columns).stack().reset_index()

    # name the columns "row", "col", and "cooc"
    cooc_long.columns = ['row', 'col', 'cooc']

    # remove self-co-occurrences, one side of the matrix
    cooc_long = cooc_long[cooc_long['row'] < cooc_long['col']]

    return cooc_long

def get_taxon_cooccurrences(wide_taxonomic_profile, db_taxonomy):
    """
    Take a horizontally merged taxonomic profile matrix
    calculate the co-occurrence of clades for each taxonomic level
    and return a dataframe with co-occurrence counts for each sample pair and taxonomic level
    i.e. shared clade, shared species, shared genera, shared families, etc.
    this does not return self-co-occurrences and only one side of the matrix
    """
    # join dataframe with taxonomy by 'clade'
    df_full = pd.merge(wide_taxonomic_profile, db_taxonomy, on='clade', how='left')

    # get the available taxonomic levels from db_taxonomy header
    tax_levels = db_taxonomy.columns
    sample_cols = [col for col in df_full.columns if col not in tax_levels]
    coocs = []

    for tax_level in tax_levels:

        # groupsum the dataframe to the specified taxonomic level
        df = df_full[[tax_level] + sample_cols].groupby(tax_level).sum()
        cooc = get_taxon_cooccurrence(df)

        # renaming the cooc column to include the taxonomic level
        cooc.rename(columns={'cooc': 'shared_' + tax_level}, inplace=True)

        # get the cooccurrence
        coocs += [cooc]

    # merge the cooccurrence dataframes by row, col,
    # and fill in the missing values with 0
    cooc = reduce(lambda left, right: pd.merge(left, right, on=[
                  'row', 'col'], how='outer'), coocs).fillna(0)

    return cooc


def get_taxon_counts(wide_taxonomic_profile, db_taxonomy):
    """
    Take a horizontally merged taxonomic profile matrix
    sum the number of taxa for each sample and taxonomic level
    and return them in a dataframe
    """
    # join dataframe with taxonomy by 'clade'
    df_full = pd.merge(wide_taxonomic_profile, db_taxonomy, on='clade', how='left')

    # get the available taxonomic levels from db_taxonomy header
    tax_levels = db_taxonomy.columns
    sample_cols = [col for col in df_full.columns if col not in tax_levels]
    counts = []

    for tax_level in tax_levels:

        # groupsum the dataframe to the specified taxonomic level
        df = df_full[[tax_level] + sample_cols].groupby(tax_level).sum()

        # get the number of taxa for each sample
        count = binarize_taxonomic_profile(df).sum(axis=0)
        count.rename('count_' + tax_level, inplace=True)
        counts += [count]

    # merge the cooccurrence dataframes by row, col,
    # and fill in the missing values with 0
    counts = pd.concat(counts, axis=1)

    # reset index and name col 'Name'
    counts.reset_index(inplace=True)
    counts.rename(columns={'index': 'Name'}, inplace=True)

    return counts
