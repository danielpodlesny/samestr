import pandas as pd
import numpy as np
from functools import reduce
from os.path import basename
import re

import logging

LOG = logging.getLogger(__name__)


def get_metaphlan_taxonomy_levels():
    """
    Return a dict of the taxonomy levels that are used in MetaPhlAn.
    """
    return(
        {'kingdom':  {'prefix': 'k__', 'index': 0},
            'phylum':   {'prefix': 'p__', 'index': 1},
            'class':    {'prefix': 'c__', 'index': 2},
            'order':    {'prefix': 'o__', 'index': 3},
            'family':   {'prefix': 'f__', 'index': 4},
            'genus':    {'prefix': 'g__', 'index': 5},
            'species':  {'prefix': 's__', 'index': 6},
            'strain':   {'prefix': 't__', 'index': 7}
         }
    )


def read_metaphlan_profile(file_fn):
    """
    Read and parse a table to a dataframe, dropping the 'NCBI_tax_id' column and
    keeping only the 'clade_name' and 'relative_abundance' columns, and adding the
    'file_fn' column to the dataframe.
    """
    with open(file_fn, 'r') as f:
        lines = f.readlines()

    # Determine the format of the input file
    if lines[0].startswith('#SampleID\tMetaphlan_Analysis'):
        data_lines = lines[1:]
        col_names = ['clade_name', 'relative_abundance']
    else:
        data_lines = [line for line in lines if not line.startswith('#')]
        col_names = ['clade_name', 'NCBI_tax_id', 'relative_abundance']

    # Parse the data
    # TODO: what is the \t#?
    data = [line.strip().split('\t#')[0].split('\t')[:3] for line in data_lines] 
    df = pd.DataFrame(data, columns=col_names)
    df['file_fn'] = basename(file_fn)

    # Drop the 'NCBI_tax_id' column if present
    if 'NCBI_tax_id' in df.columns:
        df = df.drop('NCBI_tax_id', axis=1)
    return df

def read_motus_profile(file_fn):
    """
    Read and parse a table to a dataframe, dropping the header and
    keeping only the 'clade_name' and 'relative_abundance', and adding the
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
        exit(0)

    # header check
    if len(header) == 2:
        col_names = ['clade_name', 'relative_abundance']
    elif len(header) == 3:
        col_names = ['clade_name', 'NCBI_tax_id', 'relative_abundance']
    else:
        LOG.error('Failed parsing taxonomy file: %s. Too few/many columns.' % file_fn)
        exit(0)

    # Parse the data
    data_lines = lines[1:]
    data = [line.strip().split('\t')[:3] for line in data_lines]
    df = pd.DataFrame(data, columns=col_names)

    # Reduce clade name to mOTUs ID or keep the entire string if no ID is found
    def extract_or_keep_original(clade_name):
        match = re.search(r'\[([^\]]+)\]$', clade_name)
        return match.group(1) if match else clade_name
    df['clade_name'] = df['clade_name'].apply(extract_or_keep_original)

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
    dfs = [df.set_index('clade_name') for df in dfs]

    # rename columns to file name
    dfs = [df.rename(columns={'relative_abundance': df['file_fn'][0]}).drop(
        columns=['file_fn']) for df in dfs]

    # join the dataframes
    dfs = pd.concat(dfs, axis=1, join='outer').fillna(0)

    # bring the index to a 'clade_name' column
    dfs = dfs.reset_index().rename(columns={'index': 'clade_name'})

    return dfs


def filter_metaphlan_profile(df, tax_level):
    """
    Filter a metaphlan profile dataframe to only show rows where the lowest taxonomic
    level is the one specified.
    """
    # Get the index of the taxonomic level
    tax_levels = get_metaphlan_taxonomy_levels()
    level_index = tax_levels[tax_level]['index']

    # return rows where the count of pipes in the first column is level_index
    return df[df[df.columns[0]].apply(lambda x: x.count('|')) == level_index]


def separate_metaphlan_taxonomy(df, lineage_column='clade_name', all_levels=False):
    """
    Split the 'clade_name' column into separate columns for each taxonomic level,
    using the taxonomic levels as column names and returning only the lowest taxonomic
    level if all_levels is False.
    """
    # get taxonomy level dict
    tax_levels = get_metaphlan_taxonomy_levels()

    # separate the clade column into separate columns within the df
    taxonomy_df = df[lineage_column].str.split('\\|', expand=True)

    # name the columns with the taxonomic levels
    for idx in taxonomy_df.columns:
        for tax_level, tdata in tax_levels.items():
            if tdata['index'] == idx:
                taxonomy_df.rename(columns={idx: tax_level}, inplace=True)
    if not all_levels:
        # keep only the lowest taxonomic level
        taxonomy_df.rename(columns={tax_level: lineage_column}, inplace=True)
        taxonomy_df = taxonomy_df.iloc[:, -1:]

    # Concatenate the split dataframe with the original dataframe,
    # dropping the clade column
    separated_df = pd.concat(
        [taxonomy_df, df.drop(columns=[lineage_column])], axis=1)

    return separated_df


def get_metaphlan_level(file_fns: list = None, tax_level=str, combine=['horizontal', 'vertical'], output_file=None, all_levels=True):
    """
    Read, concat, filter, and split multiple metaphlan profile files into a single dataframe,
    then save the result to a file if specified and return the result as a dataframe.
    """
    # Combine the metaphlan profile files
    if combine == 'horizontal':
        df = merge_horizontal_taxonomic_profile(read_function=read_metaphlan_profile, file_fns=file_fns)
    elif combine == 'vertical':
        df = concat_vertical_taxonomic_profile(read_function=read_metaphlan_profile, file_fns=file_fns)
    else:
        raise ValueError('Combine must be either "horizontal" or "vertical"')

    # Filter the dataframe to only show rows where the lowest taxonomic level is the one specified
    df = filter_metaphlan_profile(df, tax_level)

    # Split the 'clade_name' column into separate columns for each
    df = separate_metaphlan_taxonomy(
        df, lineage_column='clade_name', all_levels=all_levels)

    # save to file if specified
    if output_file:
        df.to_csv(output_file, sep="\t", index=False)

    return df


def get_clade_profile(ifn, db_source=['MetaPhlAn', 'mOTUs']):
    """
    Read taxonomic profile from `ifn`,
    return all clades and relative abundances as dataframe
    """

    if db_source == 'MetaPhlAn':
        profile = get_metaphlan_level(file_fns=[ifn],
                                  tax_level='strain',
                                  combine='horizontal',
                                  all_levels=False)
    elif db_source == 'mOTUs':
        profile = merge_horizontal_taxonomic_profile(read_function=read_motus_profile, file_fns=[ifn])
    return profile

def get_clade_profile_dict(ifn, db_source=['MetaPhlAn', 'mOTUs']):
    """
    Read taxonomic profile from `ifn`,
    return all clades and relative abundances as dict
    """
    profile = get_clade_profile(ifn, db_source=db_source)
    return profile.set_index('clade_name')[profile.columns[1]].to_dict()


def binarize_taxonomic_profile(df):
    """
    Take a taxonomic profile dataframe and return a binary dataframe
    """
    # index on first column
    df = df.set_index(df.columns[0])

    # convert to binary > numeric
    df = df.applymap(lambda x: 1 if x != 0 else 0)

    return df

def get_taxon_cooccurrence(wide_metaphlan_profile):
    """
    Take a horizontally merged, filtered, and (for metaphlan) tax separated taxonomic profile matrix
    and return a dataframe with the co-occurrence of clades 
    for each sample pair and taxonomic level
    this does not return self-co-occurrences and only one side of the matrix
    """
    # convert to binary > numeric, index on first column
    df = binarize_taxonomic_profile(wide_metaphlan_profile)

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


def get_taxon_cooccurrences(wide_metaphlan_profile):
    """
    Take a horizontally merged metaphlan profile matrix
    calculate the co-occurrence of clades for each taxonomic level
    and return a dataframe with co-occurrence counts for each sample pair and taxonomic level
    i.e. shared SGB, shared species, shared genera, shared families, etc.
    this does not return self-co-occurrences and only one side of the matrix
    """
    # get the taxonomic levels
    tax_levels = get_metaphlan_taxonomy_levels()
    coocs = []

    # filter the dataframe to only show rows where the lowest taxonomic level is the one specified
    for tax_level, tdata in tax_levels.items():

        # filter the dataframe to only show rows where the lowest taxonomic level is the one specified
        df = filter_metaphlan_profile(wide_metaphlan_profile, tax_level)
        df = separate_metaphlan_taxonomy(df, all_levels=False)
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


def get_taxon_counts(wide_metaphlan_profile):
    """
    Take a horizontally merged metaphlan profile matrix
    sum the number of taxa for each sample and taxonomic level
    and return them in a dataframe
    """
    # get the taxonomic levels
    tax_levels = get_metaphlan_taxonomy_levels()
    counts = []

    # filter the dataframe to only show rows where the lowest taxonomic level is the one specified
    for tax_level, tdata in tax_levels.items():

        if not tax_level == 'strain':
            # filter the dataframe to only show rows where the lowest taxonomic level is the one specified
            df = filter_metaphlan_profile(wide_metaphlan_profile, tax_level)

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
