from os import listdir
import pandas as pd
import numpy as np
import logging

LOG = logging.getLogger(__name__)


def read_distmat(directory, species, suffix):
    """
    Function reads distance matrix from file
    returns distance matrix as data frame.
    """
    fn = f"{directory}/{species}{suffix}"
    distmat = pd.read_csv(fn, sep="\t", index_col=0)
    return distmat


def distmat_nonfinite_to_NA(distmat):
    """
    Function takes distance matrix as data frame as input
    transforms non-finite values to NA, 
    returns distance matrix as data frame.
    """
    distmat = distmat.mask(~np.isfinite(distmat))
    return distmat


def distmat_rm_upper_tri(distmat, rm_diag=False):
    """
    Function takes distance matrix as data frame as input
    and returns distance matrix as data frame without upper triangle.
    It also removes diagonal if rm_diag is True.
    """
    # Remove upper triangle
    distmat = distmat.where(np.tril(np.ones(distmat.shape)).astype(bool))

    # Remove diagonal
    if rm_diag:
        distmat = distmat.where(~np.eye(distmat.shape[0], dtype=bool))
    return distmat


def distmat_to_long(distmat, value_name):
    """
    Function takes distance matrix as data frame as input
    and returns a data frame in long format,
    containing all values from the distance matrix.
    """
    # Convert to long format
    distmat = distmat.stack().reset_index()
    distmat.columns = ['row', 'col', value_name]
    return distmat


def get_merged_distmats(directory, clade):
    """
    Takes directory and clade as input
    and returns a data.frame where file pairs for each clade are merged.

    Requires:
    Input dir must contain samestr compare output files with
    extensions: .fraction.txt and .overlap.txt
    """
    # Read in the fraction and overlap files as dataframes
    distmat = read_distmat(directory, clade, '.fraction.txt')
    overlap = read_distmat(directory, clade, '.overlap.txt')

    # Convert non-finite values to NA
    distmat = distmat_nonfinite_to_NA(distmat)
    overlap = distmat_nonfinite_to_NA(overlap)

    # Convert wide to long format
    distmat_long = distmat_to_long(distmat, value_name='similarity')
    overlap_long = distmat_to_long(overlap, value_name='overlap')
    if distmat_long.shape[0] > 0 and overlap_long.shape[0] > 0:
        # Merge the two dataframes
        merged_distmat_overlap = pd.merge(distmat_long, overlap_long, on=[
                                          'row', 'col'], how='outer')

        # Add a column for clade
        merged_distmat_overlap['clade'] = clade

        # remove self-matches, one side of the matrix
        merged_distmat_overlap = merged_distmat_overlap[merged_distmat_overlap['row']
                                                        < merged_distmat_overlap['col']]
        return merged_distmat_overlap


def read_samestr_data(sstr_dir):
    """
    Takes samestr compare output directory as input and runs get_merged_distmats
    """
    distance_clade_list = [s.replace('.fraction.txt', '') for s in listdir(
        sstr_dir) if s.endswith('.fraction.txt')]
    overlap_clade_list = [s.replace('.overlap.txt', '') for s in listdir(
        sstr_dir) if s.endswith('.overlap.txt')]
    clade_intersect = set(overlap_clade_list).intersection(
        set(distance_clade_list))

    LOG.info(f'Merging matrix for {len(clade_intersect)} clades.')
    sstr_data = pd.DataFrame()
    for clade in clade_intersect:
        sstr_data = pd.concat(
            [sstr_data, get_merged_distmats(directory=sstr_dir, clade=clade)])

    return sstr_data


def analyze_strain_events(sstr_data, similarity_threshold=0.999, overlap_threshold=5000):
    """
    Takes samestr data as input and annotates it with analysis level, shared, and event.
    Also returns a data frame with the sum of strain co-occurrences
    across all clades for each file pair, including zero counts.
    """
    sstr_data['analysis_level'] = sstr_data.apply(
        lambda x: 'strain-level' if x['overlap'] > overlap_threshold else 'clade-level', axis=1)
    sstr_data['analyzed_strain'] = sstr_data['analysis_level'] == 'strain-level'
    sstr_data['shared_strain'] = (sstr_data['similarity'] >= similarity_threshold) & (
        sstr_data['analyzed_strain'] == True)
    sstr_data['event'] = sstr_data.apply(lambda x: 'shared_strain' if x['shared_strain']
                                         else 'other_strain' if x['overlap'] > overlap_threshold else 'same_clade', axis=1)
    sstr_data = sstr_data[['row', 'col', 'clade', 'similarity', 'overlap',
                           'analysis_level', 'analyzed_strain', 'shared_strain', 'event']]

    # get the sum of shared strains and analyzed strains for each file pair
    strain_cooccurrences = sstr_data.groupby(
        ['row', 'col'])[['shared_strain', 'analyzed_strain']].sum().reset_index()
    
    sstr_data = sstr_data.sort_values(["row", "col", "clade"])
        
    return sstr_data, strain_cooccurrences
