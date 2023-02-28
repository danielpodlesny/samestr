from os.path import exists
from os import makedirs
import logging
import glob

from samestr.summarize.read_metaphlan_data import merge_horizontal_metaphlan_tables, get_taxon_counts, get_taxon_cooccurrences, filter_metaphlan_profile, separate_metaphlan_taxonomy
from samestr.summarize.read_samestr_data import read_samestr_data, analyze_strain_events
import pandas as pd
import numpy as np

LOG = logging.getLogger(__name__)


def summarize(args):
    """
    Summarize SameStr results, including:
    - MetaPhlAn profiles
        - Taxonomic cooccurrences
    - SameStr results
        - Strain events
        - Strain cooccurrences
    """

    # create output directory
    if not exists(args['output_dir']):
        makedirs(args['output_dir'])

    # list metaphlan profiles
    mp_file_list = glob.glob(
        args['mp_profiles_dir'] + '/*' + args['mp_profiles_extension'])
    if not mp_file_list:
        raise ValueError('Could not find MetaPhlAn input files at: ' +
                         args['mp_profiles_dir'] + '*' + args['mp_profiles_extension'])

    # read metaphlan profiles and remove profile extension from file names
    LOG.info('Reading {} MetaPhlAn profiles from {}'.format(
        len(mp_file_list), args['mp_profiles_dir']))
    mp_profiles = merge_horizontal_metaphlan_tables(mp_file_list)
    mp_profiles.columns = mp_profiles.columns.str.replace(
        args['mp_profiles_extension'], '', regex=False)

    # get kingdom:species counts and cooccurrences
    mp_counts = get_taxon_counts(mp_profiles)
    mp_cooc = get_taxon_cooccurrences(mp_profiles)

    # filter only species level, expand taxonomy, and remove duplicates
    mp_species = filter_metaphlan_profile(mp_profiles, tax_level='species')
    mp_species = separate_metaphlan_taxonomy(mp_species, all_levels=True)
    mp_taxonomy = mp_species[['kingdom', 'phylum', 'class',
                              'order', 'family', 'genus', 'species']].drop_duplicates()

    # read samestr results
    LOG.info('Reading SameStr data from {}'.format(args['input_dir']))
    sstr_data = read_samestr_data(args['input_dir'])

    # add annotation for strain events and get possible and detected strain cooccurrences
    sstr_data, strain_cooc = analyze_strain_events(sstr_data,
                                                   overlap_threshold=args['aln_pair_min_overlap'],
                                                   similarity_threshold=args['aln_pair_min_similarity'])

    # merge all data, replace NaN with 0 only in the shared_strains and analyzed_strain column
    cooc_data = pd.merge(mp_cooc, strain_cooc, on=['row', 'col'], how='outer')
    cooc_data['analyzed_strain'] = cooc_data['analyzed_strain'].fillna(0)
    cooc_data['shared_strain'] = cooc_data['shared_strain'].fillna(0)
    cooc_data.iloc[:, 2:] = cooc_data.iloc[:, 2:].mask(
        cooc_data.iloc[:, 2:].isna(), other=np.nan).astype('Int64')

    # write output
    mp_taxonomy.to_csv(args['output_dir'] +
                       '/mp_taxonomy.tsv', sep='\t', index=False)
    mp_species.drop(['kingdom', 'phylum', 'class', 'order', 'family', 'genus'], axis=1).to_csv(
        args['output_dir'] + '/mp_species.tsv', sep='\t', index=False)
    mp_counts.to_csv(args['output_dir'] +
                     '/mp_counts.tsv', sep='\t', index=False)
    cooc_data.to_csv(args['output_dir'] +
                     '/sstr_cooccurrences.tsv', sep='\t', index=False)
    sstr_data.to_csv(args['output_dir'] +
                     '/sstr_strain_events.tsv', sep='\t', index=False)

    return args
