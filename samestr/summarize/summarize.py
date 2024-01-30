from os.path import exists
from os import makedirs
import logging
import glob

from samestr.summarize.read_taxonomic_profiles import get_clade_profile, get_taxon_counts, get_taxon_cooccurrences
from samestr.summarize.read_samestr_data import read_samestr_data, analyze_strain_events
import pandas as pd
import numpy as np

LOG = logging.getLogger(__name__)


def summarize(args):
    """
    Summarize SameStr results, including:
    - MetaPhlAn or mOTUs profiles
        - Taxonomic cooccurrences
    - SameStr results
        - Strain events
        - Strain cooccurrences
    """

    # pull out manifest data
    db_manifest = args['samestr_db']['db_manifest']
    db_taxonomy = args['samestr_db']['db_taxonomy']['records']

    # list taxonomic profiles
    tax_profile_list = glob.glob(
        args['tax_profiles_dir'] + '/*' + args['tax_profiles_extension'])
    if not tax_profile_list:
        raise ValueError('Could not find taxonomic profiles at: ' +
                         args['tax_profiles_dir'] + '*' + args['tax_profiles_extension'])

    # read taxonomic profiles and remove profile extension from file names
    LOG.info('Reading {} taxonomic profiles from {}'.format(
        len(tax_profile_list), args['tax_profiles_dir']))
    tax_profiles = get_clade_profile(tax_profile_list, db_source=db_manifest['database']['name'])
    tax_profiles.columns = tax_profiles.columns.str.replace(
        args['tax_profiles_extension'], '', regex=False)
    
    # get kingdom:species counts and cooccurrences
    taxon_counts = get_taxon_counts(tax_profiles, db_taxonomy)
    taxon_cooc = get_taxon_cooccurrences(tax_profiles, db_taxonomy)

    # read samestr results
    LOG.info('Reading SameStr data from {}'.format(args['input_dir']))
    sstr_data = read_samestr_data(args['input_dir'])

    # add annotation for strain events and get possible and detected strain cooccurrences
    sstr_data, strain_cooc = analyze_strain_events(sstr_data,
                                                   overlap_threshold=args['aln_pair_min_overlap'],
                                                   similarity_threshold=args['aln_pair_min_similarity'])

    # merge all data, replace NaN with 0 only in the shared_strains and analyzed_strain column
    cooc_data = pd.merge(taxon_cooc, strain_cooc, on=['row', 'col'], how='outer')
    cooc_data['analyzed_strain'] = cooc_data['analyzed_strain'].fillna(0)
    cooc_data['shared_strain'] = cooc_data['shared_strain'].fillna(0)
    cooc_data.iloc[:, 2:] = cooc_data.iloc[:, 2:].mask(
        cooc_data.iloc[:, 2:].isna(), other=np.nan).astype('Int64')

    # write output
    taxon_counts.to_csv(os.path.join(args['output_dir'],
                     'taxon_counts.tsv', sep='\t', index=False)
    cooc_data.to_csv(args['output_dir'] +
                     '/sstr_cooccurrences.tsv', sep='\t', index=False)
    sstr_data.to_csv(args['output_dir'] +
                     '/sstr_strain_events.tsv', sep='\t', index=False)

    return args
