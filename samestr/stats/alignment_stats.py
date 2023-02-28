
from os.path import basename, join, exists
import logging
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from samestr.filter import consensus

LOG = logging.getLogger(__name__)


def coverage(x):
    # calculate coverage for each sample
    # returns MxN numpy array
    return x.sum(axis=2)


def aln2stats(args):

    # if exists, skip
    output_name = join(args['output_dir'], basename(args['input_file']))
    if exists(output_name):
        LOG.info('Skipping %s. Output file exists.' % args['species'])
        return True

    # load sample order
    with open(args['input_name'], 'r') as file:
        samples = file.read().strip().split('\n')

    LOG.info('Gathering stats for %s found in %s samples.' %
             (args['species'], len(samples)))

    # load freqs
    x = np.load(args['input_file'], allow_pickle=True)
    total_species_markers_size = x.shape[1]
    np.seterr(divide='ignore', invalid='ignore')

    # conversion arrays
    acgt = '-NACGT'
    n_freq = [
        [0, 0, 0, 0],  # -
        [0, 0, 0, 0],  # N
        [1, 0, 0, 0],  # A
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
        [0, 0, 0, 1]  # T
    ]
    null_array = np.array([0, 0, 0, 0])

    # get dominant variants
    d = consensus(x)

    if args['dominant_variants']:
        # analyze only dominant variants
        x = d

    # suppress np warnings All-NaN Slice and Mean of empty slice
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', '', RuntimeWarning)

        # stats: vertical coverage
        cov = coverage(x)
        cov[cov == 0] = np.nan
        mean_cov = np.nanmean(cov, axis=1)
        mean_cov[np.where(np.isnan(mean_cov))] = 0
        median_cov = np.nanmedian(cov, axis=1)
        median_cov[np.where(np.isnan(median_cov))] = 0

        # stats: horizontal coverage
        n_sites = np.repeat(x.shape[1], x.shape[0])
        n_gaps = np.isnan(cov).sum(axis=1)
        n_covered = n_sites - n_gaps

        # stats: n of variant sites, monomorphic, .., polymorphic
        p_mono = ((x > 0).sum(axis=2) == 1)
        n_mono = p_mono.sum(axis=1)
        n_duo = ((x > 0).sum(axis=2) == 2).sum(axis=1)
        n_tri = ((x > 0).sum(axis=2) == 3).sum(axis=1)
        n_quat = ((x > 0).sum(axis=2) == 4).sum(axis=1)
        n_poly = ((x > 0).sum(axis=2) > 1).sum(axis=1)

        # polymorphic sites as per binomial cum. dist. func.
        illumina_error_rate = 0.3 / 100  # Q25+
        segata_error_rate = 1 / 100  # Q20
        p_value = 0.05
        n_binom = np.nansum(stats.binom.cdf(
            x.max(axis=2), cov, 1.0 - illumina_error_rate) < p_value,
            axis=1)
        f_binom = n_binom / n_covered
        n_binom_segata = np.nansum(stats.binom.cdf(
            x.max(axis=2), cov, 1.0 - segata_error_rate) < p_value,
            axis=1)
        f_binom_segata = n_binom_segata / n_covered

        # stats: fraction of covered sites,
        # stats: fraction of covered sites with variant, monomorphic, .., polymorphic
        f_covered = n_covered / n_sites
        f_mono = n_mono / n_covered
        f_duo = n_duo / n_covered
        f_tri = n_tri / n_covered
        f_quat = n_quat / n_covered
        f_poly = n_poly / n_covered

        # stats: vertical coverage of dominant variants
        # at all sites
        dom_cov = coverage(d)
        dom_cov[dom_cov == 0] = np.nan
        f_dom_cov = dom_cov / cov
        mean_dom_cov = np.nanmean(dom_cov, axis=1)
        mean_f_dom_cov = np.nanmean(f_dom_cov, axis=1)
        median_f_dom_cov = np.nanmedian(f_dom_cov, axis=1)

        # at polymorphic sites
        dom_cov[p_mono] = np.nan
        f_dom_cov = dom_cov / cov
        mean_dom_cov_polysites = np.nanmean(dom_cov, axis=1)
        median_dom_cov_polysites = np.nanmedian(dom_cov, axis=1)
        mean_f_dom_cov_polysites = np.nanmean(f_dom_cov, axis=1)
        median_f_dom_cov_polysites = np.nanmedian(f_dom_cov, axis=1)

        # mean coverage at polymorphic sites
        p_mono = np.where(np.isnan(dom_cov))
        cov[p_mono] = np.nan
        mean_cov_polysites = np.nanmean(cov, axis=1)

    # convert to pandas df
    df = pd.DataFrame(data=[
        np.array(samples), mean_cov, median_cov, n_sites, n_gaps, n_covered,
        n_mono, n_duo, n_tri, n_quat, n_poly, f_covered, f_mono, f_duo, f_tri,
        f_quat, f_poly, mean_dom_cov, mean_f_dom_cov, median_f_dom_cov,
        mean_dom_cov_polysites, median_dom_cov_polysites,
        mean_f_dom_cov_polysites, median_f_dom_cov_polysites,
        mean_cov_polysites, n_binom, f_binom, n_binom_segata, f_binom_segata
    ])
    df = df.T
    df.columns = [
        'Sample', 'mean_cov', 'median_cov', 'n_sites', 'n_gaps', 'n_covered',
        'n_mono', 'n_duo', 'n_tri', 'n_quat', 'n_poly', 'f_covered', 'f_mono',
        'f_duo', 'f_tri', 'f_quat', 'f_poly', 'mean_dom_cov', 'mean_f_dom_cov',
        'median_f_dom_cov', 'mean_dom_cov_polysites',
        'median_dom_cov_polysites', 'mean_f_dom_cov_polysites',
        'median_f_dom_cov_polysites', 'mean_cov_polysites', 'n_binom',
        'f_binom', 'n_binom_segata', 'f_binom_segata'
    ]

    # write df to file
    ofn = '%s/%s.aln_stats.txt' % (args['output_dir'], args['species'])
    df.to_csv(ofn, sep='\t', index_label=False, index=False)
