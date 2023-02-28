
from os.path import basename, join, exists
import logging
import warnings
import numpy as np

LOG = logging.getLogger(__name__)


def read_marker_positions(species_marker_file):
    """Reads species marker file.

    Returns: Dict with start, end positions and len for each marker.
    """
    marker_positions = {}
    with open(species_marker_file, 'r') as marker_pos:
        header = next(marker_pos)
        for line in marker_pos:
            line = line.rstrip().split()
            marker = str(line[0]).strip()
            marker_start = int(line[1])
            marker_len = int(line[2])
            marker_end = marker_start + marker_len
            marker_positions[marker] = [marker_start, marker_end, marker_len]
    return marker_positions


def read_marker_list(marker_list_file):
    """Reads species marker list file.

    Returns: Set of marker names.
    """
    marker_list = set()
    with open(marker_list_file) as f:
        for l in f.readlines():
            marker_list.add(l.strip())
    return marker_list


def trunc_marker_ends(marker_pos, trunc_len):
    """Sets positions for truncation of marker ends.

    Returns: List of positions which are supposed to be
             truncated given `marker_pos` and `trunc_len`.
    """
    trunc_pos = set()
    for marker, (marker_start, marker_end, marker_len) in list(marker_pos.items()):

        # drop entirely if trunc_len * 2 longer than marker:
        if trunc_len * 2 > marker_len:
            LOG.info('Marker %s shorter than total `trunc_len`' % marker)
            trunc_pos = trunc_pos.union(list(range(marker_start, marker_end)))
            continue

        # set trunc positions
        trunc_pos_start = list(range(marker_start, marker_start + trunc_len))
        trunc_pos_end = list(range(marker_end - trunc_len, marker_end))
        trunc_pos = trunc_pos.union(trunc_pos_start, trunc_pos_end)
    return trunc_pos


def remove_marker_pos(marker_pos, remove_markers=None):
    """Sets positions for marker removal.

    Returns: Set of positions in `freqs` which are supposed to be
             removed given `remove_markers`.
    """
    remove_pos = set()
    for marker, (marker_start, marker_end, marker_len) in list(marker_pos.items()):
        if marker in remove_markers:
            remove_pos = remove_pos.union(
                list(range(marker_start, marker_end)))
    return remove_pos


def coverage(x, zero_nan=False):
    """Returns MxN numpy array with coverage for each sample."""
    cov = x.sum(axis=2)
    if zero_nan:
        cov[cov == 0] = np.nan
    return cov


def z_coverage(x, covered_pos_only=False):
    """Returns MxN numpy array with standardized coverage for each sample."""
    cov = coverage(x, zero_nan=covered_pos_only)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        zcov = (cov - np.nanmean(cov, axis=1)[:, np.newaxis]) / np.nanstd(
            cov, axis=1)[:, np.newaxis]
    return zcov


def species_min_samples(args, n_samples):
    if args['species_min_samples']:
        if n_samples < args['species_min_samples']:
            LOG.debug('Skipping %s since there are fewer than %s samples.' %
                      (args['species'], args['species_min_samples']))
            return False
    return True


def filter_freqs(args):

    # if exists, skip
    output_name = join(args['output_dir'], basename(args['input_file']))
    if exists(output_name):
        LOG.info('Skipping %s. Output file exists.' % args['species'])
        return True

    # load sample order
    with open(args['input_name'], 'r') as file:
        samples = file.read().strip().split('\n')

    # load sample selection
    if args['input_select'] is not None:
        with open(args['input_select'], 'r') as file:
            input_select = file.read().strip().split('\n')

    # skip if fewer than args['species_min_samples'] samples
    if not species_min_samples(args, len(samples)):
        return False

    # load freqs
    x = np.load(args['input_file'], allow_pickle=True)
    total_species_markers_size = x.shape[1]
    np.seterr(divide='ignore', invalid='ignore')
    removed_pos = set()

    # get marker metadata
    marker_metadata_file = join(args['marker_dir'],
                                '%s.positions.txt' % args['species'])
    marker_positions = read_marker_positions(marker_metadata_file)

    # 0 Remove Samples
    if args['input_select'] is not None:
        keep_samples = [n for n in range(
            0, len(samples)) if samples[n] in input_select]
        remove_n = len(samples) - len(keep_samples)
        samples = [samples[n] for n in keep_samples]
        LOG.info('Keeping %s out of %s selected sample(s). Removing: %s.' %
                 (len(keep_samples), len(input_select), remove_n))
        x = x[keep_samples, :, :]

    LOG.info('Filtering %s found in %s samples.' %
             (args['species'], len(samples)))

    # 1 Filter Markers

    # 1.1 Truncate ends
    if args['marker_trunc_len']:
        trunc_pos = trunc_marker_ends(marker_positions,
                                      args['marker_trunc_len'])

        LOG.info('Truncating marker ends at %s [%s%%] global positions.' %
                 (len(trunc_pos),
                  round(len(trunc_pos) / total_species_markers_size, 1)))

        # x has to be float (extract -> float not int)
        x[:, list(trunc_pos), :] = np.nan

        removed_pos = removed_pos.union(trunc_pos)

    # 1.2 Keep OR Remove markers
    if args['marker_keep'] or args['marker_remove']:

        if args['marker_keep']:
            # available marker set difference w/ markers
            keep_markers = read_marker_list(args['marker_keep']) & set(
                marker_positions.keys())
            rm_marker_set = set(keep_markers) ^ set(marker_positions.keys())

        elif args['marker_remove']:
            # available marker set intersect w/ markers
            rm_marker_set = read_marker_list(args['marker_remove']) & set(
                marker_positions.keys())

        remove_pos = remove_marker_pos(marker_positions,
                                       remove_markers=rm_marker_set)

        LOG.info('Zeroing %s markers at %s [%s%%] global positions.' %
                 (len(rm_marker_set), len(remove_pos),
                  round(len(remove_pos) / total_species_markers_size, 1)))
        x[:, list(remove_pos), :] = np.nan
        removed_pos = removed_pos.union(remove_pos)

    # 2 Filter Sample Variants [per Position]

    # 2.1 min-ncov
    if args['sample_var_min_n_vcov']:
        remove_pos = x < args['sample_var_min_n_vcov']
        x[remove_pos] = 0

    # 2.2 min-fcov
    if args['sample_var_min_f_vcov']:
        remove_pos = np.true_divide(
            x,
            coverage(x)[:, :, np.newaxis]) < args['sample_var_min_f_vcov']
        x[remove_pos] = 0

    # 3 Filter Positions [per Sample]

    # 3.1 min-vcov
    if args['sample_pos_min_n_vcov']:
        remove_pos = coverage(x) < args['sample_pos_min_n_vcov']
        x[remove_pos, :] = np.nan

    # 3.2 sd-vcov [of covered positions only]
    if args['sample_pos_min_sd_vcov']:
        remove_pos = abs(z_coverage(
            x, covered_pos_only=True)) > args['sample_pos_min_sd_vcov']
        x[remove_pos, :] = np.nan

    # 4 Filter Samples

    # 4.1 min-n or f horizontal coverage
    if args['samples_min_n_hcov']:
        covered = np.sum(coverage(x, zero_nan=True) > 0, axis=1)
        keep_samples = covered > args['samples_min_n_hcov']
        keep_samples_names = set(
            [samples[i] for i, b in enumerate(keep_samples) if b])
        remove_samples_names = keep_samples_names.symmetric_difference(samples)
        LOG.info('Removing %s sample(s) due to insufficient horizontal coverage: %s' %
                 (len(remove_samples_names), ', '.join(remove_samples_names)))
        samples = [s for s in samples if s not in remove_samples_names]
        x = x[keep_samples, :, :]

    if args['samples_min_f_hcov']:
        covered = np.sum(coverage(x, zero_nan=True) > 0, axis=1)
        keep_samples = np.true_divide(
            covered, x.shape[1]) > args['samples_min_f_hcov']
        keep_samples_names = set(
            [samples[i] for i, b in enumerate(keep_samples) if b])
        remove_samples_names = keep_samples_names.symmetric_difference(samples)
        LOG.info('Removing %s sample(s) due to insufficient horizontal coverage: %s' %
                 (len(remove_samples_names), ', '.join(remove_samples_names)))
        samples = [s for s in samples if s not in remove_samples_names]
        x = x[keep_samples, :, :]

    # 4.2 mean-vcov [of covered positions only]
    if args['samples_min_m_vcov']:
        keep_samples = np.nanmean(coverage(x, zero_nan=True),
                                  axis=1) > args['samples_min_m_vcov']
        keep_samples_names = set(
            [samples[i] for i, b in enumerate(keep_samples) if b])
        remove_samples_names = keep_samples_names.symmetric_difference(samples)
        LOG.info('Removing %s sample(s) due to insufficient vertical coverage: %s' %
                 (len(remove_samples_names), ', '.join(remove_samples_names)))
        samples = [s for s in samples if s not in remove_samples_names]
        x = x[keep_samples, :, :]

    # 4.4 Sample minimum
    if not species_min_samples(args, x.shape[0]):
        return False

    # 5 Filter Positions [Global]

    # 5.1 min [n or fraction] samples with coverage at position
    if args['global_pos_min_n_vcov'] or args['global_pos_min_f_vcov']:

        if args['global_pos_min_n_vcov']:
            remove_pos = np.sum(coverage(x) > 0,
                                axis=0) < args['global_pos_min_n_vcov']
        elif args['global_pos_min_f_vcov']:
            remove_pos = (np.sum(coverage(x) > 0, axis=0) /
                          x.shape[0]) < args['global_pos_min_f_vcov']

        LOG.info(
            'Zeroing %s global positions [%s%%] covered by too few samples.' %
            (sum(remove_pos),
             round(sum(remove_pos) / total_species_markers_size, 1)))
        x[:, remove_pos, :] = np.nan
        removed_pos = set([pos for pos, b in enumerate(remove_pos)
                           if b]).union(removed_pos)

    # 6 Filter Mono- and Polymorphic Positions
    if args['keep_poly'] or args['keep_mono']:

        if args['keep_poly']:
            remove_pos = ((x > 0).sum(axis=2) > 1).sum(axis=0) == 0
            LOG.info('Zeroing %s [%s%%] global monomorphic positions.' %
                     (sum(remove_pos),
                      round(sum(remove_pos) / total_species_markers_size, 1)))
        elif args['keep_mono']:
            remove_pos = ((x > 0).sum(axis=2) > 1).sum(axis=0) > 0
            LOG.info('Zeroing %s [%s%%] global polymorphic positions.' %
                     (sum(remove_pos),
                      round(sum(remove_pos) / total_species_markers_size, 1)))

        x[:, remove_pos, :] = np.nan
        removed_pos = set([pos for pos, b in enumerate(remove_pos)
                           if b]).union(removed_pos)

    # Delete positions from array
    if args['delete_pos']:
        x = np.delete(x, list(removed_pos), axis=1)

    LOG.info('Remaining positions: %s.' % x.shape[1])

    # 7 Clean up
    # 7.1 Filter samples with no coverage
    remove_samples = np.sum(coverage(x) > 0, axis=1) == 0
    samples = [s for i, s in enumerate(samples) if not remove_samples[i]]
    x = x[~remove_samples, :, :]
    LOG.info('Removing {} sample(s) with no coverage.'.format(sum(remove_samples)))

    # 7.2 Sample minimum
    if not species_min_samples(args, x.shape[0]):
        return False

    # Save retained samples to file
    samples_file = join(args['output_dir'], basename(args['input_name']))
    with open(samples_file, 'w') as file:
        txt = '\n'.join(samples)
        file.write(txt)

    # Save Array to file
    np.save(output_name, x, allow_pickle=True)
