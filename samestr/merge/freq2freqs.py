from os.path import isfile
import logging
import numpy as np

from samestr.utils import load_numpy_file

LOG = logging.getLogger(__name__)


def merge_freqs(freqs, freq, freq_name):
    """Appends two numpy arrays."""

    if freqs.shape[1] == freq.shape[1]:
        _freqs = np.append(freqs, freq, axis=0)
        if _freqs.shape[0] == freqs.shape[0] + freq.shape[0]:
            return _freqs, True
        else:
            LOG.error('Could not add array: %s' % freq_name)
    else:
        LOG.error('Dimensions of arrays do not match: %s (%s|%s)' %
                  (freq_name, freqs.shape[1], freq.shape[1]))

    return freqs, False


def freq2freqs(args):
    """Merges nucleotide frequencies for individual species."""

    species_freqs = None
    samples = []

    # skip if species file exists
    output_file = '%s/%s' % (args['output_dir'], args['species'])
    if isfile('%s.npz' % output_file):
        LOG.info('Skipping: %s. Output file existed.' % args['species'])
        return True

    # iterate samples of species
    for sample, file_path in args['input_files']:

        # initialize species freqs arrays
        if species_freqs is None:
            species_freqs = load_numpy_file(file_path)
            if isinstance(sample, list):
                samples += sample
            else:
                samples.append(sample)

            # if more than one sample, continue samples loop
            if len(args['input_files']) > 1:
                LOG.info('Merging %s inputs: %s' %
                         (len(args['input_files']), args['species']))
                continue

        # if more than one sample, add sample freq to freqs
        if len(args['input_files']) > 1:
            species_freqs, success = merge_freqs(
                species_freqs, load_numpy_file(file_path), sample)
            if success:
                if isinstance(sample, list):
                    samples += sample
                else:
                    samples.append(sample)

    # per species, save freqs and names to file
    np.savez_compressed(output_file + '.npz', species_freqs, allow_pickle=True)

    with open(output_file + '.names.txt', 'w') as file:
        file.write('\n'.join(samples) + '\n')
