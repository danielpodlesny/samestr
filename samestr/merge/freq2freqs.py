from os.path import isfile
import logging
import numpy as np

from samestr.utils.utilities import load_numpy_file

LOG = logging.getLogger(__name__)


# def merge_freqs(freqs, freq, freq_name):
#     """Appends two numpy arrays."""

#     if freqs.shape[1] == freq.shape[1]:
#         _freqs = np.append(freqs, freq, axis=0)
#         if _freqs.shape[0] == freqs.shape[0] + freq.shape[0]:
#             return _freqs, True
#         else:
#             LOG.error('Could not add array: %s' % freq_name)
#     else:
#         LOG.error('Dimensions of arrays do not match: %s (%s|%s)' %
#                   (freq_name, freqs.shape[1], freq.shape[1]))

#     return freqs, False

def freq2freqs(args):
    """Merges nucleotide frequencies for individual clades."""

    clade_freqs = None
    clade_samples = []

    # skip if clade file exists
    output_file = '%s/%s' % (args['output_dir'], args['clade'])
    if isfile('%s.npz' % output_file):
        LOG.info('Skipping: %s. Output file existed.' % args['clade'])
        return True
    
    if len(args['input_files']) > 1:
        LOG.info(
            'Merging %s inputs: %s' % (len(args['input_files']), args['clade'])
        )

    # iterate samples of clade
    for i, (sample, file_path) in enumerate(args['input_files']):

        sample_clade_freqs = load_numpy_file(file_path)

        if clade_freqs is None:
            # final clade_freqs array will have dimensions (n_samples x clade_positions x 4)
            clade_freqs = np.zeros(
                (
                    len(args['input_files']),      # n_samples
                    sample_clade_freqs.shape[1],   # n_positions
                    sample_clade_freqs.shape[2],   # n_variants   
                )
            )

        try:
            clade_freqs[i] = sample_clade_freqs
        except ValueError:
            LOG.error(
                'Dimensions of arrays do not match: %s (%s|%s)' %
                (sample, sample_clade_freqs.shape, clade_freqs.shape)
            )
            return False

        clade_samples += (sample if isinstance(sample, list) else (sample,))

    # per clade, save freqs and names to file
    np.savez_compressed(output_file + '.npz', clade_freqs, allow_pickle=True)

    with open(output_file + '.names.txt', 'w') as _out:
        print(*clade_samples, sep="\n", file=_out)
