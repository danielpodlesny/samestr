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

    # we can have merged input from multiple samples
    # need to count the samples in order to know the
    # required size of the target array
    # keep this as list, so it can be used to avoid
    # additional instance checks downstream
    input_sizes = [
        (len(sample) if isinstance(sample, list) else 1)
        for sample, _ in args['input_files']
    ]

    n_inputs = len(args['input_files'])

    # iterate over samples of clade
    for i, ((sample, file_path), isize) in enumerate(zip(args['input_files'], input_sizes)):
        LOG.debug('Merging input file %s (%s/%s)...', sample, i + 1, n_inputs)
        j = len(clade_samples)

        sample_clade_freqs = load_numpy_file(file_path)

        if clade_freqs is None:
            # final clade_freqs array will have dimensions (n_samples x clade_positions x 4)
            clade_freqs = np.zeros(
                (
                    sum(input_sizes),                     # n_samples
                    sample_clade_freqs.shape[1],   # n_positions
                    sample_clade_freqs.shape[2],   # n_variants   
                )
            )

        try:
            clade_freqs[j:j + isize, :, :] = sample_clade_freqs
        except ValueError:
            LOG.error(
                'Dimensions of arrays do not match: %s (%s|%s)' %
                (sample, sample_clade_freqs.shape, clade_freqs.shape)
            )
            return False

        # clade_samples += (sample if isinstance(sample, list) else (sample,))
        clade_samples += (sample if isize > 1 else (sample,))


    # per clade, save freqs and names to file
    np.savez_compressed(output_file + '.npz', clade_freqs, allow_pickle=True)

    with open(output_file + '.names.txt', 'w') as _out:
        print(*clade_samples, sep="\n", file=_out)
