
from os.path import basename, join, exists
import logging
import numpy as np

from Bio import Seq, SeqRecord, AlignIO
from Bio.Align import MultipleSeqAlignment

from samestr.filter import consensus

LOG = logging.getLogger(__name__)


def compare(args):

    # if exists, skip
    output_name = join(args['output_dir'], basename(args['input_file']))
    if exists(output_name):
        LOG.info('Skipping %s. Output file exists.' % args['species'])
        return True

    # load sample order
    with open(args['input_name'], 'r') as file:
        samples = file.read().strip().split('\n')

    # skip if fewer than args['samples_min_n'] samples
    if not len(samples) > 1:
        return False

    LOG.info('Comparing %s found in %s samples.' %
             (args['species'], len(samples)))

    # load freqs
    x = np.load(args['input_file'], allow_pickle=True)
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

    # rename samples to samples.dom
    if args['dominant_variants'] or args['dominant_variants_added']:
        dom_samples = [s + '.dom' for s in samples]
        d = consensus(x)

    # analyze only dominant variants
    if args['dominant_variants']:
        x = d
        samples = dom_samples

    # add dominant variants separately as .dom samples
    elif args['dominant_variants_added']:
        x = np.append(x, d, axis=0)
        samples += dom_samples

    # set matrix colnames
    columns = 'Sample\t' + '\t'.join(samples)
    shortest_distance_matrix = [columns]
    shared_overlap_matrix = [columns]
    fraction_phenotype_matrix = [columns]

    # generate matrix: minimum variant similarity, overlap, fraction mvs
    for i in range(x.shape[0]):
        sample = samples[i]
        shortest_distance = (((x[i, :, :] > 0) *
                              (x > 0)).sum(axis=2) > 0).sum(axis=1)
        shared_overlap = ((x[i, :, :].sum(axis=1) > 0) *
                          (x.sum(axis=2) > 0)).sum(axis=1)
        fraction_phenotype = np.nan_to_num(shortest_distance / shared_overlap)
        shortest_distance_matrix.append(
            sample + '\t' +
            '\t'.join([str(n) for n in np.asarray(shortest_distance)]))
        shared_overlap_matrix.append(
            sample + '\t' +
            '\t'.join([str(n) for n in np.asarray(shared_overlap)]))
        fraction_phenotype_matrix.append(
            sample + '\t' +
            '\t'.join([str(n) for n in np.asarray(fraction_phenotype)]))

    # write matrix files
    with open('%s/%s.closest.txt' % (args['output_dir'], args['species']),
              'w') as out:
        out.write('\n'.join(shortest_distance_matrix))

    with open('%s/%s.overlap.txt' % (args['output_dir'], args['species']),
              'w') as out:
        out.write('\n'.join(shared_overlap_matrix))

    with open('%s/%s.fraction.txt' % (args['output_dir'], args['species']),
              'w') as out:
        out.write('\n'.join(fraction_phenotype_matrix))

    # output dominant variants as msa
    if 'dominant_variants_msa' in args and args['dominant_variants_msa']:

        if args['dominant_variants'] or args['dominant_variants_added']:
            samples = dom_samples
        else:
            d = consensus(x)

        # convert seqs to Bio SeqIO alignment object
        seqs_list = []
        for idx, freqs in enumerate((d > 0).astype(int).tolist()):
            seq = []
            for freq in freqs:
                try:
                    seq.append(acgt[n_freq.index(freq)])
                except ValueError:
                    max_idx = np.array([
                        i for i in range(0, len(freq)) if freq[i] == max(freq)
                    ])
                    random_max = np.random.choice(
                        max_idx, 1)[0] + 2  # for - and N in n_freq
                    seq.append(acgt[random_max])
            seqs_list.append(
                SeqRecord.SeqRecord(id=samples[idx],
                                    description=samples[idx],
                                    seq=Seq.Seq(''.join(seq))))
        seqs_msa = MultipleSeqAlignment(seqs_list)

        # write alignment fasta
        msa_filename = args['output_dir'] + '/' + args['species'] + '.msa.fa'
        with open(msa_filename, 'w') as out:
            AlignIO.write(seqs_msa, out, 'fasta')
