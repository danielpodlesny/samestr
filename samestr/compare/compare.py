
import logging
import os
from os.path import basename, exists

import numpy as np
from Bio import Seq, SeqRecord, AlignIO
from Bio.Align import MultipleSeqAlignment

from samestr.utils.utilities import load_numpy_file
from samestr.filter import consensus


LOG = logging.getLogger(__name__)


def compare(args):

    # if exists, skip
    output_name = os.path.join(args['output_dir'], basename(args['input_file']))
    if exists(output_name):
        LOG.info('Skipping %s. Output file exists.' % args['clade'])
        return True

    # load sample order
    with open(args['input_name'], 'r') as file:
        samples = file.read().strip().split('\n')

    # skip if fewer than args['samples_min_n'] samples
    if len(samples) < 2:
        return None

    LOG.info('Comparing %s found in %s samples.' %
             (args['clade'], len(samples)))

    # load freqs
    x = load_numpy_file(args['input_file'])
    np.seterr(divide='ignore', invalid='ignore')

    # conversion arrays
    acgt = '-NACGT'
    n_freq = [
        [0, 0, 0, 0],  # -
        [0, 0, 0, 0],  # N
        [1, 0, 0, 0],  # A
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
        [0, 0, 0, 1],  # T
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

    closest_out = open('%s/%s.closest.txt' % (args['output_dir'], args['clade']), 'w')
    overlap_out = open('%s/%s.overlap.txt' % (args['output_dir'], args['clade']), 'w')
    fraction_out = open('%s/%s.fraction.txt' % (args['output_dir'], args['clade']), 'w')

    with closest_out, overlap_out, fraction_out:
        print("Sample", *samples, sep="\t", file=closest_out)
        print("Sample", *samples, sep="\t", file=overlap_out)
        print("Sample", *samples, sep="\t", file=fraction_out)

        non_null_global = (x > 0)
        non_null_positions = (x.sum(axis=2) > 0)

        for i, (sample, sample_matrix) in enumerate(zip(samples, x)):
            LOG.debug("Processing sample %s (%s/%s)...", sample, i + 1, len(samples))
            shortest_distance = (((sample_matrix > 0) *
                                  non_null_global).sum(axis=2) > 0).sum(axis=1)
            shared_overlap = ((sample_matrix.sum(axis=1) > 0) *
                              non_null_positions).sum(axis=1)
            fraction_phenotype = np.nan_to_num(shortest_distance / shared_overlap)

            print(sample, *shortest_distance, sep="\t", file=closest_out)
            print(sample, *shared_overlap, sep="\t", file=overlap_out)
            print(sample, *fraction_phenotype, sep="\t", file=fraction_out)


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
        msa_filename = os.path.join(args['output_dir'], args['clade'] + '.msa.fa')
        with open(msa_filename, 'w') as out:
            AlignIO.write(seqs_msa, out, 'fasta')
