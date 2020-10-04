from os.path import isdir, basename, isfile, join, exists, dirname
import logging

from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def summarize(args):
    """
        summarizes pairwise alignment comparison data
    """
    oosp = ooSubprocess.ooSubprocess()

    LOG.debug('Summarizing: %s' % (args['input_dir']))
    oosp.ex('Rscript',
            args=[
                '--vanilla',
                '%s/summarize.R' % dirname(__file__),
                '--input-dir',
                '%s' % args['input_dir'],
                '--mp-profiles-dir',
                '%s' % args['mp_profiles_dir'],
                '--mp-profiles-extension',
                '%s' % args['mp_profiles_extension'],
                '--output-dir',
                '%s/' % args['output_dir'],
                '--aln-pair-min-overlap',
                '%s' % args['aln_pair_min_overlap'],
                '--aln-pair-min-similarity',
                '%s' % args['aln_pair_min_similarity'],
            ],
            verbose=False)

    LOG.debug('Finished: %s' % args['input_dir'])
    return args
