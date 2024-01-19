from os.path import isfile
from shutil import copyfile

import logging

from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def sam2bam(arg):
    """
        formats mapping files
        from sequence alignment format (sam)
        by filtering, converting, sorting, 
        to binary sequence alignment format (bam), 
        unless the bam already exists.
        Finally, bam files are indexed.
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=arg['tmp_dir'])

    if not isfile(arg['sorted_bam']):

        # metaphlan
        if isfile(arg['sam']) and not isfile(arg['bam']):
            error_pipe = None
            LOG.debug('Converting: %s' % (arg['sam']))
            # decompress
            p1 = oosp.chain('dump_file.py',
                            args=['--input_file', arg['sam']],
                            stderr=error_pipe)

            # filter len
            p2 = oosp.chain(
                'filter_sam.py',
                args=[
                    str(int(arg['min_aln_identity'] * 100)),  # float
                    str(arg['min_aln_len'])
                ],
                in_pipe=p1,
                stderr=error_pipe)

            # sam to bam
            p3 = oosp.chain('samtools',
                            args=['view', '-bS', '-'],
                            in_pipe=p2,
                            stderr=error_pipe)

            # sort
            oosp.chain('samtools',
                    args=['sort', '-', '-o', arg['bname']],
                    in_pipe=p3,
                    stderr=error_pipe,
                    out_fn=arg['sorted_bam'],
                    stop=True)
        
        # motus
        elif not isfile(arg['sam']) and isfile(arg['bam']):
            # copy existing bam to sorted_bam output position
            copyfile(arg['bam'], arg['sorted_bam'])

        else:
            LOG.error('Error parsing input sam/bam files: [%s/%s]' % (arg['sam'], arg['bam']))
            exit(0)

        # index
        if not isfile(arg['sorted_bam'] + '.bai'):
            oosp.ex('samtools', args=['index', arg['sorted_bam']], verbose=False)

    LOG.debug('Finished: %s' % arg['sorted_bam'])
    return arg
