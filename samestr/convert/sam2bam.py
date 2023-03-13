from os.path import isfile

import logging

from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def sam2bam(arg):
    """
        formats mapping files
        from sequence alignment format (sam)
        by filtering, converting, sorting, indexing
        to binary sequence alignment format (bam).
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=arg['tmp_dir'])
    output_exists = isfile(arg['bam'])

    if not output_exists:
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
                   out_fn=arg['bam'],
                   stop=True)

        # index
        oosp.ex('samtools', args=['index', arg['bam']], verbose=False)

    LOG.debug('Finished: %s' % arg['bam'])
    return arg
