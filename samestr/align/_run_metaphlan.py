from os.path import isfile
import logging

from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def run_metaphlan(arg):
    """MetaPhlAn2 function for ooSubprocess."""

    oosp = ooSubprocess.ooSubprocess()
    output_exists = all([
        isfile(arg['mp_profile']),
        isfile(arg['sam']),
        isfile(arg['bowtie2out'])
    ])

    # TODO(dp): if not mp_profile or sam but bowtie2out, rerun fast

    if not output_exists:

        if not isfile(arg['fastq_clean']):
            fastq_clean_zipped = '%s.gz' % arg['fastq_clean']
            if isfile(fastq_clean_zipped):
                arg['fastq_clean'] = fastq_clean_zipped

        LOG.debug('Running: %s' % arg['fastq_clean'])
        oosp.ex('/opt/metaphlan2/metaphlan2.py',
                args=[
                    arg['fastq_clean'], '--bowtie2db', arg['mpa'], '--mpa_pkl',
                    arg['mpa_pkl'], '--input_type', 'fastq', '--nproc',
                    str(arg['nprocs']), '--bowtie2out', arg['bowtie2out'],
                    '-o', arg['mp_profile'], '-s', arg['sam']
                ],
                verbose=False)

    LOG.debug('Finished: %s, %s' % (arg['mp_profile'], arg['sam']))
    return arg
