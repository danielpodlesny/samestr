from os.path import isfile
from os import remove
import logging

from glob import glob
from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def bam2freq(arg):
    """
        converts mapping files
        from binary sequence alignment (bam) format
        to SNV profiles (numpy arrays).
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=arg['tmp_dir'])
    if not arg:
        LOG.error('Empty input')
        return False

    if not isfile(arg['kp']):
        LOG.debug('Piling: %s' % arg['bam'])
        oosp.ex(
            'kpileup.py',
            args=[
                arg['bname'],
                arg['bam'],
                arg['gene_file'],
                str(arg['min_base_qual']),
                str(arg['min_aln_qual']),
                str(arg['min_vcov'])
            ],
            out_fn=arg['kp'],
            verbose=False)

    if not glob('%s/*.np.cPickle' % arg['np']) or not glob(
            '%s/*.npy' % arg['np']):
        LOG.debug('Converting: %s' % arg['kp'])
        oosp.ex('kp2np.py',
                args=[
                    '--kp', arg['kp'], '--map', arg['contig_map'], '--sample',
                    arg['bname'], '--gene-file', arg['gene_file'],
                    '--output-dir', arg['np']
                ],
                verbose=False)
    else:
        LOG.warning('Skipping: %s. '
                    'Directory already contained np files: %s' %
                    (arg['bname'], arg['np']))
        
    # clean up
    if not arg['keep_tmp_files']:
        remove(arg['gene_file'])
        remove(arg['contig_map'])
        remove(arg['kp'])
        remove(arg['bam'])
        remove(arg['bam'] + '.bai')

    LOG.debug('Finished: %s' % arg['np'])
    return arg
