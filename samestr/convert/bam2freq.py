import logging
import os
import subprocess

from glob import glob
from os import remove

from samestr.convert.convert import kp2np, pileup


LOG = logging.getLogger(__name__)


def results_incomplete(outdir):
    return not glob('%s/*.np.cPickle' % outdir) or \
       not glob('%s/*.npy' % outdir) or \
       not glob('%s/*.npz' % outdir) or \
       not glob('%s/*.npy.gz' % outdir)


def bam2freq(arg):

    if not arg:
        LOG.error('Empty input')
        return None

    if not os.path.isfile(arg["kp"]):
        LOG.debug('Piling: %s' % arg['sorted_bam'])
        cmd = ["samtools", "view", arg['sorted_bam']]
        with open(arg["kp"], "wt") as kpileup_out, subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True) as bam_process:
            position_stream = pileup(
                # arg['bname'],
                # arg['sorted_bam'],
                bam_process.stdout,
                arg['gene_file'],
                arg['min_base_qual'],
                arg['min_aln_qual'],
                arg['min_vcov'],
                outstream=kpileup_out,
            )
            if results_incomplete(arg['np']):            
                LOG.debug('Converting: %s' % arg['kp'])
                kp2np(
                    position_stream,
                    arg['contig_map'],
                    arg['bname'],
                    arg['gene_file'],
                    arg['np'],
                )
            else:
                LOG.warning(
                    'Skipping: %s. '
                    'Directory already contained np files: %s' %
                    (arg['bname'], arg['np'])
                )

            if not arg['keep_tmp_files']:
                remove(arg['gene_file'])
                remove(arg['contig_map'])
                remove(arg['kp'])
                remove(arg['sorted_bam'])
                remove(arg['sorted_bam'] + '.bai')

            LOG.debug('Finished: %s' % arg['np'])
        
    return arg
