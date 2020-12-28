from os import remove
from os.path import isfile, dirname
from tempfile import mkdtemp
from shutil import copy, rmtree

import logging

from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def run_kneaddata(arg):
    """kneaddata function for ooSubprocess."""

    oosp = ooSubprocess.ooSubprocess()
    output_exists = isfile(
        arg['fastq_clean']) or isfile(arg['fastq_clean'] + '.gz')

    if not output_exists:
        LOG.debug('Running: %s,%s' % (arg['1%s' % arg['input_extension']],
                                      arg['2%s' % arg['input_extension']]))

        cmd = []
        if arg['input_sequence_type'] == 'single-end':
            cmd += ['--input', arg[arg['input_extension']]]
        else:
            cmd += [
                '--input', arg['1%s' % arg['input_extension']], '--input',
                arg['2%s' % arg['input_extension']]
            ]
        cmd = cmd + \
            [
                '--cat-final-output',
                '-o', dirname(arg['fastq_clean']),
                '--output-prefix', arg['bname'],
                '-db', arg['host_bowtie2db'],
                '-t', str(arg['nprocs']),
                '--max-memory', arg['max_ram']
            ]
        oosp.ex(arg['kneaddata_exe'], args=cmd, verbose=False)

    if not arg['keep_intermediate_fastq']:
        LOG.debug('Removing temporary files')
        # fix for tmp file removal warning
        # v0.7.2 got stuck @trimmomatic, downgraded, manual intermediate removal
        rm_names = [
            '%s.trimmed.%s.fastq', '%s.trimmed.single.%s.fastq',
            '%s_unmatched_%s.fastq', '%s_paired_%s.fastq',
            '%s_paired_Homo_sapiens_bowtie2_contam_%s.fastq',
            '%s_Homo_sapiens_bowtie2_paired_contam_%s.fastq',
            '%s_Homo_sapiens_bowtie2_unmatched_%s_contam.fastq',
            '%s_unmatched_%s_Homo_sapiens_bowtie2_contam.fastq'
        ]

        for file in rm_names:
            for read in [1, 2]:
                rm_file = dirname(
                    arg['fastq_clean']) + '/' + file % (arg['bname'], read)

                if isfile(rm_file):
                    try:
                        remove(rm_file)
                        LOG.debug('Deleted: %s' % rm_file)
                    except OSError:
                        LOG.debug('File not removed: %s' % rm_file)

    if not isfile(arg['qc_stats']):
        # kneaddata needs dir as input to create qc tables,
        # create tmp dir to feed single files
        LOG.info('Writing kneaddata statistics')
        tmp_dir = mkdtemp(dir=arg['output_dir'])
        copy(arg['qc_log'], tmp_dir)

        oosp.ex('%s/kneaddata_read_count_table' % dirname(arg['kneaddata_exe']),
                args=['--input', tmp_dir, '--output', arg['qc_stats']],
                verbose=False)

        rmtree(tmp_dir)

    LOG.info('Finished: %s' % arg['fastq_clean'])
    return arg
