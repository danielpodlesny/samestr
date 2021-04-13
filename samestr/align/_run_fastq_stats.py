from os.path import isfile
import logging

from samestr.utils import ooSubprocess

LOG = logging.getLogger(__name__)


def parse_fastq_stats(fastq_stats_report):
    fastq_stats = {}

    save_metrics = [
        'reads', 'len', 'len_mean', 'len_stdev', 'qual_mean', 'qual_stdev',
        'total_bases'
    ]
    for line in fastq_stats_report.splitlines():
        line = line.split('\t')
        metric = line[0].strip().replace(' ', '_')
        if metric.startswith('No_reads_in_'):
            LOG.error('Not enough reads')
            exit(0)
        elif metric in save_metrics:
            fastq_stats[metric] = line[1].strip()
    values = [fastq_stats[m] for m in save_metrics]
    fastq_stats_text = '\t'.join(save_metrics) + '\n' + \
                        '\t'.join(values) + '\n'
    return fastq_stats, fastq_stats_text


def run_fastq_stats(arg):
    """fastq-stats function for ooSubprocess."""
    LOG.debug('Started fastq-stats')

    oosp = ooSubprocess.ooSubprocess()
    output_exists = all([isfile(arg['fastq_summary'])])
    if not output_exists:

        LOG.debug('Running: %s' % arg['fastq_clean'])

        if not isfile(arg['fastq_clean']):
            fastq_clean_zipped = '%s.gz' % arg['fastq_clean']
            if isfile(fastq_clean_zipped):
                arg['fastq_clean'] = fastq_clean_zipped

        fastq_stats_report = oosp.ex(arg['fastq_stats_exe'],
                                     args=[arg['fastq_clean']],
                                     get_output=True,
                                     verbose=False)

        _, fastq_stats_text = parse_fastq_stats(fastq_stats_report)
        with open(arg['fastq_summary'], 'wb') as file:
            file.write(fastq_stats_text)

    LOG.debug('Finished: %s' % (arg['fastq_summary']))
    return arg
