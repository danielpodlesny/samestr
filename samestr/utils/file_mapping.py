from os.path import isfile, basename, join, dirname, abspath, exists
from samestr.utils import ooSubprocess
from samestr.utils.utilities import list_str_all_endswith

import logging
LOG = logging.getLogger(__name__)


def set_output_structure(args):
    """Sets predefined output filenames.

    Additionally creates predefined output dirs if they don't exist.
    """

    # dest
    out_dir = abspath(args[0]['output_dir']) + '/'
    # in_dir = abspath(args[0]['input_dir']) + '/'
    cmd = args[0]['command']

    # dir names
    kneaddata = 'kneaddata/'
    fastqstats = 'fastq_stats/'
    metaphlan = 'metaphlan/'
    bam = 'bam/'
    freq = 'samestr/'

    ooSubprocess.mkdir(out_dir)

    # dir and file paths
    if cmd == 'align':
        kneaddata_dir, \
        fastqstats_dir, \
        metaphlan_dir, \
        bam_dir, \
        freq_dir = \
            out_dir + kneaddata, \
            out_dir + fastqstats, \
            out_dir + metaphlan, \
            out_dir + bam, \
            out_dir + freq

        # make base dirs
        _ = [
            ooSubprocess.mkdir(d)
            for d in [kneaddata_dir, fastqstats_dir, metaphlan_dir]
        ]

        for arg in args:
            n = arg['bname']

            # kneaddata
            arg['fastq_clean'] = kneaddata_dir + n + '.fastq'
            arg['qc_stats'] = kneaddata_dir + n + '.qc_stats.txt'
            arg['qc_log'] = kneaddata_dir + n + '.log'

            # fastq-stats
            arg['fastq_summary'] = fastqstats_dir + n + '.fastq_summary.txt'

            # metaphlan
            arg['sam'] = metaphlan_dir + n + '.sam.bz2'
            arg['bowtie2out'] = metaphlan_dir + n + '.bowtie2out'
            arg['mp_profile'] = metaphlan_dir + n + '.profile.txt'

    elif cmd == 'convert':

        bam_dir, \
        freq_dir = \
            out_dir + bam, \
            out_dir + freq

        # make base dirs
        _ = [ooSubprocess.mkdir(d) for d in [bam_dir, freq_dir]]

        for arg in args:
            n = arg['bname']

            # metaphlan
            arg['sam'] = arg['input_dir'] + n + '.sam.bz2'
            arg['bowtie2out'] = arg['input_dir'] + n + '.bowtie2out'

            ## metaphlan profiles
            if not arg['mp_profiles_dir']:
                arg['mp_profiles_dir'] = arg['input_dir']
            else:
                arg['mp_profiles_dir'] = abspath(arg['mp_profiles_dir']) + '/'

            ### file name
            arg['mp_profile'] = arg['mp_profiles_dir'] + n + arg[
                'mp_profiles_extension']

            ### exists
            if not exists(arg['mp_profile']):
                LOG.error('MetaPhlAn file not found: %s' % arg['mp_profile'])
                exit(1)

            # sam2bam
            arg['bam'] = bam_dir + n + '.bam'

            # bam2freq
            _freq_sample_dir = freq_dir + n + '/'
            arg['gene_file'] = _freq_sample_dir + n + '.gene_file.txt'
            arg['contig_map'] = _freq_sample_dir + n + '.contig_map.txt'
            arg['kp'] = _freq_sample_dir + n + '.kp.txt'
            arg['np'] = _freq_sample_dir

            # make sample dirs
            ooSubprocess.mkdir(_freq_sample_dir)

    return args


def spread_args_by_input_files(args):
    # spread args to list of args per sample/sample-pair
    spread_args = []
    for idx, (base_name,
              input_files) in enumerate(args['input_files'].items()):

        spread_args.append({})
        for arg in args:
            if not arg == 'input_files':
                spread_args[idx][arg] = args[arg]

        group_size = len(input_files)
        spread_args[idx]['bname'] = base_name
        spread_args[idx]['input_dir'] = dirname(abspath(input_files[0])) + '/'

        # for paired-end: sanity check for pair counts of two
        if args['input_sequence_type'] == 'paired':
            if group_size > 2:
                LOG.error('More than two samples (%s) found: %s [%s]' %
                          (group_size, base_name, ','.join(input_files)))
                exit(1)
            elif group_size < 2:
                LOG.error('Not all samples are paired: %s.'
                          'Check file extension after read index [%s]' %
                          (base_name, ','.join(input_files)))
                exit(1)
            else:
                spread_args[idx]['1%s' % args['input_extension']] = [
                    f for f in input_files if '1.fastq' in f
                ][0]
                spread_args[idx]['2%s' % args['input_extension']] = [
                    f for f in input_files if '2.fastq' in f
                ][0]

        # for single-end: sanity check 1 observation per base name
        else:
            if group_size > 1:
                LOG.error('More than one sample (%s) found: %s [%s]' %
                          (group_size, base_name, ','.join(input_files)))
                exit(1)
            else:
                spread_args[idx]['%s' %
                                 args['input_extension']] = input_files[0]

    return spread_args


def get_accepted_extension(file, accepted_extensions):
    e = '.' + '.'.join(file.rsplit('.', 1)[-1:])
    if e in accepted_extensions:
        return e
    else:
        e = '.' + '.'.join(file.rsplit('.', 2)[-2:])
        if e in accepted_extensions:
            return e
    return False


def get_uniform_extension(files, accepted_extensions):
    input_extension = get_accepted_extension(files[0], accepted_extensions)

    if not input_extension:
        LOG.error('Files must be supplied with accepted file extensions: %s' %
                  ', '.join(accepted_extensions))
        exit(1)

    if not list_str_all_endswith(files, input_extension):
        LOG.error('Not all files have the same file extension.')
        exit(1)

    return input_extension
