#!/usr/bin/env python3
"""
    ref2freq generates multiple sequence alignments from MetaPhlAn
    marker sequences and supplied reference genomes [e.g. RefSeq, NCBI-Genome]
    which are further reduced to all positions in the initial
    MetaPhlAn marker positions and concatenated.

    Adapted StrainPhlAn's code [Truong et al.; Segata Lab]
    for adding references to MSA of sample consensus sequences
    to output to compatible numpy format
"""

from os import remove, makedirs
from os.path import join, isdir
from sys import exit

import random
from shutil import copyfile
from glob import glob
from re import sub

import bz2
import gzip
from collections import defaultdict
from tempfile import SpooledTemporaryFile, NamedTemporaryFile

import logging
import logging.config

import numpy as np

from Bio import SeqIO, Seq, SeqRecord, AlignIO
from Bio.Align import MultipleSeqAlignment

from samestr.utils import ooSubprocess
from samestr.utils.ooSubprocess import trace_unhandled_exceptions
from samestr.utils.utilities import all_exe

LOG = logging.getLogger(__name__)


def truncate_markers(args, reference_markers):
    """
        For each reference and marker in `reference_markers` dict,
        truncates marker by `marker_trunc_len` bases at beginning and end of the sequence
        and returns dict with truncated MSA
    """

    for reference in reference_markers:
        initial_markers_total_length = 0
        truncated_markers_total_length = 0
        remove_markers = []

        for marker in reference_markers[reference]:
            initial_markers_total_length += len(
                reference_markers[reference][marker]['seq'])

            # if marker_trunc_len, truncate
            if args['marker_trunc_len'] > 0:
                reference_markers[reference][marker]['seq'] = \
                    reference_markers[reference][marker]['seq'][args['marker_trunc_len']
                        :-args['marker_trunc_len']]
            truncated_markers_total_length += len(
                reference_markers[reference][marker]['seq'])

            # if empty, drop
            if not reference_markers[reference][marker]['seq']:
                remove_markers.append(marker)

        for marker in remove_markers:
            LOG.debug('Removing marker: %s' % marker)
            del reference_markers[reference][marker]

    return reference_markers


def align(args, marker_fn, aln_file):
    """
        Aligning `marker_fn` with `alignment program`
        returning alignment.
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=args['tmp_dir'])

    # TODO: implement muscle 5 (double free corruption error if just one sequence)
    # return marker if marker has just one sequence
    # if len(list(SeqIO.parse(marker_fn, 'fasta'))) == 1:
    #     copyfile(marker_fn, aln_file)
    #     return

    if args['aln_program'] == 'muscle':
        alignment = oosp.ex('muscle',
                            # muscle 5 format, currently buggy
                            # args=['-quiet', '-align', marker_fn, '-output', aln_file],
                            # muscle v3.8.1551
                            args=['-quiet', '-in', marker_fn, '-out', aln_file],
                            get_output=True,  # wait for process to finish
                            verbose=False)

    elif args['aln_program'] == 'mafft':
        alignment = oosp.ex('mafft',
                            args=['--auto', '--quiet', marker_fn],
                            out_fn=aln_file,
                            verbose=False)
    return


def align_clean(args):
    """
        Creating clean alignment for marker of ooSubprocess instance using :align:
        returning multiple sequence alignment as dicts of numpy arrays:
            seq, DNA sequences: ['A', 'T', ..]
    """
    tmp_marker_file = NamedTemporaryFile(dir=args['tmp_dir'], mode='a+t', delete=False)
    tmp_marker_fn = tmp_marker_file.name

    # writing blast result of marker for each reference from dict to fasta file
    ref_count = 0
    for reference in list(args['reference_markers'].keys()):
        if args['species_marker'] in args['reference_markers'][reference]:
            ref_count += 1
            SeqIO.write(
                SeqRecord.SeqRecord(
                    id=reference,
                    description='',
                    seq=Seq.Seq(
                        args['reference_markers'][reference][args['species_marker']]['seq'])),
                tmp_marker_file, 'fasta')
    tmp_marker_file.close()

    # aligning sequences in fasta file
    tmp_alignment_file = NamedTemporaryFile(
        dir=args['tmp_dir'], mode='a+t', delete=not args['keep_tmp_files'])

    align(args, tmp_marker_fn, tmp_alignment_file.name)

    # adding alignment to out dicts
    seq = {}
    for rec in SeqIO.parse(tmp_alignment_file, 'fasta'):
        reference = rec.name
        seq[reference] = list(str(rec.seq))

    if args['alignment_fn']:
        copyfile(tmp_alignment_file.name, args['alignment_fn'])
    tmp_alignment_file.close()

    # delete temporary directory
    if not args['keep_tmp_files']:
        remove(tmp_marker_fn)

    if not seq:
        LOG.error('Fatal error in alignment step: %s' % species_marker)
        exit(1)

    return seq


def run_align_clean(args, reference_markers, species_markers_ordered):
    """
        Run :align_clean: on each marker in parallel with ooSubprocess
    """

    # parallelize align_clean
    args_list = []
    for i, _ in enumerate(species_markers_ordered):
        args_list.append({})
        args_list[i]['species_marker'] = species_markers_ordered[i]
        args_list[i]['species_marker_name'] = args['species_marker_name']
        args_list[i]['tmp_dir'] = args['tmp_dir']
        args_list[i]['output_dir'] = args['output_dir']
        args_list[i]['aln_program'] = args['aln_program']
        args_list[i]['reference_markers'] = reference_markers
        args_list[i]['keep_tmp_files'] = args['keep_tmp_files']
        if args['save_marker_aln']:
            formatted_marker_id = sub(r'[^a-zA-Z0-9 \n\.]', '_',
                                      species_markers_ordered[i])
            alignment_fn = join(args['output_dir'],
                                formatted_marker_id + '.marker_aligned.fasta')
            if args['marker_trunc_len'] > 0:
                alignment_fn = alignment_fn + '.marker_trunc_len_' + str(
                    args['marker_trunc_len'])
            args_list[i]['alignment_fn'] = alignment_fn
        else:
            args_list[i]['alignment_fn'] = False

    LOG.debug('Running align_clean for all markers')
    seqs = ooSubprocess.parallelize(align_clean,
                                    args_list,
                                    args['nprocs'],
                                    use_threads=False)

    if not seqs:
        LOG.error('No markers found!')
        exit(1)

    cat_seqs = defaultdict(list)

    # reformat dict to concatenate markers from seq[marker][ref] to seq[ref][marker]
    # remove positions where gaps were introduced into the species_marker
    # save positions where markers are concatenated, as well as marker lengths
    pos = 0
    marker_pos = []

    # seqs are ordered tuples
    for marker, _ in enumerate(seqs):
        if seqs[marker]:
            references_with_alignment = list(seqs[marker].keys())

            # index gaps in species_marker
            gap_positions = [
                idx for idx, n in enumerate(seqs[marker][
                    args['species_marker_name']]) if n == '-'
            ]
            marker_length = len(
                seqs[marker][args['species_marker_name']]) - len(gap_positions)

            if len(gap_positions):
                LOG.debug(
                    'Removing %s gap positions introduced during alignment: %s'
                    % (len(gap_positions), species_markers_ordered[marker]))

            # reduce marker to original mp_maker pos & length
            # add alignment seqs to reference dictionary
            for reference in list(reference_markers.keys()):
                if reference in references_with_alignment:
                    cat_seqs[reference] += [
                        n for idx, n in enumerate(seqs[marker][reference])
                        if idx not in gap_positions
                    ]

                # add empty arrays for references without marker in alignment
                else:
                    dummy_seqs = ['-' for p in range(marker_length)]
                    cat_seqs[reference] += dummy_seqs

            marker_pos.append(
                [species_markers_ordered[marker], pos, marker_length])
            pos += marker_length
        else:
            LOG.error('species_marker was not aligned: %s' %
                      species_markers_ordered[marker])
            exit(1)

    LOG.debug('Finished align_clean')
    return cat_seqs, marker_pos


def blast_markers_against_references(args, species_markers):
    """
        Using extract to find `species_markers` in BLASTdb generated from `input_files`
        returning search results in `reference_markers` as dict of markers:
            reference_markers[reference][marker_query]['seq']
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=args['tmp_dir'])

    unique_markers = set(species_markers.keys())
    references = sorted(list(set(args['input_files'])))
    reference_markers = defaultdict(dict)
    contigs = defaultdict(dict)

    # read references from files
    p1 = SpooledTemporaryFile(mode='a+t', dir=args['tmp_dir'])
    for ref_fn in references:

        genome = ooSubprocess.splitext(ref_fn)[0]
        LOG.debug('Adding contigs for %s' % genome)
        if ref_fn.endswith('.bz2'):
            ref = bz2.open(ref_fn, 'rt')
        elif ref_fn.endswith('.gz'):
            ref = gzip.open(ref_fn, 'rt')
        elif ref_fn.endswith('.fna') or \
                ref_fn.endswith('.fa') or \
                ref_fn.endswith('.fasta'):
            ref = open(ref_fn, 'r')
        else:
            LOG.error('Unknown file type: %s '
                      'Fasta file should be raw [.fa, .fna, .fasta] '
                      'or zipped with bzip2 [.bz2] or gzip [.gz]' % ref_fn)
            exit(1)

        # check duplicate contigs
        unique_contigs_only = True
        for rec in SeqIO.parse(ref, 'fasta'):
            if rec.name in list(contigs.keys()):
                LOG.error('Contig %s is not unique' %
                          (rec.name.split('___')[-1]))
                unique_contigs_only = False
                break

        # only add references with before unseen contigs
        if not unique_contigs_only:
            LOG.warning('Skipping %s due to redundant contigs.' % genome)
            continue

        # load contigs from reference to dict, cat to tmp file
        else:
            LOG.debug('Contigs found for %s are unique' % genome)
            ref.seek(0)
            for rec in SeqIO.parse(ref, 'fasta'):
                contigs[rec.name]['seq'] = str(rec.seq)
                contigs[rec.name]['genome'] = genome
                SeqIO.write(rec, p1, 'fasta')

        ref.close()

    if not contigs:
        LOG.error('Failed to aggregate contigs.')
        exit(0)

    blastdb_prefix = oosp.ftmp('genome_extract_db_%s' % (random.random()))
    if glob('%s*' % blastdb_prefix):
        LOG.error('BLASTdb exists! Please remove it or rerun!')
        exit(1)

    # build blastdb from contigs of references in tmp file
    p1.seek(0)
    oosp.ex('makeblastdb',
            args=[
                '-dbtype', 'nucl', '-title', 'genome_db', '-out',
                blastdb_prefix
            ],
            in_pipe=p1,
            verbose=True)

    p1 = SpooledTemporaryFile(mode='a+t', dir=args['tmp_dir'])

    for unique_marker in unique_markers:
        SeqIO.write(species_markers[unique_marker], p1, 'fasta')

    # blast markers against reference contigs
    extract_args = [
        '-db', blastdb_prefix, '-outfmt', '6', '-evalue', '1e-10',
        '-word_size', '10', '-max_target_seqs', '1000000000'
    ]
    if args['nprocs'] > 1:
        extract_args += ['-num_threads', str(args['nprocs'])]

    p1.seek(0)
    extract_output = oosp.ex('blastn',
                             args=extract_args,
                             in_pipe=p1,
                             get_out_pipe=True,
                             verbose=True)

    # process extract output
    for line in extract_output:
        if line.strip() == '':
            break
        line = line.decode().strip().split()
        marker_query = line[0]
        target_contig = line[1]
        pstart = int(line[8]) - 1
        pend = int(line[9]) - 1

        reference = contigs[target_contig]['genome']

        # extract extract positions from reference contigs
        if marker_query not in reference_markers[reference]:
            reference_markers[reference][marker_query] = {}
            if pstart < pend:
                reference_markers[reference][marker_query]['seq'] = contigs[
                    target_contig]['seq'][pstart:pend + 1].upper()
            else:
                reference_markers[reference][marker_query]['seq'] = \
                    str(Seq.Seq(contigs[target_contig]['seq']
                        [pend:pstart+1]).reverse_complement()).upper()

    # remove database
    for fn in glob('%s*' % blastdb_prefix):
        remove(fn)

    LOG.debug('Number of reference genomes added: %d' % len(reference_markers))

    return reference_markers


def seqs2freqs(args, seqs, marker_pos):
    """
        Converts seqs to freqs and saves as numpy ndarrays,
        Reduces the alignment to positions only covered by the species_marker

        Info: Frequencies will always be 1 for references, since generated
        from BLAST and depth == 1. For mapped reads from samples,
        frequencies will be equivalent to the count of nucleotide observations
        at each position of the alignment e.g. [[14, 1, 0, 0], [0, 16, 0, 0],
        hence the sum of each vector equals coverage at given position.
    """

    # reduce to species_marker positions, duplicate, removing @align_clean to update positions
    species_marker_pos = np.array([
        pos for pos, n in enumerate(seqs[args['species_marker_name']])
        if n != '-'
    ])
    for reference in list(seqs.keys()):
        seqs[reference] = np.array(seqs[reference])[species_marker_pos]
    len_after_removal = len(seqs[args['species_marker_name']])

    species_marker_length_delta = len_after_removal - args[
        'species_markers_total_length']

    if species_marker_length_delta != 0:
        LOG.error('Difference to initial species_marker length: %s' %
                  str(species_marker_length_delta))
        exit(1)

    # nucleotide conversion arrays
    n_freq = np.array(
        [
            np.array([0, 0, 0, 0]),  # -
            np.array([0, 0, 0, 0]),  # N
            np.array([1, 0, 0, 0]),  # A
            np.array([0, 1, 0, 0]),  # C
            np.array([0, 0, 1, 0]),  # G
            np.array([0, 0, 0, 1])  # T
        ],
        dtype=float)
    acgt = '-NACGT'

    # convert seqs to freqs, set dimensions to [m=sample, n=width, k=4]
    freqs = defaultdict(list)
    for reference in list(seqs.keys()):
        freqs[reference] = np.array([
            n_freq[acgt.index(n)] if n in acgt else n_freq[acgt.index('-')]
            for n in seqs[reference]
        ])
        freqs[reference] = np.expand_dims(freqs[reference], axis=0)
        seqs[reference] = np.expand_dims(seqs[reference], axis=0)

    # makedirs
    if not isdir(args['output_dir']):
        makedirs(args['output_dir'])

    def _format_reference_name(name):
        # fasta, fna, fa, + gz
        if name.endswith(('.fa.gz', '.fasta.gz', '.fna.gz', '.fa.bz2',
                          '.fasta.bz2', '.fna.bz2')):
            return '.'.join(name.split('.')[:-2])
        elif name.endswith(('.fa', '.fna', '.fasta')):
            return '.'.join(name.split('.')[:-1])
        return name

    # save freqs array to file .npy
    output_base = args['output_dir'] + '/' + args['species']
    for reference in list(freqs.keys()):
        if reference.endswith('.markers'):
            continue

        fname = output_base + '.%s.npy' % _format_reference_name(reference)
        np.save(fname, freqs[reference], allow_pickle=True)

    # convert seqs to Bio SeqIO alignment object
    seqs_list = []
    name_order = []
    for reference in list(seqs.keys()):
        name_order.append(_format_reference_name(reference))
        seqs_list.append(
            SeqRecord.SeqRecord(id=_format_reference_name(reference),
                                description='',
                                seq=Seq.Seq(''.join(seqs[reference][0]))))
    seqs_msa = MultipleSeqAlignment(seqs_list)

    # write alignment fasta
    LOG.info('Saving alignment.')
    with open(output_base + '.msa.fa', 'w') as out:
        AlignIO.write(seqs_msa, out, 'fasta')

    return True


def check_dependencies(args):
    return all_exe([args['aln_program'], 'extract', 'makeblastdb'])


@trace_unhandled_exceptions
def ref2freq(args):
    # check_dependencies(args)

    # set species, tmp_dir, ofn, check output exists
    args['species_marker_name'] = args['species'] + '.markers'
    args['species_marker'] = args['marker_dir'] + '/' + args[
        'species_marker_name'] + '.fa'

    if args['tmp_dir'] == '':
        args['tmp_dir'] = args['output_dir']

    # load contigs from metaphlan marker
    species_markers = {}
    species_markers_ordered = []
    args['species_markers_total_length'] = 0
    for rec in SeqIO.parse(open(args['species_marker'], 'r'), 'fasta'):
        species_markers[rec.id] = rec
        species_markers_ordered.append(rec.id)
        args['species_markers_total_length'] += len(rec)

    # report basic info
    txt = 'Species: %s\nMarker File: %s\nMarkers: %s\nTotal Nucleotides: %s\nReferences: %s'
    LOG.info(
        txt %
        (args['species'], args['species_marker'], len(species_markers_ordered),
         args['species_markers_total_length'], ', '.join(args['input_files'])))

    # blast metaphlan markers against provided references
    reference_markers = blast_markers_against_references(args, species_markers)

    # add species_markers to extract output dict
    for species_marker in list(species_markers.keys()):
        reference_markers[args['species_marker_name']][species_marker] = {}
        reference_markers[
            args['species_marker_name']][species_marker]['seq'] = str(
                species_markers[species_marker].seq)

    # truncate marker ends
    if args['marker_trunc_len'] > 0:
        reference_markers = truncate_markers(args, reference_markers)

    # align markers, concatenate in order of species_markers
    seqs, marker_pos = run_align_clean(args, reference_markers,
                                       species_markers_ordered)

    # convert seqs to freqs, save msa, seqs, freqs, marker_pos
    seqs2freqs(args, seqs, marker_pos)
