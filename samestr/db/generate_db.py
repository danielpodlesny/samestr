from os.path import join
import bz2
import pickle as pickle
from Bio import SeqIO
from os.path import isdir
from os import makedirs

import logging

from samestr.utils import clade_path

LOG = logging.getLogger(__name__)    

def generate_db(input_args):
    """
    Function expands and generates a database of clade markers 
    from MetaPhlAn or mOTUs database that is required for SameStr.
    """

    # set markers
    all_clades = set()
    all_markers = {}
    n_markers = 0

    # get markers fasta
    if input_args['markers_fasta'].endswith('.bz2'):
        markers_fasta = SeqIO.to_dict(
            SeqIO.parse(bz2.open(input_args['markers_fasta'], 'rt'), 'fasta'))
    else:
        markers_fasta = SeqIO.to_dict(
            SeqIO.parse(input_args['markers_fasta'], 'fasta'))

    # get markers mapping   
    if input_args['db_source'] == 'MetaPhlAn':
        markers_pickle = pickle.load(bz2.BZ2File(input_args['markers_pkl']))
    elif input_args['db_source'] == 'mOTUs':
        markers_pickle = {"markers": {}}

        # Populate marker_pickle based on db_mOTU_DB_CEN.fasta
        for marker in markers_fasta:
            clade, _ = marker.split('.', 1)
            if marker not in markers_pickle['markers']:
                if clade != 'NA':
                    markers_pickle['markers'][marker] = {'clade': clade}

    for marker in markers_pickle['markers']:

        # retrieve clade names
        clade = markers_pickle['markers'][marker]['clade']
        all_clades.add(clade)
        n_markers += 1

        # get only selected clades
        if 'clade' in input_args and input_args['clade'] is not None:
            if clade not in input_args['clade']:
                continue

        if clade not in all_markers:
            all_markers[clade] = set()
        all_markers[clade].add(marker)

    LOG.debug('Database contains %s clades, %s markers' %
              (len(all_clades), n_markers))

    # report selected clades
    if 'clade' in input_args and input_args['clade'] is not None:
        clade_intersect = set(
            input_args['clade']).intersection(all_clades)
        clade_diff = set(
            input_args['clade']).difference(clade_intersect)

        LOG.debug('Selected clades found in the database: %s/%s' %
                  (len(clade_intersect), len(input_args['clade'])))
        if len(clade_intersect) < len(input_args['clade']):
            LOG.debug('Clades not found in the database: %s' %
                      ', '.join(clade_diff))
        all_clades = clade_intersect

    for clade in all_clades:

        # set output directory and names
        output_dir = join(input_args['output_dir'], clade_path(clade))

        # Create dir if not exists
        if not isdir(output_dir):
            makedirs(output_dir)

        output_base = output_dir + clade
        marker_filename = output_base + '.markers.fa'
        pos_filename = output_base + '.positions.txt'
        gene_filename = output_base + '.gene_file.txt'
        cmap_filename = output_base + '.contig_map.txt'

        LOG.debug('Output directory: %s, %s' % (output_base, marker_filename))

        # check marker count
        clade_markers = sorted(all_markers[clade])
        n_markers = len(clade_markers)
        LOG.debug('Markers found for %s: %s' % (clade, n_markers))
        if not n_markers:
            continue

        total_marker_len = 0
        gene_data = []
        pos_data = ['%s\t%s\t%s' % ('contig', 'pos', 'len')]
        remove_markers = set()

        with open(marker_filename, 'w') as marker_file:

            for marker in clade_markers:
                if marker in markers_fasta:
                    seq = markers_fasta[marker]
                    marker_len = len(seq)

                    contig_id = marker
                    gene_id = marker
                    strand = '+'
                    beg = 1
                    end = marker_len

                    pos_data += [
                        '%s\t%s\t%s' %
                        (marker, total_marker_len, marker_len)
                    ]
                    gene_data += ['%s\t%s\t%s\t%s\t%s\t%s' %
                                  (contig_id, gene_id, beg, end, strand, seq.seq)]
                    SeqIO.write(seq, marker_file, 'fasta')

                    total_marker_len += marker_len
                else:
                    remove_markers.add(marker)

        if remove_markers:
            LOG.debug('Dropping markers for %s: %s' %
                      (clade, len(remove_markers)))
            clade_markers = [
                m for m in clade_markers if m not in remove_markers]

        with open(gene_filename, 'w') as gene_file, \
                open(cmap_filename, 'w') as cmap_file, \
                open(pos_filename, 'w') as pos_file:

            pos_file.write('\n'.join(pos_data) + '\n')
            gene_file.write('\n'.join(gene_data) + '\n')
            cmap_file.write('\n'.join(
                ['\t'.join([clade, m]) for m in clade_markers]) + '\n')
    return input_args
