from os.path import join
import bz2
import pickle as pickle
from Bio import SeqIO
from os.path import isdir
from os import makedirs

import logging

# TODO: do we still need this?
# from samestr.db import TaxClade, TaxTree
from samestr.utils import clade_path

# TODO: refactor species/clade

LOG = logging.getLogger(__name__)    


def generate_db(input_args):
    """
    Function expands and generates a database of clade markers 
    from MetaPhlAn or mOTUs database that is required for SameStr.
    """

    # set species markers
    all_clades = set()
    all_markers = {}
    n_markers = 0

    # get markers fasta
    if input_args['mpa_markers'].endswith('.bz2'):
        marker_fasta = SeqIO.to_dict(
            SeqIO.parse(bz2.open(input_args['mpa_markers'], 'rt'), 'fasta'))
    else:
        marker_fasta = SeqIO.to_dict(
            SeqIO.parse(input_args['mpa_markers'], 'fasta'))

    # get markers mapping   
    # TODO fix to real var
    db = 'MetaPhlAn'
    if db == 'MetaPhlAn':
        marker_pickle = pickle.load(bz2.BZ2File(input_args['mpa_pkl']))
    elif db == 'mOTUs':
        marker_pickle = {"markers": {}}

        # Populate marker_pickle based on db_mOTU_DB_CEN.fasta
        for marker in marker_fasta:
            clade, _ = marker.split('.', 1)
            if marker not in marker_pickle['markers']:
                if clade != 'NA':
                    marker_pickle['markers'][marker] = {'clade': clade}

    for marker in marker_pickle['markers']:

        # retrieve clade names
        clade = marker_pickle['markers'][marker]['clade']
        all_clades.add(clade)
        n_markers += 1

        # get only selected clades
        if 'species' in input_args and input_args['species'] is not None:
            if clade not in input_args['species']:
                continue

        if clade not in all_markers:
            all_markers[clade] = set()
        all_markers[clade].add(marker)

    LOG.debug('Database contains %s species, %s markers' %
              (len(all_clades), n_markers))

    # report selected species
    if 'species' in input_args and input_args['species'] is not None:
        species_intersect = set(
            input_args['species']).intersection(all_clades)
        species_diff = set(
            input_args['species']).difference(species_intersect)

        LOG.debug('Selected species found in the database: %s/%s' %
                  (len(species_intersect), len(input_args['species'])))
        if len(species_intersect) < len(input_args['species']):
            LOG.debug('Species not found in the database: %s' %
                      ', '.join(species_diff))
        all_clades = species_intersect

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
        species_markers = sorted(all_markers[clade])
        n_markers = len(species_markers)
        LOG.debug('Markers found for %s: %s' % (clade, n_markers))
        if not n_markers:
            continue

        total_marker_len = 0
        gene_data = []
        pos_data = ['%s\t%s\t%s' % ('contig', 'pos', 'len')]
        remove_markers = set()

        with open(marker_filename, 'w') as marker_file:

            for marker in species_markers:
                if marker in marker_fasta:
                    seq = marker_fasta[marker]
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
            species_markers = [
                m for m in species_markers if m not in remove_markers]

        with open(gene_filename, 'w') as gene_file, \
                open(cmap_filename, 'w') as cmap_file, \
                open(pos_filename, 'w') as pos_file:

            pos_file.write('\n'.join(pos_data) + '\n')
            gene_file.write('\n'.join(gene_data) + '\n')
            cmap_file.write('\n'.join(
                ['\t'.join([clade, m]) for m in species_markers]) + '\n')
    return input_args
