import bz2
import gzip
import pickle as pickle
from Bio import SeqIO
from os.path import isdir, isfile
from os import makedirs
import pandas as pd
import re

from samestr.utils import clade_path, write_json, write_tsv

import logging
LOG = logging.getLogger(__name__)    

def generate_db(input_args):
    """
    Function expands and generates a database of clade markers 
    from MetaPhlAn or mOTUs database that is required for SameStr.
    """
    # pull out manifest data
    db_manifest = input_args['samestr_db']['db_manifest']
    db_clades = input_args['samestr_db']['db_clades']
    db_existed = input_args['samestr_db']['db_existed']
    db_taxonomy = input_args['samestr_db']['db_taxonomy']

    # read db version file before starting extraction
    with open(input_args['db_version'], 'r') as file:
        db_version_lines = file.readlines()

    # if version file says MetaPhlAn
    if db_version_lines[0].startswith('mpa_'):
        db_tool = 'MetaPhlAn'
        db_version = db_version_lines[0].strip()

        ## and database exists and is..
        if db_existed and 'name' in db_manifest['database']:
    
            ## not MetaPhlAn or versions don't match: exit
            if db_manifest['database']['name'] != db_tool or db_manifest['database']['version'] != db_version:
                LOG.error('Database source in manifest (%s, %s) does not match '
                          'source identified in `--db-version` file (%s, %s).' % (db_manifest['database']['name'], db_manifest['database']['version'],
                                                                                         db_tool, db_version))
                exit(1)
        
        ## and pickle file not provided: exit
        if not input_args['markers_info'][0].endswith('.pkl'):
            LOG.error('Metaphlan pickle file not provided with `--markers-info`. Required when using MetaPhlAn as the source database.')
            exit(1)

    # if version file says mOTUs
    elif any([v.startswith('motus') for v in db_version_lines]):

        ## and info file likely incorrect: exit
        if not input_args['markers_info'][0].endswith('-mOTUs.tsv') or not input_args['markers_info'][1].endswith('-mOTUs.tsv'):
            LOG.error('File(s) provided with `--markers-info` seem to not have the expected file extension for mOTUs (for example: db_mOTU_taxonomy_meta-mOTUs.tsv).')
            exit(1)

        # and centroid version can't be found: exit
        if not any([v.startswith('cen') for v in db_version_lines]):
            LOG.error('Could not parse mOTUs centroid database version in `--db-version` file.')
            exit(1)

        db_tool = 'mOTUs'
        db_version = [v.strip().split('\t')[1] for v in db_version_lines if v.startswith('cen')][0]

        ## and database exists and is..
        if db_existed and 'name' in db_manifest['database']:

            ## not mOTUs or versions don't match: exit
            if db_manifest['database']['name'] != db_tool or db_manifest['database']['version'] != db_version:
                LOG.error('Database source in manifest (%s, %s) does not match '
                          'source identified in `--db-version` file (%s, %s).' % (db_manifest['database']['name'], db_manifest['database']['version'],
                                                                                         db_tool, db_version))
                exit(1)
                
    # if version not detectable: exit
    else:
        LOG.error('Could not identify database source from `--db-version` file (%s).' % (input_args['db_version']))
        exit(1)

    # print info and insert data into manifest var
    LOG.info('Source of provided data: %s, %s' % (db_tool, db_version))
    if db_existed:
        LOG.info('Collecting data to update existing database.')
    else:
        db_manifest['database']['name'] = db_tool
        db_manifest['database']['version'] = db_version 
        LOG.info('Collecting data to generate database.')

    # set markers
    orig_clades = set(db_clades['records'].keys()) # either empty or existing
    all_clades = set()
    all_markers = {}
    n_markers = 0
    all_taxonomy = []

    # get markers fasta
    if input_args['markers_fasta'].endswith('.bz2'):
        markers_fasta = SeqIO.to_dict(
            SeqIO.parse(bz2.open(input_args['markers_fasta'], 'rt'), 'fasta'))
    else:
        markers_fasta = SeqIO.to_dict(
            SeqIO.parse(input_args['markers_fasta'], 'fasta'))

    # get markers mapping   
    if db_tool == 'MetaPhlAn':
        markers_info = pickle.load(bz2.BZ2File(input_args['markers_info'][0]))
    elif db_tool == 'mOTUs':
        markers_info = {"markers": {}}
        
        # read taxonomy info files
        df1 = pd.read_csv(input_args['markers_info'][0], sep='\t')
        df2 = pd.read_csv(input_args['markers_info'][1], sep='\t')

        # Function to find the ID column
        def find_id_column(df):
            for col in df.columns:
                if re.search(r'mOTU_v\d+_ID', col):
                    return col
            return None
        
        # Rename columns from ~meta-mOTU_v2_ID, ~ref-mOTU_v3_ID etc to 'clade'
        id_col1 = find_id_column(df1)
        id_col2 = find_id_column(df2)
        if id_col1 and id_col2:
            df1 = df1.rename(columns={id_col1: 'clade', 'mOTU': 'species'})
            df2 = df2.rename(columns={id_col2: 'clade', 'mOTU': 'species'})
        else:
            raise ValueError("mOTUs ID column not found in one or both tables: %s" % (', '.join(input_args['markers_info'])))

        # Select relevant columns
        columns_to_keep = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'clade']
        df1 = df1[columns_to_keep]
        df2 = df2[columns_to_keep]

        # Concatenate the tables and merge taxonomy
        taxonomy_df = pd.concat([df1, df2], ignore_index=True)
       
        # Populate marker_info based on db_mOTU_DB_CEN.fasta
        for marker in markers_fasta:
            clade, _ = marker.split('.', 1)
            if marker not in markers_info['markers']:
                if clade != 'NA' and clade in taxonomy_df['clade'].values:
                    # TODO: this is easy but expensive, change
                    taxon = '|'.join(taxonomy_df.loc[taxonomy_df['clade'] == clade, ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']].values[0])
                    markers_info['markers'][marker] = {'clade': clade, 'taxon': taxon}

    for marker in markers_info['markers']:
        n_markers += 1

        # retrieve clade names
        clade = markers_info['markers'][marker]['clade']
        all_clades.add(clade)
        if clade in orig_clades:
            continue

        # get only selected clades
        if 'clade' in input_args and input_args['clade'] is not None:
            if clade not in input_args['clade']:
                continue

        if clade not in all_markers:
            all_markers[clade] = set()
            all_taxonomy.append({'clade': markers_info['markers'][marker]['clade'], 
                                 'taxon': markers_info['markers'][marker]['taxon']})   
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
        clade_diff = [c for c in clade_diff if c not in orig_clades]
        if len(clade_diff) > 0:
            LOG.debug('Clades not found in the database: %s' %
                      ', '.join(clade_diff))
        all_clades = clade_intersect

    added_clades = {}
    for clade in all_clades:

        if clade in orig_clades:
            LOG.debug('Clade was already present in the database: %s' % clade)
            continue

        # set output directory and names
        output_dir = input_args['output_dir'] + '/' + clade_path(clade)

        # Create dir if not exists
        if not isdir(output_dir):
            makedirs(output_dir)

        output_base = output_dir + clade
        marker_filename = output_base + '.markers.fa.gz'
        pos_filename = output_base + '.positions.txt.gz'
        gene_filename = output_base + '.gene_file.txt.gz'
        cmap_filename = output_base + '.contig_map.txt.gz'

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

        with gzip.open(marker_filename, 'wt') as marker_file:

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

        with gzip.open(gene_filename, 'wt') as gene_file, \
                gzip.open(cmap_filename, 'wt') as cmap_file, \
                gzip.open(pos_filename, 'wt') as pos_file:

            pos_file.write('\n'.join(pos_data) + '\n')
            gene_file.write('\n'.join(gene_data) + '\n')
            cmap_file.write('\n'.join(
                ['\t'.join([clade, m]) for m in clade_markers]) + '\n')

        # data for manifest        
        fcount = 0
        for fn in [marker_filename, pos_filename, cmap_filename, gene_filename]:
            if isfile(fn):
                fcount += 1
        added_clades.update({clade: {'n_files': fcount, 
                                     'n_markers': n_markers, 
                                     'n_positions': total_marker_len}})
        
    
    # update taxonomy
    added_taxonomy = pd.DataFrame([entry for entry in all_taxonomy if entry['clade'] in added_clades])
    added_taxonomy[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = added_taxonomy['taxon'].str.split('|', expand=True).iloc[:, :7]
    db_taxonomy_combined = pd.concat([db_taxonomy['records'], added_taxonomy.drop(columns=['taxon'], inplace=False)], ignore_index=True).drop_duplicates(subset='clade')
    write_tsv(db_taxonomy['fpath'], db_taxonomy_combined)

    # update clades
    db_clades['records'].update(added_clades)
    write_json(db_clades['fpath'], db_clades)

    # update manifest
    db_manifest['records']['total_n_files'] = sum(item.get('n_files', 0) for item in db_clades['records'].values())
    db_manifest['records']['total_n_clades'] = len(db_clades['records'].keys())
    db_manifest['records']['total_n_markers'] = sum(item.get('n_markers', 0) for item in db_clades['records'].values())
    db_manifest['records']['total_n_positions'] = sum(item.get('n_positions', 0) for item in db_clades['records'].values())
    write_json(db_manifest['fpath'], db_manifest)
    
    return input_args

        