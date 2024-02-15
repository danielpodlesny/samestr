from os.path import isfile
import logging

from samestr.utils import ooSubprocess, clade_path, read_json
from samestr.summarize.read_taxonomic_profiles import get_clade_profile_dict, identify_taxonomic_profile_database

LOG = logging.getLogger(__name__)


def concatenate_gene_files(arg):
    """
        concatenates gene_files and contig_maps
        for all clades of a given sample, but check that the databases match, first.
    """

    # pull out manifest data and check database version
    db_manifest = arg['samestr_db']['db_manifest']

    if not 'db_force' in arg:
        profile_tool, profile_version = identify_taxonomic_profile_database(arg['tax_profile'])

        if not (profile_tool == db_manifest['database']['name'] and \
            profile_version == db_manifest['database']['version']):

                LOG.error('Database source and version of taxonomic profile (%s, %s) '
                        'does not match to SameStr database (%s, %s).' % (profile_tool, profile_version, 
                                                                        db_manifest['database']['name'], 
                                                                        db_manifest['database']['version']))
                exit(1)

    # read taxonomic profiles with dedicated reader
    LOG.debug('Gathering: %s' % arg['tax_profile'])
    clades = get_clade_profile_dict([arg['tax_profile']], db_source=db_manifest['database']['name'])
    if not clades:
        LOG.warning('No clade classified: %s' % arg['bname'])
        return None
    LOG.debug('Detected clades: %s' % ', '.join(clades))

    # create a lists of gene and contig map files
    gene_files = []
    contig_maps = []
    for clade, _ in sorted(iter(clades.items()),
                             key=lambda k_v: (k_v[1], k_v[0]),
                             reverse=True):
        gene_fname = '%s/%s.gene_file.txt.gz' % (arg['marker_dir'], clade_path(clade, filebase = True))
        map_fname = '%s/%s.contig_map.txt.gz' % (arg['marker_dir'], clade_path(clade, filebase = True))
        LOG.debug('Gene file: %s' % gene_fname)

        if isfile(gene_fname) and isfile(map_fname):
            gene_files.append(gene_fname)
            contig_maps.append(map_fname)

    # for each sample, concatenate all lists fo gene and contig map files, respectively
    oosp = ooSubprocess.ooSubprocess(tmp_dir=arg['tmp_dir'])
    if not isfile(arg['gene_file']):
        LOG.debug('Running: cat for %s gene_files' % len(gene_files))
        oosp.ex('cat', args=gene_files, out_fn=arg['gene_file'], verbose=False)

    if not isfile(arg['contig_map']) and len(contig_maps):
        LOG.debug('Running: cat for %s contig_maps' % len(contig_maps))
        oosp.ex('cat',
                args=contig_maps,
                out_fn=arg['contig_map'],
                verbose=False)

    LOG.debug('Generated: %s, %s' % (arg['contig_map'], arg['gene_file']))
    return arg
