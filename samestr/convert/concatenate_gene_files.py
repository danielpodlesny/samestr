from os.path import isfile
import logging

from samestr.utils import ooSubprocess, clade_path
from samestr.summarize.read_taxonomic_profiles import get_clade_profile_dict

LOG = logging.getLogger(__name__)


def concatenate_gene_files(arg):
    """
        concatenates gene_files and contig_maps
        for all clades of a given sample.
    """

    oosp = ooSubprocess.ooSubprocess(tmp_dir=arg['tmp_dir'])

    LOG.debug('Gathering: %s' % arg['tax_profile'])
    clades = get_clade_profile_dict(arg['tax_profile'], db_source=arg['db_source'])
    if not clades:
        LOG.error('No clade classified: %s' % arg['bname'])
        return False

    LOG.debug('Detected clades: %s' % ', '.join(clades))
    gene_files = []
    contig_maps = []

    for clade, _ in sorted(iter(clades.items()),
                             key=lambda k_v: (k_v[1], k_v[0]),
                             reverse=True):
        gene_fname = '%s/%s.gene_file.txt' % (arg['marker_dir'], clade_path(clade, filebase = True))
        map_fname = '%s/%s.contig_map.txt' % (arg['marker_dir'], clade_path(clade, filebase = True))
        LOG.debug('Gene file: %s' % gene_fname)

        if isfile(gene_fname) and isfile(map_fname):
            gene_files.append(gene_fname)
            contig_maps.append(map_fname)

    if not isfile(arg['gene_file']):
        LOG.debug('Running: cat for %s gene_files' % len(gene_files))
        oosp.ex('cat', args=gene_files, out_fn=arg['gene_file'], verbose=False)

    if not isfile(arg['contig_map']) and len(contig_maps):
        LOG.debug('Running: cat for %s contig_maps' % len(contig_maps))
        oosp.ex('cat',
                args=contig_maps,
                out_fn=arg['contig_map'],
                verbose=False)

    LOG.debug('Finished: %s, %s' % (arg['contig_map'], arg['gene_file']))
    return arg
