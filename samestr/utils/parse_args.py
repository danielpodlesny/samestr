import argparse as ap

from sys import argv, stderr

from samestr import __author__, __version__

def read_params():
    """Read command line arguments and return them as a dictionary."""
    parser = ap.ArgumentParser(
        prog='samestr',
        description='Welcome to SameStr! SameStr identifies shared strains'
                    ' between pairs of metagenomic samples '
                    'based on the similarity of their Single Nucleotide Variant (SNV) profiles.',
        formatter_class=ap.ArgumentDefaultsHelpFormatter
    )

    # Show help by default if no arguments are passed
    if len(argv) < 2:
        parser.print_help(stderr)
        exit(1)

    # Retrieve version of the program
    parser.add_argument(
        '--version',
        action='version',
        version=f'samestr-{__version__}',
        help='Show version of the program and exit.'
    )

    # Print citation
    parser.add_argument(
        '--citation',
        required=False,
        metavar='STR',
        nargs='?',
        const='Text',
        choices=['Text', 'BibTex', 'Endnote', 'RIS', 'DOI'],
        type=str,
        help='Print citation and exit. '
             'Options: Text, BibTex, Endnote, RIS, DOI.')

    # Set communication verbosity
    parser.add_argument(
        '--verbosity',
        required=False,
        default='INFO',
        metavar='STR',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        type=str,
        help='Set the verbosity of the program. '
             'Options: DEBUG, INFO, WARNING, ERROR, CRITICAL.')

    parser.add_argument(
        '--random-seed',
        type=int,
        help='Set random seed for reproducible runs (affects compare, extract, filter.)'
    )

    # add subparsers for different commands
    subparser = parser.add_subparsers(
        title='commands',
        dest='command',
        help='Use one of the following commands for different tasks:',)
    # subparser.required = True  # http://bugs.python.org/issue9253#msg186387

    # DB
    db_parser = subparser.add_parser(
        'db',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Make database from MetaPhlAn or mOTUs markers.')

    # integrity check
    db_intergrity = db_parser.add_argument_group('Database check arguments')
    db_intergrity.add_argument(
        '--db-check',
        action='store_true',
        help='Performs just a database integrity check, if an existing SameStr database is provided. All other options will be ignored.'
    ) 
    db_intergrity.add_argument(
        '--marker-dir',
        required=False,
        metavar='DIR',
        type=str,
        help='Path to existing MetaPhlAn or mOTUs clade marker database.')

    # input
    db_input = db_parser.add_argument_group('Input arguments')
    db_input.add_argument(
        '--db-version',
        required=False,
        default='PATH',
        type=str,
        help='Path to version file of database (`mpa_latest` for MetaPhlAn, or `db_mOTU_versions` for mOTUs.)'
    )   
    db_input.add_argument(
        '--markers-fasta',
        required=False,
        metavar='FASTA',
        type=str,
        help='Markers fasta file (MetaPhlAn [mpa_vJan21_CHOCOPhlAnSGB_202103.fna or higher], or mOTUs [db_mOTU_DB_CEN.fasta]'
    )
    db_input.add_argument(
        '--markers-info',
        required=False,
        metavar='PATH',
        nargs='+',
        type=str,
        help='Path(s) to marker info files. For MetaPhlAn should be MetaPhlAn pickle file (mpa_vJan21_CHOCOPhlAnSGB_202103.pkl or higher). For mOTUs, should be both `db_mOTU_taxonomy_meta-mOTUs.tsv` and  `db_mOTU_taxonomy_ref-mOTUs.tsv`.')
    db_input.add_argument(
        '--clade',
        required=False,
        metavar='CLADE',
        nargs='+',
        type=str,
        help='Clade to process from input files. Processing all if not specified. '
             'Names must correspond to the taxonomy of the respective database '
             '[e.g. t__SGB10068 for MetaPhlAn or ref_mOTU_v3_00095 for mOTUs]')

    # output
    db_output = db_parser.add_argument_group('Output arguments')
    db_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_db/',
        type=str,
        help='Path to output directory.')

    # CONVERT
    convert_parser = subparser.add_parser(
        'convert',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Convert sequence alignments to SNV Profiles.')

    # general
    convert_general = convert_parser.add_argument_group('General arguments')
    convert_general.add_argument(
        '--nprocs',
        required=False,
        default=1,
        metavar='INT',
        type=int,
        help='The number of processing units to use.')
    convert_general.add_argument(
        '--tmp-dir',
        required=False,
        metavar='DIR',
        default='tmp/',
        type=str,
        help='Path to temporary directory')
    convert_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')
    convert_general.add_argument(
        '--db-force',
        action='store_true',
        help='Force execution, even when database version is not an exact match.'
    ) 

    # input
    convert_input = convert_parser.add_argument_group('Input arguments')
    convert_input.add_argument(
        '--input-files',
        required=True,
        metavar='SAM|SAM.BZ2|BAM',
        nargs='+',
        default=[],
        type=str,
        help='Path to input MetaPhlAn or mOTUs marker alignments (.sam, .sam.bz2, .bam).')
    
    # convert_input.add_argument(
    #     '--input-dir',
    #     required=True,
    #     metavar='INPUT_DIR',
    #     default='.',
    #     type=str,
    #     help='Path to input MetaPhlAn or mOTUs marker alignments (.sam, .sam.bz2, .bam).')
    
    convert_input.add_argument(
        '--recursive-input',
        required=False,
        default=False,
        action='store_true',
        help='Allow checking subdirectories of the input directory for input files.')

    # output
    convert_output = convert_parser.add_argument_group('Output arguments')
    convert_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_convert/',
        type=str,
        help='Path to output directory.')
    convert_output.add_argument(
        '--keep-tmp-files',
        required=False,
        default=False,
        action='store_true',
        help='Keeps intermediate files from transformation steps on disk.')

    # mapping
    convert_alignment = convert_parser.add_argument_group(
        'Alignment arguments')
    convert_alignment.add_argument(
        '--min-aln-identity',
        required=False,
        metavar='FLOAT',
        default=0.9,
        type=float,
        help='Minimum percent identity in alignment.')
    convert_alignment.add_argument(
        '--min-aln-len',
        required=False,
        metavar='INT',
        default=40,
        type=int,
        help='Minimum alignment length.')
    convert_alignment.add_argument(
        '--tax-profiles-dir',
        required=False,
        metavar='DIR',
        type=str,
        help='Path to directory with taxonomic profiles (MetaPhlAn or mOTUs, default extension: .profile.txt). '
        'When not specified, will look for profiles in `input-files` directory.'
    )
    convert_alignment.add_argument(
        '--tax-profiles-extension',
        required=False,
        metavar='EXT',
        default='.profile.txt',
        type=str,
        help='File extension of taxonomic profiles.'
    )
    convert_alignment.add_argument(
        '--min-base-qual',
        required=False,
        metavar='INT',
        default=20,
        type=int,
        help='Minimum base call quality.')
    convert_alignment.add_argument(
        '--min-aln-qual',
        required=False,
        metavar='INT',
        default=0,
        type=int,
        help='Minimum alignment quality. Increasing this threshold can '
             'drastically reduce the number of considered variants.')
    convert_alignment.add_argument(
        '--min-vcov',
        required=False,
        metavar='INT',
        default=3,
        type=int,
        help='Minimum vertical coverage.')

    # EXTRACT
    extract_parser = subparser.add_parser(
        'extract',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Extract SNV Profiles from Reference Genomes.')

    # general
    extract_general = extract_parser.add_argument_group('General arguments')
    extract_general.add_argument(
        '--nprocs',
        required=False,
        default=1,
        metavar='INT',
        type=int,
        help='The number of processing units to use.')
    extract_general.add_argument(
        '--tmp-dir',
        required=False,
        metavar='DIR',
        default='tmp/',
        type=str,
        help='Path to temporary directory')
    extract_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')

    # input
    extract_input = extract_parser.add_argument_group('Input arguments')
    extract_input.add_argument(
        '--input-files',
        required=True,
        metavar='FASTA|FNA|FA|FASTA.GZ|FNA.GZ|FA.GZ',
        nargs='+',
        default=[],
        type=str,
        help='Reference genomes in fasta format.')
    extract_input.add_argument(
        '--clade',
        required=True,
        metavar='CLADE',
        type=str,
        help='Clade to process from input files. '
             'Names must correspond to the taxonomy of the respective database '
             '[e.g. t__SGB10068 for MetaPhlAn or ref_mOTU_v3_00095 for mOTUs]')

    # output
    extract_output = extract_parser.add_argument_group('Output arguments')
    extract_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_extract/',
        type=str,
        help='Path to output directory.')
    extract_output.add_argument(
        '--keep-tmp-files',
        required=False,
        default=False,
        action='store_true',
        help='If not working from memory, '
             'keeps extracted alignments per sample on disk.')
    extract_output.add_argument(
        '--save-marker-aln',
        required=False,
        action='store_true',
        help='Keep alignment files for individual markers.')

    # alignment
    extract_alignment = extract_parser.add_argument_group(
        'Alignment arguments')
    extract_alignment.add_argument(
        '--aln-program',
        required=False,
        default='muscle',
        choices=['muscle', 'mafft'],
        type=str,
        help='Program to use for alignment of marker sequences.')
    extract_alignment.add_argument(
        '--marker-trunc-len',
        required=False,
        metavar='INT',
        default=0,
        type=int,
        help='Number of Nucleotides to be cut from each side of a marker.')

    # MERGE
    merge_parser = subparser.add_parser(
        'merge',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Merge SNV Profiles from multiple sources.')

    # general
    merge_general = merge_parser.add_argument_group('General arguments')
    merge_general.add_argument(
        '--nprocs',
        required=False,
        metavar='INT',
        default=1,
        type=int,
        help='The number of processing units to use.')
    merge_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')

    # input
    merge_input = merge_parser.add_argument_group('Input arguments')
    merge_input.add_argument(
        '--input-files',
        nargs='+',
        required=False,
        metavar='NPY',
        default=[],
        type=str,
        help='Path to input SNV profiles. Should have .npy, .npz or .npy.gz extension.')

    merge_input.add_argument(
        '--input-dir',
        required=False,
        metavar='INPUT_DIR',
        type=str,
        help='Path to input SNV profiles. Should have .npy, .npz or .npy.gz extension.')

    # output
    merge_output = merge_parser.add_argument_group('Output arguments')
    merge_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_merge/',
        type=str,
        help='Path to output directory.')

    # clades
    merge_clades = merge_parser.add_argument_group('Clade arguments')
    merge_clades.add_argument(
        '--clade',
        required=False,
        metavar='CLADE',
        nargs='+',
        type=str,
        help='Clade to process from input files. Processing all if not specified. '
             'Names must correspond to the taxonomy of the respective database '
             '[e.g. t__SGB10068 for MetaPhlAn or ref_mOTU_v3_00095 for mOTUs]')
    
    # FILTER
    filter_parser = subparser.add_parser(
        'filter',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Filter SNV Profiles.')

    # general
    filter_general = filter_parser.add_argument_group('General arguments')
    filter_general.add_argument(
        '--nprocs',
        required=False,
        default=1,
        metavar='INT',
        type=int,
        help='The number of processing units to use.')
    filter_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')

    # input
    filter_input = filter_parser.add_argument_group('Input arguments')
    filter_input.add_argument(
        '--input-files',
        required=True,
        nargs='+',
        metavar='NPY',
        default=[],
        type=str,
        help='Path to input SNV Profiles. Should have .npy, .npz or .npy.gz extension.')
    filter_input.add_argument(
        '--input-names',
        required=True,
        nargs='+',
        metavar='TXT',
        default=[],
        type=str,
        help='Path to input name files.')

    # output
    filter_output = filter_parser.add_argument_group('Output arguments')
    filter_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_filter/',
        type=str,
        help='Path to output directory.')
    filter_output.add_argument(
        '--keep-poly',
        required=False,
        action='store_true',
        help='Keep only positions that are polymorphic in at least one sample')
    filter_output.add_argument(
        '--keep-mono',
        required=False,
        action='store_true',
        help='Keep only positions that are monomorphic in all samples')
    filter_output.add_argument(
        '--delete-pos',
        required=False,
        action='store_true',
        help='Delete masked marker and global positions from array instead of np.nan'
    )

    # filter settings

    # clades
    filter_clade = filter_parser.add_argument_group(
        'Clade arguments')
    filter_clade.add_argument(
        '--clade',
        required=False,
        metavar='CLADE',
        nargs='+',
        type=str,
        help='Clade to process from input files. Processing all if not specified. '
              'Names must correspond to the taxonomy of the respective database '
              '[e.g. t__SGB10068 for MetaPhlAn or ref_mOTU_v3_00095 for mOTUs]')
    filter_clade.add_argument(
        '--clade-min-samples',
        required=False,
        metavar='INT',
        default=2,
        type=int,
        help='Skipping clades with fewer than `clade-min-samples` samples.')

    # markers
    filter_markers = filter_parser.add_argument_group(
        'Clade Marker arguments')
    filter_markers.add_argument(
        '--marker-remove',
        required=False,
        metavar='TXT',
        type=str,
        help='List of Markers to remove for selected clade. '
             'Requires `clade` to be specified.')
    filter_markers.add_argument(
        '--marker-keep',
        required=False,
        metavar='TXT',
        type=str,
        help='List of Markers to keep for selected clade. '
             'Requires `clade` to be specified. Overrides `marker-remove`.')
    filter_markers.add_argument(
        '--marker-trunc-len',
        required=False,
        metavar='INT',
        default=50,
        type=int,
        help='Number of Nucleotides to be cut from each two sides of a marker.')

    # sample variants
    filter_sample_pos = filter_parser.add_argument_group(
        'Sample Variant Filtering arguments')
    filter_sample_pos.add_argument(
        '--sample-var-min-n-vcov',
        required=False,
        metavar='INT',
        default=2,
        type=int,
        help='Remove variants with coverage below `sample-var-min-n-vcov` nucleotides.')
    filter_sample_pos.add_argument(
        '--sample-var-min-f-vcov',
        required=False,
        metavar='FLOAT',
        default=0.05,
        type=float,
        help='Remove variants with coverage below `sample-var-min-f-vcov` percent.')

    # sample positions
    filter_sample_pos = filter_parser.add_argument_group(
        'Sample Position Filtering arguments')
    filter_sample_pos.add_argument(
        '--sample-pos-min-n-vcov',
        required=False,
        metavar='INT',
        default=1,
        type=int,
        help='Remove positions with coverage below `sample-pos-min-n-vcov` nucleotides.'
    )
    filter_sample_pos.add_argument(
        '--sample-pos-min-sd-vcov',
        required=False,
        metavar='FLOAT',
        default=3.0,
        type=float,
        help='Remove positions with coverage +-`sample-pos-min-sd-vcov` from the mean.'
    )

    # samples
    filter_samples = filter_parser.add_argument_group(
        'Sample Filtering arguments')
    filter_samples.add_argument(
        '--samples-select',
        required=False,
        nargs='+',
        metavar='TXT',
        default=[],
        type=str,
        help='Path to names file with subsample of input names.')
    filter_samples.add_argument(
        '--samples-min-n-hcov',
        required=False,
        metavar='INT',
        type=int,
        default=5000,
        help='Remove samples with horizontal coverage below `samples-min-n-hcov`.')
    filter_samples.add_argument(
        '--samples-min-f-hcov',
        required=False,
        metavar='FLOAT',
        type=float,
        help='Remove samples with fraction of horizontal coverage below `samples-min-f-hcov`.')
    filter_samples.add_argument(
        '--samples-min-m-vcov',
        required=False,
        metavar='FLOAT',
        type=float,
        help='Remove samples with mean coverage below `samples-min-m-vcov`.')

    ### global positions
    filter_global_pos = filter_parser.add_argument_group(
        'Global Position Filtering arguments')
    filter_global_pos.add_argument(
        '--global-pos-min-n-vcov',
        required=False,
        metavar='INT',
        default=2,
        type=int,
        help='Remove positions covered by fewer than `global-pos-min-n-vcov` number of samples. '
        'Overrides `global-pos-min-f-vcov`.')
    filter_global_pos.add_argument(
        '--global-pos-min-f-vcov',
        required=False,
        metavar='FLOAT',
        default=False,
        type=float,
        help='Remove positions covered by fewer than `global-pos-min-f-vcov` fraction of samples.')

    # STATS
    stats_parser = subparser.add_parser(
        'stats',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Report alignment statistics.')
    # general
    stats_general = stats_parser.add_argument_group('General arguments')
    stats_general.add_argument(
        '--nprocs',
        required=False,
        default=1,
        metavar='INT',
        type=int,
        help='The number of processing units to use.')
    stats_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')

    # input
    stats_input = stats_parser.add_argument_group('Input arguments')
    stats_input.add_argument(
        '--input-files',
        required=True,
        nargs='+',
        metavar='NPY',
        default=[],
        type=str,
        help='Path to input SNV profiles. Should have .npy, .npz or .npy.gz extension.')
    stats_input.add_argument(
        '--input-names',
        required=True,
        nargs='+',
        metavar='FILEPATH',
        default=[],
        type=str,
        help='Path to input name files.')

    # output
    stats_output = stats_parser.add_argument_group('Output arguments')
    stats_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='marker_db/',
        type=str,
        help='Path to output directory.')
    stats_output.add_argument(
        '--dominant-variants',
        required=False,
        action='store_true',
        help='Report statistics only for dominant variants as obtained from consensus call.')
    
    # compare subparser
    compare_parser = subparser.add_parser(
        'compare',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Calculate pairwise sequence similarity.')
    # general
    compare_general = compare_parser.add_argument_group('General arguments')
    compare_general.add_argument(
        '--nprocs',
        required=False,
        default=1,
        metavar='INT',
        type=int,
        help='The number of processing units to use.')
    compare_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')

    # input
    compare_input = compare_parser.add_argument_group('Input arguments')
    compare_input.add_argument(
        '--input-files',
        required=True,
        nargs='+',
        metavar='NPY',
        default=[],
        type=str,
        help='Path to input SNV Profiles. Should have .npy, .npz or .npy.gz extension.')
    compare_input.add_argument(
        '--input-names',
        required=True,
        nargs='+',
        metavar='TXT',
        default=[],
        type=str,
        help='Path to input name files.')

    # output
    compare_output = compare_parser.add_argument_group('Output arguments')
    compare_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_compare/',
        type=str,
        help='Path to output directory.')

    compare_output.add_argument(
        '--dominant-variants',
        required=False,
        action='store_true',
        help='Compare only dominant variants as obtained from consensus call.')
    compare_output.add_argument(
        '--dominant-variants-added',
        required=False,
        action='store_true',
        help='Add dominant variants as additional entries to data.')
    compare_output.add_argument(
        '--dominant-variants-msa',
        required=False,
        action='store_true',
        help='Output alignment of dominant variants as fasta.')

    # SUMMARIZE
    summarize_parser = subparser.add_parser(
        'summarize',
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
        help='Summarize Taxonomic Co-Occurrence.')
    
    # general
    summarize_general = summarize_parser.add_argument_group('General arguments')
    summarize_general.add_argument(
        '--marker-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to MetaPhlAn or mOTUs clade marker database.')

    # input
    summarize_input = summarize_parser.add_argument_group('Input arguments')
    summarize_input.add_argument(
        '--input-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to `samestr compare` output directory. Must contain '
             'pairwise comparison of alignment files (.fraction.txt, .overlap.txt)'
    )
    summarize_input.add_argument(
        '--tax-profiles-dir',
        required=True,
        metavar='DIR',
        type=str,
        help='Path to directory with taxonomic profiles (MetaPhlAn or mOTUs, default extension: .profile.txt).'
    )
    summarize_input.add_argument(
        '--tax-profiles-extension',
        required=False,
        metavar='EXT',
        default='.profile.txt',
        type=str,
        help='File extension of taxonomic profiles.'
    )

    # output
    summarize_output = summarize_parser.add_argument_group('Output arguments')
    summarize_output.add_argument(
        '--output-dir',
        required=False,
        metavar='DIR',
        default='out_summarize/',
        type=str,
        help='Path to output directory.')

    # samples
    summarize_thresholds = summarize_parser.add_argument_group(
        'Summary threshold arguments')
    summarize_thresholds.add_argument(
        '--aln-pair-min-overlap',
        required=False,
        metavar='INT',
        default=5000,
        type=int,
        help='Minimum number of overlapping positions which have to be covered in both '
        'alignments in order to evaluate alignment similarity.'
    )
    summarize_thresholds.add_argument(
        '--aln-pair-min-similarity',
        required=False,
        metavar='FLOAT',
        default=0.999,
        type=float,
        help='Minimum pairwise alignment similarity to call shared strains.'
    )

    return vars(parser.parse_args())
