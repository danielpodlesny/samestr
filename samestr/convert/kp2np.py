#!/usr/bin/env python3
import argparse
import numpy as np
from os.path import isdir, basename
from os import makedirs

# Input arguments
# ---------------

parser = argparse.ArgumentParser()
parser.add_argument('--kp', help='Kpileup alignments (.kp.txt)')
parser.add_argument('--map', help='Map of genomes to contigs (tab-delimited)')
parser.add_argument('--sample', help='Sample Name', type=str, required=True)
parser.add_argument('--gene-file', help='kpileup gene file')
parser.add_argument('--output-dir', help='Output dir', default='./')
args = parser.parse_args()

# Read data
# ---------

# Numpy alignments
# x[contig] = numpy alignment
# y[genome] = concatenated numpy alignments
# Merge kpileups from multiple samples. Write dictionary of (M, N, 4) numpy arrays where:
# M = samples
# N = alignment sites
# 4 = nucleotides (ACGT)
# Entry (i,j,k) of this array corresponds to the count of nucleotide k at position j of sample i

# Initialize data
nts = 'ACGT'
M = 1

# Initialize numpy arrays for each contig
x = {}
for line in open(args.gene_file):
    line = line.rstrip().split()
    contig = line[0]
    beg = int(line[2])
    end = int(line[3])
    x[contig] = np.zeros([M, end, 4])

# Add kpileup results to numpy arrays
with open(args.kp, 'r') as f:
    for line in f.readlines():
        line = line.rstrip().split()
        if len(line) == 10 and line[0] != 'Sample':
            sample = line[0]
            i = 0
            contig = line[1]
            j = int(line[2])
            nt = line[7]
            k = nts.index(nt)
            count = int(line[8])
            x[contig][i, j - 1, k] = count

# Sample list
# M = [sample1]
M = np.array([args.sample])

# Contig map
# cmap[genome] = [contig1, contig2, ...]
cmap = {}
for line in open(args.map):
    line = line.strip().split()
    genome = line[0]
    contig = line[1]
    if genome not in cmap:
        cmap[genome] = []
    cmap[genome].append(contig)

# Create dir if not exists
if not isdir(args.output_dir):
    makedirs(args.output_dir)

# Mapping stats
cols = '\t'.join([
    'alignment', 'genome', 'mean_cov', 'median_cov', 'n_sites', 'n_gaps',
    'n_covered', 'n_mono', 'n_duo', 'n_tri', 'n_quat', 'n_poly', 'f_covered',
    'f_mono', 'f_duo', 'f_tri', 'f_quat', 'f_poly'
])
stats = [cols]

# Concatenate contigs
# -------------------
y = {}
for genome in cmap:
    contigs = cmap[genome]

    # Initialize array
    m = len(M)
    n = sum([np.shape(x[c])[1] for c in contigs])
    k = 4
    y[genome] = np.zeros([m, n, k])

    # Add alignment data
    beg = 0
    end = 0
    for contig in contigs:
        end += np.shape(x[contig])[1]
        y[genome][0, beg:end, :] = x[contig]
        beg = end

    np_filepath = '%s/%s.%s.npy' % (args.output_dir, genome, args.sample)
    np.save(np_filepath, y[genome], allow_pickle=True)

    species = y[genome]
    cov = species.sum(axis=2)

    # coverage [depth]
    mean_cov = round(np.mean(cov), 4)
    median_cov = round(np.median(cov), 4)

    # coverage [width]
    n_sites = species.shape[1]
    n_gaps = (cov == 0).sum()
    n_covered = n_sites - n_gaps

    # n of variant sites, monomorphic, .., polymorphic
    n_mono = ((species > 0).sum(axis=2) == 1).sum()
    n_duo = ((species > 0).sum(axis=2) == 2).sum()
    n_tri = ((species > 0).sum(axis=2) == 3).sum()
    n_quat = ((species > 0).sum(axis=2) == 4).sum()
    n_poly = ((species > 0).sum(axis=2) > 1).sum()

    # fraction of covered sites,
    # fraction of covered sites with variant, monomorphic, .., polymorphic
    if not n_covered == 0:
        f_covered = round(n_covered / n_sites, 4)
        f_mono = round(n_mono / n_covered, 4)
        f_duo = round(n_duo / n_covered, 4)
        f_tri = round(n_tri / n_covered, 4)
        f_quat = round(n_quat / n_covered, 4)
        f_poly = round(n_poly / n_covered, 4)
    else:
        f_covered, f_mono, f_duo, \
            f_tri, f_quat, f_poly = 0, 0, 0, 0, 0, 0

    stat = [
        basename(np_filepath), genome, mean_cov, median_cov, n_sites, n_gaps,
        n_covered, n_mono, n_duo, n_tri, n_quat, n_poly, f_covered, f_mono,
        f_duo, f_tri, f_quat, f_poly
    ]

    stat = [str(s) for s in stat]
    stats.append('\t'.join(stat))

with open('%s/%s.aln_stats.txt' % (args.output_dir, args.sample), 'w') as file:
    file.write('\n'.join(stats))
