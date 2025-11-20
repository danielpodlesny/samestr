#!/usr/bin/env python3
import logging
import pathlib
import re
import sys

from collections import defaultdict

import numpy as np

from samestr.convert.buffered_reader import stream_file


LOG = logging.getLogger(__name__)
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP])')


def parse_positions(f):
    positions = {}
    for line in stream_file(f):
        contig_id, _, start, end, _, _ = line.rstrip().split("\t")
        positions[contig_id] = int(start), int(end)
    return positions

def decode_cigar(cigar):
    for n, c in CIGAR_RE.findall(cigar):
        if c != "H":
            yield int(n), c

def convert_qual(qual_string):
    for q in qual_string:
        yield ord(q) - 33


def pileup(bam_stream, gene_file, min_bq, min_mq, min_depth, outstream=sys.stdout):
    bases = parse_positions(gene_file)

    f_table = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    cur_rname = None

    contig_bases = None
    for line in bam_stream:
        _, _, rname, begin, mapq, cigar, _, _, _, seq, qual, *_ = line.strip().split("\t")

        if int(mapq) < min_mq:
            continue

        if rname != cur_rname:
            cur_rname = rname            
            contig_bases = bases.get(rname)

        if not contig_bases:
            continue
        start, end = contig_bases

        p = int(begin)
        
        aln_string = decode_cigar(cigar)
        
        bases_and_quals = iter(zip(seq, convert_qual(qual)))
        
        for oplen, cigar_op in aln_string:
            if cigar_op == "H":
                continue
            elif cigar_op == "D":
                base = "-"                
            else:
                is_insertion = cigar_op in ("I", "S")
                for pp in range(p, p + oplen):
                    cur_base, cur_qual = next(bases_and_quals)
                    if not is_insertion:
                        base = (cur_base, "-")[(cur_qual < min_bq)]
                        if base != "-" and start <= pp <= end:
                            f_table[rname][pp][base] += 1
                
                if is_insertion:
                    continue
            
            p += oplen

    for c, positions in f_table.items():
        for pos, nucs in sorted(positions.items(), key=lambda x:x[0]):
            for nuc, counts in sorted(nucs.items(), key=lambda x:x[0]):
                if counts >= min_depth:
                    print(c, pos, nuc, counts, sep="\t", file=outstream)
                    yield c, pos, nuc, counts


BASES = {b: i for b, i in zip('ACGT', range(4))}

def add_pileups(positions, contigs):
    for contig, pos, allele, count in positions:
        if allele != "N":
            contigs[contig][0, pos - 1, BASES[allele]] = count


def read_contig_map(contig_map):
    # Contig map
    # cmap[genome] = [contig1, contig2, ...]
    cmap = {}
    for line in stream_file(contig_map):
        line = line.strip().split()
        genome = line[0]
        contig = line[1]
        cmap.setdefault(genome, []).append(contig)
    return cmap

def initialise_contigs(gene_file):
    # Initialize numpy arrays for each contig
    contigs = {}
    for line in stream_file(gene_file):
        line = line.rstrip().split()
        contig = line[0]
        end = int(line[3])
        contigs[contig] = np.zeros([1, end, 4])
    return contigs



def kp2np(kpileups, contig_map, sample, gene_file, output_dir):

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    x = initialise_contigs(gene_file)
    add_pileups(kpileups, x)
    cmap = read_contig_map(contig_map)

    # Concatenate contigs
    # -------------------
    m, k = 1, 4
    y = {}
    for genome, contigs in cmap.items():
        n = sum(np.shape(x[c])[1] for c in contigs)
        y[genome] = np.zeros([m, n, k])

        # Add alignment data
        beg = 0
        end = 0
        for contig in contigs:
            end += np.shape(x[contig])[1]
            y[genome][0, beg:end, :] = x[contig]
            beg = end

        species = y[genome]
        cov = species.sum(axis=2)

        n_sites = species.shape[1]
        n_gaps = (cov == 0).sum()
        n_covered = n_sites - n_gaps

        # only write to numpy file if there is coverage left after convert criteria
        if n_covered:
            np.savez_compressed(
                output_dir / f"{genome}.{sample}",
                y[genome],
                allow_pickle=True
            )