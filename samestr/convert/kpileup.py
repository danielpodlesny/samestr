#!/usr/bin/env python3
import sys
import subprocess
import argparse
import re
from collections import defaultdict

parser = argparse.ArgumentParser(
    description='This script is for parsing the BAM file and look for reads overlapping with the target genes and report the pileup.')
parser.add_argument('sample_id', help='sample ID')
parser.add_argument('bam_file', help='an aligned bam file')
parser.add_argument(
    'gene_file', help='a tab-delimited file of six columns in this order: contigId, geneId, begin, end, strand, DNA transcript seq. (Note: begin<end)')
parser.add_argument('min_bq', type=int,
                    help='the minimum base quality score of a sequenced base')
parser.add_argument('min_mq', type=int,
                    help='the minimum MQ mapping score of the aligned reads')
parser.add_argument('min_d', type=int, help='the minimum depth, an integer.')
args = parser.parse_args()

usage = f"""{sys.argv[0]} <sampleId> <bam file> <gene file> <minimum BQ> <minimum MQ> <minimum D> <maximum D>

{args.sample_id}: sample ID
{args.bam_file}: an aligned bam file
{args.gene_file}: a tab-delimited file of six columns in this order: contigId, geneId, begin, end, strand, DNA transcript seq. (Note: begin<end)
{args.min_bq}: the minimum base quality score of a sequenced base
{args.min_mq}: the minimum MQ mapping score of the aligned reads
{args.min_d}: the minimum depth, an integer.

This script is for parsing the BAM file and look for reads overlapping with the target genes and report the pileup.
Reads must pass the minimum MQ threshold.
For each position, an allele must have a depth >= <minimum D>.
The frequency of each allele is calculated after the filterings.

Note: the quality score must be Sanger Phred score (ascii 33).

Dependencies:
Use Samtools

CIGAR parser: the embedded cigar parser is not versatile, please review it to make sure that it is handling the cigar code appropriately.


The script generates a tab-delimited table directly to STDOUT.
Each row is a base of the queried genes.
A base could be reported repeatedly up to four times as four rows when polymorphisms are observed: A,T, G or C.
The columns are:
sample ID
contig ID
This Base Position
Gene ID
Ref allele of this Base
Condon Position (if coding gene)
Observed Consensus Allele of this Base (in this BAM)
Observed Allele
Coverage Depth of the Observed Allele
Allele Frequency of the Observed Allele

"""

def parse_gene(file):
    """
    Parse the input gene file

    Args:
        gene_file (str): the gene file name

    Returns:
        dict: a dictionary containing the gene name as key and the contig, start, end, strand, and sequence as values
    """
    data = {}
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            contig_id, gene_id, begin, end, strand, seq = line.split("\t")
            data[gene_id] = {
                'contig': contig_id,
                'begin': int(begin),
                'end': int(end),
                'strand': strand,
                'seq': seq
            }
    return data


def parse_bases(genes):
    """
    Go through each gene and add the nucleotide positions to a dictionary indexed by contigs

    Args:
        genes (dict): a dictionary containing gene data

    Returns:
        dict: a dictionary indexed by contigs and containing the gene name, reference nucleotide, and codon position as values
    """
    nuc = defaultdict(dict)
    for g, gene_data in genes.items():
        begin = gene_data['begin']
        end = gene_data['end']
        c = gene_data['contig']
        strand = gene_data['strand']
        temp = list(gene_data['seq'])

        for i in range(begin, end + 1):
            codon_pos = (i - begin + 1) % 3
            if strand == '-' and codon_pos != 2:
                codon_pos = 1 if codon_pos == 0 else 0
            codon_pos = 3 if codon_pos == 0 else codon_pos
            nuc[c][i] = f"{g}\t{temp[i-begin]}\t{codon_pos}"

    return nuc

def decode_cigar(cigar):
    """
    Decode the cigar string

    Args:
        cigar (str): the cigar string

    Returns:
        str: the decoded cigar string
    """
    cigar_parts = re.findall(r'(\d+)([MIDNSHP])', cigar)
    new_string = ''.join(c * int(n) for n, c in cigar_parts)
    return new_string


def convert_qual(qual_string):
    """
    Convert the quality string to a list of quality scores

    Args:
        qual_string (str): the quality string

    Returns:
        list: a list of quality scores
    """
    scores = [ord(q) - 33 for q in qual_string]
    return scores


def pileup(sample_id, bam_file, gene_file, min_bq, min_mq, min_depth):
    """
    Parse the BAM file and look for reads overlapping with the target genes and report the pileup

    Args:
        sample_id (str): the sample ID
        bam_file (str): the BAM file name
        gene_file (str): the gene file name
        min_bq (int): the minimum base quality score of a sequenced base
        min_mq (int): the minimum MQ mapping score of the aligned reads
        min_depth (int): the minimum depth, an integer.

    Returns:
        None
    """
    genes = parse_gene(gene_file)
    bases = parse_bases(genes)

    f_table = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with subprocess.Popen(["samtools", "view", bam_file], stdout=subprocess.PIPE, universal_newlines=True) as bam_process:
        for line in bam_process.stdout:
            qname, flag, rname, begin, mapq, cigar, mrnm, mpos, isize, seq, qual, *info = line.strip().split('\t')

            if rname == "*" or int(mapq) < min_mq or rname not in bases:
                continue

            begin = int(begin)
            end = begin + len(seq) - 1
            qual_scores = convert_qual(qual)

            s = decode_cigar(cigar)
            b = list(seq)
            ci = list(s)

            new = []
            read_i = 0

            for cigar_i in range(len(ci)):
                base = "-"
                if ci[cigar_i] == "D":
                    base = "-"
                elif ci[cigar_i] == "H":
                    continue
                elif ci[cigar_i] in ["I", "S"]:
                    read_i += 1
                    continue
                elif qual_scores[read_i] < min_bq:
                    base = "-"
                    read_i += 1
                else:
                    base = b[read_i]
                    read_i += 1

                new.append(base)

            b = new
            for i in range(begin, begin + len(b)):
                nuc = b[i - begin]

                if bases[rname].get(i):
                    if nuc != "-":
                        f_table[rname][i][nuc] += 1

    for c in f_table:
        print(f"{c}===")
    print("Sample\tContig\tPosition\tGene\tRef\tCodon\tConsensus\tAllele\tCounts\tFrequency")

    for c in f_table:
        for pos in sorted(f_table[c]):
            total = 0
            major = ""

            for nuc in sorted(f_table[c][pos]):
                counts = f_table[c][pos][nuc]
                if counts < min_depth:
                    continue
                total += counts
                if major == "":
                    major = nuc
                elif f_table[c][pos][major] < counts:
                    major = nuc

            for nuc in sorted(f_table[c][pos]):
                counts = f_table[c][pos][nuc]
                if counts < min_depth:
                    continue
                percent = 100 * counts / total

                print(f"{sample_id}\t{c}\t{pos}\t{bases[c][pos]}\t{major}\t{nuc}\t{counts}\t{percent:.0f}")


def main():

    pileup(args.sample_id, 
           args.bam_file, 
           args.gene_file, 
           args.min_bq, 
           args.min_mq, 
           args.min_d)

if __name__ == "__main__":
    main()
