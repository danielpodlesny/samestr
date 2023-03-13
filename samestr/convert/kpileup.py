#!/usr/bin/env python3
import argparse
import pysam
import sys
from collections import Counter

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

usage = f"""{sys.argv[0]} <sampleId> <bam file> <gene file> <minimum BQ> <minimum MQ> <minimum D>

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


def parse_gene(gene_file):
    """
    Parse the input gene file

    Args:
        gene_file (str): the gene file name

    Returns:
        dict: a dictionary containing the gene name as key and the contig, start, end, strand, and sequence as values
    """
    data = {}
    with open(gene_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            gene = fields[1]
            data[gene] = {
                'contig': fields[0],
                'begin': int(fields[2]),
                'end': int(fields[3]),
                'strand': fields[4],
                'seq': fields[5],
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
    nuc = {}
    for gene, values in genes.items():
        begin = values['begin']
        end = values['end']
        c = values['contig']
        strand = values['strand']
        seq = values['seq']
        for i in range(begin, end + 1):
            codon_pos = (i - begin + 1) % 3
            if strand == '-' and codon_pos != 2:
                codon_pos = 1 if codon_pos == 0 else 0
            codon_pos = 3 if codon_pos == 0 else codon_pos
            nuc[c, i] = (gene, seq[i - begin], codon_pos)
    return nuc


def pileup(sample_id, bam_file, gene_file, min_bq, min_mq, min_d):
    """
    Parse the BAM file and look for reads overlapping with the target genes and report the pileup.
    Reads must pass the minimum MQ threshold.
    For each position, an allele must have a depth >= min_d.
    The frequency of each allele is calculated after filtering.

    Args:
        sample_id (str): the sample ID
        bam_file (str): the BAM file name
        gene_file (str): the gene file name
        min_bq (int): the minimum base quality score of a sequenced base
        min_mq (int): the minimum MQ mapping score of the aligned reads
        min_d (int): the minimum depth, an integer.

    Returns:
        None
    """
    genes = parse_gene(gene_file)
    for key in genes.keys():
        print(key + '===')
    bases = parse_bases(genes)

    print('Sample\tContig\tPosition\tGene\tRef\tCodon\tConsensus\tAllele\tCounts\tFrequency')
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for genename, gene in genes.items():
        for pileupcolumn in samfile.pileup(
                contig=gene['contig'],
                start=gene['begin'],
                stop=gene['end'],
                min_base_quality=min_bq,
                min_mapping_quality=min_mq,
                truncate=True, max_depth=100000):
            if pileupcolumn.n < min_d:
                continue
            # 1 for the reference base
            ref_base = bases[genename, pileupcolumn.pos + 1][1]
            # 2 for the codon position
            codon = bases[genename, pileupcolumn.pos + 1][2]
            counts = Counter()
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    query_base = pileupread.alignment.query_sequence[pileupread.query_position]
                    counts[query_base] += 1

            # each base must have a depth >= min_d
            counts = Counter({k: v for k, v in counts.items() if v >= min_d})
            total_count = sum(counts.values())
            if total_count < min_d:
                continue
            consensus = max(counts, key=counts.get)
            for base, count in counts.items():
                freq = int(count/total_count*100)
                print(
                    f"{sample_id}\t{genename}\t{pileupcolumn.pos+1}\t{gene['contig']}\t{ref_base}\t{codon}\t{consensus}\t{base}\t{count}\t{freq}")
    samfile.close()


def main():
    if len(sys.argv) != 7:
        print(
            f"Usage: {sys.argv[0]} <sampleId> <bam file> <gene file> <minimum BQ> <minimum MQ> <minimum D>")
        sys.exit(1)

    sample_id = sys.argv[1]
    bam_file = sys.argv[2]
    gene_file = sys.argv[3]
    min_bq = int(sys.argv[4])
    min_mq = int(sys.argv[5])
    min_d = int(sys.argv[6])

    pileup(sample_id, bam_file, gene_file, min_bq, min_mq, min_d)


if __name__ == "__main__":
    main()
