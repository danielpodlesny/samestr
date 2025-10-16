#!/usr/bin/env python3

""" Filter SAM file by percent identity and length """

import re
import sys

from samestr.convert.buffered_reader import stream_file


REVCOMP_TRANSLATION = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP])')
MISMATCH_RE = re.compile(r'[NX]M:i:(\d+)')


def decode_cigar(cigar):
    start, end, alen, tlen = 0, 0, 0, 0

    for i, cigar_op in enumerate(CIGAR_RE.findall(cigar)):
        try:            
            oplen, op = int(cigar_op[0]), cigar_op[1]
        except IndexError:
            raise IndexError(f"Failed to parse cigar='{cigar}'.")
        if op == "H":
            if i == 0:
                start += oplen
            else:
                end += oplen  # is that actually correct? should not be 'ends here'?
        elif op == "M":
            alen += oplen
        if op not in ("H", "D"):
            tlen += oplen

    return start, end, alen, tlen


def reverse_complement(seq):
    return seq.translate(REVCOMP_TRANSLATION)[::-1]


def filter_alignments(stream, seqid_threshold, minlen, with_header=True,):
    # initialize variables
    cquery = ''
    cseq = ''
    cqual = ''
    cstrand = ''

    for line in stream:

        line = line.rstrip()

        # print header
        if line[0] == "@":
            if with_header:
                yield line, False
            continue

        # get fields
        sline = line.split('\t')
        query = sline[0]
        flag = int(sline[1])
        strand = flag & 0x10
        ref = sline[2]            
        cigar = sline[5]
        seq = sline[9]
        qual = sline[10]

        # skip empty hits
        if ref == '*' or cigar == '*':
            continue

        # make sure read is mapped
        if flag & 0x4:
            raise ValueError(f"Read is unmapped: {query}")

        # calculate edit distance, total length
        start, end, alen, tlen = decode_cigar(cigar)
        mismatch = int(MISMATCH_RE.search(line).group(1))
        match = alen - mismatch

        # handle SAM asterisks
        if seq == '*' and qual == '*':
            # use last seq/qual
            if cquery != query:
                raise ValueError(f"{cquery=} != {query=} in secondary alignment.")
        else:
            # update seq/qual
            cquery = query
            cseq = seq
            cqual = qual
            cstrand = strand

        # filter by percent identity and minimum length
        if match / tlen < seqid_threshold or tlen < minlen:
            continue

        # always set the seq/qual columns
        if strand == cstrand:
            sline[9] = cseq
            sline[10] = cqual
        else:
            sline[9] = reverse_complement(cseq)
            sline[10] = cqual[::-1]

        # ensure that the cigar matches the sequence
        if tlen != (len(cseq) - start - end):
            raise ValueError(f"Calculated {tlen=} does not match sequence length {len(cseq)} (with {cigar=}).")

        yield sline, True


def main():

    # read command line arguments and/or stdin
    if len(sys.argv[1:]) == 2:
        stream = sys.stdin
        pctid, minlen = sys.argv[1:3]
    else:
        stream = stream_file(sys.argv[1])
        pctid, minlen = sys.argv[2:4]

    pctid, minlen = float(pctid), int(minlen)
    id_treshold = pctid / 100.0

    filtered = filter_alignments(stream, id_treshold, minlen, with_header=True,)

    # two steps as we don't know where the header ends
    #�break after first alignment to avoid the extra checks
    for line, is_alignment in filtered:
        if is_alignment:
            print(*line, sep="\t")
            break
        else:
            print(line)
    
    for line, _ in filtered:
        print(*line, sep="\t")

        

    



if __name__ == "__main__":
    main()