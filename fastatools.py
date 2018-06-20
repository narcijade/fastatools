#!/usr/bin/env python
# ejr - 20161029
# Rewrite of fastatools in python
# This runs much faster with pypy
# INPUT: FASTA file on STDIN or with --file
# OUTPUT: varies, depending on options.  FASTA, table or histogram

###############################################################################
### IMPORT ####################################################################
###############################################################################
import sys
import argparse
import re

###############################################################################
### MAIN ######################################################################
###############################################################################
def main():

    args = get_args()
    fasta = read_fasta(args.file)

    # output fasta stats
    if args.stats:
        print_fasta_stats(fasta)
    # output histogram of sequence gc contents
    elif args.gc_histogram:
        print_gc_histogram(fasta)
    # output histogram of sequence lengths
    elif args.len_histogram:
        print_length_histogram(fasta)
    # output table of per sequence gc content
    elif args.gc_table:
        print_gc_table(fasta)
    # output table of per sequence length
    elif args.len_table:
        print_lengths_table(fasta)
    # output table of sequences
    elif args.seq_table:
        print_seq_table(fasta)
    # filter FASTA file and print FASTA
    elif args.filter:
        fasta_out = filter_fasta(fasta, args.len_min, args.len_max, args.gc_min, args.gc_max)
        print_fasta(fasta_out)
    else:
        print "Something went wrong, this shouldn't happen"

###############################################################################
### SUBROUTINES ###############################################################
###############################################################################

###############################################################################
# Read in FASTA file from STDIN
###############################################################################
def read_fasta(filename):
    header = ""
    fasta = {}

    for line in filename:
        line = line.rstrip()
        if (line[0] == ">"):
            header = line
            fasta[header] = ''
        else:
            fasta[header] = fasta[header] + line

    return fasta

###############################################################################
# Add newlines every 80 characters 
###############################################################################
def insert_newlines(string, every=80):
    lines = []

    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])

    return '\n'.join(lines)

###############################################################################
# Print FASTA to STDOUT
###############################################################################
def print_fasta(fasta):

    for header in fasta:
        print header
        print insert_newlines(fasta[header])
    
###############################################################################
# Get command-line arguments using argparse
###############################################################################
def get_args():
    parser = argparse.ArgumentParser(description="General Purpose FASTA Manipulation")
    # file defaults to stdin 
    parser.add_argument('--file', type = argparse.FileType('r'), default = sys.stdin, help = 'Input FASTA - defaults to STDIN')
    parser.add_argument('--len_min', type = int, default = '0', help = 'Minimum sequence length to retain')
    parser.add_argument('--len_max', type = int, default = '10000000', help = 'Maximum sequence length to retain')
    parser.add_argument('--gc_min', type = int, default = '0', help = 'Minimum sequence GC percentage to retain (20 = .20)')
    parser.add_argument('--gc_max', type = int, default = '100', help = 'Maximum sequence GC percentage to retain (60 = .60)')
    # only allow one of the following option to be used
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--stats', help = 'Print FASTA Statistics', action="store_true")
    group.add_argument('--gc_histogram', help = 'Output Histogram of per sequence GC content', action="store_true")
    group.add_argument('--len_histogram', help = 'Output Histogram of per sequence length', action="store_true")
    group.add_argument('--gc_table', help = 'Output Table of GC Content', action="store_true")
    group.add_argument('--len_table', help = 'Output Table of Sequence Length', action="store_true")
    group.add_argument('--seq_table', help = 'Output Table of Sequences', action="store_true")
    group.add_argument('--filter', help = 'Filter FASTA file', action="store_true")
    args = parser.parse_args()

    return args

###############################################################################
# Filter FASTA
###############################################################################
def filter_fasta(fasta, len_min , len_max , gc_min, gc_max):
    fasta_out = {}

    for header in fasta:
        seq = fasta[header].upper()
        slen = float(len(seq))
        g = float(seq.count('G'))
        c = float(seq.count('C'))
        gc = (g + c) / slen * 100

        if (slen >= len_min and 
            slen <= len_max and 
            gc >= gc_min and 
            gc <= gc_max):
            fasta_out[header] = fasta[header]
        
    return fasta_out

###############################################################################
# Print out FASTA statistics
###############################################################################
def print_fasta_stats(fasta):
    num_bases = 0
    num_seqs = 0 
    total_gc = 0
    total_at = 0
    total_n = 0
    longest_seq = 0
    shortest_seq = 1000000000
    seq_lengths = []

    for header in fasta:
        seq = fasta[header].upper()
        slen = len(seq)
        seq_lengths.append(slen)
        g = seq.count('G')
        c = seq.count('C')
        a = seq.count('A')
        t = seq.count('T')
        n = seq.count('N')
        gc = g + c
        at = a + t
        # add lengths together for number of bases
        num_bases += slen
        # count number of sequences
        num_seqs += 1
        # track total number of bases of each type
        total_gc += gc
        total_at += at
        total_n += n
        # track longest and shortest sequence
        if (slen < shortest_seq):
            shortest_seq = slen
        if (slen > longest_seq):
            longest_seq = slen

    average_length = float(num_bases) / float(num_seqs)
    percent_gc = float(total_gc) / (float(total_gc) + float(total_at))
    percent_n = float(total_n) / float(num_bases)
    sorted_lengths = sorted(seq_lengths, reverse=True)
    fifty = float(num_bases * .5)
    ninety = float(num_bases *.9)

    running_total = 0
    n50 = 0
    n90 = 0
    for length in sorted_lengths:
        running_total += length
        if (running_total < fifty):
            n50 = length
        if (running_total < ninety):
            n90 = length

    # Average is a floating point number, convert to integer
    average_length = int(round(average_length, 0))

    print """
|                      |            |
|:-------------------- |-----------:|
| Number of Sequences: | {: >10,} |
| Total Length:        | {: >10,} |
| Average Length:      | {: >10,} |
| Longest Sequence:    | {: >10,} |
| Shortest Sequence:   | {: >10,} |
| Percent GC:          | {: >10.0%} |
| Percent N:           | {: >10.0%} |
| N50:                 | {: >10,} |
| N90:                 | {: >10,} |
""".format(
    num_seqs,
    num_bases,
    average_length,
    longest_seq,
    shortest_seq,
    percent_gc,
    percent_n,
    n50,
    n90)

###############################################################################
# Print histogram of sequence lengths: bins of 100, max 10,000
###############################################################################
def print_length_histogram(fasta):
    len_hist = [0] * 101
    print "bin\tnum.seqs"

    for header in fasta:
        seq = fasta[header].upper()
        slen = len(seq)
        if (slen % 100 != 0): 
             slen = slen - slen % 100
        slen = slen / 100
        if (slen > 100):
            slen = 100
        len_hist[slen] += 1

    i = 0
    for count in len_hist:
        print i * 100, "\t", count
        i += 1

###############################################################################
# Print histogram of GC content
###############################################################################
def print_gc_histogram(fasta):
    gc_hist = [0] * 101
    print "bin\tnum.seqs"

    for header in fasta:
        seq = fasta[header].upper()
        slen = float(len(seq))
        g = float(seq.count('G'))
        c = float(seq.count('C'))
        gc = (g + c) / slen * 100
        gc = int(round(gc, 0))
        gc_hist[gc] += 1

    i = 1
    for count in gc_hist[1:]:
        print i, "\t", count
        i += 1

###############################################################################
# Print table of sequences
###############################################################################
def print_seq_table(fasta):
    print "seq.id\tseq"

    for header in fasta:
        seq = fasta[header].upper()
        header = header[1:]
        print header, "\t", seq

###############################################################################
# Print table of sequence lengths
###############################################################################
def print_lengths_table(fasta):
    print "seq.id\tseq.length"

    for header in fasta:
        seq = fasta[header].upper()
        slen = len(seq)
        header = header[1:]
        print header, "\t", slen

###############################################################################
# Print table of GC content
###############################################################################
def print_gc_table(fasta):
    print "seq.id\tseq.gc"

    for header in fasta:
        seq = fasta[header].upper()
        slen = float(len(seq))
        g = float(seq.count('G'))
        c = float(seq.count('C'))
        gc = (g + c) / slen * 100
        header = header[1:]
        print header, "\t%0.2f" % gc

###############################################################################
### RUN MAIN ##################################################################
###############################################################################
if __name__ == "__main__":
    main()
