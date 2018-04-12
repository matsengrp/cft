#!/usr/bin/env python

import argparse
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqs')
    return parser.parse_args()


def main():
    args = get_args()
    seqs = SeqIO.parse(args.seqs, 'fasta')
    # Can't exceed 32000 for memory without having to request a large node
    n_seqs = 0
    for seq in seqs:
        n_seqs += 1
    if n_seqs > 8000:
        print '--exclusive ',
    # add baseline of 1/2 a gig
    mem_needed = 500 + int(n_seqs * 1.6)
    if mem_needed > 32000:
        print '--partition=largenode ',
    print '--mem={} '.format(mem_needed),


if __name__ == '__main__':
    main()
    
