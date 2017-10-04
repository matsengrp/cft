#!/usr/bin/env python

import argparse
import csv
from Bio import SeqIO


def inseqs_arg(filename):
    with open(filename, 'r') as fh:
        return list(SeqIO.parse(fh, 'fasta'))

def csargs(x):
    return set(x.split(','))


def get_args():
    parser = argparse.ArgumentParser(
            description="""Safely translates sequences into phylip format, renaming sequences to not exceed 10 chars""")
    parser.add_argument('inseqs', type=inseqs_arg)
    parser.add_argument('--dont-rename', type=csargs,
            help="Comma separated list of names to not shorten/translate")
    parser.add_argument('outseqs', type=argparse.FileType('w'))
    parser.add_argument('seqname_mapping', type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    mapping = []
    i = 0
    # Go through and change each sequence id to something short and unique
    for seqrecord in args.inseqs:
        if seqrecord.id in args.dont_rename:
            mapping.append({'original_id': seqrecord.id, 'new_id': seqrecord.id})
        else:
            new_id = "seq-{}".format(i)
            i += 1
            # keep a mapping of the changes
            mapping.append({'original_id': seqrecord.id, 'new_id': new_id})
            seqrecord.id, seqrecord.name = new_id, new_id
    # write output as phylip and a csv
    SeqIO.write(args.inseqs, args.outseqs, 'phylip')
    args.outseqs.close()
    writer = csv.DictWriter(args.seqname_mapping, fieldnames = ['original_id', 'new_id'])
    writer.writeheader()
    writer.writerows(mapping)
    args.seqname_mapping.close()



if __name__ == '__main__':
    main()



