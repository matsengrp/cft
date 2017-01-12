#!/usr/bin/env python

import argparse
from Bio import SeqIO, SeqRecord
import json


def get_frame(metadata_handle):
    metadata = json.load(metadata_handle)
    cdr3_start = int(metadata[0]["cdr3_start"])
    return cdr3_start % 3


def translate_seqrecord(seqrecord, frame):
    trans_seq = seqrecord.seq[frame:].translate()
    trans_seqrecord = SeqRecord.SeqRecord(trans_seq, id=seqrecord.id, name=seqrecord.name, description=seqrecord.description)
    return trans_seqrecord

def translate_seqrecords(seqrecords, frame):
    return (translate_seqrecord(sr, frame) for sr in seqrecords)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata', type=argparse.FileType('r'))
    parser.add_argument('inseqs', type=lambda x: SeqIO.parse(x, "fasta"))
    parser.add_argument('outseqs', type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    frame = get_frame(args.metadata)
    translated_seqs = translate_seqrecords(args.inseqs, frame)
    SeqIO.write(translated_seqs, args.outseqs, "fasta")
    args.outseqs.close()

if __name__ == '__main__':
    main()


