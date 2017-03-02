#!/usr/bin/env python

import argparse
from Bio import SeqIO, SeqRecord, Seq
import json


def get_frame(metadata_handle):
    metadata = json.load(metadata_handle)
    cdr3_start = int(metadata[0]["cdr3_start"])
    return cdr3_start % 3

def trim_end(seqrecord, frame):
    return ((len(seqrecord.seq) - frame) / 3) * 3 + frame

def translate_seqrecord(seqrecord, frame):
    end = trim_end(seqrecord, frame)
    trans_seq = Seq.Seq(str(seqrecord.seq[frame:end]).replace('-', '')).translate()
    trans_seqrecord = SeqRecord.SeqRecord(trans_seq, id=seqrecord.id, name=seqrecord.name, description=seqrecord.description)
    return trans_seqrecord

def translate_seqrecords(seqrecords, frame):
    return (translate_seqrecord(sr, frame) for sr in seqrecords)

def trim_seqrecord(seqrecord, frame):
    end = trim_end(seqrecord, frame)
    trimmed_seq = seqrecord.seq[frame:end]
    trimmed_seqrecord = SeqRecord.SeqRecord(trimmed_seq, id=seqrecord.id, name=seqrecord.name, description=seqrecord.description)
    return trimmed_seqrecord

def trim_seqrecords(seqrecords, frame):
    return (trim_seqrecord(sr, frame) for sr in seqrecords)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata', type=argparse.FileType('r'))
    parser.add_argument('inseqs', type=lambda x: list(SeqIO.parse(x, "fasta")))
    parser.add_argument('outseqs', type=argparse.FileType('w'))
    parser.add_argument('-t', '--trimmed-inseqs', type=argparse.FileType('w'),
        help="If any inseqs lengths aren't divisible by three, trim excess and write to this file.")
    return parser.parse_args()


def main():
    args = get_args()
    frame = get_frame(args.metadata)
    translated_seqs = translate_seqrecords(args.inseqs, frame)
    SeqIO.write(translated_seqs, args.outseqs, "fasta")
    args.outseqs.close()
    if args.trimmed_inseqs:
        trimmed_inseqs = trim_seqrecords(args.inseqs, frame)
        SeqIO.write(trimmed_inseqs, args.trimmed_inseqs, "fasta")
        args.trimmed_inseqs.close()


if __name__ == '__main__':
    main()


