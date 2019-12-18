#!/usr/bin/env python

import argparse
import csv
from Bio import SeqIO


def inseqs_arg(filename):
    with open(filename, "r") as fh:
        return list(SeqIO.parse(fh, "fasta"))


def csargs(x):
    return set(x.split(","))


def get_args():
    parser = argparse.ArgumentParser(
        description="""Safely translates sequences into phylip format, renaming sequences to not exceed 10 chars"""
    )
    parser.add_argument("inseqs", type=inseqs_arg)
    parser.add_argument("--inferred-naive-name", default="inferred_naive")
    parser.add_argument("outseqs", type=argparse.FileType("w"))
    parser.add_argument("seqname_mapping", type=argparse.FileType("w"))
    return parser.parse_args()


def main():
    args = get_args()
    if len(args.inseqs) < 2:
        raise Exception(
            "too few sequences (%d) passed to make_phylip.py (need at least two): %s"
            % (len(args.inseqs), [sr.id for sr in args.inseqs])
        )
    mapping = []
    i = 0
    # Go through and change each sequence id to something short and unique, keeping a mapping of the changes
    for seqrecord in args.inseqs:
        # Here we hardcod `naive` as the translation for the inferred_naive_name, and this is now assumed in `bin/mkconfig.py`
        new_id = (
            "naive" if seqrecord.id == args.inferred_naive_name else "seq-{}".format(i)
        )
        mapping.append({"original_id": seqrecord.id, "new_id": new_id})
        seqrecord.id, seqrecord.name = new_id, new_id
        i += 1
    # write output as phylip and a csv
    SeqIO.write(args.inseqs, args.outseqs, "phylip")
    args.outseqs.close()
    writer = csv.DictWriter(args.seqname_mapping, fieldnames=["original_id", "new_id"])
    writer.writeheader()
    writer.writerows(mapping)
    args.seqname_mapping.close()


if __name__ == "__main__":
    main()
