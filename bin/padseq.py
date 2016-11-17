#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Some collections of sequences coming out of `paris` have sequences of different lengths. 
In order to make them look like a multiple sequence alignment, we padd the sequences on the left with 'N's.

This is based on similar functionality inplemented with awk by Will 

"""

import os
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


# given a single SeqIO.SeqRecord, extend the sequence on the left with 'N' up to `maxlen`
# The supplied `rec` is modified in-place, and returned
def pad_inplace(rec, maxlen):
    need = maxlen - len(rec.seq)
    if need < 0:
        raise RuntimeError(
                      "Can only extend sequences - not truncate.\n" + 
                      ">{} current len: {}   desired len: {}".format(rec.id, len(rec.seq), maxlen))
    if need > 0:
        rec.seq = Seq("N"*need + str(rec.seq), rec.seq.alphabet)
    return rec

# pad left of shorter sequences with 'N's
def pad_fasta(file):
    records = list(SeqIO.parse(file, "fasta"))
    # the longest sequence length
    maxlen = max([len(r.seq) for r in records])
    # left pad shorter sequences with Ns
    return [pad_inplace(r, maxlen) for r in records]

def main():
    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('fasta', nargs=1, help='FASTA input', type=existing_file)
    a = parser.parse_args()
    
    for file in a.fasta:
        SeqIO.write(pad_fasta(file), sys.stdout, "fasta")

if __name__ == "__main__":
   main()
