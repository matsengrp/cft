#!/usr/bin/env python

import argparse
from Bio import SeqIO, SeqRecord, Seq


def string_get(s, i):
    try:
        return s[i]
    except IndexError:
        pass

def fix_alnd_seq(alnd, orig):
    new_seq = ""
    i = 0    # index to position in orig seq
    for x in alnd.seq:
        # We don't increment i if we have a gap, since only inferred in alnd
        if x == '-':
            pass
        # If we have a X in alnd and * in orig, switch char to write to *
        elif x == 'X' and string_get(orig.seq, i) == '*':
            i += 1
            x = '*'
        # Just truck along with index
        else:
            i += 1
        # Add given aa to seq
        new_seq += x
    return SeqRecord.SeqRecord(Seq.Seq(new_seq, Seq.Alphabet.ProteinAlphabet()), id=alnd.id, name=alnd.name, description=alnd.description)


def fix_alignment(alnd_seqrecords, orig_seqrecords):
    for id, alnd_seqrecord in alnd_seqrecords.iteritems():
        orig_seqrecord = orig_seqrecords[id]
        yield fix_alnd_seq(alnd_seqrecord, orig_seqrecord)


def seqrecord_dict_arg(x):
    return SeqIO.to_dict(SeqIO.parse(x, "fasta"))

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('aligned_translation', type=seqrecord_dict_arg)
    parser.add_argument('orig_translation', type=seqrecord_dict_arg)
    parser.add_argument('output', type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    fixed_seqs = fix_alignment(args.aligned_translation, args.orig_translation)
    SeqIO.write(fixed_seqs, args.output, "fasta")
    args.output.close()


if __name__ == '__main__':
    main()


