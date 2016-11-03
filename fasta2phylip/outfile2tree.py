#! /bin/env python

import argparse

def main():
    parser = argparse.ArgumentParser(description='give me a phylip dnaml outfile and I''ll give you ancestral sequences and a tree')
    parser.add_argument('--outfile', type=str, default='outfile', help='dnaml outfile (verbose output with inferred ancestral sequences, option 5)')
    args = parser.parse_args()

    # parse all sequences from phylip outfile
    outfiledat = [block.split('\n\n\n')[0].split('\n\n')[1:] for block in open(args.outfile, 'r').read().split('Probable sequences at interior nodes:\n')[1:]]

    # only one tree in outfile
    assert len(outfiledat) == 1

    tree_sequence_dict = {}
    for j, block in enumerate(outfiledat[0]):
        if j == 0:
            for line in block.split('\n'):
                fields = line.split()
                if len(fields) == 0:
                    continue
                name = fields[0]
                seq = ''.join(fields[1:])
                tree_sequence_dict[name] = seq
        else:
            for line in block.split('\n'):
                fields = line.split()
                if len(fields) == 0: continue
                name = fields[0]
                seq = ''.join(fields[1:])
                tree_sequence_dict[name] += seq

    for name in tree_sequence_dict:
        print '> '+name
        print tree_sequence_dict[name]

if __name__ == "__main__":
    main()

