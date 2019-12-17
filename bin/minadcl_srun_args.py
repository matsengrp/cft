#!/usr/bin/env python

import argparse
from Bio import Phylo

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tree')
    return parser.parse_args()

def main():
    args = get_args()
    tree = Phylo.parse(args.tree, 'newick').next()
    # Can't exceed 32000 for memory without having to request a large node
    n_tips = len(tree.get_terminals())
    if n_tips > 8000:
        print '--exclusive ',
    # add baseline of 1/2 a gig
    mem_needed = 500 + int(n_tips * 1.8) if n_tips < 10000 else 20000 + int((n_tips - 10000) * 6)
    if mem_needed > 32000:
        print '--partition=largenode ',
    print '--mem={} '.format(mem_needed),

if __name__ == '__main__':
    main()
