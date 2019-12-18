#!/usr/bin/env python

import argparse
from Bio import Phylo


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("tree")
    return parser.parse_args()


def main():
    args = get_args()
    tree = Phylo.parse(args.tree, "newick").next()
    # Can't exceed 32000 for memory without having to request a large node
    n_terminals = len(tree.get_terminals())
    if n_terminals > 4000:
        print "--exclusive ",
    mem_needed = n_terminals * 5
    if mem_needed > 32000:
        print "--partition=largenode ",
    print "--mem={} ".format(mem_needed),


if __name__ == "__main__":
    main()
