#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
prune a newick tree based on distnace from root-->seed lineage
print taxa from pruned tree
"""

from ete3 import Tree
from outfile2tree import find_node, reroot_tree

import subprocess
import tempfile
import argparse


def lineage_selection(args):
    tree = args.tree
    naive_node = find_node(tree, args.naive)
    tree.set_outgroup(naive_node)
    tree = reroot_tree(tree, args.naive)

    # ids of nodes on the seed lineage
    seed_node = find_node(tree, args.seed)
    if args.seed is 'seed' and seed_node is None:
        seed_node = tree.get_farthest_leaf()[0]
    seed_root_lineage = set(seed_node.get_ancestors()+[seed_node])
    # iterate over taxa and find n closest to seed lineage
    leaf_distances = []
    for leaf in tree.iter_leaves():
        distance = 0
        node = leaf
        while node not in seed_root_lineage:
            distance += node.dist
            node = node.up
        leaf_distances.append((leaf.name, distance))

    # gotta print this, since it's not a leaf in the rerooted tree
    yield naive_node.name
    for leaf, distance in sorted(leaf_distances, key=lambda x: x[1])[:(args.n_keep + 1)]:
        yield leaf


def with_temporary_handle(lines):
    """Given some lines of text to be placed in a temporary file, returns a decorator function that immediately
    calls the decorated function with as sole argument a temporary file containing the given lines of text.
    Closes the named temporary file when done."""
    def deco(f):
        handle = tempfile.NamedTemporaryFile("w+")
        handle.writelines([line + "\n" for line in lines])
        handle.flush()
        f(handle)
        handle.close()
        return f
    return deco


def min_adcl_selection(args):
    tipnames = [t.name for t in args.tree.iter_leaves()]
    if len(tipnames) <= args.n_keep:
        return tipnames
    else:
        results = []
        @with_temporary_handle([args.seed, args.naive])
        def always_include(ai_handle):
            command = "rppr min_adcl_tree --algorithm pam --leaves".split(" ") \
                + [args.n_keep, "--always-include", ai_handle.name, args.tree_file]
            command = map(str, command)
            results.append(subprocess.check_output(command))
        cut_names = [n for n in results[0].split("\n")]
        return (x for x in tipnames if x not in cut_names)


def tree_arg(tree_arg_value):
    return Tree(tree_arg_value, format=1)

def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'tree_file',
        help='input newick tree file')
    parser.add_argument(
        '--strategy', type=str, choices=["min_adcl", "seed_lineage"], default="seed_lineage")
    parser.add_argument(
        '--naive', type=str, help='id of root [default \'naive0\']', default='naive0')
    parser.add_argument(
        '--seed', type=str, help='id of leaf [default \'seed\']', default='seed')
    parser.add_argument(
        '-n', '--n-keep', type=int,
        help='number of sequences to keep [default: 100]', default=100)
    args = parser.parse_args()
    args.tree = tree_arg(args.tree_file)
    return args

def main(args):
    # read newick file into ete tree, set naive as outgroup, and reroot
    selection_fn = lineage_selection if args.strategy == "seed_lineage" else min_adcl_selection
    for name in selection_fn(args):
        # Writes to stdout
        print name
    

if __name__ == "__main__":
    main(get_args())

