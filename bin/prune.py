#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Return a reduced set of taxa based on either seed lineage selection ("pruning")
or min ADCL selection ("trimming").
"""

from ete3 import Tree
from process_asr import find_node, reroot_tree

import subprocess
import tempfile
import argparse
import copy


def seed_lineage_selection(args):
    """
    Take the closest taxa on subtrees branching off the root (naive) to seed lineage.
    E.g. given (((seed:1,(l3:1,l4:2)s2:1,(l5:3,l6:4)s3:1)b:1,(l1:1,l2:2)s1:1),naive0:1);

	     /-seed
	    |
	    |   /-l3
	  /-|--|
	 |  |   \-l4
	 |  |
	 |  |   /-l5
       /-|   \-|
      |  |      \-l6
      |  |
    --|  |   /-l1
      |   \-|
      |      \-l2
      |
       \-naive0

    return (naive0, l3, l1, l5, seed)

    We cycle through these subtrees from root to seed, always taking the
    closest taxon to the root to seed lineage, until we have gotten the required
    number of taxa.
    """
    tree = args.tree
    naive_node = find_node(tree, args.naive)
    tree.set_outgroup(naive_node)
    tree = reroot_tree(tree, args.naive)

    # Collect ids of nodes on the seed lineage.
    seed_node = find_node(tree, args.seed)
    if args.seed is 'seed' and seed_node is None:
        seed_node = tree.get_farthest_leaf()[0]

    # Set args.n_keep to the number of leaves minus 1 (the naive counts as a
    # leaf) if args.n_keep is smaller than this quantity.
    args.n_keep = min(args.n_keep, len(tree) - 1)

    distance_dict = {}
    def distances(lineage_node, leaf):
        if leaf not in distance_dict:
            distance_dict[leaf] = lineage_node.get_distance(leaf)
        return distance_dict[leaf]

    # Iterate over seed lineage and find closest taxon from each branch, and
    # repeat until we have args.n_keep leaf sequences.
    leaves_to_keep = set(find_node(tree, n) for n in args.always_include if tree.get_leaves_by_name(n))
    # Extract the sequence of nodes on lineage from root to seed.
    # Note: ete doc suggests get_ancestors() would include the seed node, but it doesn't.
    #       "Returns the list of all ancestor nodes from current node to the current tree root"
    # Note this slicing (::-1) reverses the order, so we do indeed go from root to seed.
    seed_lineage = seed_node.get_ancestors()[::-1] + [seed_node]
    # We'll build a list of the subtrees off the seed lineage.
    subtrees = []
    for i, lineage_node in enumerate(seed_lineage[:-1]):
        for subtree in lineage_node.children:
            # The subtree that continues down the seed lineage doesn't count.
            if subtree != seed_lineage[i+1]:
                subtrees.append(subtree)
    # Repeatedly pass through this list of subtrees, grabbing the one
    # closest leaf (to the seed lineage) from each, until we get how many we
    # need.
    while len(leaves_to_keep) < (args.n_keep - 1): # -1 because we naive gets rerooted out, and have to manually yield as below
        for subtree in subtrees:
            # Obtain all the leaves in this subtree that aren't already in leaves_to_keep.
            leaves = [leaf for leaf in subtree.iter_leaves() if leaf not in leaves_to_keep]
            if leaves:
                # Add the leaf that's the closest to the seed lineages.
                # Use subtree.up so that we are getting the distance from the
                # node on the lineage from root to seed (rather than the root
                # of the subtree).
                leaves_to_keep.add(min(leaves, key=lambda leaf: distances(subtree.up, leaf)))
                if len(leaves_to_keep) == args.n_keep - 1:
                    break

    # Yield all the selected leaves (including naive and seed)
    yield args.naive
    for leaf in leaves_to_keep:
        yield leaf.name


def with_temporary_handle(lines):
    """
    Given some lines of text to be placed in a temporary file, returns a decorator function that immediately
    calls the decorated function with as sole argument a temporary file containing the given lines of text.
    Closes the named temporary file when done.
    """
    def deco(f):
        handle = tempfile.NamedTemporaryFile("w+")
        handle.writelines([line + "\n" for line in lines])
        handle.flush()
        f(handle)
        handle.close()
        return f
    return deco


def min_adcl_selection(args):
    """
    Minimize ADCL for a tree using pplacer suite.
    """
    tipnames = [t.name for t in args.tree.iter_leaves()]
    if len(tipnames) <= args.n_keep:
        return tipnames
    else:
        results = []
        @with_temporary_handle(args.always_include)
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
        '--strategy', choices=["min_adcl", "seed_lineage"], default="seed_lineage")
    parser.add_argument(
        '--naive', help='id of root [default \'naive0\']', default='naive0')
    parser.add_argument(
        '--seed', help='id of leaf [default \'seed\']', default='seed')
    parser.add_argument(
        '--always-include', type=lambda x: set(x.split(',')), help='comma separated list of ids to keep',
        default=set())
    parser.add_argument(
        '-n', '--n-keep', type=int,
        help='number of sequences to keep [default: 100]', default=100)
    args = parser.parse_args()
    args.tree = tree_arg(args.tree_file)
    args.always_include = set([args.seed, args.naive] \
            + [leaf_name for leaf_name in args.always_include if leaf_name in args.tree.get_leaf_names()])
    return args


def main(args):
    selection_fn = seed_lineage_selection if args.strategy == "seed_lineage" else min_adcl_selection
    for name in selection_fn(args):
        # Writes to stdout
        print name


if __name__ == "__main__":
    main(get_args())
