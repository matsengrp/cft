#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
prune a newick tree based on distnace from root-->seed lineage
print taxa from pruned tree
"""
import argparse, sys
from ete3 import Tree
from dnaml2tree import find_node, reroot_tree
        
def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'newick', type=str, default='outfile',
        help='input newick tree file')
    parser.add_argument(
        '--naive', type=str, help='id of root [default \'naive\']', default='naive')
    parser.add_argument(
        '--seed', type=str, help='id of leaf [default \'seed\']', default='seed')
    parser.add_argument(
        '--n', type=int,
        help='number of sequences to keep [default: 100]', default=100)
    args = parser.parse_args()

    # read newick file into ete tree, set naive as outgroup, and reroot
    tree = Tree(args.newick, format=1)
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
    print naive_node.name
    for leaf, distance in sorted(leaf_distances, key=lambda x: x[1])[:(args.n+1)]:
        print leaf

if __name__ == "__main__":
    main()
