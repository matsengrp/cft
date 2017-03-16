#!/usr/bin/env python

import argparse
import copy
import numpy as np
import json
from Bio import SeqIO, Phylo


# Functional helpers
# ------------------

def compose(*functions):
    def compose2(f, g):
        return lambda x: f(g(x))
    return reduce(compose2, functions, lambda x: x)



# Annotating trees
# ================

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents

def with_layout(tree):
    """Returns a copy of tree with clade, xvalue, and yvalue attributes on all nodes"""
    tree = copy.deepcopy(tree)
    parents = all_parents(tree)
    clade_n = 0
    yvalue = tree.count_terminals()
    for node in tree.find_clades(order="preorder"):
        node.clade = clade_n
        clade_n += 1
        parent = parents.get(node)
        if parent is not None:
            node.xvalue = parent.xvalue + node.branch_length
        else:
            node.xvalue = 0
            #node.tvalue = 0
        # This is a terrible terrible hack to get some more interesting data to look at; really should just add timepoint
        #node.tvalue = 0
        node.tvalue = int(2000 + node.xvalue * 100)
        # Assing yvalues for terminal nodes
        if node.is_terminal():
            node.yvalue = yvalue
            yvalue -= 1

    # Clever way of assinging yvalues to the internal nodes
    for node in tree.get_nonterminals(order="postorder"):
        node.yvalue = np.mean([x.yvalue for x in node.clades])

    return tree


def annotated_tree(tree):
    annotate = compose(with_layout)
    return annotate(tree)



# Constructing dictionary representations
# =======================================

# Most of our work here is in just constructing dictionary representations of our tree and sequences.

repr_attrs = ['xvalue', 'yvalue', 'tvalue']
repr_attr_trans = {'strain': 'name'}


def clade_repr(clade):
    rep = dict()
    clade_dict = clade.__dict__
    for x, y in repr_attr_trans.iteritems():
        rep[x] = clade_dict[y]
    for x in repr_attrs:
        rep[x] = clade_dict[x]
    if clade.clades:
        rep['children'] = map(clade_repr, clade.clades)
    return rep


def tree_repr(tree):
    tree = annotated_tree(tree)
    return clade_repr(tree.root)



# Set up argument parser
# ======================

def tree_arg(filename):
    return Phylo.parse(filename, "newick").next()

def seqs_arg(filename):
    return SeqIO.parse(filename, "fasta")

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_tree', type=tree_arg)
    parser.add_argument('input_seqs', type=seqs_arg)
    parser.add_argument('json_tree', type=argparse.FileType('w'))
    parser.add_argument('json_seqs', type=argparse.FileType('w'))
    return parser.parse_args()



# Main function tying everything together
# =======================================

def main():
    args = get_args()
    json.dump(tree_repr(args.input_tree), args.json_tree)
    args.json_tree.close()
    args.json_seqs.close()


if __name__ == '__main__':
    main()


