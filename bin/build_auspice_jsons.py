#!/usr/bin/env python

import argparse
import functools as fun
import copy
import numpy as np
import json
import csv
from Bio import SeqIO, Phylo


# Overall strategy
# ================

# We're going to frame this as annotating the tree with data in one pass, then translating those anntoations
# into dict representations of the final json data in a second pass.


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


def seq_muts(naive, seq):
    base_pairs = zip(range(len(naive)), naive, seq)
    mut_pairs = filter(lambda x: x[1] != x[2], base_pairs)
    return [x[1] + str(x[0]) + x[2] for x in mut_pairs]

def annotate_muts(seqs, tree, attr='muts'):
    naive = seqs['naive0'].seq
    tree = copy.deepcopy(tree)
    for node in tree.find_clades(order="preorder"):
        node_seq = seqs[node.name].seq
        node.__dict__[attr] = seq_muts(naive, node_seq)
    return tree

def annotate_aa_muts(seqs, tree):
    seqs = copy.deepcopy(seqs)
    for _, seqrecord in seqs.iteritems():
        seqrecord.seq = seqrecord.seq.translate()
    return annotate_muts(seqs, tree, attr='aa_muts')

def annotate_timepoints(seqmeta, tree):
    tree = copy.deepcopy(tree)
    for node in tree.find_clades(order='preorder'):
        if node.name:
            node.timepoint = seqmeta[node.name]['timepoint']
    return tree


def annotated_tree(tree, seqs, seqmeta):
    annotate = compose(
            with_layout,
            fun.partial(annotate_muts, seqs),
            fun.partial(annotate_aa_muts, seqs),
            fun.partial(annotate_timepoints, seqmeta))
    return annotate(tree)



# Constructing dictionary representations
# =======================================

# Once we've annotated the tree nodes with all of the relevant information, most of our work in creating the
# representation is just indicating which attributes we want to map to the final dict represetnation (and how)

# These things get mapped directly
repr_attrs = ['xvalue', 'yvalue', 'tvalue', 'muts', 'aa_muts']
# These things get renamed
repr_attrs_trans = {'strain': 'name',
                    'clade': 'name'}
# These things get mapped and renamed under the "attr" attr
repr_attr_attrs = {'div': 'xvalue',
                   'num_date': 'tvalue'}


def clade_repr(clade):
    rep = dict()
    clade_dict = clade.__dict__
    for x, y in repr_attrs_trans.iteritems():
        rep[x] = clade_dict[y]
    for x in repr_attrs:
        rep[x] = clade_dict[x]
    rep['attr'] = {}
    for x, y in repr_attr_attrs:
        rep['attr'][x] = clade_dict.get(y)
    if clade.clades:
        rep['children'] = map(clade_rep, clade.clades)
    # haven't remapped date yet, but not sure if I need to
    #rep['attr'] = {
                   #'date': str(rep['tvalue']) + '-02-13',
                   #}
    return rep


def tree_repr(tree, seqs, seqmeta):
    "Annotates and translates the tree into a dict representation using seqs and seqmeta"
    tree = annotated_tree(tree, seqs, seqmeta)
    return clade_repr(tree.root, seqs, seqmeta)



# Sequence representation
# =======================

def seq_repr(seqrecord):
    return {'nuc': str(seqrecord.seq)}

def seqs_repr(seqrecords):
    return {id: seq_repr(seqrecord) for id, seqrecord in seqrecords.iteritems()}



# Set up argument parser
# ======================

def tree_arg(filename):
    return Phylo.parse(filename, "newick").next()

def seqs_arg(filename):
    return SeqIO.to_dict(SeqIO.parse(filename, "fasta"))

def csv_arg(filename):
    with open(filename) as fh:
        return {row['sequence']: row for row in csv.DictReader(fh)}


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_tree', type=tree_arg)
    parser.add_argument('input_seqs', type=seqs_arg)
    parser.add_argument('input_seqmeta', type=csv_arg)
    parser.add_argument('json_tree', type=argparse.FileType('w'))
    parser.add_argument('json_seqs', type=argparse.FileType('w'))
    return parser.parse_args()



# Main function tying everything together
# =======================================

def main():
    args = get_args()
    json.dump(tree_repr(args.input_tree, args.input_seqs, args.input_seqmeta), args.json_tree)
    json.dump(seqs_repr(args.input_seqs), args.json_seqs)
    args.json_tree.close()
    args.json_seqs.close()


if __name__ == '__main__':
    main()


