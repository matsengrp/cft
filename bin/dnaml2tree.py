#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
	Given an outputfile from the PHYLIP tool `dnaml`, produce and an alignment (including ancestral sequences) and a tree.
"""
from __future__ import print_function

import argparse
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node
import sys
import re
import os
from warnings import warn
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


# iterate over recognized sections in the phylip output file.
def sections(fh):
    patterns = {
        "sequences": "\s+node\s+Reconstructed\s+sequence",
        "parents": "\s+Between\s+And\s+Length\s+Approx.\s+Confidence Limits"
    }
    patterns = {k: re.compile(v) for (k,v) in patterns.items()}
    for line in fh:
        for k, pat in patterns.items():
            if pat.match(line):
                yield k
                break

            
# iterate over entries in the sequences section
def parse_seqdict(fh):
    #  152        sssssssssG AGGTGCAGCT GTTGGAGTCT GGGGGAGGCT TGGTACAGCC TGGGGGGTCC
    seqs = defaultdict(str)
    pat = re.compile("^\s*(?P<id>[a-zA-Z0-9>_.-]*)\s+(?P<seq>[a-zA-Z ]+)")
    fh.next()
    for line in fh:
        m = pat.match(line)
        if m:
            seqs[m.group("id")] += m.group("seq").replace(" ", "")
        elif line.rstrip() == '':
            continue
        else:
            break
    return seqs


# iterate over entries in the distance section
def iter_parents(fh):
    parents = {}
    #  152          >naive2           0.01208     (     zero,     0.02525) **
    pat = re.compile("^\s*(?P<parent>[0-9]+)\s+(?P<child>[a-zA-z0-9>_.-]+)\s+(?P<distance>[0-9]+(\.[0-9]+))")
    fh.next()
    fh.next()
    for line in fh:
        # print("\"{}\"".format(line))
        m = pat.match(line)
        #print(line)
        #print(m)
        if m:
            yield (m.group("child"), (m.group("parent"), float(m.group("distance"))))
        else:
            break


# parse the dnaml output file and return data strictures containing a
# list biopython.SeqRecords and a dict containing adjacency
# relationships and distances between nodes.
def outfile2seqs(outfile='outfile'):
    sequences = []
    parents = {}
    with open(outfile, 'rU') as fh:
        for sect in sections(fh):
            if sect == 'sequences':
                d = parse_seqdict(fh)
                sequences = [ SeqRecord(Seq(v), id=k, description="") for (k,v) in d.items()]
            elif sect == 'parents':
                parents = { k: v for k, v in iter_parents(fh) }
            else:
                raise RuntimeError("unrecognized phylip setion = {}".format(sect))
    # a necessary, but not sufficient, condition for this to be a valid tree is
    # that exactly one node is parentless
    if not len(parents) == len(sequences) - 1:
        #print(len(parents), len(sequences))
        print(parents)
        raise RuntimeError('invalid results attempting to parse {}: there are {} parentless sequences'.format(outfile, len(sequences) - len(parents)))
    return sequences, parents


# build a tree from a set of sequences and an adjacency dict.
def build_tree(sequences, parents):
    # build an ete tree
    # first a dictionary of disconnected nodes
    def mknode(r):
        n = Tree(
            name=r.id,
            dist=parents[r.id][1] if r.id in parents else None,
            format=1)
        n.add_features(seq=r)
        return n
    
    nodes = { r.id: mknode(r) for r in sequences }

    # connect the nodes using the parent data
    orphan_nodes = 0
    for r in sequences:
        if r.id in parents:
            nodes[parents[r.id][0]].add_child(nodes[r.id])
        else:
            # node without parent becomes root
            tree = nodes[r.id]
            orphan_nodes += 1
    # there can only be one root
    if orphan_nodes != 1:
        raise RuntimeError("The tree is not properly rooted; expected a single root but there are {}.".format(orphan_nodes))
    return tree

def find_node(tree, pattern):
    regex = re.compile(pattern).search
    nodes =  [ node for node in tree.traverse() for m in [regex(node.name)] if m]
    if not nodes:
        warn("Cannot find matching node; looking for name matching '{}'".format(pattern))
        return
    else:
        if len(nodes) > 1:
            warn("multiple nodes found; using first one.\nfound: {}".format([n.name for n in nodes]))
        return nodes[0]


# reroot the tree on node matching regex pattern.
# Usually this is used to root on the naive germline sequence with a name matching '.*naive.*'
def reroot_tree(tree, pattern='.*naive.*'):
    # find all nodes matching pattern
    node = find_node(tree, pattern)
    if tree != node:
        tree.remove_child(node)
        node.add_child(tree)
        tree.dist = node.dist
        node.dist = 0
        tree = node
    return tree


# iterate up a lineage toward the root.
# lineage starts at node whose name matches pattern.
def iter_lineage(tree, pattern):
    node = find_node(tree, pattern)
    while node:
        yield node
        node = node.up

        
# highlight the lineage (branches) from a node to the root.
# target node name is identified by a regular expression.
# The default pattern is chosen to highlight lineage from a seed node to the root.
def highlight_lineage(tree, pattern='seed.*'):
    # find all nodes matching pattern 'seed.*'
    for node in iter_lineage(tree, pattern):
        nstyle = NodeStyle()
        nstyle['fgcolor'] = 'red'
        nstyle['size'] = 10
        node.set_style(nstyle)
    return tree

        
def main():
    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--dnaml', type=existing_file, default='outfile',
        help='dnaml outfile (verbose output with inferred ancestral sequences, option 5). [ default \'%(default)s\' ]')
    parser.add_argument(
        '--outdir', default='.',
        help='output directory where results are left.  [ default \'%(default)s\' ]')
    parser.add_argument(
        '--basename', help='basename of output files.   [ default \'basename(DNAML)\' ]')
    args = parser.parse_args()

    # basename of outputfiles is taken from --basename if supplied, otherwise basename of DNAML.
    basename = args.basename if args.basename else os.path.basename(args.dnaml)
    outbase = os.path.join(args.outdir, os.path.splitext(basename)[0])

    sequences, parents = outfile2seqs(args.dnaml)

    if not sequences or not parents:
        raise RuntimeError("No sequences were available; are you sure this is a dnaml output file?")

    # write alignment to fasta
    fname = outbase + '.fa'
    with open(fname, "w") as fh:
        SeqIO.write(sequences, fh, "fasta")

    tree = build_tree(sequences, parents)
    
    # reroot on germline outgroup, if available.
    # [csw] don't reroot for now while we experiment with ascii-art trees.
    # [wsd] ^ note that this breaks lineage iteration (won't include naive), issue #73
    # tree = reroot_tree(tree, '.*naive.*')

    # highlight lineage from seed to root.
    highlight_lineage(tree, 'seed.*')

    # write sequences along lineage from seed to root.
    fname = outbase + '.seedLineage.fa'
    with open(fname, 'w') as f:
        seq2write = [node.seq for node in iter_lineage(tree, 'seed.*')]
        # if we didn't reroot on naive, we need a special line of code to also print that one
        if 'naive' not in tree.name:
            seq2write.append(find_node(tree, '.*naive.*').seq)
        SeqIO.write(seq2write, f, 'fasta')

    # write newick file
    fname = outbase + '.newick'
    tree.write(format=1, format_root_node=True, outfile=fname)

    # render tree SVG
    ts = TreeStyle()
    ts.show_leaf_name = False

    def my_layout(node):
        F = TextFace(node.name)
        add_face_to_node(F, node, column=0, position='branch-right')

    ts.layout_fn = my_layout
    fname = outbase + '.svg'
    tree.render(fname, tree_style=ts)


if __name__ == "__main__":
    main()
