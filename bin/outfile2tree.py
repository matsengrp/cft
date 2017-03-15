#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
	Given an outputfile from the PHYLIP tool `dnaml`, produce and an alignment (including ancestral sequences) and a tree.
"""
from __future__ import print_function

import argparse
import ete3
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node
import csv
import re
import os
from warnings import warn
from collections import defaultdict
import colorbrewer

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# iterate over recognized sections in the phylip output file.
def sections(fh):
    patterns = {
        "parents": "\s+between\s+and\s+length",
        ("sequences", 'dnaml'): "\s*node\s+reconstructed\s+sequence",
        ('sequences', 'dnapars'): "from\s+to\s+any steps"}
    patterns = {k: re.compile(v, re.IGNORECASE) for (k,v) in patterns.items()}
    for line in fh:
        for k, pat in patterns.items():
            if pat.match(line):
                yield k
                break


# iterate over entries in the distance section
def iter_edges(fh):
    #  152          >naive2           0.01208     (     zero,     0.02525) **
    pat = re.compile("\s*(?P<parent>\w+)\s+(?P<child>[\w>_.-]+)\s+(?P<distance>\d+\.\d+)")
    # drop the header underline
    fh.next()
    matches = 0
    for line in fh:
        m = pat.match(line)
        if m:
            matches += 1
            yield (m.group("child"), m.group("parent"), float(m.group("distance")))
        # We only want to break at the very end of the block of matches; dnaml has an extra blank line between
        # header and rows that dnapars doesn't
        elif not matches:
            continue
        else:
            break


# iterate over entries in the sequences section
def parse_seqdict(fh, mode='dnaml'):
    #  152        sssssssssG AGGTGCAGCT GTTGGAGTCT GGGGGAGGCT TGGTACAGCC TGGGGGGTCC
    seqs = defaultdict(str)
    patterns = {
        'dnaml': re.compile("^\s*(?P<id>[a-zA-Z0-9>_.-]*)\s+(?P<seq>[a-zA-Z \-]+)"),
        'dnapars': re.compile("^\s*\S+\s+(?P<id>[a-zA-Z0-9>_.-]*)\s+(yes\s+|no\s+)?(?P<seq>[a-zA-Z \-]+)")}
    fh.next()
    for line in fh:
        m = patterns[mode].match(line)
        if m:
            seqs[m.group("id")] += m.group("seq").replace(" ", "")
        elif line.rstrip() == '':
            continue
        else:
            break
    return seqs


# parse the dnaml output file and return data strictures containing a
# list biopython.SeqRecords and a dict containing adjacency
# relationships and distances between nodes.
def parse_outfile(outfile, seqmeta):
    '''parse phylip outfile'''
    name_map = {seq_id[0:10]: seq_id for seq_id in seqmeta}
    if len(name_map) != len(seqmeta):
        raise RuntimeError("Conflicting shortened names!")
    def full_name(x):
        return name_map.get(x, x)
    with open(outfile, 'rU') as fh:
        for sect in sections(fh):
            if sect == 'parents':
                parents = { full_name(parent): (full_name(child), dist) for parent, child, dist in iter_edges(fh) }
            elif sect[0] == 'sequences':
                d = parse_seqdict(fh, sect[1])
                sequences = [ SeqRecord(Seq(seq), id=full_name(seq_id), description="")
                                for seq_id, seq in d.items() ]
            else:
                raise RuntimeError("unrecognized phylip setion = {}".format(sect))
    # sanity check;  a valid tree should have exactly one node that is parentless
    if not len(parents) == len(sequences) - 1:
        raise RuntimeError('invalid results attempting to parse {}: there are {} parentless sequences'.format(outfile, len(sequences) - len(parents)))

    return sequences, parents


# build a tree from a set of sequences and an adjacency dict.
def build_tree(sequences, parents):
    # build an ete tree
    # first a dictionary of disconnected nodes
    def mknode(r):
        n = Tree(
            name=r.id,
            dist=parents[r.id][1] if r.id in parents else 0,
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
        nstyle['size'] = 7
        node.set_style(nstyle)
    return tree


def timepoint_colors(annotations):
    # Should really improve this sort so that we're fully chronological
    timepoints = sorted(set(x['timepoint'] for x in annotations.values()))
    colors = ['#%02x%02x%02x' % c for c in colorbrewer.RdYlBu[max(len(timepoints), 3)]]
    # note:      ^^^ this bit    converts into a hex
    return dict(zip(timepoints, colors))


def timepoint_legend(ts, tp_colors):
    for timepoint, color in tp_colors.iteritems():
        ts.legend.add_face(ete3.faces.CircleFace(12, color, "circle"), 0)
        ts.legend.add_face(ete3.faces.TextFace(timepoint), 1)


def render_tree(fname, tree, annotations, highlight_node):
    "render tree SVG"
    ts = TreeStyle()
    ts.show_leaf_name = False
    tp_colors = timepoint_colors(annotations)

    def my_layout(node):
        name = node.name
        seqmeta = annotations.get(node.name)
        if seqmeta and not re.compile(".*naive.*").match(node.name):
            # Add the node name for tips
            name = node.name + " (mf={}) ".format(round(float(seqmeta['mut_freqs']), 3))
            F = TextFace(name)
            add_face_to_node(F, node, column=0, position='branch-right')
            # Style the node with color corresponding to timepoint
            nstyle = NodeStyle()
            nstyle['fgcolor'] = tp_colors[seqmeta['timepoint']]
            nstyle['size'] = 15
            node.set_style(nstyle)

    ts.layout_fn = my_layout
    timepoint_legend(ts, tp_colors)
    # whether or not we had rerooted on naive before, we want to do so for the SVG tree
    if 'naive' not in tree.name:
        tree = reroot_tree(tree, '.*naive.*')
    # highlight lineage from seed to root.
    tree = highlight_lineage(tree, highlight_node+'.*')
    tree.render(fname, tree_style=ts)



def get_args():
    def seqmeta_input(fname):
        with open(fname, 'r') as fhandle:
            return dict((row['unique_ids'], row) for row in csv.DictReader(fhandle))

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'phylip_outfile', type=existing_file,
        help='dnaml outfile (verbose output with inferred ancestral sequences, option 5).')
    parser.add_argument(
        'seqmeta', type=seqmeta_input,
        help="Sequence metadata file for annotating mut_freqs")
    parser.add_argument(
        '--outdir', default='.',
        help="output directory where results are left.  [default '%(default)s']")
    parser.add_argument(
        '--basename', help="basename of output files.   [default 'basename(DNAML)']")
    parser.add_argument(
        '--seed', type=str, help="id of leaf [default 'seed']", default='seed')
    return parser.parse_args()


def main():
    args = get_args()

    # basename of outputfiles is taken from --basename if supplied, otherwise basename of DNAML.
    basename = args.basename if args.basename else os.path.basename(args.phylip_outfile)
    outbase = os.path.join(args.outdir, os.path.splitext(basename)[0])

    sequences, parents = parse_outfile(args.phylip_outfile, args.seqmeta)

    if not sequences or not parents:
        raise RuntimeError("No sequences were available; are you sure this is a dnaml output file?")

    # write alignment to fasta
    fname = outbase + '.fa'
    with open(fname, "w") as fh:
        SeqIO.write(sequences, fh, "fasta")

    tree = build_tree(sequences, parents)
    tree = reroot_tree(tree, 'naive.*')

    # write newick file
    fname = outbase + '.newick'
    tree.write(format=1, format_root_node=True, outfile=fname)

    render_tree(outbase+'.svg', tree, args.seqmeta, args.seed)


if __name__ == "__main__":
    main()

