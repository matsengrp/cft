#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given an outputfile from one of the PHYLIP tools - `dnaml` or `dnapars` - produce an alignment (including
ancestral sequences), a newick tree (with matching internal node lables), and an svg rendering of said tree.
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
        'dnapars': re.compile("^\s*\S+\s+(?P<id>[a-zA-Z0-9>_.-]*)\s+(yes\s+|no\s+|maybe\s+)?(?P<seq>[a-zA-Z \-]+)")}
    fh.next()
    for line in fh:
        m = patterns[mode].match(line)
        if m and m.group("id") is not '':
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
    trees = []
    # Ugg... for compilation need to let python know that these will definely both be defined :-/
    sequences, parents = {}, {}
    with open(outfile, 'rU') as fh:
        for sect in sections(fh):
            if sect == 'parents':
                parents = { full_name(child):(full_name(parent), distance) for child, parent, distance in iter_edges(fh) }
            elif sect[0] == 'sequences':
                sequences = parse_seqdict(fh, sect[1])
                # sanity check;  a valid tree should have exactly one node that is parentless
                if not len(parents) == len(sequences) - 1:
                    raise RuntimeError('invalid results attempting to parse {}: there are {} parentless sequences'.format(outfile, len(sequences) - len(parents)))
                trees.append(build_tree(sequences, parents))
            else:
                raise RuntimeError("unrecognized phylip section = {}".format(sect))
    return trees

# build a tree from a set of sequences and an adjacency dict.
def build_tree(sequences, parents, naive='.*naive.*'):
    # build an ete tree
    # first a dictionary of disconnected nodes
    nodes = {}
    for name in sequences:
        node = Tree()
        node.name = name
        node.add_feature('sequence', sequences[node.name])
        node.dist = parents[name][1] if name in parents else 0
        nodes[name] = node
    for name in sequences:
        if name in parents:
            nodes[parents[name][0]].add_child(nodes[name])
        else:
            tree = nodes[name]
    # reroot on naive
    tree = reroot_tree(tree, pattern=naive):

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
    if len(node.children) != 0:
        raise ValueError('specified naive node is not a leaf')
    if node not in tree.children:
        raise ValueError('specified naive node is not an outgroup')
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
    ts.legend.add_face(ete3.faces.CircleFace(8, 'brown', "circle"), 0)
    ts.legend.add_face(ete3.faces.TextFace("seedlineage"), 1)


def render_tree(fname, tree, annotations, highlight_node):
    "render tree SVG"
    ts = TreeStyle()
    ts.show_leaf_name = False
    tp_colors = timepoint_colors(annotations)

    seed_lineage = [n for n in iter_lineage(tree, highlight_node)]

    def my_layout(node):
        name = node.name
        seqmeta = annotations.get(node.name)
        # handle leaves
        if seqmeta and not re.compile(".*naive.*").match(node.name):
            # Add the node name for tips
            name = node.name + " (mf={}) ".format(round(float(seqmeta['mut_freqs']), 3))
            F = TextFace(name)
            add_face_to_node(F, node, column=0, position='branch-right')
            # Style the node with color corresponding to timepoint
            nstyle = NodeStyle()
            if node.name == highlight_node:
                nstyle['fgcolor'] = 'brown'
            else:
                nstyle['fgcolor'] = tp_colors[seqmeta['timepoint']]
            nstyle['size'] = 14
            node.set_style(nstyle)
        # Deal with naive and true internal nodes
        else:
            # Have to set all node sizes to 0 or we end up distorting the x-axis
            nstyle = NodeStyle()
            nstyle['size'] = 0
            node.set_style(nstyle)
            # Highlight just those nodes in the seedlineage which have nonzero branch lengths, or bad things
            if node in seed_lineage and node.dist > 0:
                position = "float"
                # Note; we use circleface instead of node style to avoid borking the x-axis
                circle_face = ete3.CircleFace(radius=5, color='brown', style="circle")
                add_face_to_node(circle_face, node, column=0, position=position)

    ts.layout_fn = my_layout
    timepoint_legend(ts, tp_colors)
    # whether or not we had rerooted on naive before, we want to do so for the SVG tree
    if 'naive' not in tree.name:
        tree = reroot_tree(tree, '.*naive.*')
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

    trees = parse_outfile(args.phylip_outfile, args.seqmeta)

    # take the first one, arbitrarily
    tree = trees[0]

    sequences =[SeqRecord(Seq(node.sequence), id=full_name(node.name), description="") for node in tree.traverse()]

    # write alignment to fasta
    fname = outbase + '.fa'
    with open(fname, "w") as fh:
        SeqIO.write(sequences, fh, "fasta")

    # write newick file
    fname = outbase + '.nwk'
    tree.write(format=1, format_root_node=True, outfile=fname)

    render_tree(outbase+'.svg', tree, args.seqmeta, args.seed)

    # all the trees
    # write newick file
    fname = outbase + '.alltrees.nwk'
    for tree in trees:
        tree.write(format=1, format_root_node=True, outfile=fname)

if __name__ == "__main__":
    main()
