#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given an outputfile from one of the PHYLIP tools - `dnaml` or `dnapars` - produce an alignment (including
ancestral sequences), a newick tree (with matching internal node lables), and an svg rendering of said tree.
"""
from __future__ import print_function

import argparse
import ete3
import csv
import re
import os
from warnings import warn
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import plot_tree

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
    # EH: this I'm guessing is what the lines being parsed by the regular expression will look like?
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
    # EH: this I'm guessing is what the lines being parsed by the regular expression will look like?
    #  152        sssssssssG AGGTGCAGCT GTTGGAGTCT GGGGGAGGCT TGGTACAGCC TGGGGGGTCC
    seqs = defaultdict(str)
    patterns = {
        'dnaml': re.compile("^\s*(?P<id>[a-zA-Z0-9>_.-]*)\s+(?P<seq>[a-zA-Z \-]+)"),
        'dnapars': re.compile("^\s*\S+\s+(?P<id>[a-zA-Z0-9>_.-]*)\s+(yes\s+|no\s+|maybe\s+)?(?P<seq>[a-zA-Z \-]+)")}
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
def parse_outfile(outfile, seqmeta, inferred_naive_name, seqname_mapping=None):
    '''parse phylip outfile'''
    name_map = seqname_mapping or {seq_id[0:10]: seq_id for seq_id in seqmeta}
    if len(name_map) != len(seqmeta) and not seqname_mapping:
        raise RuntimeError("Conflicting shortened names!")

    # Keeping track of which names match up for validation and error handling
    goods = set()
    bads = set()
    def full_name(x):
        name = name_map.get(x)
        if not name and not re.match('\d*$', x):
            bads.add(x)
        elif name:
            goods.add(x)
        return name or x
    sequences, parents = [], {}
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

    if inferred_naive_name in bads:
        warn("Naive sequence unique ID not found!")
        bads.remove(inferred_naive_name)
    if bads:
        print("found sequences:", sorted(goods))
        print("not found sequences:", sorted(bads))
        print("name_map misses:")
        for k, v in sorted(name_map.items()):
            if k not in goods:
                print("  ", k, ":", v)
        raise RuntimeError("mismatched sequence/node names")

    # sanity check;  a valid tree should have exactly one node that is parentless
    if not len(parents) == len(sequences) - 1:
        raise RuntimeError('invalid results attempting to parse {}: there are {} parentless sequences'.format(outfile, len(sequences) - len(parents)))

    return sequences, parents

# build a tree from a set of sequences and an adjacency dict.
def build_tree(sequences, parents):
    # build an ete tree
    # first a dictionary of disconnected nodes
    def mknode(r):
        n = ete3.Tree(
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
# Usually this is used to root on the naive germline sequence
# NOTE duplicates fcn in plot_tree.py
def reroot_tree(tree, pattern):
    # find all nodes matching pattern
    node = find_node(tree, pattern)
    if tree != node:
        # In general this would be necessary, but we are actually assuming that naive has been set as an
        # outgroup in dnaml, and if it hasn't, we want to raise an error, as below
        #tree.set_outgroup(node)
        # Raise error if naive isn't in root's children
        if node not in tree.children:
            raise ValueError("Naive node not set as outgroup; Check dnaml/asr run to make sure this is the case")
        # This actually assumes the `not in` condition above, but we check as above for clarity
        tree.remove_child(node)
        node.add_child(tree)
        tree.dist = node.dist
        node.dist = 0
        tree = node
    return tree

def seqname_mapping_arg(filename):
    with open(filename, 'r') as fh:
        return {row['new_id']: row['original_id'] for row in csv.DictReader(fh)}

def get_args():
    def seqmeta_input(fname):
        print(fname)
        with open(fname, 'r') as fhandle:
            return dict((row['sequence'], row) for row in csv.DictReader(fhandle))

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
        '--seqname-mapping',
        type=seqname_mapping_arg,
        help="optional csv file translating original_id to new_id for phylip 10 char cludgery")
    parser.add_argument(
        '--outdir', default='.',
        help="output directory where results are left.  [default '%(default)s']")
    parser.add_argument(
        '--basename', help="basename of output files.   [default 'basename(DNAML)']")
    parser.add_argument(
        '--inferred-naive-name', type=str, required=True)
    parser.add_argument(
        '--seed', type=str, help="id of leaf [default 'seed']", default='seed')
    return parser.parse_args()

def main():
    args = get_args()

    # basename of outputfiles is taken from --basename if supplied, otherwise basename of DNAML.
    basename = args.basename if args.basename else os.path.basename(args.phylip_outfile)
    outbase = os.path.join(args.outdir, os.path.splitext(basename)[0])

    sequences, parents = parse_outfile(args.phylip_outfile, args.seqmeta, args.inferred_naive_name, args.seqname_mapping)

    if not sequences or not parents:
        raise RuntimeError("No sequences were available; are you sure this is a dnaml output file?")

    # write alignment to fasta
    fname = outbase + '.fa'
    with open(fname, "w") as fh:
        SeqIO.write(sequences, fh, "fasta")

    # write just inferred ancestors, inferred naive, and seed to a fasta; we later check to see if we sampled any of these sequences (checking against sampled sequences in the cluster) or any sequences that are close to them.
    tip_names = set(args.seqname_mapping.values())
    all_names = set([record.id for record in sequences])
    ancestors = all_names - tip_names
    ancestors_naive_and_seed = ancestors.union(set([args.inferred_naive_name, args.seed]))
    fname = outbase + '.ancestors_naive_and_seed.fa'
    with open(fname, "w") as fh:
        SeqIO.write([record for record in sequences if record.id in ancestors_naive_and_seed], fh, "fasta")

    tree = build_tree(sequences, parents)
    tree = reroot_tree(tree, args.inferred_naive_name)

    # write newick file
    fname = outbase + '.nwk'
    tree.write(format=1, format_root_node=True, outfile=fname)
    
    plot_tree.render_tree(outbase+'.svg', tree, args.seqmeta, args.seed, args.inferred_naive_name)

if __name__ == "__main__":
    main()
