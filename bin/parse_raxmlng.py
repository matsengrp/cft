#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
Parse the output ancestral state reconstruction from RAxML-NG together with the topology
to create an ete3 tree with the ancestral sequences.
'''

import re
from warnings import warn
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from ete3 import Tree
import argparse

# TODO remove unused code, and run black

class TreeFileParsingError(Exception):
    '''When ete3 fails to read the input tree.'''

def ASR_parser(args):
    try:
        tree = Tree(args.tree, format=1)
    except Exception as e:
        print(e)
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    # Write complete fasta
    write_tree_fasta(args.asr_seq, args.input_seq, args.outbase + '.fa')

    # Reroot to make the naive sequence the real root instead of just an outgroup:
    tree = reroot_tree(tree, pattern=args.naive)

    # Dump tree as newick:
    tree.write(outfile='{}.nwk'.format(args.outbase), format_root_node=True, format=1)
    print('Done parsing RAxML-NG tree')

def parse_raxmlng_ancestral_state(line):
    asr_seqid, asr_seq = line.strip().split()
    return SeqRecord(Seq(asr_seq), id=asr_seqid, name='', description='')

def write_tree_fasta(asr_seqs_fname, input_seqs_fname, outfname):
    '''Combines RAxML-NG asr states and input sequence records.'''
    input_records = list(SeqIO.parse(input_seqs_fname, "fasta"))
    with open(asr_seqs_fname) as fh:
        asr_records = [parse_raxmlng_ancestral_state(l) for l in fh]
    # Check that ASR lengths are same as input lengths
    assert({len(str(s.seq)) for s in asr_records} == {len(str(s.seq)) for s in input_records})
    SeqIO.write(input_records + asr_records, outfname, "fasta")

def find_node(tree, pattern):
    # TODO same as reroot todo
    regex = re.compile(pattern).search
    nodes =  [node for node in tree.traverse() for m in [regex(node.name)] if m]
    if not nodes:
        warn("Cannot find matching node; looking for name matching '{}'".format(pattern))
        return
    else:
        if len(nodes) > 1:
            warn("multiple nodes found; using first one.\nfound: {}".format([n.name for n in nodes]))
        return nodes[0]


# reroot the tree on node matching regex pattern.
# Usually this is used to root on the naive germline sequence with a name matching '.*naive.*'
def reroot_tree(tree, pattern='.*naive.*', outgroup=False):
    # TODO replace with reroot from CFT when putting into CFT
    # find all nodes matching pattern
    node = find_node(tree, pattern)
    tree.set_outgroup(node)
    if tree != node and outgroup:
        s = node.get_sisters()  # KBH - want the outgroup with a branch length of (near) zero
        s[0].dist = node.dist * 2
        node.dist = 0.0000001   # KBH - actual zero length branches cause problems
        tree.swap_children()    # KBH - want root to be at the last taxon in the newick file.
    elif tree != node:
        tree.remove_child(node)
        node.add_child(tree.get_children()[0])
        tree.dist = node.dist
        node.dist = 0
        tree = node
    return tree


def main():
    parser = argparse.ArgumentParser(description='Tools for running ASR with RAxML-NG.')
    parser.add_argument('--tree', required=True, metavar='NEWICK TREE', help='Input tree used for topology.')
    parser.add_argument('--asr-seq', required=True, help='Input ancestral sequences file.')
    parser.add_argument('--input-seq', required=True, help='Phylip file with input sequences.')
    parser.add_argument('--outbase', required=True, metavar='FILENAME', help='Filename for the output ASR tree.')
    parser.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    args = parser.parse_args()
    ASR_parser(args)

if __name__ == '__main__':
    main()
