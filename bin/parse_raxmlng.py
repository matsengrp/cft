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
import process_asr
# TODO (do this last): run black

class TreeFileParsingError(Exception):
    '''When ete3 fails to read the input tree.'''

def ASR_parser(args):
    try:
        tree = Tree(args.tree, format=1)
    except Exception as e:
        print(e)
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    # Write complete fasta
    write_tree_fastas(args.asr_seq, args.input_seq, args.inferred_naive_name, args.seed, args.outbase)

    # Reroot to make the naive sequence the real root instead of just an outgroup:
    tree = reroot_tree(tree, pattern=args.inferred_naive_name)

    # Dump tree as newick:
    tree.write(outfile='{}.nwk'.format(args.outbase), format_root_node=True, format=1)
    print('Done parsing RAxML-NG tree')

def write_ancestors_naive_and_seed(input_seqs, ancestors, inferred_naive_name, seed, outbase):
    """
    Write a fasta containing inferred ancestors, inferred naive, and the seed sequence,
    all of which are used as query sequences input to BLAST to search among sampled sequences
    in order to validate inference.
    """
    naive_and_seed_records = [r for r in input_seqs if r.id in {inferred_naive_name, seed}]
    SeqIO.write(ancestors + naive_and_seed_records, outbase + ".ancestors_naive_and_seed.fa", "fasta")

def parse_raxmlng_ancestral_state(line):
    asr_seqid, asr_seq = line.strip().split()
    return SeqRecord(Seq(asr_seq), id=asr_seqid, name='', description='')

def write_tree_fastas(asr_seqs_fname, input_seqs_fname, inferred_naive_name, seed, outbase):
    '''Combines RAxML-NG asr states and input sequence records.'''
    input_records = list(SeqIO.parse(input_seqs_fname, "fasta"))
    with open(asr_seqs_fname) as fh:
        asr_records = [parse_raxmlng_ancestral_state(l) for l in fh]
    # Check that ASR lengths are same as input lengths
    assert({len(str(s.seq)) for s in asr_records} == {len(str(s.seq)) for s in input_records})
    SeqIO.write(input_records + asr_records, outbase + ".fa", "fasta")
    write_ancestors_naive_and_seed(input_records, asr_records, inferred_naive_name, seed, outbase)

# reroot the tree on node matching regex pattern.
# Usually this is used to root on the naive germline sequence with a name matching '.*naive.*'
def reroot_tree(tree, pattern, outgroup=False):
    # TODO do we need "outgroup" option/case?
    # find all nodes matching pattern
    node = process_asr.find_node(tree, pattern)
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
    parser.add_argument("--inferred-naive-name", type=str, required=True, default='inferred_naive', help='naive sequence id')
    parser.add_argument(
        "--seed", type=str, help="id of seed sequence, default: 'seed'", default="seed"
    )
    args = parser.parse_args()
    ASR_parser(args)

if __name__ == '__main__':
    main()
