#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Parse the output ancestral state reconstruction from RAxML-NG together with the topology
to create an ete3 tree with the ancestral sequences.
"""

import re
from warnings import warn
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from ete3 import Tree
import argparse

# TODO (do this last): run black


class TreeFileParsingError(Exception):
    """When ete3 fails to read the input tree."""


def ASR_parser(args):
    try:
        tree = Tree(args.tree, format=1)
    except Exception as e:
        print(e)
        raise TreeFileParsingError(
            "Could not read the input tree. Is this really newick format?"
        )

    # Write complete fasta
    write_tree_fastas(
        args.asr_seq, args.input_seq, args.inferred_naive_name, args.seed, args.outbase
    )

    # Root on the naive sequence:
    rooted_tree = root_tree(tree, args.inferred_naive_name)

    # Dump tree as newick:
    rooted_tree.write(outfile="{}.nwk".format(args.outbase), format_root_node=True, format=1)
    print("Done parsing RAxML-NG tree")


def write_ancestors_naive_and_seed(
    input_seqs, ancestors, inferred_naive_name, seed, outbase
):
    """
    Write a fasta containing inferred ancestors, inferred naive, and the seed sequence,
    all of which are used as query sequences input to BLAST to search among sampled sequences
    in order to validate inference.
    """
    naive_and_seed_records = [
        r for r in input_seqs if r.id in {inferred_naive_name, seed}
    ]
    SeqIO.write(
        ancestors + naive_and_seed_records,
        outbase + ".ancestors_naive_and_seed.fa",
        "fasta",
    )


def parse_raxmlng_ancestral_state(line):
    asr_seqid, asr_seq = line.strip().split()
    return SeqRecord(Seq(asr_seq), id=asr_seqid, name="", description="")


def write_tree_fastas(
    asr_seqs_fname, input_seqs_fname, inferred_naive_name, seed, outbase
):
    """Combines RAxML-NG asr states and input sequence records."""
    input_records = list(SeqIO.parse(input_seqs_fname, "fasta"))
    with open(asr_seqs_fname) as fh:
        asr_records = [parse_raxmlng_ancestral_state(l) for l in fh]
    # Check that ASR lengths are same as input lengths
    assert {len(str(s.seq)) for s in asr_records} == {
        len(str(s.seq)) for s in input_records
    }
    SeqIO.write(input_records + asr_records, outbase + ".fa", "fasta")
    write_ancestors_naive_and_seed(
        input_records, asr_records, inferred_naive_name, seed, outbase
    )


def root_tree(tree, desired_root_name):
    """
    Root an unrooted tree on a given node. Unrooted trees are represented by ete3 as a trifurcation at the root,
    which means we end up with an extra empty-string-named node when we set an outgroup. See below for how this is handled
    in order to yield a rooted tree on the given node, without the empty node.
    """
    node = tree.search_nodes(name=desired_root_name)[0]
    if tree != node:
        tree.set_outgroup(node)
        #outgrouping an unrooted tree causes an empty string named node, whose actual name goes to the root for some reason. Swap them back:
        tree.search_nodes(name="")[0].name = tree.name
        tree.name = ""
        #Since we've outgrouped "node" the root should have two children: "node" and it's "sister" (aka the rest of the tree). We actually need "node" to be the root. First take off "node":
        tree.remove_child(node)
        #"sister" remains now as the other child of the tree
        sister = tree.get_children()[0]
        # add to the "sister" branch the branch length we lost when we removed "node" (should be equal branch length so could multiply sister.dist by 2 instead of doing it this way but this reads more clearly).
        sister.dist = sister.dist + node.dist
        node.dist = 0
        # attach "sister" (and the rest of the tree below it) to "node" to yield a tree rooted on "node"
        node.add_child(sister)
        #return "node" below since "node" is now the root
    return node

def main():
    parser = argparse.ArgumentParser(description="Tools for running ASR with RAxML-NG.")
    parser.add_argument(
        "--tree",
        required=True,
        metavar="NEWICK TREE",
        help="Input tree used for topology.",
    )
    parser.add_argument(
        "--asr-seq", required=True, help="Input ancestral sequences file."
    )
    parser.add_argument(
        "--input-seq", required=True, help="Phylip file with input sequences."
    )
    parser.add_argument(
        "--outbase",
        required=True,
        metavar="FILENAME",
        help="Filename for the output ASR tree.",
    )
    parser.add_argument(
        "--inferred-naive-name",
        type=str,
        required=True,
        default="inferred_naive",
        help="naive sequence id",
    )
    parser.add_argument(
        "--seed", type=str, help="id of seed sequence, default: 'seed'", default="seed"
    )
    args = parser.parse_args()
    ASR_parser(args)


if __name__ == "__main__":
    main()
