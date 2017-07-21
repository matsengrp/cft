#!/usr/bin/env python

# A little thing for running prank
# Ideally, this will eventually let you reconstruct according to the original multifurcated tree as well

import argparse
import subprocess
#import warnings
from Bio import SeqIO
import ete3
import uuid


def random_tmpfilename():
    return '/tmp/joker.py-' + str(uuid.uuid1())

def get_args():
    parser = argparse.ArgumentParser(description="ancestral state reconstruction using prank")
    parser.add_argument('-t', '--input-tree', help="option guide tree (will generate if not specified)")
    parser.add_argument('seqs')
    # Final outputs
    parser.add_argument('asr_tree', help="authoritative output tree, with node names and/or splits matching those of asr_seqs")
    parser.add_argument('asr_seqs', help="ancestral sequence reconstruction")
    # These might be useful options later...
    parser.add_argument('-s', '--splits', action='store_true', help='name internal nodes using splits')
    parser.add_argument('--internal_node_prefix', default='in-')
    parser.add_argument('-S', '--splits-separator', default=',',
            help="Seprator to use in splits output; default , in seqs and _ in tree file")
    parser.add_argument('--no-reroot', action='store_true',
            help="Don't reroot the tree prior to bifurcation; warning this may fail")
    parser.add_argument('--preserve-internal-nodenames', action='store_true', 
            help="Try to preserve the internal node names in rerooting and bifurcation whenever possible")
    parser.add_argument('--keep-gaps', action='store_true', help="trigger -keep gaps option of PRANK")
    parser.add_argument('-v', '--verbose', action='store_true')
    # These are more or less for tinkering with the internals and won't generally be called
    parser.add_argument('--rerooted-tree', help="optionally, output the rerooted input tree to this file, when applicable")
    parser.add_argument('--bifurcated-tree', help="optionally, output the bifurcated input tree to this file, when applicable")
    parser.add_argument('--final-input-tree', help="optionally, output the final input tree to this file, when applicable")
    parser.add_argument('--prank-output-tree', help="optionally, pranks actual output tree")
    args = parser.parse_args()
    args.outbase = '.'.join(args.asr_seqs.split('.')[:-1]) if '.' in args.asr_seqs else args.asr_seqs # dragons
    args.rerooted_tree = args.rerooted_tree or random_tmpfilename()
    args.bifurcated_tree = args.bifurcated_tree or random_tmpfilename()
    args.final_input_tree = args.final_input_tree or random_tmpfilename()
    args.prank_output_tree = args.prank_output_tree or random_tmpfilename()
    return args



def cp_file(f1, f2):
    return subprocess.check_call(['cp', f1, f2])


# Define functions for running prank, preprocessing, etc

def reroot_tree(args):
    command = ['nw_reroot', args.input_tree]
    tree_str = subprocess.check_output(command)
    with open(args.rerooted_tree, 'w') as fh:
        fh.write(tree_str)

def bifurcate_tree(args):
    with open(args.rerooted_tree) as rerooted_tree_handle:
        tree = ete3.Tree(rerooted_tree_handle.read(), format=1)
        tree.resolve_polytomy()
        tree.write(outfile=args.bifurcated_tree, format=1)

# TODO I think we may want -keep as well; erm... maybe not? need to look back...
prank_options = '-quiet -showtree -showanc -showevents -DNA -f=fasta'.split(' ')

def run_prank(args):
    outbase = '.'.join(args.asr_seqs.split('.')[:-1]) # dragons
    command = ['prank', '-d='+args.seqs, '-o='+outbase] + prank_options
    if args.input_tree:
        command += ['-t='+args.final_input_tree, '-once']
    if args.keep_gaps:
        command.append('-keep')
    if args.verbose:
        print "\nRunning prank command:"
        print " ".join(command)
    prank_out = subprocess.check_output(command)
    if args.verbose:
        print "\nPrank stdout:"
        print prank_out
    return prank_out


def process_prank_results(args):
    outfile = args.asr_seqs
    infile = args.outbase+'.best.anc.fas'
    # TODO; should be final_input_tree?
    # TODO; should also check that the toplogies match her, as in the test script
    prank_output_tree = args.outbase + '.best.dnd'
    if args.prank_output_tree:
        cp_file(prank_output_tree, args.prank_output_tree)
    treefile = args.final_input_tree if args.input_tree else prank_output_tree
    tree = ete3.Tree(treefile, format=1)
    seqs = list(SeqIO.parse(infile, 'fasta'))
    for i, seq in enumerate(seqs):
        # Should validate no # in input data; TODO
        if i % 2:
            # handle as internal node
            leaves = [tree.get_leaves_by_name(seqs[j].id)[0] for j in [i-1, i+1]]
            mrca = leaves[0].get_common_ancestor(leaves[1])
            if args.splits:
                new_name = args.splits_separator.join(n.name for n in mrca.get_leaves())
            else:
                if mrca.name and args.preserve_internal_nodenames:
                    new_name = mrca.name
                else:
                    new_name = args.internal_node_prefix + seq.id.replace('#', '')
            # Make sure that all agree
            seq.id = new_name
            seq.name = new_name
            seq.description = new_name
            mrca.name = new_name
            
        else:
            # handle as leaf node (don't need to do anything here...)
            pass
    with open(args.asr_seqs, 'w') as fh:
        SeqIO.write(seqs, fh, 'fasta')
    tree.write(outfile=args.asr_tree, format=1)
    return outfile

def main():
    args = get_args()
    if args.input_tree:
        if not args.no_reroot:
            reroot_tree(args)
        else:
            cp_file(args.input_tree, args.rerooted_tree)
        bifurcate_tree(args)
        cp_file(args.bifurcated_tree, args.final_input_tree)
    run_prank(args)
    process_prank_results(args)
    if args.verbose:
        print "Have a nice day"


if __name__ == '__main__':
    main()


