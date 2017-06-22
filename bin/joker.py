#!/usr/bin/env python

# A little thing for running prank
# Ideally, this will eventually let you reconstruct according to the original multifurcated tree as well

import argparse
import subprocess
import ete3
import uuid


def random_tmpfilename():
    return '/tmp/joker.py-' + str(uuid.uuid1())

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tree')
    parser.add_argument('seqs')
    #parser.add_argument('--outgroup')
    parser.add_argument('--rerooted-tree')
    parser.add_argument('--bifurcated-tree')
    parser.add_argument('asr_seqs')
    args = parser.parse_args()
    args.rerooted_tree = args.rerooted_tree or random_tmpfilename()
    args.bifurcated_tree = args.bifurcated_tree or random_tmpfilename()
    return args


# Define functions for running prank, preprocessing, etc

def reroot_tree(args):
    command = ['nw_reroot', args.tree]
    tree_str = subprocess.check_output(command)
    with open(args.rerooted_tree, 'w') as fh:
        fh.write(tree_str)

def bifurcate_tree(args):
    with open(args.rerooted_tree) as rerooted_tree_handle:
        tree = ete3.Tree(rerooted_tree_handle.read(), format=1)
        tree.resolve_polytomy()
        tree.write(outfile=args.bifurcated_tree, format=1)

prank_options = '-once -quiet -showanc -showevents -DNA -f=fasta'.split(' ')

def run_prank(args):
    outfile = args.asr_seqs
    outbase = '.'.join(args.asr_seqs.split('.')[:-1]) # dragons
    command = ['prank', '-d='+args.seqs, '-t='+args.bifurcated_tree, '-o='+outbase] + prank_options
    print "\nRunning prank command:"
    print " ".join(command)
    prank_out = subprocess.check_output(command)
    print "\nPrank stdout:"
    print prank_out
    prank_out = subprocess.check_output(['seqmagick', 'convert', '--pattern-replace', '#', '', outbase+'.best.anc.fas', outfile])
    print "Joker has copied the 'best.anc.fas' file to", outfile
    print "Have a nice day"
    return prank_out


def main():
    args = get_args()
    reroot_tree(args)
    bifurcate_tree(args)
    run_prank(args)


if __name__ == '__main__':
    main()



# Tried and failed still...

# 28027  prank -once -d=pruned.fa -t=asr.rerooted.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fastas
# 28028  prank -once -d=pruned.fa -t=asr_rrtd.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fastas

# 28029  less pruned.fa
# 28030  less asr_rrtd.nwk


# Realized it was the multifurcation so copying over ml results

# 28031  cp output/kate-qrs-v10-dnaml/QA255.006-VL/QA255-l-IgL/run-viterbi-best-plus-0/asr.nwk asr2.nwk
# 28032  cp output/kate-qrs-v10-dnaml/QA255.006-VL/QA255-l-IgL/run-viterbi-best-plus-0/pruned.fa pruned2.fa
# ---------------------------------------------------------------------------------------------------------


# Now trying prank

# typos..
# 28033  prank -once -d=pruned2.fa -t=asr2.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fastas
# 28034  prank -once -d=pruned2.fa -t=asr2.nwk -o=test.fasta -once -quiet -showanc -showevents -DNA -f=fastas
# 28035  prank -once -d=pruned2.fa -t=asr2.nwk -o=test.fast -once -quiet -showanc -showevents -DNA -f=fastas

# There we go; this is the one:
# 28036  prank -once -d=pruned2.fa -t=asr2.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fasta

# reroot!
# 28037  nw_reroot asr2.nwk > asr2.rrtd.nwk

# try running again
# 28038  prank -once -d=pruned2.fa -t=asr2.rrtd.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fasta

# still not working? trying condense wtf?
# 28039  nw_condense -h
# oh god, was this where I tried manually editing?
# 28040  vim asr2.nwk
# 28041  vim pruned2.fa
# 28042  prank -once -d=pruned2.fa -t=asr2.rrtd.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fasta
# 28043  vim asr2.rrtd.nwk
# 28044  prank -once -d=pruned2.fa -t=asr2.rrtd.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fasta
# 28045  j prank
# 28046  seqids pruned2.fa
# 28047  nw_labels asr2.rrtd.nwk
# 28048  prank -once -d=pruned2.fa -t=asr2.rrtd.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fasta
# 28049  prank -version
# 28050  j encap
# 28051  mkdir prank
# 28052  ls
# 28053  rmdir prank


# Oh right... tried for all hell long, and then realized it all had to do with the version of prank that was
# installed. So now updating encap...
# 28054  wget http://wasabiapp.org/download/prank/prank.linux64.150803.tgz
# 28055  tar -gxf prank.linux64.150803.tgz
# 28056  tar -zxf prank.linux64.150803.tgz
# 28057  ls prank
# 28058  mv prank prank.linux64.150803
# 28059  less prank.linux64.150803/prank.1
# 28060  ln -s ../encap/prank.linux64.150803/bin/prank ../bin/prank
# 28061  ls ../bin/prank
# 28062  l ../bin/prank
# 28063  ll ../bin/prank
# 28064  mv ../bin/prank ../bin/prank-old
# 28065  ln -s ../encap/prank.linux64.150803/bin/prank ../bin/prank
# 28066  prank -version
# 28067  l  ../bin/prank
# 28068  ./prank.linux64.150803/bin/pra
# 28069  ./prank.linux64.150803/bin/prank
# 28070  l ./prank.linux64.150803/bin/prank
# 28071  ls ./prank.linux64.150803/bin/prank
# 28072  ls ./prank.linux64.150803/bin/
# 28073  ls prank.linux64.150803
# 28074  ls
# 28075  ls astyle-r409
# 28076  l prank.linux64.150803
# 28077  tar -zxf prank.linux64.150803.tgz
# 28078  ls prank
# 28079  l prank
# 28080  chmod +x prank/bin/prank
# 28081  chmod +r prank/bin/
# 28082  chmod +r prank/
# 28083  mv prank prank.linux64.150803
# 28084  mv prank.linux64.150803 prank.linux64.150803-brokerz
# 28085  mv prank prank.linux64.150803
# 28086  l ../bin/prank
# 28087  which prank
# 28088  prank -version

# Trying to run again?
# 28089  prank -once -d=pruned2.fa -t=asr2.rrtd.nwk -o=test.fa -once -quiet -showanc -showevents -DNA -f=fasta

# And success
# Checking out results
# 28090  less test.fa.best.events
# 28091  less test.fa.best.
# 28092  less test.fa.best.anc.fas
# 28093  less test.fa.best.fas
# 28094  less test.fa.best.anc.dnd
