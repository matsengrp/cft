#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Carry partis output through prepation for display via cftweb application.

This SConstruct does the following:

* Runs process_partis.py for each seed/sequencing-run/parition combination, producing a metadata file of information
* For each such parition, the cluster with the seed is analyzed:
    * Tree construction using `FastTree`
    * Ancestral state reconstruction using `dnaml`
    * Some SVG visualization
* Results of all these analyses are then pointed to by a merged metadata file, which is then consumed by cftweb

Clusters with only two sequences cannot be analyzed by dnaml or FastTree, and so these are skipped atuomatically.
Additionally, some clusters end up with bad trees after a rerooting step, and also have trouble in our dnaml2tree step.
These are left out of the final `metadata.json` file.
Eventually, all of these cases should be caught and tombstone information should be left in the metadata file capturing as much infomration as possible.

See README for typical environment setup and usage.
'''


from __future__ import print_function
import os
import re
import sys
import csv
import subprocess
import glob
import shutil
import functools
import sconsutils
import datetime
import getpass

from os import path
from warnings import warn
from nestly import Nest
from nestly.scons import SConsWrap
from SCons.Script import Environment, AddOption

import json

# Need this in order to read csv files with sequences in the fields
csv.field_size_limit(sys.maxsize)

# No-op; Prevents analysis warnings
sconsutils


# Set up SCons environment

environ = os.environ.copy()

# If the PARTIS env var isn't already set, default to $PWD/partis (where we have a git
# submodule checkout; this is needed for bin/process_partis.py)
environ['PARTIS'] = environ.get('PARTIS', path.join(os.getcwd(), 'partis'))

env = Environment(ENV=environ)

# Add stuff to PATH
env.PrependENVPath('PATH', 'bin')
env.PrependENVPath('PATH', 'post_partis/scripts')
env.PrependENVPath('PATH', 'tree')



# Some environment sanity checks to make sure we have all prerequisits

def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# only check for dependencies if we aren't in a dry-run.
if not env.GetOption('no_exec'):
    msg = ""
    if not cmd_exists('dnaml'):
        msg += '''
           Required dependency command, 
           `dnaml` not found on PATH
           Consider using,
                $ module use ~matsengrp/modules
                $ module load phylip
            '''

    if not cmd_exists('seqmagick'):
        msg += '''
           Required dependency command, 
           `seqmagick` not found on PATH
           Consider using,
                $ module load seqmagick
            '''
        
    if not cmd_exists('FastTree'):
        msg += '''
           Required dependency command, 
           `FastTree` not found on PATH
           Consider using,
                $ module load FastTree
            '''
    # if we are missing any prerequisites, print a message and exit
    if len(msg):
        warn(msg)
        sys.exit(1)



# Setting up command line arguments/options

AddOption('--datapath',
              dest='datapath',
              type='string',
              nargs=1,
              action='store',
              metavar='DIR',
              default="/fh/fast/matsen_e/processed-data/partis/kate-qrs-2016-09-09/new",
              help='partis output directory')

AddOption('--outdir',
        dest='outdir',
        type='string',
        nargs=1,
        action='store',
        metavar='DIR',
        default="output",
        help="directory in which to output results")

AddOption('--test',
        dest='test_run',
        action='store_true',
        default=False,
        help="Setting this flag does a test run for just a couple of seeds")


datapath = env.GetOption('datapath')
outdir_base = env.GetOption('outdir') # we call this outdir_base in order to not conflict with nestly fns outdir arg
test_run = env.GetOption("test_run")

print("datapath = {}".format(datapath))
print("outdir = {}".format(outdir_base))
print("test_run = {}".format(test_run))



# Initialize nestly!
# ==================

# This lets us create parameter nests and run the given pipeline for each combination of nesting parameters.
# It also lets us "pop" out of nesting levels in order to aggregate on nested results.

# Our nesting is going to more or less model the nesting of things in our datapath directory.
# seed > parameter_set > partition

# Here we initialize nestly, and create an scons wrapper for it

nest = Nest()
w = SConsWrap(nest, outdir_base, alias_environment=env)


# Here we add some aggregate targets which we'll be build upon in order to spit out our final metadata.json file.

w.add_aggregate('metadata', list)
w.add_aggregate('svgfiles', list)



# Initialize seed nest
# --------------------

# The very first nesting is on the seed id.

# For the sake of testing, we allow for switching between the full set of seeds, as returned by `seeds_fn`, and a small test set (`test_seeds`).
# This is controllable via the `--test` cli flag

test_seeds = ["QB850.424-Vk", "QB850.043-Vk"]

def seeds_fn(datapath):
    return os.listdir(path.join(datapath, 'seeds'))

seeds = test_seeds if test_run else seeds_fn(datapath)

# Initialize our first nest level
w.add('seed', seeds)


# Initialize parameter set nest
# -----------------------------

# Next we nest on parameter sets; I believe these correspond to sequencing runs.
# They're the things that look like "Hs-LN2-5RACE-IgG", and are nested within each seed's output directory.

w.add('parameter_set', lambda c: map(lambda x: path.basename(x), glob.glob(path.join(datapath, 'seeds', c['seed'], "*"))))


# Some helpers at the seed level

def input_dir(c):
    "Return the `seed > parameter_set` directory given the closed over datapath and a nestly control dictionary"
    # This 'seeds' thing here is potentially not a very general assumption; not sure how variable that might be upstream
    return path.join(datapath, 'seeds', c['seed'], c['parameter_set'])

def path_base_root(full_path):
    return path.splitext(path.basename(full_path))[0]

def valid_partition(fname):
    with open(fname, 'r') as partition_handle:
        partition = csv.DictReader(partition_handle).next()
        unique_ids = partition['unique_ids'].split(':')
        return len(unique_ids) >= 2 # Do we want >2 or >=2?


def annotations(c):
    """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
    for actual analysis."""
    return map(path_base_root,
               # We might eventually want to handle small partitions more manually, so there's data on the other end, but for now, filter is easiest
               filter(valid_partition,
                      glob.glob(path.join(input_dir(c), "*-plus-*.csv"))))

def partitions(c):
    "Returns the `partition.csv` file path with all the partition information for every partition output by partis."
    return path.join(input_dir(c), 'partition.csv')

def parameter_dir(c):
    "The input parameter directory (as in input to partis); see --parameter-dir flag to various partis cli calls."
    return path.join(datapath, c['parameter_set'])



# Each c["partition"] value actually points to the annotations for that partition... a little weird but...
w.add('partition', annotations)


def annotation(c):
    return path.join(input_dir(c), c['partition'] + '.csv')

@w.add_target()
def process_partis(outdir, c):
    return env.Command(path.join(outdir, 'metadata.json'),
            [partitions(c), annotation(c)],
            'process_partis.py ' +
                '--partition ${SOURCES[0]} ' +
                '--annotations ${SOURCES[1]} ' +
                #'--param_dir ' + parameter_dir(c) + ' '
                '--partis_log ' + path.join(input_dir(c), 'partition.log ') +
                '--cluster_base cluster ' +
                '--output_dir ' + outdir)

# To get the input sequences we need to read in the file output by process_partis to figure out what the input fasta was,
# and then for convenience copy it over locally

def extract_inseqs(outdir, target, source, env):
    with open(str(source[0]), "r") as source_handle:
        metadata = json.load(source_handle)[0]
        shutil.copy(path.join(outdir, metadata['file']), str(target[0]))

@w.add_target()
def inseqs(outdir, c):
    return env.Command(
        path.join(outdir, "inseqs.fasta"),
        c['process_partis'],
        functools.partial(extract_inseqs, outdir))

# Forget why we do this; Chris W. seems to think it may not have been as necessary as original thought.
# Regardless of what the deal was with this, I think we can take it out once we add an alignment for inseqs
#@w.add_target()
#def padded_seqs(outdir, c):
    #return env.Command(
        #path.join(outdir, "padded_seqs.fa"),
        #c['inseqs'],
        #"python bin/padseq.py $SOURCE > $TARGET")



# Sequence Alignment
# ==================

# What follows is a rather convoluted backtranslation strategy for sequence alignment.
# The basic idea is that we translate our sequences, align them (taking advantage of coding information),
# and from this infer the a nucleotide alignment using `seqmagick backtrans-align`.
# The problem is that `seqmagick` rather finicky here about sequence ordering and nucleotide seqs being in
# lengths of multiples of three.
# So there's a bunch of housekeeping here in making sure all of this is the case for our final backtrans-align
# step.

@w.add_target()
def translated_inseqs_(outdir, c):
    return env.Command(
        #path.join(outdir, "translated_inseqs.fasta"),
        #c['inseqs'],
        #'seqmagick convert --translate dna2protein $SOURCE $TARGET')
        [path.join(outdir, "translated_inseqs.fasta"), path.join(outdir, 'inseqs_trimmed.fasta')],
        [c['process_partis'], c['inseqs']],
        'translate_seqs.py $SOURCES $TARGET -t ${TARGETS[1]}')

@w.add_target()
def translated_inseqs(outdir, c):
    return c['translated_inseqs_'][0]

@w.add_target()
def trimmed_inseqs(outdir, c):
    return c['translated_inseqs_'][1]

@w.add_target()
def aligned_translated_inseqs(outdir, c):
    return env.SRun(
        path.join(outdir, "aligned_translated_inseqs.fasta"),
        c['translated_inseqs'],
        # Replace stop codons with X, or muscle inserts gaps, which messed up seqmagick backtrans-align below
        # Note that things will break down at backtrans-align if any seq ids have * in them...
        'sed \'s/\*/X/g\' $SOURCE | muscle -in /dev/stdin -out $TARGET')

@w.add_target()
def sorted_inseqs(outdir, c):
    return env.Command(
        path.join(outdir, "sorted_inseqs.fasta"),
        c['trimmed_inseqs'],
        'seqmagick convert --sort name-asc $SOURCE $TARGET')

@w.add_target()
def sorted_aligned_translated_inseqs(outdir, c):
    return env.Command(
        path.join(outdir, "sorted_aligned_translated_inseqs.fasta"),
        c['aligned_translated_inseqs'],
        'seqmagick convert --sort name-asc $SOURCE $TARGET')

# Now, finally, we can backtranslate

@w.add_target()
def aligned_inseqs(outdir, c):
    return env.Command(
        path.join(outdir, 'aligned_inseqs.fasta'),
        [c['sorted_aligned_translated_inseqs'], c['sorted_inseqs']],
        'seqmagick backtrans-align -a warn $SOURCES -o $TARGET')



# On with trees and other things...
# =================================

# use fasttree to make newick tree from sequences
@w.add_target()
def fasttree(outdir, c):
    return env.SRun(
        path.join(outdir, "fasttree.nwk"),
        c['aligned_inseqs'],
        "FastTree -nt -quiet $SOURCE > $TARGET")

# calculate list of sequences to be pruned
@w.add_target()
def pruned_ids(outdir, c):
    return env.Command(
        path.join(outdir, "pruned_ids.txt"),
        c['fasttree'],
        "python bin/prune.py --seed " + c['seed'] + " $SOURCE > $TARGET")

# prune out sequences to reduce taxa
@w.add_target()
def pruned_seqs(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.fa"),
        [c['pruned_ids'], c['aligned_inseqs']],
        "seqmagick convert --include-from-file $SOURCES $TARGET")

# Convert to phylip for dnaml
@w.add_target()
def phy(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.phy"),
        c['pruned_seqs'],
        "seqmagick convert $SOURCE $TARGET")

# Create a "config file" for dnaml (a persnickety old program with interactive menues...)
@w.add_target()
def dnaml_config(outdir, c):
    return env.Command(
        path.join(outdir, "dnaml.cfg"),
        c['phy'],
        "python bin/mkconfig.py $SOURCE > $TARGET")

# Run dnaml by passing in the "config file" as stdin hoping the menues all stay sane
# (Aside: There is a program for programatically responding to interactive menus if this gets any hairier)
@w.add_target()
def dnaml(outdir, c):
    "run dnaml (from phylip package) to create tree with inferred sequences at internal nodes"
    tgt = env.SRun(
        map(lambda x: path.join(outdir, x), ["outtree", "outfile", "dnaml.log"]),
        c['dnaml_config'],
        'cd ' + outdir + ' && dnaml < ${SOURCE.file} > ${TARGETS[2].file}')
    # Manually depend on phy so that we rerun dnaml if the input sequences change (without this, dnaml will
    # only get rerun if one of the targets are removed or if the iput dnaml_config file is changed).
    env.Depends(tgt, c['phy'])
    return tgt


@w.add_target()
def dnaml_tree(outdir, c):
    """parse dnaml output into fasta and newick files, and make SVG format tree with ETE package.
    xvfb-run is needed because of issue https://github.com/etetoolkit/ete/issues/101"""
    tgt = env.Command(
            map(lambda x: path.join(outdir, x),
                ["dnaml.svg", "dnaml.fa", "dnaml.seedLineage.fa", "dnaml.newick"]),
            c['dnaml'],
            # Note: the `-` prefix here tells scons to keep going if this command fails.
            "- xvfb-run -a bin/dnaml2tree.py --seed " + c['seed'] + " --dnaml ${SOURCES[1]} --outdir ${TARGETS[0].dir} --basename dnaml")
    # Manually depend on dnaml2tree.py script, since it doesn't fall in the first position within the command
    # string.
    env.Depends(tgt, 'bin/dnaml2tree.py')
    # Do aggregate work
    c['svgfiles'].append(tgt[0])
    return tgt

# Might want to switch back to this approach...
#def extended_metadata(metadata, dnaml_tree_tgt):
    #"Add dnaml_tree target(s) as metadata to the given metadata dict; used prior to metadata write."
    #with open(
    #m = copy.copy(metadata)
    #m['svg'] = path.relpath(str(dnaml_tree_tgt[0]), outdir_base)
    #m['fasta'] = path.relpath(str(dnaml_tree_tgt[1]), outdir_base)
    #m['seedlineage'] = path.relpath(str(dnaml_tree_tgt[2]), outdir_base)
    #m['newick'] = path.relpath(str(dnaml_tree_tgt[3]), outdir_base)
    #del m['file']
    #return m


@w.add_target()
def cluster_aa(outdir, c):
    return env.Command(
        path.join(outdir, 'cluster_aa.fa'),
        [c['process_partis'], c['dnaml_tree'][1]],
        'translate_seqs.py $SOURCES $TARGET')

# TODO I think we can take this out now that we compute lineages dynamically in cftweb
@w.add_target()
def seedlineage_aa(outdir, c):
    return env.Command(
        path.join(outdir, 'seedlineage_aa.fa'),
        [c['process_partis'], c['dnaml_tree'][2]],
        'translate_seqs.py $SOURCES $TARGET')


@w.add_target()
def cluster_metadata(outdir, c):
    # Note: We have to do all this crazy relpath crap because the assumption of cftweb is that all paths
    # specified in the input metadata.json file are relative to _its_ location, in our context, always base_outdir
    def relpath(tgt_path):
        return path.relpath(tgt_path, outdir_base)
    def dnaml_tgt_relpath(i):
        return relpath(str(c['dnaml_tree'][i]))
    n = re.compile('run-viterbi-best-plus-(?P<step>.*)').match(c['partition']).group('step')
    tgt = env.Command(
            path.join(outdir, 'extended_metadata.json'),
            c['process_partis'] + [partitions(c)] + c['dnaml_tree'],
            # Don't like the order assumptions on dnaml_tgts here...
            'assoc_logprob.py $SOURCE ${SOURCES[1]} -n ' + n + ' /dev/stdout | ' +
                'json_assoc.py /dev/stdin $TARGET ' +
                # Not 100% on this; Pick optimal attr/key name
                #'best_partition ' + c['partition'] + ' ' +
                'clustering_step ' + n + ' ' +
                'svg ' + dnaml_tgt_relpath(0) + ' ' +
                'fasta ' + dnaml_tgt_relpath(1) + ' ' +
                'seedlineage ' + dnaml_tgt_relpath(2) + ' ' +
                'cluster_aa ' + relpath(str(c['cluster_aa'][0])) + ' ' +
                'seedlineage_aa ' + relpath(str(c['seedlineage_aa'][0])) + ' ' +
                'newick ' + dnaml_tgt_relpath(3) + ' ')
    # Note; we used to delete the 'file' attribute as well; not sure why or if that's necessary
    c['metadata'].append(tgt)
    return tgt


# Popping out
# -----------

# Here we pop out to the "seed" level so we can aggregate our metadata

w.pop('seed')


# Filtering out bad clusters:

# As mentioned above, dnaml2tree fails for some clusters, so we'd like to filter these clusters out of the final metadata results.
# We do this based on whether the svg targets of the dnaml2tree command point to actual files or not.
# Eventually it would be nice to let them pass through with information indicating there was a failure, and handle appropriately in cftweb.
# We may also want to handle clusters that are too small similarly, but for now we're filtering them out at the very beginning of the pipeline.

# For now, to the stated end, we have this filter predicate
def tgt_exists(tgt):
    p = str(tgt)
    exists = path.exists(p)
    if not exists:
        print("Path doesn't exist:", p)
    return exists

def in_pairs(xs):
    it = iter(xs)
    for x in it:
        yield x, next(it)

def node_metadata(node):
    with open(str(node), 'r') as handle:
        return json.load(handle)

def git(*args):
    return subprocess.check_output(['git'] + list(args))

def write_metadata(target, source, env):
    # Here's where we filter out the seeds with clusters that didn't compute through dnaml2tree.py
    good_clusters = map(lambda x: node_metadata(x[1]), filter(lambda x: tgt_exists(x[0]), in_pairs(source)))
    metadata = {'clusters': good_clusters,
                'build_info': {'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                               'command': " ".join(sys.argv),
                               'workdir': os.getcwd(),
                               'datapath': datapath,
                               'user': getpass.getuser(),
                               'commit': git('rev-parse', 'HEAD'),
                               'status': git('status', '--porcelain')}}
    with open(str(target[0]), "w") as fh:
        json.dump(metadata, fh, sort_keys=True,
                       indent=4, separators=(',', ': '))

@w.add_target()
def metadata(outdir, c):
    tgt = env.Command(
        path.join(outdir, "metadata.json"),
        zip(c['svgfiles'], c['metadata']),
        write_metadata)
    env.AlwaysBuild(tgt)
    return tgt

import textwrap

def print_hints(target, source, env):
    msg = """\
		hint: to run the cft web interface,
			$ cd cftweb && python -m cftweb --file {}
	""".format(os.path.abspath(str(source[0])))
    print(textwrap.dedent(msg))

@w.add_target()
def hints(outdir, c):
    hints = env.Command(
        None,
        c['metadata'],
        print_hints)
    env.AlwaysBuild(hints)
    return hints


