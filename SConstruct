#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Carry partis output through prepation for display via cftweb application.

This SConstruct does the following:

* Runs process_partis.py for each seed/sequencing-run/parition combination, producing a metadata file of information
* For each such parition, the cluster with the seed is analyzed:
    * Tree construction using `FastTree`
    * Ancestral state reconstruction using `dnapars` or `dnaml`
    * Some SVG visualization
* Results of all these analyses are then pointed to by a merged metadata file, which is then consumed by cftweb

Clusters with only two sequences cannot be analyzed by dnapars, dnaml, or FastTree, and so these are skipped atuomatically.
Additionally, some clusters end up with bad trees after a rerooting step, and also have trouble in our dnaml2tree.py/dnapars.py step.
These are left out of the final `metadata.json` file.
Eventually, all of these cases should be caught and tombstone information should be left in the metadata file capturing as much infomration as possible.

See README for typical environment setup and usage.
'''


from __future__ import print_function
import os
import re
import sys
import csv
import time
import subprocess
import glob
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

# Setting up command line arguments/options

AddOption('--datapath',
              dest='datapath',
              type='string',
              nargs=1,
              action='store',
              metavar='DIR',
              default="/fh/fast/matsen_e/processed-data/partis/kate-qrs-2016-09-09/latest",
              help='partis output directory')

AddOption('--outdir',
        dest='outdir',
        type='string',
        nargs=1,
        action='store',
        metavar='DIR',
        help="directory in which to output results; defaults to `output/<dataset-tag>-<datapath-basename>-<asr-prog>")

AddOption('--asr-prog',
        dest='asr_prog',
        type='string',
        nargs=1,
        action='store',
        metavar='dnapars/dnaml',
        default='dnapars',
        help="dnapars (default) or dnaml tree program")

AddOption('--test',
        dest='test_run',
        action='store_true',
        default=False,
        help="Setting this flag does a test run for just a couple of seeds")

AddOption('--dataset-id',
        dest='dataset_id',
        metavar='IDENTIFIER',
        help="""Assign a dataset identitifier for distinguishing between multiple datasets loaded into cftweb;
                Defaults to [<dataset-tag>-]<datapath-basename>-<asr-prog>-<date>""")

AddOption('--dataset-tag',
        dest='dataset_tag',
        metavar='TAG',
        help="Adds a tag to the beginning of the automatically generated dataset id")


# prefer realpath so that running latest vs explicit vN doesn't require rerun; also need for defaults below
datapath = path.realpath(env.GetOption('datapath'))
asr_prog = env.GetOption('asr_prog')
test_run = env.GetOption("test_run")
dataset_tag = env.GetOption('dataset_tag') or path.basename(path.dirname(datapath))[:-11] + '-' + path.basename(datapath)
dataset_base = dataset_tag + '-' + asr_prog
# we call this outdir_base in order to not conflict with nestly fns outdir arg
outdir_base = env.GetOption('outdir') or path.join('output', dataset_base)
dataset_id = env.GetOption('dataset_id') or dataset_base + '-' + time.strftime('%Y.%m.%d')


# pruning is dependent on which program we use, dnapars seems to handle bigger trees more quickly
if asr_prog == 'dnapars':
    prune_n = 300
elif asr_prog == 'dnaml':
    prune_n = 100

print("datapath = {}".format(datapath))
print("outdir = {}".format(outdir_base))
print("asr_prog = {}".format(asr_prog))
print("test_run = {}".format(test_run))
print('pruned size = {}'.format(prune_n))
print('dataset id = {}'.format(dataset_id))


# Some environment sanity checks to make sure we have all prerequisits

def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# only check for dependencies if we aren't in a dry-run.
if not env.GetOption('no_exec'):
    msg = ""
    if not cmd_exists(asr_prog):
        msg += '''
           Required dependency command,
           `'''+asr_prog+'''` not found on PATH
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

#test_seeds = ["QB850.049-Vh", "QB850.043-Vh"]
test_seeds = ["QB850.405-Vh", "QB850.430-Vh", "QB850.091-Vh"]

def seeds_fn(datapath):
    return os.listdir(path.join(datapath, 'seeds'))

seeds = test_seeds if test_run else seeds_fn(datapath)

# Initialize our first nest level
w.add('seed', seeds)


# Initialize parameter set nest
# -----------------------------

# Next we nest on parameter sets; I believe these correspond to sequencing runs.
# They're the things that look like "Hs-LN2-5RACE-IgG", and are nested within each seed's output directory.

merged_timepoints_glob = "Q*-*-Ig*"
w.add('parameter_set',
      lambda c: map(lambda x: path.basename(x),
                    glob.glob(path.join(datapath, 'seeds', c['seed'], merged_timepoints_glob))))


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

def partis_log(c):
    return path.join(input_dir(c), 'partition.log')

@w.add_target()
def process_partis_(outdir, c):
    # Should get this to explicitly depend on cluster0.fa
    return env.Command(
            [path.join(outdir, x) for x in ['metadata.json', 'cluster0.fa', 'base_seqmeta.csv']],
            [partitions(c), annotation(c)],
            'process_partis.py -F ' +
                '--partition ${SOURCES[0]} ' +
                '--annotations ${SOURCES[1]} ' +
                #'--param_dir ' + parameter_dir(c) + ' '
                '--partis_log ' + partis_log(c) + ' ' +
                '--cluster_base cluster ' +
                '--melted_base base_seqmeta ' +
                '--output_dir ' + outdir + ' ' +
                '--paths-relative-to ' + outdir_base)

@w.add_target()
def base_metadata(outdir, c):
    return c['process_partis_'][0]

@w.add_target()
def inseqs(outdir, c):
    return c['process_partis_'][1]

@w.add_target()
def base_seqmeta(outdir, c):
    return c['process_partis_'][2]




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
        #path.join(outdir, "translated_inseqs.fa"),
        #c['inseqs'],
        #'seqmagick convert --translate dna2protein $SOURCE $TARGET')
        [path.join(outdir, "translated_inseqs.fa"), path.join(outdir, 'inseqs_trimmed.fa')],
        [c['base_metadata'], c['inseqs']],
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
        path.join(outdir, "aligned_translated_inseqs.fa"),
        c['translated_inseqs'],
        # Replace stop codons with X, or muscle inserts gaps, which messed up seqmagick backtrans-align below
        # Note that things will break down at backtrans-align if any seq ids have * in them...
        'sed \'s/\*/X/g\' $SOURCE | muscle -in /dev/stdin -out $TARGET')

# Not sure why we're not using this any more; I think we want to be, but it must have been causing some
# issue... need to look at this XXX
# This script will replace the X characters in our alignment with stop codons where they were taken out in the
# sed step prior to alignment above. We should try to use these when possible in alignments we display.
# However, this may be less useful here than it is further down where we translate our ancestral state
# inferences, since those are the seqs we use elsewhere.

#@w.add_target()
#def fixed_aligned_translated_inseqs(outdir, c):
    #return env.Command(
        #path.join(outdir, 'fixed_aligned_translated_inseqs.fa'),
        #[c['aligned_translated_inseqs'], c['translated_inseqs']],
        #'fix_stop_deletions.py $SOURCES $TARGET')

# Sort the sequences to have the same order, so that seqmagick doesn't freak out

@w.add_target()
def sorted_inseqs(outdir, c):
    return env.Command(
        path.join(outdir, "sorted_inseqs.fa"),
        c['trimmed_inseqs'],
        'seqmagick convert --sort name-asc $SOURCE $TARGET')

@w.add_target()
def sorted_aligned_translated_inseqs(outdir, c):
    return env.Command(
        path.join(outdir, "sorted_aligned_translated_inseqs.fa"),
        #c['fixed_aligned_translated_inseqs'],
        c['aligned_translated_inseqs'],
        'seqmagick convert --sort name-asc $SOURCE $TARGET')

# Now, finally, we can backtranslate

@w.add_target()
def aligned_inseqs(outdir, c):
    return env.Command(
        path.join(outdir, 'aligned_inseqs.fa'),
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
        "python bin/prune.py --n " + str(prune_n) + " --seed " + c['seed'] + " $SOURCE > $TARGET")

# prune out sequences to reduce taxa, making sure to cut out columns in the alignment that are now entirely
# gaps from insertions in sequences that have been pruned out.
@w.add_target()
def pruned_seqs(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.fa"),
        [c['pruned_ids'], c['aligned_inseqs']],
        "seqmagick convert --include-from-file $SOURCES - | " +
        "seqmagick convert --squeeze - $TARGET")

infname_regex = re.compile('.*--infname\s+(?P<infname_base>\S+)\.fa')
@w.add_target()
def merge_translations(outdir, c):
    "Returns the filename of the translation file for the merge operation, which includes timepoint info"
    with open(partis_log(c), 'r') as fh:
        partis_command = fh.readline()
    infname_base = infname_regex.match(partis_command).group('infname_base')
    return infname_base + '-translations.csv'

@w.add_target()
def timepoint_mapping(outdir, c):
    "Subset the merge_translations mapping to just the sequence name and the timepoint for the pruned seqs"
    return env.Command(
        path.join(outdir, "timepoint_mapping.csv"),
        [c['pruned_ids'], c['merge_translations']],
        "csvgrep -c new -f $SOURCES | sed '1 s/new/sequence/' | csvcut -c sequence,timepoint > $TARGET")

@w.add_target()
def seqmeta(outdir, c):
    "Merge the base_seqmeta from process_partis with the timepoint mappings"
    return env.Command(
        path.join(outdir, 'seqmeta.csv'),
        [c['base_seqmeta'], c['timepoint_mapping']],
        'csvjoin -c unique_ids,sequence $SOURCES > $TARGET')

# Convert to phylip format for dnapars/dnaml
@w.add_target()
def phy(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.phy"),
        c['pruned_seqs'],
        "seqmagick convert $SOURCE $TARGET")

# Create a config file for dnapars/dnaml (a persnickety old program with interactive menues...)
@w.add_target()
def asr_config(outdir, c):
    return env.Command(
        path.join(outdir, asr_prog + ".cfg"),
        c['phy'],
        'python bin/mkconfig.py $SOURCE ' + asr_prog + ' > $TARGET')

# Run dnapars/dnaml by passing in the "config file" as stdin hoping the menues all stay sane
# (Aside: If gets any messier can look at Expect; https://en.wikipedia.org/wiki/Expect)
@w.add_target()
def asr(outdir, c):
    "run dnapars/dnaml (from phylip package) to create tree with inferred sequences at internal nodes"
    tgt = env.SRun(
        path.join(outdir, "outfile"),
        c['asr_config'],
        'cd ' + outdir + ' && rm -f outtree && ' + asr_prog + ' < $SOURCE.file > ' + asr_prog + '.log')
    # Manually depend on phy so that we rerun dnapars/dnaml if the input sequences change (without this, dnapars/dnaml will
    # only get rerun if one of the targets are removed or if the iput asr_config file is changed).
    env.Depends(tgt, c['phy'])
    return tgt


@w.add_target()
def processed_asr(outdir, c):
    """parse dnapars/dnaml output into fasta and newick files, and make SVG format tree with ETE package.
    xvfb-run is needed because of issue https://github.com/etetoolkit/ete/issues/101"""
    basename = 'asr'
    tgt = env.Command(
            [path.join(outdir, basename + '.' + ext) for ext in ['nwk', 'svg', 'fa']],
            [c['asr'], c['seqmeta']],
            "xvfb-run -a bin/process_asr.py --seed " + c['seed'] + " --outdir " + outdir + 
                " --basename " + basename + " $SOURCES")
    # Manually depend on dnaml2tree.py/dnapars.py script, since it doesn't fall in the first position within the command
    # string.
    env.Depends(tgt, 'bin/process_asr.py')
    # Do aggregate work
    c['svgfiles'].append(tgt[0])
    return tgt

@w.add_target()
def asr_tree(outdir, c):
    return c['processed_asr'][0]

@w.add_target()
def asr_tree_svg(outdir, c):
    return c['processed_asr'][1]

@w.add_target()
def asr_seqs(outdir, c):
    return c['processed_asr'][2]


# Might want to switch back to this approach...
#def extended_metadata(metadata, dnapars_tree_tgt):
    #"Add dnapars_tree target(s) as metadata to the given metadata dict; used prior to metadata write."
    #with open(
    #m = copy.copy(metadata)
    #m['svg'] = path.relpath(str(dnapars_tree_tgt[0]), outdir_base)
    #m['fasta'] = path.relpath(str(dnapars_tree_tgt[1]), outdir_base)
    #m['newick'] = path.relpath(str(dnapars_tree_tgt[2]), outdir_base)
    #del m['file']
    #return m

@w.add_target()
def cluster_aa(outdir, c):
    return env.Command(
        path.join(outdir, 'cluster_aa.fa'),
        c['asr_seqs'],
        "sed -i 's/\?/N/g' $SOURCE && seqmagick convert --translate dna2protein $SOURCE $TARGET")

@w.add_target()
def cluster_metadata(outdir, c):
    # Note: We have to do all this crazy relpath crap because the assumption of cftweb is that all paths
    # specified in the input metadata.json file are relative to _its_ location, in our context, always base_outdir
    def relpath(tgt_key):
        tgt = c[tgt_key]
        if type(tgt) == list:
            tgt = tgt[0]
        return path.relpath(str(tgt), outdir_base)
    n = re.compile('run-viterbi-best-plus-(?P<step>.*)').match(c['partition']).group('step')
    tgt = env.Command(
            path.join(outdir, 'extended_metadata.json'),
            [c['base_metadata'], partitions(c)] + c['processed_asr'],
            # Don't like the order assumptions on dnaml_tgts here...
            'assoc_logprob.py $SOURCE ${SOURCES[1]} -n ' + n + ' /dev/stdout | ' +
                'json_assoc.py /dev/stdin $TARGET ' +
                # Not 100% on this; Pick optimal attr/key name
                #'best_partition ' + c['partition'] + ' ' +
                ' dataset_id ' + dataset_id +
                ' clustering_step ' + n +
                ' svg ' + relpath('asr_tree_svg') +
                ' fasta ' + relpath('asr_seqs') +
                ' cluster_aa ' + relpath('cluster_aa') +
                ' newick ' + relpath('asr_tree'))
    # Note; we used to delete the 'file' attribute as well; not sure why or if that's necessary
    c['metadata'].append(tgt)
    return tgt


# Popping out
# -----------

# Here we pop out to the "seed" level so we can aggregate our metadata

w.pop('seed')


# Filtering out bad clusters:

# As mentioned above, dnaml2tree.py/dnapars.py fails for some clusters, so we'd like to filter these clusters out of the final metadata results.
# We do this based on whether the svg targets of the dnaml2tree.py/dnapars.py command point to actual files or not.
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
    # Here's where we filter out the seeds with clusters that didn't compute through dnaml2tree.py/dnapars.py
    good_clusters = map(lambda x: node_metadata(x[1]), filter(lambda x: tgt_exists(x[0]), in_pairs(source)))
    metadata = {'clusters': good_clusters,
                'dataset_id': dataset_id,
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
			$ cd cftweb && python -m cftweb {}
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
