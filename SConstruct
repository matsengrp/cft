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
import datetime
import getpass
import glob
import sconsutils
import itertools
import copy
#import functools as fun


from os import path
#from warnings import warn


#from nestly import Nest
#from nestly.scons import SConsWrap
from nestly.nestly import Nest
from nestly.nestly.scons import SConsWrap
from datascripts import heads


#from SCons import Node
from SCons.Script import Environment, AddOption

#import json

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

AddOption('--datapaths',
        dest='datapaths',
        metavar='DIR_LIST',
        default="laura-mb/latest:kate-qrs/latest",
        help="""Specify ':' separated list of partis output directories to process on; if full path not specified,
        assumed to be in --base-datapath. Dataset names will be assined in relation to this --base-datapath
        if present. Note: symlinks may not work properly here unless they point to things in base-datapath also.""")

AddOption('--base-datapath',
        dest='base_datapath',
        metavar='DIR',
        default="/fh/fast/matsen_e/processed-data/partis/",
        help="""Location in which to find the --datapaths directories. Defauilts to %default%""")

AddOption('--asr-progs',
        dest='asr_progs',
        metavar='LIST',
        #default='dnaml:dnapars:raxml',
        default='dnaml',
        help="""Specify ':' separated list of ancestral state reconstruction programs to run. Options are `dnaml` and
        `dnapars`. Defaults to running dnaml.""")

AddOption('--prune-strategies',
        dest='prune_strategies',
        metavar='LIST',
        default='min_adcl:seed_lineage',
        help="""Specify ':' separated list of pruning strategies. Options are 'min_adcl' and 'seed_lineage'.
        Defaults to both.""")

AddOption('--separate-timepoints',
        dest='separate_timepoints',
        action='store_true',
        help='Include separate timepoints (vs just running on unmerged)?')

AddOption('--test',
        dest='test_run',
        action='store_true',
        default=False,
        help="Setting this flag does a test run for just a couple of seeds")

AddOption('--dataset-tag',
        dest='dataset_tag',
        metavar='TAG',
        help="Adds a tag to the beginning of the automatically generated dataset ids")

AddOption('--outdir`',
        dest='outdir',
        metavar='DIR',
        default='output',
        help="Directory in which to output results; defaults to `output`")


# prefer realpath so that running latest vs explicit vN doesn't require rerun; also need for defaults below
base_datapath = env.GetOption('base_datapath')
datapaths = map(path.realpath,
                [path.join(base_datapath, x) for x in env.GetOption('datapaths').split(':')])
asr_progs = env.GetOption('asr_progs').split(':')
prune_strategies = env.GetOption('prune_strategies').split(':')
test_run = env.GetOption("test_run")
separate_timepoints = env.GetOption("separate_timepoints")
dataset_tag = env.GetOption('dataset_tag') or ('test' if test_run else None)
outdir_base = env.GetOption('outdir')


print("outdir = {}".format(outdir_base))
print("test_run = {}".format(test_run))



# Utility functions
# -----------------

def merge_dicts(d1, d2):
    d = d1.deepcopy(d1)
    d.update(d2)
    return d


# Initialize nestly!
# ==================

# This lets us create parameter nests and run the given pipeline for each combination of nesting parameters.
# It also lets us "pop" out of nesting levels in order to aggregate on nested results.

# Our nesting is going to more or less model the nesting of things in our datapath directory.
# seed > sample > partition

# Here we initialize nestly, and create an scons wrapper for it

build_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def git(*args):
    return subprocess.check_output(['git'] + list(args))

import tripl.tripl.nestly as nestly_tripl

nest = Nest()
w = SConsWrap(nest, outdir_base, alias_environment=env)
w = nestly_tripl.NestWrap(w,
        name='build',
        # Need to base hashing off of this for optimal incrementalization
        metadata={'id': 'cft-build-' + build_time.replace(' ', '-'),
                  'time': build_time,
                  'command': " ".join(sys.argv),
                  'workdir': os.getcwd(),
                  'base_datapath': base_datapath,
                  'user': getpass.getuser(),
                  'commit': git('rev-parse', 'HEAD'),
                  # This will be really cool :-)
                  'diff': git('diff'),
                  'status': git('status', '--porcelain')},
        base_namespace='cft',
        id_attrs=['cft.dataset:id', 'cft.build:id'])


software = {
    'dnaml': None,
    'muscle': 'muscle -version',
    'seqmagick': 'seqmagick --version',
    'FastTree': None,
    'prank': 'prank -v'
    }


def software_info(prog):
    version_command = software[prog]
    return {'cft.software:name': prog,
            'cft.software:version': subprocess.check_output(version_command.split()) if version_command else None,
            'cft.software:which': subprocess.check_output(['which', prog])}


@w.add_target('cft.build:software')
def software(outdir, c):
    return [software_info(prog) for prog in software]


# Could do the metadata as a separate target?
#@w.add_target()
#def build_data(outdir, c):



# Dataset nest level
# =================

def _dataset_outdir(dataset_params):
    "Outputs the basename of the path to which this dataset's output lives; e.g. 'kate-qrs-v10-dnaml'"
    base = dataset_tag + '-' if dataset_tag else ''
    base += path.relpath(dataset_params['datapath'], base_datapath).replace('/', '-')
    if dataset_params['prune_strategy'] == 'min_adcl':
        base += '-' + 'minadcl'
    if dataset_params['separate_timepoints']:
        base += '-' + 'septmpts'
    return base + '-' + dataset_params['asr_prog']


def _dataset_id(outdir):
    """The return dataset_id corresponding to the control dictionary; This is the key under which datasets are
    organized in CFTWeb."""
    return outdir + '-' + time.strftime('%Y.%m.%d')


def with_data_id_and_outdir(dataset_params):
    outdir = _dataset_outdir(dataset_params)
    d = {'outdir': outdir,
         'id': _dataset_id(outdir)}
    d.update(dataset_params)
    return d


print("Running for")
print("  datapaths:", datapaths)
print("  asr_progs:", asr_progs)
print("  prune_strategies:", prune_strategies)

# A collection of datasets, where the `datapath` key is the full realpath to a leaf node input directory
datasets = [with_data_id_and_outdir({
                 'datapath': datapath,
                 'study': datapath.split('/')[-2],
                 'version': datapath.split('/')[-1],
                 'asr_prog': asr_prog,
                 'prune_strategy': prune_strategy,
                 'separate_timepoints': separate_timepoints,
                 })
             for datapath, asr_prog, prune_strategy
             in itertools.product(datapaths, asr_progs, prune_strategies)]

# This is a little idiosynchratic the way we're doing things here, because of how we wanted to think about
# uniqueness before; This will probably get scrapped and we'll have everything build off of attributes of the
# cluster, and let tripl take care of identity and uniquness. We may not be able to preserve uniquness though,
# unless we specify a route mapping via a pull expression.

w.add('dataset', datasets,
        full_dump=True,
        label_func=lambda d: d['outdir'],
        metadata=lambda c, d: d,
        id_attrs=['cft.subject:id', 'cft.sample:id'])

# Helpers for accessing info about the dataset

def dataset_outdir(c):
    "Returns _dataset_outdir of `c['dataset']`, for easier access via `c` below."
    return c['dataset']['outdir']

@w.add_target('cft.dataset:seqmeta')
def dataset_seqmeta(outdir, c):
    "The sequence metadata file path for the entire dataset."
    return env.Command(
        path.join(outdir, 'seqmeta.csv'),
        glob.glob(path.join(base_datapath, c['dataset']['datapath'], '*-translations.csv')),
        'csvstack $SOURCES > $TARGET')

def dataset_id(c):
    return c['dataset']['id']

def asr_prog(c):
    "Retrieve the asr_prog from the dataset map"
    return c['dataset']['asr_prog']

def prune_strategy(c):
    "Retrieve the asr_prog from the dataset map"
    return c['dataset']['prune_strategy']

def datapath(c):
    "Returns the input realpath datapath value of the dataset"
    return c['dataset']['datapath']



def wrap_test_run(take_n=2):
    def deco(nestables_fn):
        def f(c):
            nestables = nestables_fn(c)
            nestables = nestables[:take_n]
            return nestables
        f.__name__ = nestables_fn.__name__
        return f if test_run else nestables_fn
    return deco


@w.add_nest()
@wrap_test_run(take_n=2)
def subject(c):
    study = c['dataset']['study']
    return list(set(d['subject'] for d in heads.read_metadata(study).values()))


# Initialize seed nest
# --------------------

# The very first nesting is on the seed id.

# For the sake of testing, we allow for switching between the full set of seeds, as returned by `seeds_fn`, or
# just a small subsampling thereof via execution with the `--test` cli flag

@w.add_target()
def _seeds(outdir, c):
    study = c['dataset']['study']
    subject = c['subject']
    seeds = heads.get_seeds(study, subject)
    # If we don't have a directory for a seed, it might not be done running yet, so just run on all complete
    seeds = filter(lambda x: path.isdir(path.join(datapath(c), 'seeds', x)), seeds)
    return seeds

# Initialize our first sub dataset nest level
@w.add_nest('seed')
# would like to have a lower number here but sometimes we get no good clusters for the first two seeds?
# (on laura-mb for example).
@wrap_test_run(take_n=4)
def seeds_fn(c):
    return c['_seeds']


# Initialize parameter set nest
# -----------------------------

# Next we nest on parameter sets; I believe these correspond to sequencing runs.
# They're the things that look like "Hs-LN2-5RACE-IgG", and are nested within each seed's output directory.
#
# What this should eventually look like is that we have already specified the relationships between the data
# and the directories (timepoints etc) upstream, so we don't have to muck around with this.

is_merged = re.compile('[A-Z]+\d+-\w-Ig[A-Z]').match
is_unmerged = re.compile('Hs-(LN-?\w+)-.*').match

def timepoint(c, sample_filename):
    study = c['dataset']['study']
    d = heads.read_metadata(study)[sample_filename]
    return d['timepoint']


def sample_metadata(c, filename): # control dict as well?
    study = c['dataset']['study']
    d = copy.deepcopy(heads.read_metadata(study)[filename])
    d.update(
        {'id': filename,
         'timepoints': [{'cft.timepoint:id': timepoint(c, filename)}]})
    return d

#dataset,shorthand,species,timepoint,subject,locus
@w.add_nest(metadata=sample_metadata)
@wrap_test_run()
def sample(c):
    def keep(filename):
        return is_unmerged(filename) if separate_timepoints else is_merged(filename)
    return filter(keep,
                  os.listdir(path.join(datapath(c), 'seeds', c['seed'])))


# Some helpers at the seed level

def input_dir(c):
    "Return the `seed > sample` directory given the closed over datapath and a nestly control dictionary"
    # This 'seeds' thing here is potentially not a very general assumption; not sure how variable that might be upstream
    return path.join(datapath(c), 'seeds', c['seed'], c['sample'])

def path_base_root(full_path):
    return path.splitext(path.basename(full_path))[0]

# If there are duplicates 
cluster_step_re = re.compile('.*-plus-(?P<cluster_step>\d+)').match
def cluster_step(partition_filename):
    return int(cluster_step_re(partition_filename).group('cluster_step'))

def partition_size(fname):
    with open(fname, 'r') as partition_handle:
        partition = csv.DictReader(partition_handle).next()
        unique_ids = partition['unique_ids'].split(':')
        # add 1 for naive
        return len(unique_ids) + 1


# Setting up the partition nest level

def partitions(c):
    "Returns the `partition.csv` file path with all the partition information for every partition output by partis."
    return path.join(input_dir(c), 'partition.csv')

def parameter_dir(c):
    "The input parameter directory (as in input to partis); see --parameter-dir flag to various partis cli calls."
    return path.join(datapath(c), c['sample'])

def annotation(c):
    return path.join(input_dir(c), c['partition'] + '.csv')

def partis_log(c):
    return path.join(input_dir(c), 'partition.log')

def translations(c):
    return path.join(input_dir(c), 'partition.log')

def partition_file_metadata(partition_handle, cluster_step):
    parts = list(csv.DictReader(partition_handle))
    for i, x in enumerate(parts):
        x['index'] = i
        x['n_clusters'] = int(x['n_clusters'])
        x['logprob'] = float(x['logprob'])
    best_part = max(parts, key=lambda x: x['logprob'])
    best_i = best_part['index']
    cluster_i = best_i + cluster_step
    return parts[cluster_i]

def partition_metadata(c, filename):
    step = cluster_step(filename)
    metadata = partition_file_metadata(file(partitions(c)), step)
    return {'cluster_step': step,
            'logprob': metadata['logprob'],
            'n_clusters': metadata['n_clusters']}

# Each c["partition"] value actually points to the annotations for that partition... a little weird but...
@w.add_nest(metadata=partition_metadata)
@wrap_test_run()
def partition(c):
    """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
    for actual analysis."""
    partitions = sorted(glob.glob(path.join(input_dir(c), "*-plus-*.csv")), key=cluster_step)
    keep_partitions = []
    for partition in partitions:
        size = partition_size(partition)
        # We only add clusters bigger than two, since we can only make trees if we have hits
        if size > 2:
            keep_partitions.append(partition)
        # Once we get a cluster of size 50, we don't need later cluster steps
        if size >= 50:
            break
    return map(path_base_root, keep_partitions)

# This is a little silly, but gives us the right semantics for partitions > clusters
#w.add('cluster', ['cluster0'], metadata=lambda _, cluster_id: {'id': cluster_id}) # set true
w.add('cluster', ['cluster0'])

@w.add_metadata()
def _process_partis(outdir, c):
    # Should get this to explicitly depend on cluster0.fa
    return env.Command(
            [path.join(outdir, x) for x in ['partis_metadata.json', 'cluster0.fa', 'partis_seqmeta.csv']],
            [partitions(c), annotation(c)],
            'process_partis.py -F ' +
                '--partition ${SOURCES[0]} ' +
                '--annotations ${SOURCES[1]} ' +
                #'--param_dir ' + parameter_dir(c) + ' '
                '--remove-frameshifts ' +
                '--partis-log ' + partis_log(c) + ' ' +
                '--cluster-base cluster ' +
                '--melted-base partis_seqmeta ' +
                '--output-dir ' + outdir + ' ' +
                '--paths-relative-to ' + dataset_outdir(c))

@w.add_target(ingest=True)
def partis_metadata(outdir, c):
    return c['_process_partis'][0]

@w.add_target()
def inseqs(outdir, c):
    return c['_process_partis'][1]

@w.add_target()
def partis_seqmeta(outdir, c):
    return c['_process_partis'][2]




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
        [c['partis_metadata'], c['inseqs']],
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
        'sed \'s/\*/X/g\' $SOURCE | muscle -in /dev/stdin -out $TARGET 2> $TARGET-.log')

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
        "FastTree -nt -quiet $SOURCE > $TARGET 2> $TARGET-.log")

# pruning is dependent on which program we use, dnapars seems to handle bigger trees more quickly
def prune_n(c):
    prune_n_ = {'dnapars': 300, 'dnaml': 100, 'raxml': 100}
    return prune_n_[asr_prog(c)]

# calculate list of sequences to be pruned
@w.add_target()
def pruned_ids(outdir, c):
    strategy = prune_strategy(c)
    return env.Command(
        path.join(outdir, "pruned_ids.txt"),
        c['fasttree'],
        "prune.py -n " + str(prune_n(c)) + " --always-include " + ','.join(c['_seeds']) + " --strategy " + strategy + " --naive naive0 --seed " + c['seed'] + " $SOURCE > $TARGET")


@w.add_target()
def cluster_mapping(outdir, c):
    if prune_strategy(c) == 'min_adcl':
        return env.Command(
            path.join(outdir, 'cluster_mapping.csv'),
            [c['fasttree'], c['pruned_ids']],
            'minadcl_clusters.py $SOURCES $TARGET')



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
def infname_base(outdir, c):
    "Returns the filename of the translation file for the merge operation, which includes timepoint info"
    with open(partis_log(c), 'r') as fh:
        partis_command = fh.readline()
    # Hacky temp fix for the dataset rename 2017/04/03
    partis_command = partis_command.replace('kate-qrs-2016-09-09', 'kate-qrs')
    partis_command = partis_command.replace('laura-mb-2016-12-22', 'laura-mb')
    infname_base = infname_regex.match(partis_command).group('infname_base')
    return infname_base

@w.add_target(ingest=True, attr_map={'bio.seq:id': 'sequence', 'cft.timepoint:id': 'timepoint',
    'cft.seq:duplicity': 'duplicity'})
    #'cft.seq:duplicity': 'duplicity', 'cft.seq:cluster_duplicity': 'cluster_duplicity'})
def seqmeta(outdir, c):
    """The merge of process_partis output with pre sequence metadata spit out by datascripts containing
    timepoint mappings. Base input duplicity is coded into the original input sequence names from vlad as N-M,
    where N is the ranking of vlads untrimmed deduplication, and M is the duplicity of said deduplication."""
    sources = [c['partis_seqmeta'], c['cft.dataset:seqmeta']]
    base_call =  'merge_timepoints_and_duplicity.py '
    # This option controls which sequences get joined on in the merge for the partis_seqmeta file, which has
    # orig/new names, joined on sequence from the other file
    if separate_timepoints:
        base_call += '--timepoint ' + timepoint(c, c['sample']) + ' '
    if prune_strategy(c) == 'min_adcl':
        sources = [c['cluster_mapping']] + sources
        base_call += '--cluster-mapping '
    return env.Command(path.join(outdir, 'seqmeta.csv'), sources, base_call + '$SOURCES $TARGET')

@w.add_target()
def full_orig_seqs(outdir, c):
    if c['dataset']['study'] == 'kate-qrs':
        return c['infname_base'] + '.fa'
    elif c['dataset']['study'] == 'laura-mb':
        return c['infname_base'] + '.fasta'

@w.add_target(ingest=True)
def orig_seqs(outdir, c):
    "The original input seqs to partis, including the non-VDJ, trimmed out regions"
    return env.Command(
        path.join(outdir, 'orig_seqs.fa'),
        [c['pruned_ids'], c['full_orig_seqs']],
        "seqmagick convert --include-from-file $SOURCES $TARGET")

# Convert to phylip format for dnapars/dnaml
@w.add_target()
def phy(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.phy"),
        c['pruned_seqs'],
        "seqmagick convert $SOURCE $TARGET")


# Run dnapars/dnaml by passing in the "config file" as stdin hoping the menues all stay sane
# (Aside: If gets any messier can look at Expect; https://en.wikipedia.org/wiki/Expect)
@w.add_target()
def _asr_tree(outdir, c):
    "run dnapars/dnaml (from phylip package) to create tree with inferred sequences at internal nodes"
    if asr_prog(c) in {'dnapars', 'dnaml'}:
        config = env.Command(
            path.join(outdir, asr_prog(c) + ".cfg"),
            c['phy'],
            'python bin/mkconfig.py $SOURCE ' + asr_prog(c) + ' > $TARGET')
        phylip_out = env.SRun(
            path.join(outdir, "outfile"),
            config,
            'cd ' + outdir + ' && rm -f outtree && ' + asr_prog(c) + ' < $SOURCE.file > ' + asr_prog(c) + '.log',
            ignore_errors=True)
        # Manually depend on phy so that we rerun dnapars/dnaml if the input sequences change (without this, dnapars/dnaml will
        # only get rerun if one of the targets are removed or if the iput asr_config file is changed). IMPORTANT!
        env.Depends(phylip_out, c['phy'])
        # Now process the phylip output into something that isn't shit
        basename = 'asr_input'
        tgt = env.Command(
                [path.join(outdir, basename + '.' + ext) for ext in ['nwk', 'svg', 'fa']],
                [phylip_out, c['seqmeta']],
                # Note that `-` at the beggining lets things keep running if there's an error here; This is
                # protecting us at the moment from clusters with 2 seqs. We should be catching this further
                # upstream and handling more appropriately, but for now this is an easy stopgap...
                "xvfb-run -a bin/process_asr.py --seed " + c['seed'] + " --outdir " + outdir + 
                    " --basename " + basename + " $SOURCES")
        asr_tree, asr_tree_svg, asr_seqs = tgt
        # manually depnd on this because the script isn't in first position
        env.Depends(tgt, 'bin/process_asr.py')
        env.Depends(tgt, 'bin/plot_tree.py')
        return [asr_tree, asr_tree_svg, asr_seqs]
    elif asr_prog(c) == 'raxml':
        asr_supports_tree = env.SRun(
            path.join(outdir, 'asr.sup.nwk'),
            c['pruned_seqs'],
            # Question should use -T for threads? how many?
            # Don't know if the reroot will really do what we want here
            'raxml.py --rapid-bootstrap 30 -x 3243 -o naive0 $SOURCE $TARGET')
        asr_tree_svg = env.Command(
            path.join(outdir, 'asr.svg'),
            [asr_supports_tree, c['seqmeta']],
            'xvfb-run -a bin/plot_tree.py $SOURCES $TARGET --supports --seed ' + c['seed'])
        asr_tree = env.Command(
            path.join(outdir, 'asr.nwk'),
            asr_supports_tree,
            'name_internal_nodes.py $SOURCE $TARGET')
        return [asr_tree, asr_tree_svg, asr_supports_tree]
    else:
        print("something has gone terribly wrong")



#@w.add_target(ingest=True)
#def asr_input_tree(outdir, c):
    #return c['_asr_tree'][0]

@w.add_target(ingest=True)
def asr_tree_svg(outdir, c):
    return c['_asr_tree'][1]

#@w.add_target()
#def asr_supports_tree(outdir, c):
    #vals = c['_asr_tree']
    #if len(vals) > 2:
        #return c['_asr_tree'][2]

#@w.add_target(ingest=True)
#def _asr(outdir, c):
    #return env.Command(
        #[path.join(outdir, 'asr_tree.nwk'), path.join(outdir, 'asr_seqs.fa')],
        #[c['asr_input_tree'], c['pruned_seqs']],
        #'joker.py ' + ('--no-reroot ' if asr_prog(c) == 'raxml' else '') + '-t $SOURCES $TARGETS')

#@w.add_target(ingest=True)
#def asr_tree(outdir, c):
    #return c['_asr'][0]

#@w.add_target(ingest=True)
#def asr_seqs(outdir, c):
    #return c['_asr'][1]

@w.add_target(ingest=True)
def asr_tree(outdir, c):
    return c['_asr_tree'][0]

@w.add_target(ingest=True)
def asr_seqs(outdir, c):
    return c['_asr_tree'][2]


@w.add_target(ingest=True)
def cluster_aa(outdir, c):
    return env.Command(
        path.join(outdir, 'cluster_aa.fa'),
        c['asr_seqs'],
        "sed 's/\?/N/g' $SOURCE | seqmagick convert --translate dna2protein - $TARGET")



# Popping out
# -----------

# Here we pop out to the "seed" level so we can aggregate our metadata

# I don't know if we have to do this anymore
w.pop('seed')


# Go back to the base nest level, forcing a metadata write, and print help for cftweb

w.pop('dataset')


