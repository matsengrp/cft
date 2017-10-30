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


# Basic imports
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
#import json
#import functools as fun

from os import path
#from warnings import warn


# Nestly things
# this is temporary
#import nestly
#import nestly_scons
from nestly import nestly
from nestly.nestly import scons as nestly_scons

# Partis and datascripts things

# If the PARTIS env var isn't already set, default to $PWD/partis (where we have a git
# submodule checkout; this is needed for bin/process_partis.py)
default_partis_path = path.join(os.getcwd(), 'partis')
partis_path = os.environ.get('PARTIS', default_partis_path)
sys.path.append(path.join(partis_path, 'python'))
import clusterpath
from datascripts import heads

# Scons requirements
from SCons.Script import Environment, AddOption


# Build modules (in site_scons):
import backtrans_align




# Need this in order to read csv files with sequences in the fields
csv.field_size_limit(sys.maxsize)

# No-op; Prevents analysis warnings
sconsutils


# Set up SCons environment

environ = os.environ.copy()

# install partis path as env var if not already set
environ['PARTIS'] = partis_path

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

AddOption('--outdir',
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

nest = nestly.Nest()
w = nestly_scons.SConsWrap(nest, outdir_base, alias_environment=env)
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


# Recording software versions
# ---------------------------

# First import some libs we'll need versions for
import ete3, Bio, dendropy #, pandas
# tripl version is a little messy because sometimes we load from local checkout
def tripl_version():
    try:
        import tripl
        return tripl.__version__
    except:
        from tripl import tripl
        return tripl.__version__

# The contract here is that a string val mapped to here is a command string to get a software version. A
# function value is called to get a version. And as long as it's not a function value, then it's assumed it's
# assumed the program name is something that `which` can be called on. If function value, assumed to be a lib
# and which is not called. For now... this is a little arbitrary and specific to our use case here.
software = {
    'dnaml': None,
    'muscle': 'muscle -version',
    'seqmagick': 'seqmagick --version',
    'FastTree': None,
    'prank': 'prank -v',
    'tripl': tripl_version,
    'nestly': lambda: nestly.__version__,
    'ete3': lambda: ete3.__version__,
    'biopython': lambda: Bio.__version__,
    'scons': 'scons -v',
    'dendropy': lambda: dendropy.__version__,
    # For some reason pandas won't load here... so no version tracking for now
    #'pandas': lambda: pandas.__version__,
    # For minadcl
    'rppr': 'rppr --version'
    }


def software_info(prog):
    version_command = software[prog]
    return {'cft.software:name': prog,
            'cft.software:version': version_command() if callable(version_command) else (
                 subprocess.check_output(version_command.split()) if version_command else None),
            'cft.software:which': subprocess.check_output(['which', prog]) if not callable(version_command) else None}


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
    if dataset_params['separate_timepoints']:
        base += '-' + 'septmpts'
    return base


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

# A collection of datasets, where the `datapath` key is the full realpath to a leaf node input directory
datasets = [with_data_id_and_outdir({
                 'datapath': datapath,
                 'study': datapath.split('/')[-2],
                 'version': datapath.split('/')[-1],
                 'separate_timepoints': separate_timepoints,
                 })
             for datapath in datapaths]

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


def is_merged(c, x):
    if c['dataset']['study'] in {'kate-qrs', 'laura-mb'}:
        return re.compile('[A-Z]+\d+-\w-Ig[A-Z]').match(x)
    elif c['dataset']['study'] in {'laura-mb-2'}:
        return re.compile('\w+-\w+-merged').match(x)

def is_unmerged(c, x):
    if c['dataset']['study'] in {'kate-qrs', 'laura-mb'}:
        #return re.compile('Hs-(LN-?\w+)-.*').match(x)
        return re.compile('Hs-(LN-?\w+)(?:-5RACE)?-Ig\w').match(x)
    elif c['dataset']['study'] in {'laura-mb-2'}:
        return re.compile('\w+-\w+-(?!merged)\w+').match(x)


def canonical_sample_name(c, sample_filename):
    study_metadata = heads.read_metadata(c['dataset']['study'])
    sample_name = [k for k in study_metadata if k in sample_filename][0]
    return sample_name

def raw_sample_metadata(c, sample_filename):
    study_metadata = heads.read_metadata(c['dataset']['study'])
    return study_metadata[canonical_sample_name(c, sample_filename)]

def timepoint(c, sample_filename):
    return raw_sample_metadata(c, sample_filename)['timepoint']

def sample_metadata(c, filename): # control dict as well?
    d = copy.deepcopy(raw_sample_metadata(c, filename))
    d.update(
        {'id': filename,
         'timepoints': [{'cft.timepoint:id': timepoint(c, filename)}]})
    return d


# Initialize sample nest
# -----------------------------

# Next we nest on sample.
# These things look like "Hs-LN2-5RACE-IgG", and are nested within each seed's output directory.

#dataset,shorthand,species,timepoint,subject,locus
@w.add_nest(metadata=sample_metadata)
@wrap_test_run()
def sample(c):
    def keep(filename):
        return is_unmerged(c, filename) if separate_timepoints else is_merged(c, filename)
    results = filter(keep,
                  os.listdir(path.join(datapath(c), 'seeds', c['seed'])))
    #print('results', results)
    return results


# Some helpers at the seed level


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


def input_dir(c):
    """If seeded, return the `seed > sample` directory given the closed over datapath and a nestly control
    dictionary. Else return the `sample` directory given the closed over datapath and a nestly control dictionary."""
    # This 'seeds' thing here is potentially not a very general assumption; not sure how variable that might be upstream
    if 'seed' in c:
        return path.join(datapath(c), 'seeds', c['seed'], c['sample'])
    else:
        return path.join(datapath(c), 'partitions', c['sample'])


# Should we add target so it's just in the c dictionary?
#@w.add_target()
def path_base_root(full_path):
    return path.splitext(path.basename(full_path))[0]

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


# These functions define the "interface" of the nesting for a cluster analysis, if you will

def partis_log(c):
    return path.join(input_dir(c), 'partition.log')

def annotation(c):
    if 'seed' in c:
        return path.join(input_dir(c), c['partition'] + '.csv')
    else:
        return path.join(input_dir(c), 'partition-cluster-annotations.csv')

def partitions(c):
    "Returns the `partition.csv` file path with all the partition information for every partition output by partis."
    return path.join(input_dir(c), 'partition.csv')

# This one we're not using anymore; should delete, but for now.
#def parameter_dir(c):
    #"The input parameter directory (as in input to partis); see --parameter-dir flag to various partis cli calls."
    #return path.join(datapath(c), c['sample'])


def add_cluster_analysis(w):
    @w.add_metadata()
    def _process_partis(outdir, c):
        # Should get this to explicitly depend on cluster0.fa
        return env.Command(
                [path.join(outdir, x) for x in ['partis_metadata.json', 'cluster0.fa', 'partis_seqmeta.csv']],
                [partitions(c), annotation(c)],
                'process_partis.py -F' +
                    ' --partition ${SOURCES[0]}' +
                    ' --annotations ${SOURCES[1]}' +
                    #' --param_dir ' + parameter_dir(c) +
                    ' --remove-frameshifts' +
                    ' --partis-log ' + partis_log(c) +
                    ' --cluster-base cluster' +
                    ' --melted-base partis_seqmeta' +
                    ' --output-dir ' + outdir +
                    ' --paths-relative-to ' + dataset_outdir(c) +
                    ('' if c.get('seed') else ' --unique-ids ' + ':'.join(c['cluster']['unique_ids'])))

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
    # ------------------

    backtrans_align.add(env, w)


    # On with trees and other things...
    # ---------------------------------

    # use fasttree to make newick tree from sequences
    @w.add_target()
    def fasttree(outdir, c):
        return env.SRun(
            path.join(outdir, "fasttree.nwk"),
            c['aligned_inseqs'],
            "FastTree -nt -quiet $SOURCE > $TARGET 2> $TARGET-.log")

    @w.add_nest('reconstruction',
        label_func=lambda d: d['id'],
        metadata=lambda c, d: d)
    def reconstruction(c):
        return [{'id': prune_strategy + '-' + asr_prog,
                 'prune_strategy': prune_strategy,
                 'asr_prog': asr_prog,
                 # Just 100 for everyone now
                 'prune_count': 100}
                 for prune_strategy, asr_prog
                 in itertools.product(
                     ['min_adcl', 'seed_lineage'] if 'seed' in c else ['min_adcl'],
                     #['dnaml', 'ecgtheow']]
                     # ^ in the future?
                     ['dnaml'])]


    # calculate list of sequences to be pruned
    @w.add_target()
    def pruned_ids(outdir, c):
        recon = c['reconstruction']
        builder = env.SRun if recon['prune_strategy'] == 'min_adcl' else env.Command
        return builder(
            path.join(outdir, "pruned_ids.txt"),
            c['fasttree'],
            "prune.py -n " + str(recon['prune_count'])
                + ((" --always-include " + ','.join(c['_seeds'])) if c['_seeds'] else '')
                + " --strategy " + recon['prune_strategy']
                + " --naive naive0"
                + (" --seed " + c['seed'] if 'seed' in c else '')
                + " $SOURCE $TARGET")


    @w.add_target()
    def cluster_mapping(outdir, c):
        if c['reconstruction']['prune_strategy'] == 'min_adcl':
            return env.SRun(
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
        #'cft.seq:multiplicity': 'multiplicity'})
        'cft.seq:multiplicity': 'multiplicity', 'cft.seq:cluster_multiplicity': 'cluster_multiplicity'})
    def seqmeta(outdir, c):
        """The merge of process_partis output with pre sequence metadata spit out by datascripts containing
        timepoint mappings. Base input multiplicity is coded into the original input sequence names from vlad as N-M,
        where N is the ranking of vlads untrimmed deduplication, and M is the multiplicity of said deduplication."""
        sources = [c['partis_seqmeta'], c['cft.dataset:seqmeta']]
        base_call =  'merge_timepoints_and_multiplicity.py '
        # This option controls which sequences get joined on in the merge for the partis_seqmeta file, which has
        # orig/new names, joined on sequence from the other file
        if separate_timepoints:
            base_call += '--timepoint ' + timepoint(c, c['sample']) + ' '
        if c['reconstruction']['prune_strategy'] == 'min_adcl':
            sources = [c['cluster_mapping']] + sources
            base_call += '--cluster-mapping '
        return env.Command(path.join(outdir, 'seqmeta.csv'), sources, base_call + '$SOURCES $TARGET')

    @w.add_target()
    def _phy(outdir, c):
        "Save seqs in phylip format for dnaml, renaming sequences as necessary to avoid seqname limits"
        return env.Command(
            [path.join(outdir, x) for x in ('pruned.phy', 'seqname_mapping.csv')],
            c['pruned_seqs'],
            'make_phylip.py $SOURCE $TARGETS --dont-rename naive0')


    @w.add_target()
    def phy(outdir, c):
        return c['_phy'][0]

    @w.add_target()
    def seqname_mapping(outdir, c):
        """Seqname translations for reinterpretting dnaml output in terms of original seqnames, due to phylip name
        length constraints."""
        return c['_phy'][1]


    # Run dnapars/dnaml by passing in the "config file" as stdin hoping the menues all stay sane
    # (Aside: If gets any messier can look at Expect; https://en.wikipedia.org/wiki/Expect)
    @w.add_target()
    def _asr_tree(outdir, c):
        "run dnapars/dnaml (from phylip package) to create tree with inferred sequences at internal nodes"
        asr_prog = c['reconstruction']['asr_prog']
        if asr_prog in {'dnapars', 'dnaml'}:
            config = env.Command(
                path.join(outdir, asr_prog + ".cfg"),
                c['phy'],
                'python bin/mkconfig.py $SOURCE ' + asr_prog + ' > $TARGET')
            phylip_out = env.SRun(
                path.join(outdir, "outfile"),
                config,
                'cd ' + outdir + ' && rm -f outtree && ' + asr_prog + ' < $SOURCE.file > ' + asr_prog + '.log',
                ignore_errors=True)
            # Manually depend on phy so that we rerun dnapars/dnaml if the input sequences change (without this, dnapars/dnaml will
            # only get rerun if one of the targets are removed or if the iput asr_config file is changed). IMPORTANT!
            env.Depends(phylip_out, c['phy'])
            # Now process the phylip output into something that isn't shit
            basename = 'asr_input'
            tgt = env.Command(
                    [path.join(outdir, basename + '.' + ext) for ext in ['nwk', 'svg', 'fa']],
                    [c['seqname_mapping'], phylip_out, c['seqmeta']],
                    # Note that adding `-` at the beginning of this command string can keep things running if
                    # there's an error trying to build a tree from 2 sequences; should be filtered prior now...
                    "xvfb-run -a bin/process_asr.py"
                        + (" --seed " + c['seed'] if 'seed' in c else '')
                        + " --outdir " + outdir
                        + " --basename " + basename
                        + " --seqname-mapping $SOURCES")
            asr_tree, asr_tree_svg, asr_seqs = tgt
            # manually depnd on this because the script isn't in first position
            env.Depends(tgt, 'bin/process_asr.py')
            env.Depends(tgt, 'bin/plot_tree.py')
            return [asr_tree, asr_tree_svg, asr_seqs]
        elif asr_prog == 'raxml':
            asr_supports_tree = env.SRun(
                path.join(outdir, 'asr.sup.nwk'),
                c['pruned_seqs'],
                # Question should use -T for threads? how many?
                # Don't know if the reroot will really do what we want here
                'raxml.py --rapid-bootstrap 30 -x 3243 -o naive0 $SOURCE $TARGET')
            asr_tree_svg = env.Command(
                path.join(outdir, 'asr.svg'),
                [asr_supports_tree, c['seqmeta']],
                'xvfb-run -a bin/plot_tree.py $SOURCES $TARGET --supports'
                    + (' --seed ' + c['seed'] if 'seed' in c else ''))
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



# Temporarily turn off to debug above
# Now actually add all of these build targets
add_cluster_analysis(w)



# Popping out
# -----------

# Here we pop off the "seed" nest level and everything down from it.
# We're doing this so we can add_cluster_analysis for the unseeded runs.

w.pop('seed')





# Unseeded cluster analysis metadata:
# ===================================

def add_unseeded_analysis(w):

    # Extract metadata per sample.
    def sample_metadata(c, filename): # control dict as well?
        study = c['dataset']['study']
        matcher = is_unmerged if separate_timepoints else is_merged
        sample_name = matcher(c, filename).group(0)
        meta = heads.read_metadata(study)
        #print('meta keys', sample_name)
        #print('meta keys', meta.keys())
        d = copy.deepcopy(meta[sample_name])
        d.update(
            {'id': filename,
             'timepoints': [{'cft.timepoint:id': d['timepoint']}]})
        return d

    #dataset,shorthand,species,timepoint,subject,locus
    @w.add_nest(metadata=sample_metadata)
    @wrap_test_run()
    def sample(c):
        def keep(filename):
            return is_unmerged(c, filename) if separate_timepoints else is_merged(c, filename)
        return filter(keep,
                      os.listdir(path.join(datapath(c), 'partitions')))

    # Setting up the partition nest level

    def path_base_root(full_path):
        return path.splitext(path.basename(full_path))[0]

    # Each c["partition"] value actually points to the annotations for that partition... a little weird but...
    @w.add_nest(
            label_func=lambda d: d['base_root'],
            #metadata=lambda c, d: {'clusters': map("".join, d['clusters'])})
            metadata=lambda c, d: {'clusters': 'elided'})
    @wrap_test_run()
    def partition(c):
        """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
        for actual analysis."""
        partition_filename = partitions(c)
        cp = clusterpath.ClusterPath()
        cp.readfile(partition_filename)
        clusters = cp.partitions[cp.i_best]
        return [{'clusters': clusters,
                 # Should cluster step be partition step?
                 'cluster_step': 0,
                 'n_clusters': len(clusters),
                 'logprob': cp.logprobs[cp.i_best],
                 'base_root': path_base_root(partition_filename)}]


    def has_seeds(cluster, c):
        """Manual selection of seeds to check for from Laura; If this becomes generally useful can put in a
        data hook."""
        return any((('BF520.1-ig' + x) in cluster) for x in ['h', 'k']) and len(cluster) > 2

    # Add cluster nesting level

    @w.add_nest(label_func=lambda d: d['id'])
    def cluster(c):
        return [{'id': 'cluster' + str(i),
                 'unique_ids': clust}
                for i, clust
                # Sort by len (dec) and apply index i
                in enumerate(sorted(c['partition']['clusters'], key=len, reverse=True))
                # Select top 5 or any matching seeds of interest
                if i < 5 or has_seeds(clust, c)]


    # Finally call out to the separate cluster analyses as defined above
    add_cluster_analysis(w)


add_unseeded_analysis(w)

# Next we recreate our whole nesting thing.


# Go back to the base (build) nest level, forcing a metadata write, and print help for cftweb

w.pop('dataset')


