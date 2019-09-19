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
import sys
import csv
import time
import subprocess
import datetime
import getpass
import itertools
import yaml
import json
import re
import functools as fun
import traceback
import string

from bin import process_partis, translate_seqs

from os import path
from warnings import warn


# Nestly things
# this temporarily switches between a local checkout and whatever is installed
# Uncomment this line for a local checkout
#sys.path.append(path.join(os.getcwd(), 'nestly'))
import nestly
from nestly import scons as nestly_scons

# Tripl data modelling
# Uncomment this line for a local checkout
sys.path = [path.join(os.getcwd(), 'deps', 'tripl')] + sys.path
from tripl import nestly as nestly_tripl


# Partis and datascripts things

# If the PARTIS env var isn't already set, default to $PWD/partis (where we have a git
# submodule checkout; this is needed for bin/process_partis.py)
default_partis_path = path.join(os.getcwd(), 'partis')
partis_path = os.environ.get('PARTIS', default_partis_path)
sys.path.append(path.join(partis_path, 'python'))
import utils as partisutils
import glutils
#print("partis utils loading from:", partisutils.__file__)

# Scons requirements
from SCons.Script import Environment


# Build modules (in site_scons):
import sconsutils
import backtrans_align
import options
import software_versions


# Need this in order to read csv files with sequences in the fields
csv.field_size_limit(sys.maxsize)

# No-op; Prevents analysis warnings
sconsutils #lint


# Set up SCons environment

environ = os.environ.copy()

# install partis path as env var if not already set
environ['PARTIS'] = partis_path

env = Environment(ENV=environ)

# Add stuff to PATH
env.PrependENVPath('PATH', 'bin')
env.PrependENVPath('PATH', 'post_partis/scripts')
env.PrependENVPath('PATH', 'tree')

# Setting up command line arguments/options. See `site_scons/options.py` to see the option parsing setup.
options = options.get_options(env)




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


print("\nscons build command:", " ".join(sys.argv))

nest = nestly.Nest()
w = nestly_scons.SConsWrap(nest, options['outdir_base'], alias_environment=env)
w = nestly_tripl.NestWrap(w,
        name='build',
        # Need to base hashing off of this for optimal incrementalization
        metadata={'id': 'cft-build-' + build_time.replace(' ', '-'),
                  'time': build_time,
                  'command': " ".join(sys.argv),
                  'workdir': os.getcwd(),
                  'user': getpass.getuser(),
                  'commit': git('rev-parse', 'HEAD'),
                  # This will be really cool :-)
                  'diff': git('diff'),
                  'status': git('status', '--porcelain')},
        always_build_metadata=options['always_build_metadata'],
        base_namespace='cft',
        id_attrs=['cft.dataset:id', 'cft.build:id'])


# Recording software versions
# ---------------------------

software_versions.add_software_versions(w)



# Dataset nest level
# =================


# A dataset is a collection of data pointed to by one of the infiles.

def dataset_metadata(infile):
    with open(infile) as fp:
        if re.match('.*\.json$', infile):
            d = json.load(fp)
        else:
            d = yaml.load(fp)
    label = (options['dataset_tag'] + '-' if options['dataset_tag'] else '') + d['id']
    outdir = path.join(options['outdir_base'], label)
    return sconsutils.merge_dicts(d, {'id': label + '-' + time.strftime('%Y.%m.%d'), 'label': label, 'outdir': outdir})

@w.add_nest(full_dump= True,
        label_func= lambda d: d['label'],
        metadata= lambda c, d: {'samples': None},
        id_attrs= ['cft.subject:id', 'cft.sample:id'])
def dataset(c):
    return map(dataset_metadata, options['infiles'])


# Helpers for accessing info about the dataset

def dataset_outdir(c):
    "Returns _dataset_outdir of `c['dataset']`, for easier access via `c` below."
    return c['dataset']['outdir']

def dataset_id(c):
    return c['dataset']['id']


# Helper for running test runs on a subset of the data, togglable via the `--test` cli flag

def wrap_test_run(take_n=2):
    def deco(nestables_fn):
        def f(c):
            nestables = nestables_fn(c)
            nestables = nestables[:take_n]
            return nestables
        f.__name__ = nestables_fn.__name__
        return f if options['test_run'] else nestables_fn
    return deco


# Subject nest level
# ------------------

# We don't really do much here other than nest through the subjects in the input file for the sake of
# establishing the metadata heirarchy


def keep_sample(sample):
    return (sample.get('partition-file') and options['only_seeds'] is None) or \
           [seed for seed_id, seed in sample.get('seeds', {}).items() if seed.get('partition-file') and (options['only_seeds'] is None or seed_id in options['only_seeds'])]

def samples(c):
    return {sample_id: sample
            for sample_id, sample in c['dataset']['samples'].items()
            if keep_sample(sample)}


@w.add_nest(label_func=str)
@wrap_test_run(take_n=2)
def subject(c):
    return list(set(sconsutils.get_in(sample, ['meta', 'subject'])
                    for sample_id, sample in samples(c).items()))



# Initialize sample nest
# -----------------------------

# Samples can either point to a specific timepoint through the yaml meta, or can have a "merged" attribute
# value there, if it is a sample composed of many timepoints. These metadata will be processed accordingly.

# Samples can have partitions, and they can also have other-partitions, and seeded partitions.
# These are handled in separate nest loops below, with a pop in between.

# There may eventually be some required arguments here as this is where we get our locus and isotype and such
@w.add_nest(metadata=lambda c, d: sconsutils.merge_dicts(d.get('meta', {}),
                                                    {'id': d['id'], 'seeds': None, 'meta': None}))
@wrap_test_run(take_n=2)
def sample(c):
    # Make sure to add timepoints here as necessary
    return [sconsutils.merge_dicts(sample, {'id': sample_id})
            for sample_id, sample in samples(c).items()
            if sconsutils.get_in(sample, ['meta', 'subject']) == c['subject']]

def locus(c):
    sample = c['sample']
    locus = sample.get('locus') or sample.get('meta').get('locus')
    return locus



# Initialize seed nest
# --------------------

# We start with the seed nest level.
# Herein we'll loop over partitions and fetch the seeded clusters from the partitions of interest, as defined below.
# Eventually, we'll pop off this seed nest level so we can renest these partitions and clusters directly from the sample nest level.

# Initialize our first sub dataset nest level
@w.add_nest(metadata=lambda c, d: sconsutils.merge_dicts(d.get('meta', {}), {'id': d['id']}))
# would like to have a lower number here but sometimes we get no good clusters for the first two seeds?
# (on laura-mb for example).
@wrap_test_run(take_n=3)
def seed(c):
    return [sconsutils.merge_dicts(seed, {'id': seed_id})
            for seed_id, seed in samples(c)[c['sample']['id']].get('seeds', {}).items()
            if seed.get('partition-file') and (options['only_seeds'] is None or seed_id in options['only_seeds'])]

# Some accessor helpers

def timepoint(c):
    return sconsutils.get_in(c, ['sample', 'meta', 'timepoint'])

def is_merged(c):
    return timepoint(c) == 'merged'

def is_unmerged(c):
    return not is_merged(c)




# Seeded partitions nest level
# ---------------------

# For seeded clusters, we pick the "best" logprob partition, and if it doesn't have at least 50 seqs, we keep
# going through partitions until we find a seeded cluster that does.
# In general though, we'll end up with one partition per seed; the "best" according to the logprob.

def seed_cluster(cp, i_part_step, seed_id):
    for cluster in cp.partitions[i_part_step]:
        if seed_id in cluster:
            return cluster
    warn("unable to find seed cluster in partition")

def seed_cluster_size(cp, i_part_step, seed_id):
    return len(seed_cluster(cp, i_part_step, seed_id))

def get_alt_naive_probabilities(annotation):
    alternatives = annotation.get('alternative-annotations', {})
    naive_probabilities = alternatives.get('naive-seqs') if alternatives else None
    return naive_probabilities if naive_probabilities and len(naive_probabilities) > 0 else None

def partition_steps(cp):
    return range(len(cp.partitions)) if options['process_all_partis_partition_steps'] else [cp.i_best]

def partition_metadata(part, annotation_list, cp, i_step, seed=None, other_id=None):
    clusters = cp.partitions[i_step]
    seed_cluster_annotation = None
    if seed:
        seed_cluster_annotation = process_partis.choose_cluster(part['partition-file'], annotation_list, cp, i_step)
    
    meta = {'id': ('seed-' if seed else 'unseeded-') + (other_id + '-' if other_id else '') + 'part-' + str(i_step),
            'clusters': clusters,
            # Should cluster step be partition step?
            'step': i_step,
            'n_clusters': len(clusters),
            'largest_cluster_size': max(map(len, clusters)),
            'logprob': cp.logprobs[i_step],
            'partition-file': part['partition-file'],
            'seed_cluster_annotation': seed_cluster_annotation
            }
    if seed:
        meta['seed_cluster_size'] = seed_cluster_size(cp, i_step, seed)
    return sconsutils.merge_dicts(meta, part.get('meta') or {})

# Whenever we iterate over partitions, we always want to assume there could be an `other-partitions` mapping,
# and iterate over all these things, while tracking their other_id keys in the `other-partitions` dict.
def with_other_partitions(node):
    parts = []
    if node.get('partition-file'):
        parts.append(node)
    if node.get('other-partitions'):
        parts += [sconsutils.merge_dicts(part, {'other_id': other_id})
                  for other_id, part in node['other-partitions'].items()
                  if part.get('partition-file')]
    return parts

def meets_min_cluster_size(size, seed_cluster=False):
    """By default just checks for >= 3 sequences for seed clusters (otherwise, we can't build a tree downstream) and >= 6 for unseeded (somewhat arbitrary, though we often dont see smaller especially without processing all partition steps). This simple function exists just to track different min cluster sizes in one place"""
    return size >= (3 if seed_cluster else 6)

def valid_cluster(annotation_list, part, clust):
    """Reads the corresponding cluster annotation and return True iff after applying our health metric filters
    we still have greater than 2 sequences (otherwise, we can't build a tree downstream)."""
    for line in annotation_list:
        if line.get('unique_ids') == clust:
            func_list = [partisutils.is_functional(line, iseq) for iseq in range(len(line['unique_ids']))]
            n_good_seqs = func_list.count(True)
            return meets_min_cluster_size(n_good_seqs, seed_cluster=True)
    raise Exception('couldn\'t find requested uids %s in %s' % (clust, part['partition-file']))

def valid_seed_partition(annotation_list, cp, part, i_step, seed_id, seed_cluster_size, max_size_to_check=10):
    """If seed cluster size is less than max_size_to_check, read the corresponding cluster annotation and return True iff after applying our health metric filters
    we still have greater than 2 sequences (otherwise, we can't build a tree downstream)."""
    if seed_cluster_size > max_size_to_check:
        return True
    elif meets_min_cluster_size(seed_cluster_size, seed_cluster=True):  
        clust = seed_cluster(cp, i_step, seed_id)
        return valid_cluster(annotation_list, part, clust)
    return False

# The actual nest construction for this

# Try to read partition file; If fails, it is possibly because it's empty. Catch that case and warn
def read_partition_file(part, c):
    try:
        glfo, annotation_list, cpath = process_partis.read_partis_output(part['partition-file'], c['sample']['glfo-dir'], locus(c))
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        print(''.join((' ' * 8) + line for line in lines))
        warn("Unable to parse partition file (see error above, ommitting from results): {}".format(part))
        return []
    return annotation_list, cpath

# note we elide the nested partitions > clusters lists (as well as the seed cluster annotation)
# so as not to kill tripl when it tries to load them as a value and can't hash 
@w.add_nest(metadata=lambda c, d: {'clusters': 'elided', 'seed_cluster_annotation': 'elided'})
def partition(c):
    """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
    for actual analysis."""
    keep_partitions = []
    seed_id = c['seed']['id']
    if options['only_seeds'] is not None and seed_id not in options['only_seeds']:
        return [] 
    for part in with_other_partitions(c['seed']):
        annotation_list, cp = read_partition_file(part, c)
        if cp:
            for i_step in partition_steps(cp):
                meta = partition_metadata(part, annotation_list, cp, i_step, seed=seed_id, other_id=part.get('other_id'))
                if valid_seed_partition(annotation_list, cp, part, i_step, seed_id, meta['seed_cluster_size']):
                    keep_partitions.append(meta)
    return keep_partitions


# The cluster level
# -----------------

# For seeded clusters we only process the seed containing cluster.


# This is a little silly, but gives us the right semantics for partitions > clusters
#w.add('cluster', ['cluster'], metadata=lambda _, cluster_id: {'id': cluster_id}) # set true
@w.add_nest(label_func=lambda d: d['id'], metadata=lambda c, d: {'annotation': 'elided', 'naive_probabilities': 'elided'})
def cluster(c):
    part = c['partition']
    seed_cluster_annotation = part['seed_cluster_annotation']
    unique_ids = ':'.join(seed_cluster_annotation['unique_ids'])
    naive_probabilities = get_alt_naive_probabilities(seed_cluster_annotation)
    return [{'id': 'seed-cluster',
             'seed_name': c['seed']['id'],
             'size': part['seed_cluster_size'],
             'annotation': seed_cluster_annotation,
             'unique_ids': unique_ids,
             'naive_probabilities': naive_probabilities}]


def add_cluster_analysis(w):

    @w.add_target(name='path')
    def path_fn(outdir, c):
        return outdir

    @w.add_metadata()
    def _process_partis(outdir, c):
        # Should get this to explicitly depend on cluster0.fa
        sources = [c['partition']['partition-file']]
        perseq_metafile = c['sample'].get('per-sequence-meta-file')
        if perseq_metafile:
            sources.append(perseq_metafile)
        cluster_seqs_fname = 'cluster_seqs.fa'
        if options['match_indels_in_uid']:
            cluster_seqs_fname = '{}_indel_filtered_cluster_seqs.fa'.format(options['match_indels_in_uid'])
        return env.Command(
                [path.join(outdir, x) for x in ['partis_metadata.json', cluster_seqs_fname, 'partis_seqmeta.csv']],
                sources,
                'process_partis.py' +
                    ' --remove-stops --remove-frameshifts --remove-mutated-invariants' +
                   (' --indel-reversed-seqs ' if not options['preserve_indels'] else '') +
                    ' --partition-file ${SOURCES[0]}' +
                    ' --partition {}'.format(c['partition']['step']) +
                   (' --upstream-seqmeta ${SOURCES[1]}' if perseq_metafile else '') +
                   (' --glfo-dir ' + c['sample']['glfo-dir'] if partisutils.getsuffix(c['partition']['partition-file']) == '.csv' else '') +
                    ' --locus ' + locus(c) +
                    ' --max-sequences 10000' +
                    ' --paths-relative-to ' + dataset_outdir(c) +
                    ' --namespace cft.cluster' +
                    ' --inferred-naive-name ' + options['inferred_naive_name'] +
                  ((' --match-indels-in-uid ' + options['match_indels_in_uid']) if options['match_indels_in_uid'] is not None else '')  +
                   (' --ignore-seed-indels' if options['ignore_seed_indels'] else '') +
                  ((" --always-include " + ','.join(c['sample']['seeds'])) if c['sample'].get('seeds') else '') +
                   (' --cluster {}'.format(c['cluster']['sorted_index']) if not c.get('seed') else '') +
                    ' --cluster-meta-out ${TARGETS[0]}' +
                    ' --seqs-out ${TARGETS[1]}' +
                    ' --seqmeta-out ${TARGETS[2]}')

    @w.add_target(ingest=True)
    def partis_metadata(outdir, c):
        return c['_process_partis'][0]

    @w.add_target()
    def inseqs(outdir, c):
        return c['_process_partis'][1]

    @w.add_target()
    def partis_seqmeta(outdir, c):
        return c['_process_partis'][2]

    # Partis alternative naives
    # -------------------------

    @w.add_target()
    def alternative_naive_probabilities(outdir, c):
        '''
        Write partis alternative naives to a fasta in order of probability
        '''
        if c['cluster']['naive_probabilities'] is not None:
            cluster_name = c['cluster'].get('seed_name', c['cluster']['id'])

            naives_sorted_by_prob =  list(sorted(c['cluster']['naive_probabilities'], key=lambda x: x[1], reverse=True))
            def write_naive_fastas(target, source, env):
                '''
                This is the action for this target. Because it is a function, not a file being executed,
                the rules SCons follows for determining whether to rebuild this target are less well defined.
                '''
                targets = [str(fname) for fname in target]
                with open(targets[0], 'w') as ranked_fasta, open(targets[1], 'w') as aa_ranked_fasta:
                    for rank, (naive_seq, probability) in enumerate(naives_sorted_by_prob):
                        aa_seq = translate_seqs.translate(naive_seq)
                        ranked_fasta.write('>%s\n%s\n' % ('naive_{}_probability_{}'.format(rank, probability), naive_seq))
                        aa_ranked_fasta.write('>%s\n%s\n' % ('naive_{}_probability_{}'.format(rank, probability), aa_seq))
            
            naive_probs_fname = path.join(outdir, 'ranked_naive_probabilities_%s.fasta' % cluster_name)
            aa_naive_probs_fname = path.join(outdir, 'ranked_aa_naive_probabilities_%s.fasta' % cluster_name)

            return env.Command([naive_probs_fname, aa_naive_probs_fname], c['partition']['partition-file'], write_naive_fastas)
    
    @w.add_target()
    def alternative_naive_logo_plots(outdir, c):
        '''
        Create logo plot according to probabilities
        '''
        if c['cluster']['naive_probabilities'] is not None:
            
            annotation = c['cluster']['annotation']
            cluster_name = c['cluster'].get('seed_name', c['cluster']['id'])  

            aa_input_fasta_path = str(c['alternative_naive_probabilities'][1])
            aa_cdr3_start, aa_cdr3_end = int(annotation['codon_positions']['v']/3), int((annotation['codon_positions']['j'] + 3)/3)
            aa_naive_len = int(len(annotation['naive_seq'])/3)
            logo_out, cdr3_logo_out = 'naive_logo_{}.png'.format(cluster_name), 'naive_logo_cdr3_{}.png'.format(cluster_name) 
            logo_plots = env.Command( [path.join(outdir, logo_out), path.join(outdir, cdr3_logo_out)],
                                aa_input_fasta_path,
                                'python bin/create_partis_naive_logo.py $SOURCE' +
                                ' --aa-cdr3-start=%d' % aa_cdr3_start +
                                ' --aa-cdr3-end=%d' % aa_cdr3_end +
                                ' --aa-naive-len=%d' % aa_naive_len +
                                ' --logo-fname=%s' % logo_out +
                                ' --cdr3-logo-fname=%s' % cdr3_logo_out +
                                ' --outdir=%s' % outdir )
 
            env.Depends(logo_plots, 'bin/create_partis_naive_logo.py')
            return logo_plots

    # Sequence Alignment
    # ------------------

    backtrans_align.add(env, w, options)


    # On with trees and other things...
    # ---------------------------------

    # use fasttree to make newick tree from sequences
    @w.add_target()
    def fasttree(outdir, c):
        return env.SRun(
            path.join(outdir, "fasttree.nwk"),
            c['aligned_inseqs'],
            "FastTree -nt -quiet $SOURCE > $TARGET 2> $TARGET-.log")

    # @w.add_target(ingest=True)
    # def selection_metrics(baseoutdir, c):
    #     outdir = path.join(baseoutdir, 'selection-metrics')
    #     return env.Command(
    #         [path.join(outdir, "tree-stats.json")],
    #         [path.join(baseoutdir, "fasttree.nwk"), c['aligned_inseqs']],  # sources # don't use fasttree any more
    #         # TODO lonr a.t.m. is still making its own trees, which should probably change TODO
    #         "%s/bin/calculate_tree_metrics.py --seqfile ${SOURCES[1]} XXX --treefile ${SOURCES[0]} XXX --outfile ${TARGETS[0]} --naive-seq-name %s --debug" % (partis_path, options['inferred_naive_name'])
    #     )

    @w.add_nest(metadata=lambda c, d: d)
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
        tgt = path.join(outdir, "pruned_ids.txt")
        # This whole thing is a safety mechanism to prevent prune files with 0 sequences from hanging around,
        # which can happen from failed builds
        try:
            remove = False
            with open(tgt, 'r') as fh:
                remove = len(fh.readlines()) == 0
            if remove:
                os.remove(tgt)
        except:
            pass
        recon = c['reconstruction']
        builder = fun.partial(env.SRun, srun_args='`minadcl_srun_args.py $SOURCE`') \
                if recon['prune_strategy'] == 'min_adcl' \
                else env.Command
        return builder(tgt,
            c['fasttree'],
            "prune.py -n " + str(recon['prune_count'])
                + ((" --always-include " + ','.join(c['sample']['seeds'])) if c['sample'].get('seeds') else '')
                + " --strategy " + recon['prune_strategy']
                + " --naive %s" % options['inferred_naive_name']
                + (" --seed " + c['seed']['id'] if 'seed' in c else '')
                + " $SOURCE $TARGET")


    if options['fasttree_png']:
        # create png showing included seqs (kept in pruning) as red
        @w.add_target()
        def pruned_cluster_fasttree_png(outdir, c):
            cluster = c['cluster']
            max_cluster_size = 4500
            if cluster.get('size') < max_cluster_size:
                pruned_cluster_fasttree_png = env.Command(
                    path.join(outdir, "pruned_cluster_fasttree.png"),
                    [c["fasttree"], c["pruned_ids"]],
                    # The `-` at the start here tells scons to ignore if it doesn't build; this may occasionally be
                    # the case for large clusters. Also, redirect stdin/out to dev/null because the errors messages
                    # here can be pretty noisy.
                    "- xvfb-run -a bin/annotate_fasttree_tree.py $SOURCES " + " --naive %s" % options['inferred_naive_name'] +
                        (" --seed " + c['seed']['id'] if 'seed' in c else '')  +
                        " --output-path $TARGET &>> /dev/null" )
                env.Depends(pruned_cluster_fasttree_png, "bin/annotate_fasttree_tree.py")
                return pruned_cluster_fasttree_png

    @w.add_target()
    def cluster_mapping(outdir, c):
        if c['reconstruction']['prune_strategy'] == 'min_adcl':
            return env.SRun(
                path.join(outdir, 'cluster_mapping.csv'),
                [c['fasttree'], c['pruned_ids']],
                'minadcl_clusters.py $SOURCES $TARGET',
                srun_args='`minadcl_clusters_srun_args.py $SOURCE`')


    # prune out sequences to reduce taxa, making sure to cut out columns in the alignment that are now entirely
    # gaps from insertions in sequences that have been pruned out.
    @w.add_target()
    def pruned_seqs(outdir, c):
        return env.Command(
            path.join(outdir, "pruned.fa"),
            [c['pruned_ids'], c['aligned_inseqs']],
            "seqmagick convert --include-from-file $SOURCES - | " +
            "seqmagick convert --squeeze - $TARGET")

    if options['write_linearham_yaml_input']:
        @w.add_target()
        def pruned_partis_outfile(outdir, c):
            if 'seed' in c:
                #these are not the unique ids we want to use when subsetting the cluster, they are just a way to identify the cluster we want to subset
                clust_ids_string = c['cluster']['unique_ids']
                yaml_format = partisutils.getsuffix(c['partition']['partition-file']) == '.yaml'
                return env.Command(
                    path.join(outdir, 'pruned_partis_output.yaml'),
                    [c['partition']['partition-file'], c['pruned_ids']],
                    'python bin/write_subset_partis_outfile.py $SOURCES $TARGET' +
                    ' --partition-step={}'.format(c['partition']['step']) +
                    ' --sw-cache={}'.format(c['sample']['sw-cache']) + 
                    (' --glfo-dir={}'.format(c['sample']['glfo-dir']) if not yaml_format else '') +
                    (' --locus={}'.format(locus(c)) if not yaml_format else '') )
        
        @w.add_target()
        def linearham_base_command(outdir, c):
            """ This allows us to not have to manually build the minimum necessary commmad for running this seed
            cluster through Linearham. Eventually we may want to run this command here in the CFT SCons pipeline,
            but for now we make life a little easier by just being able to copy this command to run Linearham."""
            if 'seed' in c:
                return env.Command(
                    path.join(outdir, 'linearham_base_command.txt'),
                    c['pruned_partis_outfile'], 
                    'echo "scons --run-linearham --template-path=templates/revbayes_template.rev ' +
                    ' --parameter-dir={}'.format(c['sample']['parameter-dir']) + 
                    ' --partis-yaml-file={}'.format(path.join(os.getcwd(), str(c['pruned_partis_outfile'][0]))) + #get cwd for absolute path
                    ' --seed-seq={}" > $TARGET'.format(c['cluster']['seed_name']))

    @w.add_target()
    def tip_seqmeta(outdir, c):
        """The merge of process_partis output with pre sequence metadata spit out by datascripts containing
        timepoint mappings. Base input multiplicity is coded into the original input sequence names from vlad as N-M,
        where N is the ranking of vlads untrimmed deduplication, and M is the multiplicity of said deduplication."""
        # This option controls which sequences get joined on in the merge for the partis_seqmeta file, which has
        # orig/new names, joined on sequence from the other file
        sources = {'--partis-seqmeta': c['partis_seqmeta'],
                   '--cluster-mapping': c['cluster_mapping'] if c['reconstruction']['prune_strategy'] == 'min_adcl' else None,
                   '--pruned-ids': c['pruned_ids'] if c['reconstruction']['prune_strategy'] == 'seed_lineage' else None, 
                  }
        sources = {k: v for k, v in sources.items() if v}
        base_call = 'aggregate_minadcl_cluster_multiplicities.py '
        for i, (k, v) in enumerate(sources.items()):
            base_call += k + ' ${SOURCES[' + str(i) + ']} '
        return env.Command(path.join(outdir, 'tip_seqmeta.csv'), sources.values(), base_call + '$TARGET')

    @w.add_target()
    def _phy(outdir, c):
        "Save seqs in phylip format for dnaml, renaming sequences as necessary to avoid seqname limits"
        return env.Command(
            [path.join(outdir, x) for x in ('pruned.phy', 'seqname_mapping.csv')],
            c['pruned_seqs'],
            'make_phylip.py $SOURCE $TARGETS --inferred-naive-name ' + options['inferred_naive_name'])

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
    def _asr(outdir, c):
        "run dnapars/dnaml (from phylip package) to create tree with inferred sequences at internal nodes"
        asr_prog = c['reconstruction']['asr_prog']
        if asr_prog in {'dnapars', 'dnaml'}:
            config = env.Command(
                path.join(outdir, asr_prog + ".cfg"),
                c['phy'],
                'python bin/mkconfig.py $SOURCE ' + asr_prog + '> $TARGET')
            phylip_out = env.SRun(
                path.join(outdir, "outfile"),
                config,
                'cd ' + outdir + ' && rm -f outtree && ' + asr_prog + ' < $SOURCE.file > ' + asr_prog + '.log')
                #ignore_errors=True) # before; dropping
            # Manually depend on phy so that we rerun dnapars/dnaml if the input sequences change (without this, dnapars/dnaml will
            # only get rerun if one of the targets are removed or if the iput asr_config file is changed). IMPORTANT!
            env.Depends(phylip_out, c['phy'])
            # Now process the phylip output into something that isn't shit
            basename = 'asr'
            tgt = env.Command(
                    [path.join(outdir, basename + '.' + ext) for ext in ['nwk', 'svg', 'fa', 'ancestors_naive_and_seed.fa']],
                    [c['seqname_mapping'], phylip_out, c['tip_seqmeta']],
                    # Note that adding `-` at the beginning of this command string can keep things running if
                    # there's an error trying to build a tree from 2 sequences; should be filtered prior now...
                    "xvfb-run -a bin/process_asr.py"
                        + (" --seed " + c['seed']['id'] if 'seed' in c else '')
                        + " --outdir " + outdir
                        + " --basename " + basename
                        + " --inferred-naive-name " + options['inferred_naive_name']
                        + " --seqname-mapping $SOURCES")
            asr_tree, asr_tree_svg, asr_seqs, ancestors_naive_and_seed = tgt
            # manually depnd on this because the script isn't in first position
            env.Depends(tgt, 'bin/process_asr.py')
            env.Depends(tgt, 'bin/plot_tree.py')
            return [asr_tree, asr_tree_svg, asr_seqs, ancestors_naive_and_seed]
        elif asr_prog == 'raxml':
            asr_supports_tree = env.SRun(
                path.join(outdir, 'asr.sup.nwk'),
                c['pruned_seqs'],
                # Question should use -T for threads? how many?
                # Don't know if the reroot will really do what we want here
                'raxml.py --rapid-bootstrap 30 -x 3243 -o %s $SOURCE $TARGET' % options['inferred_naive_name'])
            asr_tree_svg = env.Command(
                path.join(outdir, 'asr.svg'),
                [asr_supports_tree, c['tip_seqmeta']],
                'xvfb-run -a bin/plot_tree.py $SOURCES $TARGET --supports --inferred-naive-name ' + options['inferred_naive_name']
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
        #return c['_asr'][0]

    @w.add_target(ingest=True)
    def asr_tree_svg(outdir, c):
        return c['_asr'][1]

    #@w.add_target()
    #def asr_supports_tree(outdir, c):
        #vals = c['_asr']
        #if len(vals) > 2:
            #return c['_asr'][2]

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
        return c['_asr'][0]

    @w.add_target(ingest=True)
    def asr_seqs(outdir, c):
        return c['_asr'][2]

    @w.add_target(ingest=True)
    def ancestors_naive_and_seed(outdir, c):
        return c['_asr'][3]

    @w.add_target()
    def sampled_ancestors(outdir, c):
        sampled_ancestors_output = env.Command(
                path.join(outdir, "sampled_ancestors.json"),
                [c["aligned_inseqs"], c["ancestors_naive_and_seed"]],
                "time python bin/find_sampled_ancestors.py $SOURCES $TARGET")
        env.Depends(sampled_ancestors_output, "bin/find_sampled_ancestors.py")
        return sampled_ancestors_output

    @w.add_target()
    def selection_metrics(baseoutdir, c):
        outdir = path.join(baseoutdir, 'selection-metrics')
        def fix_file(target, source, env):
            extra_ambig_chars = 'RYSWKMBDHV-'
            target = str(target[0])
            source = str(source[0])
            ambiguous_char_translations = string.maketrans(extra_ambig_chars, ''.join('N' for _ in range(len(extra_ambig_chars))))
            seqfos = partisutils.read_fastx(source)
            with open(target, 'w') as nfile:
                for sfo in seqfos:
                    nfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq'].translate(ambiguous_char_translations)))
        new_partis_infname = path.join(baseoutdir, 'asr-with-only-Ns.fa')
        new_partis_infile = env.Command(new_partis_infname, c['asr_seqs'], fix_file)

        tree_metrics = env.Command(
            [path.join(outdir, "selection-metrics.yaml")],
            [path.join(baseoutdir, 'asr.nwk')],  # sources
            "%s/bin/get-tree-metrics.py ${SOURCES[0]} ${TARGETS[0]}" % (partis_path)
        )
        env.Depends(tree_metrics, "partis/bin/get-tree-metrics.py")
        return tree_metrics

    @w.add_target(ingest=True, attr_map={'bio.seq:id': 'sequence', 'cft.timepoint:id': 'timepoint',
        'cft.seq:timepoint': 'timepoint', 'cft.seq:timepoints': 'timepoints', 'cft.seq:cluster_timepoints': 'cluster_timepoints',
        'cft.seq:multiplicity': 'multiplicity', 'cft.seq:cluster_multiplicity': 'cluster_multiplicity',
        'cft.seq:timepoint_multiplicities': 'timepoint_multiplicities',
        'cft.seq:cluster_timepoint_multiplicities': 'cluster_timepoint_multiplicities',
        'cft.seq:affinity': 'affinity',
        'cft.tree.node:lbi': 'lbi', 'cft.tree.node:lbr': 'lbr'})
    def seqmeta(outdir, c):
        return env.Command(path.join(outdir, "seqmeta.csv"),
            [c['tip_seqmeta'], c['selection_metrics']],
            "merge_selection_metrics.py $SOURCES $TARGET")


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



# Unseeded cluster analysis
# -------------------------


# Now we define a function that builds the cluster analysis defined above for the unseeded
# partitions/clusters. We do this so that the functions names we use for building things with nestly don't
# overlap with those defined for the seeded analysis above

def add_unseeded_analysis(w):

    # Setting up the partition nest level

    # Each c["partition"] value actually points to the annotations for that partition... a little weird but...
    # We elide the annotations here because we pass them through in order to avoid re-reading them from
    # the partition file and because we don't need to write them to the metadata for the partition.
    # See https://github.com/matsengrp/cft/pull/270#discussion_r267502415 for details on why we decided to do things this way.
    @w.add_nest(metadata=lambda c, d: {'clusters': 'elided', 'cp': 'elided' })
    def partition(c):
        """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
        for actual analysis."""
        keep_partitions = []
        for partition_run in with_other_partitions(c['sample']):
            annotation_list, cp = read_partition_file(partition_run, c)
            if cp:
                keep_partitions += [partition_metadata(partition_run, annotation_list, cp, i_step, other_id=partition_run.get('other_id')) for i_step in partition_steps(cp)]
        return keep_partitions

    # Add cluster nesting level
    @w.add_nest(label_func=lambda d: d['id'], metadata=lambda c, d: {'unique_ids': 'elided', 'annotation': 'elided', 'naive_probabilities': 'elided'})
    def cluster(c):
        part = c['partition']
        clusters = []
        annotation_list = None
        # Sort by len (dec) and apply index i
        for i, clust in enumerate(sorted(part['clusters'], key=len, reverse=True)):
            # Select top N or any matching seeds of interest
            if (i < options['depth']) and meets_min_cluster_size(len(clust)):
                if annotation_list is None:
                    # Here we reread the partition file instead of caching annotations of the partition along with its metadata above in partition_metadata. This saves on memory and slows the process down, but we are restricted by memory use more than time at the moment.
                    annotation_list, cp = read_partition_file(part, c)
                if valid_cluster(annotation_list, part, clust):
                    cluster_annotation = process_partis.choose_cluster(part['partition-file'], annotation_list, cp, part['step'], i)
                    # It seems like we might only need to check that one of these clusters has alternative naive info and then  we could assume it is the case for
                    # all of them (unless --queries was set for --calculate-alternative-naive-seqs). Leaving it as is for now but may speed up the SConstruct process to do this later. (EH)
                    naive_probabilities = get_alt_naive_probabilities(cluster_annotation)       
                    cluster_meta = {'id': 'clust-' + str(i),
                                    'sorted_index': i,
                                    'unique_ids': clust,
                                    'size': len(clust),
                                    'annotation': cluster_annotation,
                                    'naive_probabilities': naive_probabilities}
                    clusters.append(cluster_meta)
        return clusters

    add_cluster_analysis(w)

    # end function


# Then we immediately call this function we define

add_unseeded_analysis(w)


# Write metadata
# --------------------

# could also do this by a final w.pop('dataset'), but we want a little more control here for snapshots

w.pop('subject')

@w.add_target()
def metadata_snapshot(outdir, c):
    return env.Command(
            path.join(outdir, time.strftime('%Y-%m-%d') + ('-test' if options['test_run'] else '') + '-metadata.json'),
            path.join(outdir, 'metadata.json'),
            'cp $SOURCE $TARGET')

w.pop('dataset')


