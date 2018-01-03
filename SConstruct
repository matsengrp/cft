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
sys.path = [path.join(os.getcwd(), 'tripl')] + sys.path
from tripl import nestly as nestly_tripl


# Partis and datascripts things

# If the PARTIS env var isn't already set, default to $PWD/partis (where we have a git
# submodule checkout; this is needed for bin/process_partis.py)
default_partis_path = path.join(os.getcwd(), 'partis')
partis_path = os.environ.get('PARTIS', default_partis_path)
sys.path.append(path.join(partis_path, 'python'))
import clusterpath

# Scons requirements
from SCons.Script import Environment


# Build modules (in site_scons):
import sconsutils
import utils
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
    return utils.merge_dicts(d, {'id': d['id'] + '-' + time.strftime('%Y.%m.%d'), 'label': label, 'outdir': outdir})

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

def datapath(c):
    "Returns the input realpath datapath value of the dataset"
    return c['dataset']['datapath']


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

@w.add_nest()
@wrap_test_run(take_n=2)
def subject(c):
    return list(set(utils.get_in(sample, ['meta', 'subject'])
                    for sample_id, sample in c['dataset']['samples'].items()))



# Initialize sample nest
# -----------------------------

# Samples an either point to a specific timepoint through the yaml meta, or can have a "merged" attribute
# value there, if it is a sample composed of many timepoints. These metadata will be processed accordingly.

# Samples can have partitions, and they can also have other-partitions, and seeded partitions.
# These are handled in separate nest loops below, with a pop in between.

# There may eventually be some required arguments here as this is where we get our locus and isotype and such
@w.add_nest(metadata=lambda c, d: utils.merge_dicts(d.get('meta', {}),
                                                    {'id': d['id'], 'seeds': None, 'meta': None}))
@wrap_test_run(take_n=2)
def sample(c):
    # Make sure to add timepoints here as necessary
    return [utils.merge_dicts(sample, {'id': sample_id})
            for sample_id, sample in c['dataset']['samples'].items()
            if utils.get_in(sample, ['meta', 'subject']) == c['subject']]



# Initialize seed nest
# --------------------

# We start with the seed nest level.
# Herein we'll loop over partitions and fetch the seeded clusters from the partitions of interest, as defined below.
# Eventually, we'll pop off this seed nest level so we can renest these partitions and clusters directly from the sample nest level.

# Initialize our first sub dataset nest level
@w.add_nest(metadata=lambda c, d: utils.merge_dicts(d.get('meta', {}), {'id': d['id']}))
# would like to have a lower number here but sometimes we get no good clusters for the first two seeds?
# (on laura-mb for example).
@wrap_test_run(take_n=4)
def seed(c):
    return [utils.merge_dicts(seed, {'id': seed_id})
            for seed_id, seed in c['dataset']['samples'][c['sample']['id']].get('seeds', {}).items()]

# Some accessor helpers

def timepoint(c):
    return utils.get_in(c, ['sample', 'meta', 'timepoint'])

def is_merged(c):
    return timepoint(c) == 'merged'

def is_unmerged(c):
    return not is_merged(c)




# Seeded partitions nest level
# ---------------------

# For seeded clusters, we pick the "best" logprob partition, and if it doesn't have at least 50 seqs, we keep
# going through partitions until we find a seeded cluster that does.
# In general though, we'll end up with one partition per seed; the "best" according to the logprob.

def seed_cluster(cp, best_plus_i, seed_id):
    i = cp.i_best + best_plus_i
    for cluster in cp.partitions[i]:
        if seed_id in cluster:
            return cluster
    warn("unable to find seed cluster in partition")

def seed_cluster_size(cp, best_plus_i, seed_id):
    return len(seed_cluster(cp, best_plus_i, seed_id))

def partition_metadata(part, cp, best_plus_i, seed=None, other_id=None):
    i = cp.i_best + best_plus_i
    clusters = cp.partitions[i]
    meta = {'id': ('seed-' if seed else 'unseeded-') + (other_id + '-' if other_id else '') + 'part-' + str(best_plus_i),
            'clusters': clusters,
            # Should cluster step be partition step?
            'step': best_plus_i,
            'n_clusters': len(clusters),
            'largest_cluster_size': max(map(len, clusters)),
            'logprob': cp.logprobs[i],
            'partition-file': part['partition-file'],
            'cluster-annotation-file': part['cluster-annotation-file']}
    if seed:
        meta['seed_cluster_size'] = seed_cluster_size(cp, best_plus_i, seed)
    return utils.merge_dicts(meta, part.get('meta') or {})

# Whenever we iterate over partitions, we always want to assume there could be an `other-partitions` mapping,
# and iterate over all these things, while tracking their other_id keys in the `other-partitions` dict.
def with_other_partitions(node):
    parts = []
    if node.get('partition-file') and node.get('cluster-annotation-file'):
        parts.append(node)
    if node.get('other-partitions'):
        parts += [utils.merge_dicts(part, {'other_id': other_id})
                  for other_id, part in node['other-partitions'].items()
                  if part.get('partition-file') and part.get('cluster-annotation-file')]
    return parts


def valid_cluster(cp, part, clust):
    """Reads the corresponding cluster annotation and return True iff after applying our health metric filters
    we still have greater than 2 sequences (otherwise, we can't build a tree downstream)."""
    clust_sig = ':'.join(clust)
    with open(part['cluster-annotation-file']) as fh:
        for cluster_row in csv.DictReader(fh):
            if cluster_row['unique_ids'] == clust_sig:
                n_good_seqs = sum(map(lambda x: x[0] and not (x[1] or x[2]),
                    zip(*(map(lambda col: [x == 'True' for x in cluster_row[col].split(':')],
                              ['in_frames', 'stops', 'mutated_invariants'])))))
                return n_good_seqs > 2

def valid_seed_partition(cp, part, best_plus_i, seed_id):
    """Reads the corresponding cluster annotation and return True iff after applying our health metric filters
    we still have greater than 2 sequences (otherwise, we can't build a tree downstream)."""
    clust = seed_cluster(cp, best_plus_i, seed_id)
    return valid_cluster(cp, part, clust)

# The actual nest construction for this

# note we elide the nested partitions > clusters lists so as not to kill tripl when it tries to load them as a
# value and can't hash
@w.add_nest(metadata=lambda c, d: {'clusters': 'elided'})
@wrap_test_run()
def partition(c):
    """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
    for actual analysis."""
    keep_partitions = []
    for part in with_other_partitions(c['seed']):
        cp = clusterpath.ClusterPath()
        cp.readfile(part['partition-file'])
        # important to start from i_best
        for best_plus_i in range(len(cp.partitions) - cp.i_best):
            meta = partition_metadata(part, cp, best_plus_i, seed=c['seed']['id'], other_id=part.get('other_id'))
            # We only add clusters bigger than two, since we can only make trees if we have hits
            if meta['seed_cluster_size'] > 10:
                # if we have 10 sequences, assume enough of them will be good
                keep_partitions.append(meta)
            elif meta['seed_cluster_size'] > 2 and valid_seed_partition(cp, part, best_plus_i, c['seed']['id']):
                # if less than 10 sequences, make sure we still have enough sequences after health filters
                keep_partitions.append(meta)
            # Once we get a cluster of size 50, we don't need later cluster steps
            #if meta['seed_cluster_size'] >= 50:
                #break
            # Ignore the break above... For now, we break as soon as we get a single
            # cluster until psathyrella is able to fix partis to export cluster annotations for all
            # partition steps to the same file.
            break
    return keep_partitions



# The cluster level
# -----------------

# For seeded clusters we only process the seed containing cluster.


# This is a little silly, but gives us the right semantics for partitions > clusters
#w.add('cluster', ['cluster'], metadata=lambda _, cluster_id: {'id': cluster_id}) # set true
w.add('cluster', ['seed-cluster'])


# This one we're not using anymore; should delete, but for now.
#def parameter_dir(c):
    #"The input parameter directory (as in input to partis); see --parameter-dir flag to various partis cli calls."
    #return path.join(datapath(c), c['sample'])


def add_cluster_analysis(w):
    @w.add_metadata()
    def _process_partis(outdir, c):
        # Should get this to explicitly depend on cluster0.fa
        sources = [c['partition']['partition-file'], c['partition']['cluster-annotation-file']]
        return env.Command(
                [path.join(outdir, x) for x in ['partis_metadata.json', 'cluster_seqs.fa', 'partis_seqmeta.csv']],
                sources,
                'process_partis.py' +
                    ' --remove-stops --remove-frameshifts --remove-mutated-invariants' +
                    ' --partition-file ${SOURCES[0]}' +
                    ' --cluster-annotation-file ${SOURCES[1]}' +
                    ' --parameter-dir ' + c['sample']['parameter-dir'] +
                    ' --locus ' + c['sample']['meta']['locus'] +
                    ' --cluster-meta-out ${TARGETS[0]}' +
                    ' --seqs-out ${TARGETS[1]}' +
                    ' --seqmeta-out ${TARGETS[2]}' +
                    ' --paths-relative-to ' + dataset_outdir(c) +
                    ' --namespace cft.cluster' +
                    (' --partition {}'.format(c['partition']['step']) if c.get('seed') else '') +
                    (' --cluster {}'.format(c['cluster']['sorted_index']) if not c.get('seed') else ''))

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
        recon = c['reconstruction']
        builder = env.SRun if recon['prune_strategy'] == 'min_adcl' else env.Command
        return builder(
            path.join(outdir, "pruned_ids.txt"),
            c['fasttree'],
            "prune.py -n " + str(recon['prune_count'])
                + ((" --always-include " + ','.join(c['sample']['seeds'])) if c['sample'].get('seeds') else '')
                + " --strategy " + recon['prune_strategy']
                + " --naive naive"
                + (" --seed " + c['seed']['id'] if 'seed' in c else '')
                + " $SOURCE $TARGET")


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


    @w.add_target(ingest=True, attr_map={'bio.seq:id': 'sequence', 'cft.timepoint:id': 'timepoint',
        #'cft.seq:multiplicity': 'multiplicity'})
        'cft.seq:multiplicity': 'multiplicity', 'cft.seq:cluster_multiplicity': 'cluster_multiplicity'})
    def seqmeta(outdir, c):
        """The merge of process_partis output with pre sequence metadata spit out by datascripts containing
        timepoint mappings. Base input multiplicity is coded into the original input sequence names from vlad as N-M,
        where N is the ranking of vlads untrimmed deduplication, and M is the multiplicity of said deduplication."""
        # This option controls which sequences get joined on in the merge for the partis_seqmeta file, which has
        # orig/new names, joined on sequence from the other file
        sources = {'--partis-seqmeta': c['partis_seqmeta'],
                   '--cluster-mapping': c['cluster_mapping'] if c['reconstruction']['prune_strategy'] == 'min_adcl' else None,
                   '--upstream-seqmeta': c['sample'].get('per-sequence-meta-file')
                   }
        sources = {k: v for k, v in sources.items() if v}
        base_call = 'merge_timepoints_and_multiplicity.py '
        for i, (k, v) in enumerate(sources.items()):
            base_call += k + ' ${SOURCES[' + str(i) + ']} '
        return env.Command(path.join(outdir, 'seqmeta.csv'), sources.values(), base_call + '$TARGET')

    @w.add_target()
    def _phy(outdir, c):
        "Save seqs in phylip format for dnaml, renaming sequences as necessary to avoid seqname limits"
        return env.Command(
            [path.join(outdir, x) for x in ('pruned.phy', 'seqname_mapping.csv')],
            c['pruned_seqs'],
            'make_phylip.py $SOURCE $TARGETS --dont-rename naive')


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
                        + (" --seed " + c['seed']['id'] if 'seed' in c else '')
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
                'raxml.py --rapid-bootstrap 30 -x 3243 -o naive $SOURCE $TARGET')
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



# Unseeded cluster analysis
# -------------------------


# Now we define a function that builds the cluster analysis defined above for the unseeded
# partitions/clusters. We do this so that the functions names we use for building things with nestly don't
# overlap with those defined for the seeded analysis above

def add_unseeded_analysis(w):

    # Setting up the partition nest level

    # Each c["partition"] value actually points to the annotations for that partition... a little weird but...
    @w.add_nest(metadata=lambda c, d: {'clusters': 'elided', 'cp': 'elided'})
    @wrap_test_run()
    def partition(c):
        """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
        for actual analysis."""
        def meta(part):
            cp = clusterpath.ClusterPath()
            cp.readfile(part['partition-file'])
            return utils.merge_dicts(partition_metadata(part, cp, 0, other_id=part.get('other_id')),
                                     {'cp': cp})
        return [meta(part)
                for part in with_other_partitions(c['sample'])]

    def has_seeds(cluster, c):
        """Manual selection of seeds to check for from Laura; If this becomes generally useful can put in a
        data hook."""
        # Is the len(clusters) check working here?
        return any((('BF520.1-ig' + x) in cluster) for x in ['h', 'k']) and len(cluster) > 2

    # Add cluster nesting level

    @w.add_nest(label_func=lambda d: d['id'])
    def cluster(c):
        part = c['partition']
        cp = part['cp']
        return [{'id': 'clust-' + str(i),
                 'sorted_index': i,
                 'unique_ids': clust,
                 'size': len(clust)}
                for i, clust
                # Sort by len (dec) and apply index i
                in enumerate(sorted(c['partition']['clusters'], key=len, reverse=True))
                # Select top 5 or any matching seeds of interest
                if (len(clust) > 5) and (i < 5 or has_seeds(clust, c)) and valid_cluster(cp, part, clust)]


    # Finally call out to the separate cluster analyses as defined above
    add_cluster_analysis(w)

    # end function


# Then we immediately call this function we define

add_unseeded_analysis(w)


# Popping off metadata
# --------------------

# Go back to the base (build) nest level, forcing a metadata write

w.pop('dataset')



