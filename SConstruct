#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This scons pipeline does the following:

* Runs process_partis.py for each seed/sequencing-run/parition combination, producing a metadata file of information
* For each such parition, the cluster with the seed is analyzed:
    * Tree construction using `FastTree`
    * Ancestral state reconstruction using `raxml-ng` (default) and/or `dnaml`
* Results of all these analyses are then pointed to by a merged metadata file, which is then consumed by Olmsted https://github.com/matsengrp/olmsted

Clusters with only two sequences cannot be analyzed by raxml-ng, dnaml, or FastTree, and so these are skipped atuomatically.
See README for typical environment setup and usage.
"""

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
# sys.path.append(path.join(os.getcwd(), 'nestly'))
import nestly
from nestly import scons as nestly_scons

# Tripl data modelling
# Uncomment this line for a local checkout
sys.path = [path.join(os.getcwd(), "deps", "tripl")] + sys.path
from tripl import nestly as nestly_tripl

# Partis and datascripts things

# If the PARTIS env var isn't already set, default to $PWD/partis (where we have a git
# submodule checkout; this is needed for bin/process_partis.py)
default_partis_path = path.join(os.getcwd(), "partis")
partis_path = os.environ.get("PARTIS", default_partis_path)
sys.path.append(path.join(partis_path, "python"))
import utils as partisutils
import glutils

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
sconsutils  # lint

# Set up SCons environment
environ = os.environ.copy()

# install partis path as env var if not already set
environ["PARTIS"] = partis_path

env = Environment(ENV=environ)

# Add stuff to PATH
env.PrependENVPath("PATH", "bin")
env.PrependENVPath("PATH", "post_partis/scripts")
env.PrependENVPath("PATH", "tree")

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
    return subprocess.check_output(["git"] + list(args))


print("\nscons build command:", " ".join(sys.argv))

nest = nestly.Nest()
w = nestly_scons.SConsWrap(nest, options["outdir_base"], alias_environment=env)
w = nestly_tripl.NestWrap(
    w,
    name="build",
    metadata={
        "id": "cft-build-" + build_time.replace(" ", "-"),
        "time": build_time,
        "command": " ".join(sys.argv),
        "workdir": os.getcwd(),
        "user": getpass.getuser(),
        "commit": git("rev-parse", "HEAD"),
        "diff": git("diff"),
        "status": git("status", "--porcelain"),
    },
    always_build_metadata=options["always_build_metadata"],
    base_namespace="cft",
    id_attrs=["cft.dataset:id", "cft.build:id"],
)

# Recording software versions
# ---------------------------
software_versions.add_software_versions(w)

# Dataset nest level
# =================

# A dataset is a collection of data pointed to by one of the infiles.


def dataset_metadata(infile):
    with open(infile) as fp:
        if re.match(".*\.json$", infile):
            d = json.load(fp)
        else:
            d = yaml.load(fp)
    label = (options["dataset_tag"] + "-" if options["dataset_tag"] else "") + d["id"]
    outdir = path.join(options["outdir_base"], label)
    return sconsutils.merge_dicts(
        d,
        {
            "id": label + "-" + time.strftime("%Y.%m.%d"),
            "label": label,
            "outdir": outdir,
        },
    )


# See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
@w.add_nest(
    full_dump=True,
    label_func=lambda d: d["label"],
    metadata=lambda c, d: {"samples": None},
    id_attrs=["cft.subject:id", "cft.sample:id"],
)
def dataset(c):
    return map(dataset_metadata, options["infiles"])


# Helpers for accessing info about the dataset


def dataset_outdir(c):
    "Returns _dataset_outdir of `c['dataset']`, for easier access via `c` below."
    return c["dataset"]["outdir"]


def dataset_id(c):
    return c["dataset"]["id"]


# Helper for running test runs on a subset of the data, togglable via the `--test` cli flag
def wrap_test_run(take_n=2):
    def deco(nestables_fn):
        def f(c):
            nestables = nestables_fn(c)
            nestables = nestables[:take_n]
            return nestables

        f.__name__ = nestables_fn.__name__
        return f if options["test_run"] else nestables_fn

    return deco


# Subject nest level
# ------------------

# Whenever we iterate over partitions, we always want to assume there could be an `other-partitions` mapping,
# and iterate over all these things, while tracking their other_id keys in the `other-partitions` dict.
def get_partitions(node):
    parts = []
    if node.get("partition-file"):
        parts.append(node)
    if node.get("other-partitions"):
        parts += [
            sconsutils.merge_dicts(part, {"other_id": other_id})
            for other_id, part in node["other-partitions"].items()
            if part.get("partition-file")
        ]
    return parts


def keep_sample(sample):
    return (len(get_partitions(sample)) > 0 and options["only_seeds"] is None) or [
        seed
        for seed_id, seed in sample.get("seeds", {}).items()
        if seed.get("partition-file")
        and (options["only_seeds"] is None or seed_id in options["only_seeds"])
    ]


def samples(c):
    return {
        sample_id: sample
        for sample_id, sample in c["dataset"]["samples"].items()
        if keep_sample(sample)
    }


# See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
@w.add_nest(label_func=str)
@wrap_test_run(take_n=2)
def subject(c):
    return list(
        set(
            sconsutils.get_in(sample, ["meta", "subject"])
            for sample_id, sample in samples(c).items()
        )
    )


# Initialize sample nest
# -----------------------------

# Samples can either point to a specific timepoint through the yaml meta, or can have a "merged" attribute
# value there, if it is a sample composed of many timepoints. These metadata will be processed accordingly.

# Samples can have partitions, and they can also have other-partitions, and seeded partitions.
# These are handled in separate nest loops below, with a pop in between.

# There may eventually be some required arguments here as this is where we get our locus and isotype and such
# See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
@w.add_nest(
    metadata=lambda c, d: sconsutils.merge_dicts(
        d.get("meta", {}), {"id": d["id"], "seeds": None, "meta": None}
    )
)
@wrap_test_run(take_n=2)
def sample(c):
    # Make sure to add timepoints here as necessary
    return [
        sconsutils.merge_dicts(sample, {"id": sample_id})
        for sample_id, sample in samples(c).items()
        if sconsutils.get_in(sample, ["meta", "subject"]) == c["subject"]
    ]


def locus(c):
    sample = c["sample"]
    locus = sample.get("locus") or sample.get("meta").get("locus")
    return locus


# Initialize seed nest
# --------------------

# We start with the seed nest level.
# Herein we'll loop over partitions and fetch the seeded clusters from the partitions of interest, as defined below.
# Eventually, we'll pop off this seed nest level so we can renest these partitions and clusters directly from the sample nest level.

# Initialize our first sub dataset nest level
# See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
@w.add_nest(
    metadata=lambda c, d: sconsutils.merge_dicts(d.get("meta", {}), {"id": d["id"]})
)
# would like to have a lower number here but sometimes we get no good clusters for the first two seeds?
@wrap_test_run(take_n=3)
def seed(c):
    return [
        sconsutils.merge_dicts(seed, {"id": seed_id})
        for seed_id, seed in samples(c)[c["sample"]["id"]].get("seeds", {}).items()
        if seed.get("partition-file")
        and (options["only_seeds"] is None or seed_id in options["only_seeds"])
    ]


# Some accessor helpers


def timepoint(c):
    return sconsutils.get_in(c, ["sample", "meta", "timepoint"])


def is_merged(c):
    return timepoint(c) == "merged"


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


def get_alt_naive_probabilities(annotation):
    alternatives = annotation.get("alternative-annotations", {})
    naive_probabilities = alternatives.get("naive-seqs") if alternatives else None
    return (
        naive_probabilities
        if naive_probabilities and len(naive_probabilities) > 0
        else None
    )


def partition_steps(cp):
    if len(cp.partitions) == 0:
        return []
    return (
        range(len(cp.partitions))
        if options["process_all_partis_partition_steps"]
        else [cp.i_best]
    )


def partition_metadata(part, annotation_list, cp, i_step, seed=None, other_id=None):
    clusters = cp.partitions[i_step]
    seed_cluster_annotation = None
    if seed:
        try:
            seed_cluster_annotation = process_partis.choose_cluster(
                part["partition-file"], annotation_list, cp, i_step
            )
        except ValueError, e:
            warn(
                "Due to the following error: {},\n no annotation was found in {} for seed cluster {}. Skipping this cluster".format(
                    " ".join(e.args), part["partition-file"], seed
                )
            )
            return None

    meta = {
        "id": ("seed-" if seed else "unseeded-")
        + (other_id + "-" if other_id else "")
        + "part-"
        + str(i_step),
        "clusters": clusters,
        "step": i_step,
        "n_clusters": len(clusters),
        "largest_cluster_size": max(map(len, clusters)),
        "logprob": cp.logprobs[i_step],
        "partition-file": part["partition-file"],
        "seed_cluster_annotation": seed_cluster_annotation,
    }
    return sconsutils.merge_dicts(meta, part.get("meta") or {})


def min_cluster_size(is_seed_cluster=False):
    """A function to track the different min cluster sizes for seed clusters vs unseeded.
    Numbers are somewhat arbitrary, though we often dont see smaller especially without processing all partition steps."""
    return 4 if is_seed_cluster else 6


def meets_cluster_size_reqs(unique_ids, is_seed_cluster=False):
    """Takes cluster's list of unique ids and returns True if cluster is within size requirements here, otherwise false.
    Upper bound is computational-resource-driven and lower bounds are practical use based."""
    size = len(unique_ids)
    if size > 10000:
        message = """
                    cluster size limit exceeded: {}
                    clusters are limited to 10,000 sequences in order to make tree building
                    possible in reasonable time and not exceed memory resources. Downsample
                    this cluster or rerun partis with --max-cluster-size.
                  """.format(
            len(unique_ids)
        )
        if not options["skip_large_clusters"]:
            message += " To skip large clusters and build the rest, run again with --skip-large-clusters."
            raise Exception(message)
        warn("Skipping a cluster!" + message)
        return False
    return size >= min_cluster_size(is_seed_cluster=is_seed_cluster)


def valid_cluster(annotation_list, part, unique_ids, is_seed_cluster=False):
    """Reads the corresponding cluster annotation and return True iff after applying our health metric filters
    we still have greater than 2 sequences (otherwise, we can't build a tree downstream)."""
    for line in annotation_list:
        if line.get("unique_ids") == unique_ids:
            functional_seqs_uids = [
                uid
                for iseq, uid in enumerate(line["unique_ids"])
                if partisutils.is_functional(line, iseq)
            ]
            return meets_cluster_size_reqs(
                functional_seqs_uids, is_seed_cluster=is_seed_cluster
            )
    raise Exception(
        "couldn't find requested uids %s in %s" % (unique_ids, part["partition-file"])
    )


def valid_seed_partition(
    annotation_list, cp, part, i_step, seed_id, max_size_to_check=10
):
    """If seed cluster size is less than max_size_to_check, read the corresponding cluster annotation and return True iff after applying our health metric filters
    we still have greater than 2 sequences (otherwise, we can't build a tree downstream)."""
    seed_cluster_unique_ids = seed_cluster(cp, i_step, seed_id)
    if seed_cluster_unique_ids is not None and meets_cluster_size_reqs(
        seed_cluster_unique_ids, is_seed_cluster=True
    ):
        if len(seed_cluster_unique_ids) > max_size_to_check:
            return True
        return valid_cluster(
            annotation_list, part, seed_cluster_unique_ids, is_seed_cluster=True
        )
    return False


# Try to read partition file; If fails, it is possibly because it's empty. Catch that case and warn
def read_partition_file(part, c):
    try:
        glfo, annotation_list, cpath = process_partis.read_partis_output(
            part["partition-file"], c["sample"]["glfo-dir"], locus(c)
        )
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        print("".join((" " * 8) + line for line in lines))
        warn(
            "Unable to parse partition file (see error above, ommitting from results): {}".format(
                part
            )
        )
        return []
    return annotation_list, cpath


# note we elide the nested partitions > clusters lists (as well as the seed cluster annotation)
# so as not to kill tripl when it tries to load them as a value and can't hash
# See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
@w.add_nest(
    metadata=lambda c, d: {"clusters": "elided", "seed_cluster_annotation": "elided"}
)
def partition(c):
    """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
    for actual analysis."""
    keep_partitions = []
    seed_id = c["seed"]["id"]
    if options["only_seeds"] is not None and seed_id not in options["only_seeds"]:
        return []
    for part in get_partitions(c["seed"]):
        annotation_list, cp = read_partition_file(part, c)
        if cp:
            for i_step in partition_steps(cp):
                meta = partition_metadata(
                    part,
                    annotation_list,
                    cp,
                    i_step,
                    seed=seed_id,
                    other_id=part.get("other_id"),
                )
                if valid_seed_partition(annotation_list, cp, part, i_step, seed_id):
                    keep_partitions.append(meta)
    return keep_partitions


# The cluster level
# -----------------

# For seeded clusters we only process the seed containing cluster.
# See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
@w.add_nest(
    label_func=lambda d: d["id"],
    metadata=lambda c, d: {"annotation": "elided", "naive_probabilities": "elided"},
)
def cluster(c):
    seed_cluster_annotation = c["partition"]["seed_cluster_annotation"]
    return [
        {
            "id": "seed-cluster",
            "seed_name": c["seed"]["id"],
            "size": len(seed_cluster_annotation["unique_ids"]),
            "unique_ids": seed_cluster_annotation["unique_ids"],
            "annotation": seed_cluster_annotation,
            "naive_probabilities": get_alt_naive_probabilities(seed_cluster_annotation),
        }
    ]


def add_cluster_analysis(w):
    @w.add_target(name="path")
    def path_fn(outdir, c):
        return outdir

    @w.add_target()
    def partis_cluster_fasta(outdir, c):
        return env.Command(
            path.join(outdir, "unfiltered_partis_cluster.fa"),
            c["partition"]["partition-file"],
            "process_partis.py"
            + " --partition-file $SOURCE"
            + " --partition {}".format(c["partition"]["step"])
            + (
                " --glfo-dir " + c["sample"]["glfo-dir"]
                if partisutils.getsuffix(c["partition"]["partition-file"]) == ".csv"
                else ""
            )
            + " --ignore-seed-indels "
            + " --locus "
            + locus(c)
            + " --paths-relative-to "
            + dataset_outdir(c)
            + " --inferred-naive-name "
            + options["inferred_naive_name"]
            + (
                " --cluster {}".format(c["cluster"]["sorted_index"])
                if not c.get("seed")
                else ""
            )
            + " --seqs-out $TARGET",
        )

    @w.add_metadata()
    def _process_partis(outdir, c):
        sources = [c["partition"]["partition-file"]]
        perseq_metafile = c["sample"].get("per-sequence-meta-file")
        if perseq_metafile:
            sources.append(perseq_metafile)
        cluster_seqs_fname = "cluster_seqs.fa"
        if options["match_indel_in_uid"]:
            cluster_seqs_fname = "{}_indel_filtered_cluster_seqs.fa".format(
                options["match_indel_in_uid"]
            )
        return env.Command(
            [
                path.join(outdir, x)
                for x in [
                    "partis_metadata.json",
                    cluster_seqs_fname,
                    "partis_seqmeta.csv",
                ]
            ],
            sources,
            "process_partis.py"
            + " --remove-stops --remove-frameshifts --remove-mutated-invariants"
            + (" --indel-reversed-seqs " if not options["preserve_indels"] else "")
            + " --partition-file ${SOURCES[0]}"
            + " --partition {}".format(c["partition"]["step"])
            + (" --upstream-seqmeta ${SOURCES[1]}" if perseq_metafile else "")
            + (
                " --glfo-dir " + c["sample"]["glfo-dir"]
                if partisutils.getsuffix(c["partition"]["partition-file"]) == ".csv"
                else ""
            )
            + " --locus "
            + locus(c)
            + " --paths-relative-to "
            + dataset_outdir(c)
            + " --namespace cft.cluster"
            + " --inferred-naive-name "
            + options["inferred_naive_name"]
            + (
                (" --show-indel-in-trees " + options["show_indel_in_trees"])
                if options["show_indel_in_trees"] is not None
                else ""
            )
            + (
                (" --match-indel-in-uid " + options["match_indel_in_uid"])
                if options["match_indel_in_uid"] is not None
                else ""
            )
            + (" --ignore-seed-indels" if options["ignore_seed_indels"] else "")
            + (
                (" --always-include " + ",".join(c["sample"]["seeds"]))
                if c["sample"].get("seeds")
                else ""
            )
            + (
                " --cluster {}".format(c["cluster"]["sorted_index"])
                if not c.get("seed")
                else ""
            )
            + " --cluster-meta-out ${TARGETS[0]}"
            + " --seqs-out ${TARGETS[1]}"
            + " --seqmeta-out ${TARGETS[2]}",
        )

    @w.add_target(ingest=True)
    def partis_metadata(outdir, c):
        return c["_process_partis"][0]

    @w.add_target()
    def inseqs(outdir, c):
        return c["_process_partis"][1]

    @w.add_target()
    def partis_seqmeta(outdir, c):
        return c["_process_partis"][2]

    # Partis alternative naives
    # -------------------------

    @w.add_target()
    def alternative_naive_probabilities(outdir, c):
        """
        Write partis alternative naives to a fasta in order of probability
        """
        if c["cluster"]["naive_probabilities"] is not None:
            cluster_name = c["cluster"].get("seed_name", c["cluster"]["id"])

            naives_sorted_by_prob = list(
                sorted(
                    c["cluster"]["naive_probabilities"],
                    key=lambda x: x[1],
                    reverse=True,
                )
            )

            def write_naive_fastas(target, source, env):
                """
                This is the action for this target. Because it is a function, not a file being executed,
                the rules SCons follows for determining whether to rebuild this target are less well defined.
                """
                targets = [str(fname) for fname in target]
                with open(targets[0], "w") as ranked_fasta, open(
                    targets[1], "w"
                ) as aa_ranked_fasta:
                    for rank, (naive_seq, probability) in enumerate(
                        naives_sorted_by_prob
                    ):
                        aa_seq = translate_seqs.translate(naive_seq)
                        ranked_fasta.write(
                            ">%s\n%s\n"
                            % (
                                "naive_{}_probability_{}".format(rank, probability),
                                naive_seq,
                            )
                        )
                        aa_ranked_fasta.write(
                            ">%s\n%s\n"
                            % (
                                "naive_{}_probability_{}".format(rank, probability),
                                aa_seq,
                            )
                        )

            naive_probs_fname = path.join(
                outdir, "ranked_naive_probabilities_%s.fasta" % cluster_name
            )
            aa_naive_probs_fname = path.join(
                outdir, "ranked_aa_naive_probabilities_%s.fasta" % cluster_name
            )

            return env.Command(
                [naive_probs_fname, aa_naive_probs_fname],
                c["partition"]["partition-file"],
                write_naive_fastas,
            )

    @w.add_target()
    def alternative_naive_logo_plots(outdir, c):
        """
        Create logo plot according to probabilities
        """
        if c["cluster"]["naive_probabilities"] is not None:

            annotation = c["cluster"]["annotation"]
            cluster_name = c["cluster"].get("seed_name", c["cluster"]["id"])

            aa_input_fasta_path = str(c["alternative_naive_probabilities"][1])
            aa_cdr3_start, aa_cdr3_end = (
                int(annotation["codon_positions"]["v"] / 3),
                int((annotation["codon_positions"]["j"] + 3) / 3),
            )
            aa_naive_len = int(len(annotation["naive_seq"]) / 3)
            logo_out, cdr3_logo_out = (
                "naive_logo_{}.png".format(cluster_name),
                "naive_logo_cdr3_{}.png".format(cluster_name),
            )
            logo_plots = env.Command(
                [path.join(outdir, logo_out), path.join(outdir, cdr3_logo_out)],
                aa_input_fasta_path,
                "python bin/create_partis_naive_logo.py $SOURCE"
                + " --aa-cdr3-start=%d" % aa_cdr3_start
                + " --aa-cdr3-end=%d" % aa_cdr3_end
                + " --aa-naive-len=%d" % aa_naive_len
                + " --logo-fname=%s" % logo_out
                + " --cdr3-logo-fname=%s" % cdr3_logo_out
                + " --outdir=%s" % outdir,
            )

            env.Depends(logo_plots, "bin/create_partis_naive_logo.py")
            return logo_plots

    # Sequence Alignment
    # ------------------

    backtrans_align.add(env, w, options)

    # Trees
    # ---------------------------------

    # use fasttree to make newick tree from sequences
    @w.add_target()
    def fasttree(outdir, c):
        return env.SRun(
            path.join(outdir, "fasttree.nwk"),
            c["aligned_inseqs"],
            "FastTree -nt -quiet $SOURCE > $TARGET 2> $TARGET-.log",
        )

    if options["show_indel_in_trees"]:

        @w.add_target()
        def indel_matching_ids(outdir, c):
            def write_indel_matching_ids(target, source, env):
                with open(str(source[0])) as seqmeta:
                    lines = list(csv.DictReader(seqmeta))
                    uids = [
                        row["unique_id"]
                        for row in lines
                        if row["indel_match"] == "True"
                    ]
                with open(str(target[0]), "w") as indel_matches_file:
                    for uid in uids:
                        indel_matches_file.write(uid + "\n")

            return env.Command(
                os.path.join(
                    outdir,
                    "uids_matching_indel_{}.txt".format(options["show_indel_in_trees"]),
                ),
                c["partis_seqmeta"],
                write_indel_matching_ids,
            )

        @w.add_target()
        def indel_fasttree_svg(outdir, c):
            """create graphic showing indel-matching (matching uid passed in --show-indel-in-trees) in the fasttree as red"""
            indel_svg = env.Command(
                path.join(
                    outdir,
                    "{}_indel_fasttree.svg".format(options["show_indel_in_trees"]),
                ),
                [c["fasttree"], c["indel_matching_ids"]],
                # The `-` at the start here tells scons to ignore if it doesn't build; this may occasionally be
                # the case for large clusters. Also, redirect stdin/out to dev/null because the errors messages
                # here can be pretty noisy.
                "- xvfb-run -a bin/annotate_tree.py $SOURCES "
                + " --naive %s" % options["inferred_naive_name"]
                + (" --seed " + c["seed"]["id"] if "seed" in c else "")
                + " --set-root"
                + " --size 100"
                + " --output-path $TARGET &>> /dev/null",
            )
            env.Depends(indel_svg, "bin/annotate_tree.py")
            return indel_svg

    # See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
    @w.add_nest(metadata=lambda c, d: d)
    def reconstruction(c):
        return [
            {
                "id": prune_strategy + "-" + asr_prog,
                "prune_strategy": prune_strategy,
                "asr_prog": asr_prog,
                "prune_count": 100,
            }
            for prune_strategy, asr_prog in itertools.product(
                ["min_adcl", "seed_lineage"] if "seed" in c else ["min_adcl"],
                ["raxml_ng"] if not options["run_dnaml"] else ["dnaml"],
            )
        ]

    # calculate list of sequences to be pruned
    @w.add_target()
    def pruned_ids(outdir, c):
        tgt = path.join(outdir, "pruned_ids.txt")
        # This whole thing is a safety mechanism to prevent prune files with 0 sequences from hanging around,
        # which can happen from failed builds
        try:
            remove = False
            with open(tgt, "r") as fh:
                remove = len(fh.readlines()) == 0
            if remove:
                os.remove(tgt)
        except:
            pass
        recon = c["reconstruction"]
        builder = (
            fun.partial(env.SRun, srun_args="`minadcl_srun_args.py $SOURCE`")
            if recon["prune_strategy"] == "min_adcl"
            else env.Command
        )
        return builder(
            tgt,
            c["fasttree"],
            "prune.py -n "
            + str(recon["prune_count"])
            + (
                (" --always-include " + ",".join(c["sample"]["seeds"]))
                if c["sample"].get("seeds")
                else ""
            )
            + " --strategy "
            + recon["prune_strategy"]
            + " --naive %s" % options["inferred_naive_name"]
            + (" --seed " + c["seed"]["id"] if "seed" in c else "")
            + " $SOURCE $TARGET",
        )

    if options["fasttree_png"]:

        @w.add_target()
        def pruned_fasttree_png(outdir, c):
            """create png showing included seqs (kept in pruning) as red"""
            if c["cluster"].get("size") < 4500:
                pruned_cluster_fasttree_png = env.Command(
                    path.join(outdir, "pruned_cluster_fasttree.png"),
                    [c["fasttree"], c["pruned_ids"]],
                    # The `-` at the start here tells scons to ignore if it doesn't build; this may occasionally be
                    # the case for large clusters. Also, redirect stdin/out to dev/null because the errors messages
                    # here can be pretty noisy.
                    "- xvfb-run -a bin/annotate_tree.py $SOURCES "
                    + " --naive %s" % options["inferred_naive_name"]
                    + (" --seed " + c["seed"]["id"] if "seed" in c else "")
                    + " --output-path $TARGET &>> /dev/null",
                )
                env.Depends(pruned_cluster_fasttree_png, "bin/annotate_tree.py")
                return pruned_cluster_fasttree_png

    @w.add_target()
    def cluster_mapping(outdir, c):
        if c["reconstruction"]["prune_strategy"] == "min_adcl":
            return env.SRun(
                path.join(outdir, "cluster_mapping.csv"),
                [c["fasttree"], c["pruned_ids"]],
                "minadcl_clusters.py $SOURCES $TARGET",
                srun_args="`minadcl_clusters_srun_args.py $SOURCE`",
            )

    # prune out sequences to reduce taxa, making sure to cut out columns in the alignment that are now entirely
    # gaps from insertions in sequences that have been pruned out.
    @w.add_target()
    def pruned_seqs(outdir, c):
        return env.Command(
            path.join(outdir, "pruned.fa"),
            [c["pruned_ids"], c["aligned_inseqs"]],
            "seqmagick convert --include-from-file $SOURCES - | "
            + "seqmagick convert --squeeze - $TARGET",
        )

    if options["write_linearham_yaml_input"]:

        @w.add_target()
        def pruned_partis_outfile(outdir, c):
            if "seed" in c:
                yaml_format = (
                    partisutils.getsuffix(c["partition"]["partition-file"]) == ".yaml"
                )
                subset_partis_outfile = env.Command(
                    path.join(outdir, "pruned_partis_output.yaml"),
                    [c["partition"]["partition-file"], c["pruned_ids"]],
                    "python bin/write_subset_partis_outfile.py $SOURCES $TARGET"
                    + " --partition-step={}".format(c["partition"]["step"])
                    + " --original-cluster-unique-ids={}".format(
                        ":".join(c["cluster"]["unique_ids"])
                    )
                    + " --sw-cache={}".format(c["sample"]["sw-cache"])
                    + (
                        " --glfo-dir={}".format(c["sample"]["glfo-dir"])
                        if not yaml_format
                        else ""
                    )
                    + (" --locus={}".format(locus(c)) if not yaml_format else ""),
                )
                env.Depends(subset_partis_outfile, "bin/write_subset_partis_outfile.py")
                return subset_partis_outfile

        @w.add_target()
        def linearham_base_command(outdir, c):
            """ This allows us to not have to manually build the minimum necessary commmad for running this seed
            cluster through Linearham. Eventually we may want to run this command here in the CFT SCons pipeline,
            but for now we make life a little easier by just being able to copy this command to run Linearham."""
            if "seed" in c:
                return env.Command(
                    path.join(outdir, "linearham_base_command.txt"),
                    c["pruned_partis_outfile"],
                    'echo "scons --run-linearham --template-path=templates/revbayes_template.rev '
                    + " --parameter-dir={}".format(c["sample"]["parameter-dir"])
                    + " --partis-yaml-file={}".format(
                        path.join(os.getcwd(), str(c["pruned_partis_outfile"][0]))
                    )
                    + ' --seed-seq={}" > $TARGET'.format(  # get cwd for absolute path
                        c["cluster"]["seed_name"]
                    ),
                )

    @w.add_target()
    def tip_seqmeta(outdir, c):
        """The merge of process_partis output with pre sequence metadata spit out by datascripts containing
        timepoint mappings. Base input multiplicity is coded into the original input sequence names from vlad as N-M,
        where N is the ranking of vlads untrimmed deduplication, and M is the multiplicity of said deduplication."""
        # This option controls which sequences get joined on in the merge for the partis_seqmeta file, which has
        # orig/new names, joined on sequence from the other file
        sources = {
            "--partis-seqmeta": c["partis_seqmeta"],
            "--cluster-mapping": c["cluster_mapping"]
            if c["reconstruction"]["prune_strategy"] == "min_adcl"
            else None,
            "--pruned-ids": c["pruned_ids"]
            if c["reconstruction"]["prune_strategy"] == "seed_lineage"
            else None,
        }
        sources = {k: v for k, v in sources.items() if v}
        base_call = "aggregate_minadcl_cluster_multiplicities.py "
        for i, (k, v) in enumerate(sources.items()):
            base_call += k + " ${SOURCES[" + str(i) + "]} "
        return env.Command(
            path.join(outdir, "tip_seqmeta.csv"),
            sources.values(),
            base_call + "$TARGET",
        )

    # Run raxml-ng/dnaml
    @w.add_target()
    def _asr(outdir, c):
        "run raxml-ng and/or dnaml(from phylip package) to create tree with inferred sequences at internal nodes"
        asr_prog = c["reconstruction"]["asr_prog"]
        if asr_prog == "raxml_ng":
            raxml_base_cmd = (
                "raxml-ng --model GTR+G --threads 2 --redo --force msa_allgaps"
                + " --msa {}".format(str(c["pruned_seqs"][0]))
            )
            # run once to infer tree
            basename = "treeInference"
            log, raxml_best_tree = env.SRun(
                [
                    path.join(outdir, basename + ".raxml." + ext)
                    for ext in ["log", "bestTree"]
                ],
                c["pruned_seqs"],
                raxml_base_cmd
                + " --prefix {}".format(path.join(outdir, basename))
                + " > ${TARGETS[0]}",
            )
            # run again to reconstruct ancestral sequences (ASR)
            basename = "ASR"
            log, raxml_asr_tree, raxml_asr_seqs = env.SRun(
                [
                    path.join(outdir, basename + ".raxml." + ext)
                    for ext in ["log", "ancestralTree", "ancestralStates"]
                ],
                [c["pruned_seqs"], raxml_best_tree],
                raxml_base_cmd
                + " --prefix {}".format(path.join(outdir, basename))
                + " --ancestral"
                + " --tree ${SOURCES[1]}"
                + " > ${TARGETS[0]}",
            )
            rooted_asr_tree, asr_seqs, ancestors_naive_and_seed = env.Command(
                [
                    path.join(outdir, basename + "." + ext)
                    for ext in ["nwk", "fa", "ancestors_naive_and_seed.fa"]
                ],
                [raxml_asr_tree, raxml_asr_seqs, c["pruned_seqs"]],
                "bin/parse_raxmlng.py"
                + " --tree ${SOURCES[0]}"
                + " --asr-seq ${SOURCES[1]}"
                + " --input-seq ${SOURCES[2]}"
                + " --outbase {}".format(path.join(outdir, basename))
                + " --inferred-naive-name {}".format(options["inferred_naive_name"])
                + (" --seed " + c["seed"]["id"] if "seed" in c else ""),
            )
            return [rooted_asr_tree, asr_seqs, ancestors_naive_and_seed]
        elif asr_prog == "dnaml":
            # Seqname translations for reinterpretting dnaml output in terms of original seqnames, due to phylip name length constraints.
            pruned_seqs_phylip, seqname_mapping = env.Command(
                [path.join(outdir, x) for x in ("pruned.phy", "seqname_mapping.csv")],
                c["pruned_seqs"],
                "make_phylip.py $SOURCE $TARGETS --inferred-naive-name "
                + options["inferred_naive_name"],
            )
            basename = "asr"
            config = env.Command(
                path.join(outdir, asr_prog + ".cfg"),
                pruned_seqs_phylip,
                "python bin/mkconfig.py $SOURCE " + asr_prog + "> $TARGET",
            )
            phylip_out = env.SRun(
                path.join(outdir, "outfile"),
                config,
                "cd "
                + outdir
                + " && rm -f outtree && "
                + asr_prog
                + " < $SOURCE.file > "
                + asr_prog
                + ".log",
            )
            # Manually depend on phylip sequences so that we rerun dnaml if the input sequences change; without this, dnaml will only get rerun if one of the targets are removed or if the iput asr_config file is changed.
            env.Depends(phylip_out, pruned_seqs_phylip)
            # process the phylip output
            tgt = env.Command(
                [
                    path.join(outdir, basename + "." + ext)
                    for ext in ["nwk", "fa", "ancestors_naive_and_seed.fa"]
                ],
                [seqname_mapping, phylip_out, c["tip_seqmeta"]],
                # Note that adding `-` at the beginning of this command string can keep things running if
                # there's an error trying to build a tree from 2 sequences; should be filtered prior now...
                "xvfb-run -a bin/process_asr.py"
                + (" --seed " + c["seed"]["id"] if "seed" in c else "")
                + " --outdir "
                + outdir
                + " --basename "
                + basename
                + " --inferred-naive-name "
                + options["inferred_naive_name"]
                + " --seqname-mapping $SOURCES",
            )
            asr_tree, asr_seqs, ancestors_naive_and_seed = tgt
            # manually depnd on this because the script isn't in first position
            env.Depends(tgt, "bin/process_asr.py")
            return [asr_tree, asr_seqs, ancestors_naive_and_seed]

    @w.add_target(ingest=True)
    def asr_tree(outdir, c):
        return c["_asr"][0]

    @w.add_target(ingest=True)
    def asr_seqs(outdir, c):
        return c["_asr"][1]

    @w.add_target(ingest=True)
    def ancestors_naive_and_seed(outdir, c):
        return c["_asr"][2]

    if options["show_indel_in_trees"]:

        @w.add_target()
        def indel_mltree_svg(outdir, c):
            """create graphic showing indel-matching (matching uid passed in --show-indel-in-trees) in the max-likelihood tree as red"""
            indel_svg = env.Command(
                path.join(
                    outdir, "{}_indel_mltree.svg".format(options["show_indel_in_trees"])
                ),
                [c["asr_tree"], c["indel_matching_ids"]],
                # The `-` at the start here tells scons to ignore if it doesn't build; this may occasionally be
                # the case for large clusters. Also, redirect stdin/out to dev/null because the errors messages
                # here can be pretty noisy.
                "- xvfb-run -a bin/annotate_tree.py $SOURCES "
                + " --naive %s" % options["inferred_naive_name"]
                + (" --seed " + c["seed"]["id"] if "seed" in c else "")
                + " --size 100"
                + " --output-path $TARGET &>> /dev/null",
            )
            env.Depends(indel_svg, "bin/annotate_tree.py")
            return indel_svg

    @w.add_target()
    def observed_ancestors(outdir, c):
        blast_db_files = ["blast_db.{}".format(ext) for ext in ("nin", "nhr", "nsq")]
        blast_results_files = [
            "blast_results.{}.tsv".format(method) for method in ("blastn")
        ]
        targets = blast_db_files + blast_results_files
        observed_ancestors_output = env.Command(
            [path.join(outdir, fname) for fname in targets],
            [c["partis_cluster_fasta"], c["ancestors_naive_and_seed"]],
            "python bin/blast.py $SOURCES --outdir {}".format(outdir),
        )
        env.Depends(observed_ancestors_output, "bin/blast.py")
        return observed_ancestors_output

    @w.add_target()
    def selection_metrics(baseoutdir, c):
        outdir = path.join(baseoutdir, "selection-metrics")

        def fix_file(target, source, env):
            extra_ambig_chars = "RYSWKMBDHV-"
            target = str(target[0])
            source = str(source[0])
            ambiguous_char_translations = string.maketrans(
                extra_ambig_chars, "".join("N" for _ in range(len(extra_ambig_chars)))
            )
            seqfos = partisutils.read_fastx(source)
            with open(target, "w") as nfile:
                for sfo in seqfos:
                    nfile.write(
                        ">%s\n%s\n"
                        % (
                            sfo["name"],
                            sfo["seq"].translate(ambiguous_char_translations),
                        )
                    )

        new_partis_infname = path.join(baseoutdir, "asr-with-only-Ns.fa")
        new_partis_infile = env.Command(new_partis_infname, c["asr_seqs"], fix_file)
        plotdir_name = "partis-plot"
        tree_metrics = env.Command(
            [
                path.join(outdir, outfile)
                for outfile in [
                    "partition.yaml",
                    "partis.std.out.log",
                    path.join(plotdir_name, "inferred-tree-metrics/overview.html"),
                ]
            ],
            [c["asr_seqs"], c["sample"]["parameter-dir"], c["asr_tree"]],  # sources
            os.path.join(partis_path, "bin/partis")
            + " annotate "
            + " --locus "
            + locus(c)
            + " --get-selection-metrics "
            + " --all-seqs-simultaneous "
            + " --min-selection-metric-cluster-size %d "
            % min_cluster_size(is_seed_cluster=("seed" in c))
            + " --no-partition-plots "
            + " --plotdir "
            + path.join(outdir, plotdir_name)
            + " --infname ${SOURCES[0]} "
            + " --parameter-dir ${SOURCES[1]} "
            + " --treefname ${SOURCES[2]} "
            + " --outfname ${TARGETS[0]} "
            + " > ${TARGETS[1]}",
        )
        return tree_metrics

    @w.add_target(
        ingest=True,
        attr_map={
            "bio.seq:id": "sequence",
            "cft.timepoint:id": "timepoint",
            "cft.seq:timepoint": "timepoint",
            "cft.seq:timepoints": "timepoints",
            "cft.seq:cluster_timepoints": "cluster_timepoints",
            "cft.seq:multiplicity": "multiplicity",
            "cft.seq:cluster_multiplicity": "cluster_multiplicity",
            "cft.seq:timepoint_multiplicities": "timepoint_multiplicities",
            "cft.seq:cluster_timepoint_multiplicities": "cluster_timepoint_multiplicities",
            "cft.seq:affinity": "affinity",
            "cft.tree.node:lbi": "lbi",
            "cft.tree.node:lbr": "lbr",
        },
    )
    def seqmeta(outdir, c):
        return env.Command(
            path.join(outdir, "seqmeta.csv"),
            [c["tip_seqmeta"], c["selection_metrics"][0]],
            "merge_selection_metrics.py $SOURCES $TARGET",
        )

    @w.add_target(ingest=True)
    def cluster_aa(outdir, c):
        return env.Command(
            path.join(outdir, "cluster_aa.fa"),
            c["asr_seqs"],
            "sed 's/\?/N/g' $SOURCE | seqmagick convert --translate dna2protein - $TARGET",
        )


# This calls the above function and runs the whole anaylsis for the seeded clusters which have been added to nestly's nest structure (stored in the variable called 'w')
add_cluster_analysis(w)

# Then we call pop() here to go up one nest level and add the unseeded clusters next to the seeded clusters in the build hierarchy (DAG) so we can run a slightly differeny analysis for unseeded clusters.
w.pop("seed")

# Unseeded cluster analysis
# -------------------------
# Now we define a function that builds the cluster analysis defined above for the clusters from unseeded partitions.
# We do this so that the functions names we use for building things with nestly don't overlap with those defined for the seeded analysis above.
# Note that the "partition" and "cluster" nest levels are defined here separately from above in order to apply different criteria for building the set of unseeded clusters.


def add_unseeded_analysis(w):

    # Setting up the partition nest level

    # Each c["partition"] value actually points to the annotations for that partition... a little weird but...
    # We elide the annotations here because we pass them through in order to avoid re-reading them from
    # the partition file and because we don't need to write them to the metadata for the partition.
    # See https://github.com/matsengrp/cft/pull/270#discussion_r267502415 for details on why we decided to do things this way.
    # See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
    @w.add_nest(metadata=lambda c, d: {"clusters": "elided", "cp": "elided"})
    def partition(c):
        """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
        for actual analysis."""
        keep_partitions = []
        for partition_run in get_partitions(c["sample"]):
            annotation_list, cp = read_partition_file(partition_run, c)
            if cp:
                keep_partitions += [
                    partition_metadata(
                        partition_run,
                        annotation_list,
                        cp,
                        i_step,
                        other_id=partition_run.get("other_id"),
                    )
                    for i_step in partition_steps(cp)
                ]
        return keep_partitions

    # Add cluster nesting level
    # See https://nestly.readthedocs.io/en/latest/index.html for a definition of add_nest and more info on the "nestly" package which governs the nesting levels of things getting built in this pipeline
    @w.add_nest(
        label_func=lambda d: d["id"],
        metadata=lambda c, d: {
            "unique_ids": "elided",
            "annotation": "elided",
            "naive_probabilities": "elided",
        },
    )
    def cluster(c):
        part = c["partition"]
        clusters = []
        annotation_list = None
        # Sort by len (dec) and apply index i
        for i, unique_ids in enumerate(sorted(part["clusters"], key=len, reverse=True)):
            # Select top N or any matching seeds of interest
            if (i < options["depth"]) and meets_cluster_size_reqs(unique_ids):
                if annotation_list is None:
                    # Here we reread the partition file instead of caching annotations of the partition along with its metadata above in partition_metadata. This saves on memory and slows the process down, but we are restricted by memory use more than time at the moment.
                    annotation_list, cp = read_partition_file(part, c)
                if valid_cluster(annotation_list, part, unique_ids):
                    cluster_annotation = process_partis.choose_cluster(
                        part["partition-file"], annotation_list, cp, part["step"], i
                    )
                    # It seems like we might only need to check that one of these clusters has alternative naive info and then  we could assume it is the case for
                    # all of them (unless --queries was set for --calculate-alternative-naive-seqs). Leaving it as is for now but may speed up the SConstruct process to do this later. (EH)
                    naive_probabilities = get_alt_naive_probabilities(
                        cluster_annotation
                    )
                    cluster_meta = {
                        "id": "clust-" + str(i),
                        "sorted_index": i,
                        "unique_ids": unique_ids,
                        "size": len(unique_ids),
                        "annotation": cluster_annotation,
                        "naive_probabilities": naive_probabilities,
                    }
                    clusters.append(cluster_meta)
        return clusters

    # do the cluster analysis defined above for all the unseeded clusters
    add_cluster_analysis(w)


add_unseeded_analysis(w)

w.pop("subject")

# Write metadata
@w.add_target()
def metadata_snapshot(outdir, c):
    return env.Command(
        path.join(
            outdir,
            time.strftime("%Y-%m-%d")
            + ("-test" if options["test_run"] else "")
            + "-metadata.json",
        ),
        path.join(outdir, "metadata.json"),
        "cp $SOURCE $TARGET",
    )


w.pop("dataset")
