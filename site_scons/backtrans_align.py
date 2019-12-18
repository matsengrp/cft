#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import options
import os
import sys

# submodule checkout; this is needed for bin/process_partis.py)
default_partis_path = os.path.join(os.getcwd(), "partis")
partis_path = os.environ.get("PARTIS", default_partis_path)
sys.path.append(os.path.join(partis_path, "python"))
import utils as partisutils

# What follows is a rather convoluted backtranslation strategy for sequence alignment.
# The basic idea is that we translate our sequences, align them (taking advantage of coding information),
# and from this infer the a nucleotide alignment using `seqmagick backtrans-align`.
# The problem is that `seqmagick` rather finicky here about sequence ordering and nucleotide seqs being in
# lengths of multiples of three.
# So there's a bunch of housekeeping here in making sure all of this is the case for our final backtrans-align
# step.


def add(env, w, options):
    if options["preserve_indels"]:

        @w.add_target()
        def translated_inseqs_(outdir, c):
            return env.Command(
                [
                    os.path.join(outdir, "translated_inseqs.fa"),
                    os.path.join(outdir, "inseqs_trimmed.fa"),
                ],
                [c["partis_metadata"], c["inseqs"]],
                "translate_seqs.py $SOURCES $TARGET -t ${TARGETS[1]}",
            )

        @w.add_target()
        def translated_inseqs(outdir, c):
            return c["translated_inseqs_"][0]

        @w.add_target()
        def trimmed_inseqs(outdir, c):
            return c["translated_inseqs_"][1]

        @w.add_target()
        def aligned_translated_inseqs(outdir, c):
            return env.SRun(
                os.path.join(outdir, "aligned_translated_inseqs.fa"),
                c["translated_inseqs"],
                # Replace stop codons with X, or muscle inserts gaps, which messed up seqmagick backtrans-align below
                # Note that things will break down at backtrans-align if any seq ids have * in them...
                "sed 's/\*/X/g' $SOURCE | muscle -in /dev/stdin -out $TARGET 2> $TARGET-.log",
                srun_args="`alignment_srun_args.py $SOURCE`",
            )

        # Sort the sequences to have the same order, so that seqmagick doesn't freak out

        @w.add_target()
        def sorted_inseqs(outdir, c):
            return env.Command(
                os.path.join(outdir, "sorted_inseqs.fa"),
                c["trimmed_inseqs"],
                "seqmagick convert --sort name-asc $SOURCE $TARGET",
            )

        @w.add_target()
        def sorted_aligned_translated_inseqs(outdir, c):
            return env.Command(
                os.path.join(outdir, "sorted_aligned_translated_inseqs.fa"),
                c["aligned_translated_inseqs"],
                "seqmagick convert --sort name-asc $SOURCE $TARGET",
            )

    # Now, finally, we can backtranslate

    @w.add_target()
    def aligned_inseqs(outdir, c):
        aligned_inseqs_fname = os.path.join(outdir, "aligned_inseqs.fa")
        if options["preserve_indels"]:
            msg = """
                    Aligning partis cluster sequences because --preserve-indels was passed
                    (either explicitly or by passing --match-indels-in-uid).

                    Check this alignment using information from partis about indels in these sequences
                    before relying on this alignment for inference and moving forward with analysis of these clonal families.
                """
            warnings.warn(partisutils.color("red", msg))
            return env.Command(
                aligned_inseqs_fname,
                [c["sorted_aligned_translated_inseqs"], c["sorted_inseqs"]],
                "seqmagick backtrans-align -a warn $SOURCES -o $TARGET",
            )
        else:
            return env.Command(aligned_inseqs_fname, c["inseqs"], "cp $SOURCE $TARGET")
