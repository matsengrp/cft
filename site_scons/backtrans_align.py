#!/usr/bin/env python
# -*- coding: utf-8 -*-


from os import path


# What follows is a rather convoluted backtranslation strategy for sequence alignment.
# The basic idea is that we translate our sequences, align them (taking advantage of coding information),
# and from this infer the a nucleotide alignment using `seqmagick backtrans-align`.
# The problem is that `seqmagick` rather finicky here about sequence ordering and nucleotide seqs being in
# lengths of multiples of three.
# So there's a bunch of housekeeping here in making sure all of this is the case for our final backtrans-align
# step.



def add(env, w):
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
            'sed \'s/\*/X/g\' $SOURCE | muscle -in /dev/stdin -out $TARGET 2> $TARGET-.log',
            srun_args='`alignment_srun_args.py $SOURCE`')

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

