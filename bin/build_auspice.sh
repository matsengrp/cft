#!/usr/bin/env sh

datadir=output/kate-qrs-v10-dnaml/QB850.405-Vk/QB850-k-IgK/run-viterbi-best-plus-0

./bin/build_auspice_jsons.py \
  $datadir/outfile2tree.newick \
  $datadir/outfile2tree.fa \
  $datadir/seqmeta.csv \
  $datadir/tree.json \
  $datadir/sequence.json

cp $datadir/tree.json auspice/data/bcell_tree.json
cp $datadir/sequence.json auspice/data/bcell_sequences.json

