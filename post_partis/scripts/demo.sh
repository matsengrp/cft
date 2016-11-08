#!/usr/bin/bash

set -eux

# Run processing on annotations file
seed1="QA255.006-Vh"
seed2="QB850.424-Vk"
seed3="QB850.043-Vk"

datapath="/fh/fast/matsen_e/dshaw/cft/data/seeds"

# Separate fasta files for each cluster
python scripts/process_partis.py \
    --incsv $datapath/$seed2/Hs-LN1-5RACE-IgK-100k-cluster-annotations.csv \
    --cluster_base cluster \
    --input_dir $datapath \
    --output_dir _output/seeds \
    --separate

# Aggregated fasta file with all seqs
python scripts/process_partis.py \
    --incsv $datapath/$seed1/Hs-LN2-5RACE-IgG-new-cluster-annotations.csv \
    --cluster_base cluster \
    --input_dir $datapath \
    --output_dir _output/seeds

