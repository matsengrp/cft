#!/usr/bin/bash

# Run processing on annotations file

# Separate fasta files for each cluster
python scripts/process_partis.py \
    --incsv data/seeds/QA255.006-Vh/Hs-LN2-5RACE-IgG-new-cluster-annotations.csv \
    --cluster_base cluster \
    --input_dir data/seeds/ \
    --output_dir _output/seeds/ \
    --separate

# Aggregated fasta file with all seqs
python scripts/process_partis.py \
    --incsv data/seeds/QB850.424-Vk/Hs-LN1-5RACE-IgK-100k-cluster-annotations.csv \
    --cluster_base cluster \
    --input_dir data/seeds/ \
    --output_dir _output/seeds/

# BASELINe-style output
python scripts/process_partis.py \
    --incsv data/seeds/QB850.043-Vk/Hs-LN1-5RACE-IgK-100k-cluster-annotations.csv \
    --cluster_base cluster \
    --input_dir data/seeds/ \
    --output_dir _output/seeds/ \
    --baseline
