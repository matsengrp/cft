#!/usr/bin/bash

# Run processing on annotations file
python scripts/process_partis.py \
    --incsv data/seeds/QA255.006-Vh/Hs-LN2-5RACE-IgG-new-cluster-annotations.csv \
    --cluster_base cluster \
    --input_dir data/seeds/ \
    --output_dir output/seeds/

## BASELINe-style output
#python scripts/process_partis.py \
#    --incsv data/seeds/QA255.006-Vh/Hs-LN2-5RACE-IgG-new-cluster-annotations.csv \
#    --cluster_base cluster \
#    --input_dir data/seeds/ \
#    --output_dir output/seeds/ \
#    --baseline
