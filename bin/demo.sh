#!/usr/bin/env bash

set -eux
shopt -s nullglob

# add the current directory to the path.
# we expect to find process_partis.py
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null ; pwd -P )"
PATH="${SCRIPTPATH}:${PATH}"
which process_partis.py

# path where we expect to find the seed directories
# containing processed partis outpout.
# note you can override this with an environment variable,
#    datapath="/my/home" demo.sh

datapath="${datapath:-/fh/fast/matsen_e/processed-data/partis/kate-qrs-2016-09-09/seeds}"
outdir="${outdir:-output}"


# Run processing on annotations file
# BUG: we should scan for these seeds or they should be provided on
# the command line.
#seed1="$datapath/QA255.006-Vh"
seed2="$datapath/QB850.424-Vk"
seed3="$datapath/QB850.043-Vk"
#seeds=( $seed1 $seed2 $seed3)
seeds=( $seed2 $seed3)

# Uncomment below to run on all seeds.
# Will pass as argument eventually, but maybe this is a job for
# python/scons?
seeds=($datapath/Q*)

# Separate fasta files for each cluster
for seeddir in "${seeds[@]}"
do
  seed="$(basename $seeddir)"
  files=($datapath/$seed/*-cluster-annotations.csv)
  for annotation in  "${files[@]}"
  do
    path="${annotation%/*}"
    basename="${annotation##*/}"
    basename="${basename%-cluster-annotations.csv}"

    partition=${path}/${basename}.csv
    logfile=${path}/${basename}.log
    process_partis.py \
           --annotations ${annotation} \
           --partition ${partition} \
           --partis_log ${logfile} \
           --cluster_base cluster \
           --output_dir ${outdir}/$seed/$basename \
           --separate
  done
done

