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

datapath="${datapath:-/fh/fast/matsen_e/processed-data/partis/kate-qrs-2016-09-09/new}"
outdir="${outdir:-output}"

# Run processing on annotations file
# BUG: we should scan for these seeds or they should be provided on
# the command line.
#seed1="$datapath/seeds/QA255.006-Vh"
seed2="$datapath/seeds/QB850.001-Vh"
seed3="$datapath/seeds/QB850.043-Vk"
#seeds=( $seed1 $seed2 $seed3)
seeds=( $seed2 $seed3)

# Uncomment below to run on all seeds.
# Will pass as argument eventually, but maybe this is a job for
# python/scons?
#seeds=($datapath/seeds/Q*)

echo "${seeds[@]}"

# Separate fasta files for each cluster
for seeddir in "${seeds[@]}"
do
    seed="$(basename $seeddir)"
    timedirs=($datapath/seeds/$seed/*)
    for annotation in  "${timedirs[@]}"
    do
	path="${annotation%/*}"
	basename="${annotation##*/}"
	timept="${basename%-100k}"
    paramdir="$datapath/$timept"
	
	partition=${path}/${basename}/partition.csv
	annotation=${path}/${basename}/partition-cluster-annotations.csv
	logfile=${path}/${basename}.log

    # we can run this with either --partis_log ${logfile} or
    # --param_dir ${paramdir}
    # the latter is probably preferred so log files do not have to
    # be processed...

	process_partis.py \
	       --annotations ${annotation} \
	       --partition ${partition} \
	       --param_dir ${paramdir} \
	       --cluster_base cluster \
	       --output_dir ${outdir}/$seed/$basename
    done
done

