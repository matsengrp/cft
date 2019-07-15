#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import os.path
import csv
import sys
import json
import textwrap
import process_partis


# Figure out where partis is so that partis utils and glutils ccan be loaded below
partis_path = os.environ.get('PARTIS')
if not partis_path or not os.path.exists(partis_path):
    msg = """\
    Must set environment variable PARTIS to point to an `partis` installation.
    You can clone the partis repo with,
        git clone --depth 1 git@github.com:psathyrella/partis.git
        export  PARTIS=$PWD/partis
    """
    print(textwrap.dedent(msg))
    sys.exit(1)

sys.path.insert(1, os.path.join(partis_path, 'python'))
import utils

def iseqs_from_uids(uids_file, cluster_annotation):
    with open(args.subset_ids_path) as f:
        ids = f.readlines()
    ids = [id.rstrip("\n") for id in ids]
    return [iseq for iseq, uid in enumerate(cluster_annotation['unique_ids']) if uid in ids]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rewrite partis output file with subset of ids.")
    parser.add_argument(
        'partition_file', type=str,
        help="Path to partis partition file")
    parser.add_argument(
        'subset_ids_path', type=str,
        help="Path to text file containing ids we want to use for the subset.")
    parser.add_argument(
        'outfname', type=str,
        help="The name of partis yaml file to write.")
    parser.add_argument(
        '--partition-step', type=str, required=True,
        help="The partition step to use.")
    parser.add_argument(
        '--original-cluster-unique-ids', type=str, required=True,
        help="Colon separated list of unique ids to identify the cluster we want to subset")
    parser.add_argument(
        '--glfo-dir',
        help='path to germline info, only necessary for deprecated .csv output files')
    parser.add_argument(
        '--locus',
        help='Sample locus, only necessary for deprecated .csv output files')

    args = parser.parse_args()
    
    #print('CLUSTER IDS {}'.format(args.original_cluster_unique_ids))
    glfo, annotation_list, cpath = process_partis.read_partis_output(args.partition_file, args.glfo_dir, args.locus)
    cluster_annotation = process_partis.choose_cluster(args.partition_file, annotation_list, cpath, args.partition_step, i_cluster=None, unique_ids=args.original_cluster_unique_ids.split(':'))
    iseqs = iseqs_from_uids(args.subset_ids_path, cluster_annotation)
    # TODO: restrict to iseqs should rewrite linearham info among other things according to the new subset
    new_annotation_list = [utils.restrict_to_iseqs(cluster_annotation, iseqs, glfo)]

    utils.write_annotations(args.outfname, glfo, new_annotation_list, utils.sw_cache_headers)
