#!/usr/bin/env python

import argparse
import json
import copy
import os
import shutil
from os import path


path_attrs_to_copy = ['cluster_aa', 'fasta', 'newick', 'seedlineage', 'seedlineage_aa', 'svg']


def copy_file(inbase, outbase, filepath):
    newpath = path.join(outbase, path.dirname(filepath))
    try:
        os.makedirs(newpath)
    except OSError:
        pass
    shutil.copy(path.join(inbase, filepath), newpath)

def copy_cluster_data(inbase, outbase, cluster_meta):
    for attr in path_attrs_to_copy:
        copy_file(inbase, outbase, cluster_meta[attr])

def copy_data(inbase, outbase, metadata):
    for cluster in metadata['clusters']:
        copy_cluster_data(inbase, outbase, cluster)

def with_dataset_id(metadata, dataset_id):
    metadata = copy.deepcopy(metadata)
    metadata['dataset_id'] = dataset_id
    for cluster in metadata['clusters']:
        cluster['dataset_id'] = dataset_id
    return metadata


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata')
    parser.add_argument('outdir')
    parser.add_argument('-d', '--dataset-id', help='specify a new dataset id for the data')
    return parser.parse_args()

def main():
    args = get_args()
    inbase = path.dirname(args.metadata)
    with open(args.metadata, 'r') as meta_handle:
        metadata = json.load(meta_handle)
    copy_data(inbase, args.outdir, metadata)
    if args.dataset_id:
        metadata = with_dataset_id(metadata, args.dataset_id)
        with open(path.join(args.outdir, 'metadata.json'), 'w') as fh:
            json.dump(metadata, fh)
    else:
        shutil.copy(args.metadata, args.outdir)


if __name__ == '__main__':
    main()

