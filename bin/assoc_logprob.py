#!/usr/bin/env python

import argparse
import csv
import json


def get_logprob(partition_handle, cluster_step):
    parts = [x for x in csv.DictReader(partition_handle)]
    for x in parts:
        x['n_clusters'] = int(x['n_clusters'])
        x['logprob'] = float(x['logprob'])
    best_part = min(parts, key=lambda x: x['n_clusters'])
    best_n = best_part['n_clusters']
    cluster_n = best_n + cluster_step
    return [x['logprob'] for x in parts if x['n_clusters'] == cluster_n][0]


def merge_logprob(in_handle, logprob, out_handle):
    doc = json.load(in_handle)
    doc[0]['logprob'] = logprob
    json.dump(doc, out_handle)
    out_handle.close()
    

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata_in', type=argparse.FileType('r'))
    parser.add_argument('partition', type=argparse.FileType('r'))
    parser.add_argument('-n', '--cluster-step', type=int, required=True)
    parser.add_argument('metadata_out', type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    logprob = get_logprob(args.partition, args.cluster_step)
    merge_logprob(args.metadata_in, logprob, args.metadata_out)
    

if __name__ == '__main__':
    main()


