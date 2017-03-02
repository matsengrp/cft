#!/usr/bin/env python

import argparse
import csv
import json



def get_logprob(partition_handle, cluster_step):
    parts = [x for x in csv.DictReader(partition_handle)]
    for i, x in enumerate(parts):
        x['index'] = i
        x['n_clusters'] = int(x['n_clusters'])
        x['logprob'] = float(x['logprob'])
    best_part = max(parts, key=lambda x: x['logprob'])
    best_i = best_part['index']
    cluster_i = best_i + cluster_step
    return parts[cluster_i]["logprob"]


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


