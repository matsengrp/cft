#!/usr/bin/env python

import argparse
import json
import csv
from tripl import tripl


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', nargs='+')
    parser.add_argument('-c', '--csv', action="store_true")
    parser.add_argument('-o', '--output')
    return parser.parse_args()

def strip_ns(a):
    return a.split(':')[-1]

def clean_record(d):
    try:
        del d['cft.cluster:unique_ids']
    except Exception:
        pass
    d = {strip_ns(k): clean_record(v) if isinstance(v, dict) else v
            for k, v in d.items()}
    return d

pull_pattern = [
    "cft.cluster:id",
    "cft.cluster:naive_seq",
    "cft.cluster:has_seed",
    "cft.cluster:n_seqs",
    "cft.cluster:n_sampled_seqs",
    "cft.cluster:size",
    "cft.cluster:sorted_index",
    "cft.cluster:v_end",
    "cft.cluster:v_start",
    "cft.cluster:v_gene",
    "cft.cluster:d_end",
    "cft.cluster:d_gene",
    "cft.cluster:d_start",
    "cft.cluster:j_end",
    "cft.cluster:j_gene",
    "cft.cluster:j_start",
    "cft.cluster:cdr3_length",
    "cft.cluster:cdr3_start",
    "cft.cluster:naive_seq",
    "cft.cluster:mean_mut_freq",
    "db:ident",
    "tripl:type",
    {"cft.cluster:sample": ["cft.sample:id", "cft.sample:timepoint"],
     "cft.cluster:dataset": ["cft.dataset:id"],
     "cft.cluster:partition": ["cft.partition:id", "cft.partition:logprob", "cft.partition:step"],
     "cft.cluster:subject": ["cft.subject:id"],
     "cft.cluster:v_per_gene_support": ["cft.gene_support:gene", "cft.gene_support:prob"],
     "cft.cluster:d_per_gene_support": ["cft.gene_support:gene", "cft.gene_support:prob"],
     "cft.cluster:j_per_gene_support": ["cft.gene_support:gene", "cft.gene_support:prob"]}]

def get_data(t):
    return map(clean_record, t.pull_many(pull_pattern, {'tripl:type': 'cft.cluster'}))

def main():
    args = get_args()
    t = tripl.TripleStore.loads(args.inputs)
    data = get_data(t)
    if args.csv:
        data = [{k: v for k, v in d.items()}
                for d in data]
        with open(args.output, 'w') as fh:
            writer = csv.DictWriter(fh, fieldnames=sorted(data[0].keys()))
            writer.writeheader()
            writer.writerows(data)
    else:
        # Then assume json
        with open(args.output, 'w') as fh:
            json.dump(data, fh, default=list,
                    indent=4)


if __name__ == '__main__':
    main()


