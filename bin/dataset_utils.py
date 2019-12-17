#!/usr/bin/env python

import argparse
import yaml
import json
import uuid
import sys
import functools as fun

def comp(f1, f2):
    def f(*args, **kw_args):
        return f1(f2(*args, **kw_args))
    return f

def merged(d1, d2):
    d = d1.copy()
    d.update(d2)
    return d

def map_vals(f, d):
    return {k: f(v) for k, v in d.items()}

def map_keys(f, d):
    return {f(k): v for k, v in d.items()}

def set_input(input_str):
    return input_str.split(":") if input_str else []

def set_input_of(parse_fn):
    return comp(fun.partial(map, parse_fn), set_input)

def dataset_infile(filename):
    split = filename.split(':')[-1]
    frmt = split[-1] if len(split) > 1 else 'yaml'
    with open(filename) as fh:
        return (json if frmt == 'json' else yaml).load(fh)

def get_args():
    parser = argparse.ArgumentParser(description="Swiss army knife for cft datasets. Filter and merge cft dataset files.")
    parser.add_argument('-i', '--input', type=set_input_of(dataset_infile), help="':'-separated set of yaml or json dataset input files")
    parser.add_argument('-o', '--output', help="dataset output yaml or json file (defaults to stdout as yaml)")
    parser.add_argument('--id', help="output dataset id")
    parser.add_argument('--samples', type=set_input, help="':'-separated set of sample ids to filter to")
    parser.add_argument('--subjects', type=set_input, help="':'-separated set of subject ids to filter to")
    parser.add_argument('--seeds', type=set_input, help="':'-separated set of seed ids to filter to")
    parser.add_argument('--unseeded', action='store_true', help="include unseeded partitions")
    parser.add_argument('--seeded', action='store_true', help="include seeded partitions")
    parser.add_argument('--sample-prefix', help="prefix all sample ids with this value")
    args = parser.parse_args()
    # apply defaults
    args.id = args.id or (args.output.replace('/', '-') if args.output else str(uuid.uuid4()))
    return args

def subset_by_keys(d, keys):
    return {k: v for k, v in d.items() if k in keys}

def apply_seeds_filter(args, sample):
    sample = sample.copy()
    if args.seeds:
        seeds = subset_by_keys(sample['seeds'], args.seeds)
    elif args.seeded or not args.unseeded:
        seeds = sample['seeds']
    else:
        seeds = {}
    sample['seeds'] = seeds
    return sample

def apply_unseeded_filter(args, sample):
    if not args.unseeded and (args.seeded or args.seeds):
        return subset_by_keys(sample, set(sample.keys()).difference(['partition-file', 'cluster-annotation-file', 'other-partitions']))
    else:
        return sample

def remove_empty_samples(args, samples):
    return {sample_id: sample
            for sample_id, sample in samples.items()
            if sample.get('seeds') or (sample.get('partition-file') and sample.get('cluster-annotation-file')) or sample.get('other-partitions')}

def prefix_sample_ids(args, samples):
    return map_keys(lambda k: args.sample_prefix + '-' + k, samples) if args.sample_prefix else samples

def process_samples(args, samples):
    if args.samples:
        samples = subset_by_keys(samples, args.samples)
    samples = map_vals(fun.partial(apply_seeds_filter, args), samples)
    samples = map_vals(fun.partial(apply_unseeded_filter, args), samples)
    samples = remove_empty_samples(args, samples)
    samples = prefix_sample_ids(args, samples)
    return samples

def process_dataset(args, dataset):
    return merged(dataset, {'samples': process_samples(args, dataset['samples'])})

def merge_datasets(args, datasets):
    dataset = reduce(merged, datasets)
    samples = reduce(merged, [d['samples'] for d in datasets])
    dataset['samples'] = samples
    dataset['subjects'] = reduce(merged, [d.get('subjects', {}) for d in datasets])
    dataset['id'] = args.id
    return dataset

def main():
    args = get_args()
    datasets = map(fun.partial(process_dataset, args), args.input)
    dataset = merge_datasets(args, datasets)
    out_frmt = args.output.split('.')[-1] if args.output else 'yaml'
    with open(args.output, 'w') if args.output else sys.stdout as fh:
        (json if out_frmt == 'json' else yaml).dump(dataset, fh)

if __name__ == '__main__':
    main()
