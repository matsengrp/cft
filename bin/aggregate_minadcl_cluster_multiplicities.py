#!/usr/bin/env python

import argparse
import csv
import collections
#import itertools



def aggregate_clusters(merge_results, cluster_mapping):
    merge_results = {row['sequence']: row for row in merge_results}
    for centroid_id, sequence_ids in cluster_mapping.items():
        centroid_data = merge_results[centroid_id]
        sequences = [merge_results[sequence_id] for sequence_id in sequence_ids]
        cluster_duplicates = reduce(set.union, (seq['duplicates'] for seq in sequences), set())
        timepoint_multiplicities = collections.defaultdict(lambda: 0)
        # Aggregate our timepoint multiplicities
        for seq in sequences:
            for timepoint, timepoint_multiplicity in zip(seq['timepoints'], seq['timepoint_multiplicities']):
                timepoint_multiplicities[timepoint] += int(timepoint_multiplicity)
        centroid_data.update({
            'cluster_duplicates': cluster_duplicates,
            'cluster_multiplicity': sum(int(seq['multiplicity']) for seq in sequences),
            'cluster_timepoints': timepoint_multiplicities.keys(),
            'cluster_timepoint_multiplicities': timepoint_multiplicities.values(),
            })
        yield centroid_data


def format_list(row, key):
    if key in row:
        row[key] = ':'.join(map(str, row[key]))


def format_results(results):
    for row in results:
        format_list(row, 'timepoints')
        format_list(row, 'timepoint_multiplicities')
        format_list(row, 'duplicates')
        format_list(row, 'cluster_duplicates')
        format_list(row, 'cluster_timepoints')
        format_list(row, 'cluster_timepoint_multiplicities')
        yield row


# Wrapping together the CLI

def csv_reader(index=None, filter_by=None):
    def f(filename):
        with open(filename) as fh:
            reader = csv.DictReader(fh)
            if filter_by:
                reader = filter(filter_by, reader)
            reader = list(reader)
            if index:
                return {r[index]: r for r in reader}
            else:
                return list(reader)
    return f

def cluster_reader(filename):
    data = csv_reader()(filename)
    #return {k: [v['sequence'] for v in vs] for k, vs in itertools.groupby(data, lambda d: d['centroid'])}
    result = collections.defaultdict(set)
    for d in data:
        result[d['centroid']].add(d['sequence'])
    return result


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster-mapping', type=cluster_reader)
    parser.add_argument('--partis-seqmeta', type=csv_reader())
    parser.add_argument('output', type=argparse.FileType('w'))
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    fieldnames = ['sequence', 'unique_id', 'orig_seqid', 'timepoint', 'mut_freq', 'is_seed', 'stops',
            'mutated_invariants', 'frameshifted', 'multiplicity', 'timepoints',
            'timepoint_multiplicities', 'duplicates', ]
    if args.cluster_mapping:
        fieldnames += ['cluster_multiplicity', 'cluster_timepoints', 'cluster_timepoint_multiplicities', 'cluster_duplicates']

    out_writer = csv.DictWriter(args.output, fieldnames=fieldnames, extrasaction='ignore')
    out_writer.writeheader()
    results = args.partis_seqmeta
    if args.cluster_mapping:
        results = aggregate_clusters(results, args.cluster_mapping)
    out_writer.writerows(format_results(results))
    args.output.close()


if __name__ == '__main__':
    main()

