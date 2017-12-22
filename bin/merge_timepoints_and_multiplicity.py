#!/usr/bin/env python

import argparse
import csv
import collections
#import itertools


# The core merge/processing logic here

def get_upstream_row(args, seqid):
    row = args.upstream_seqmeta.get(seqid, {'timepoint': '', 'multiplicity': 1})
    row['sequence'] = seqid
    return row

def merge(args):
    """"Initial merge of upstream (pre-partis) metadata, including timepoint and multiplicity info coded in orig
    seqids, with the metadata output of process_partis (partis_seqmeta)."""
    # For each row of our partis sequence metadata (as output from process_partis)
    for seqid, row in args.partis_seqmeta.items():
        # We get the corresponding row of the upstream metadata (which contains our full lenght sequence multiplicities...)
        upstream_row = get_upstream_row(args, seqid)
        duplicates = filter(lambda x: x, row['duplicates'].split(':'))
        seqids = [seqid] + duplicates
        timepoints_dict = collections.defaultdict(lambda: 0)
        for dup_seqid in seqids:
            dup_upstream_row = get_upstream_row(args, seqid)
            # pre-partis filtering multiplicity
            dup_multiplicity = int(dup_upstream_row.get('multiplicity', 1))
            timepoints_dict[dup_upstream_row.get('timepoint')] += dup_multiplicity
        timepoints = sorted(timepoints_dict.items())
        multiplicity = sum(t[1] for t in timepoints)
        result_row = {
                'sequence': seqid,
                'timepoint': upstream_row.get('timepoint'),
                'duplicates': seqids,
                'mut_freq': row['mut_freq'],
                'is_seed': row.get('is_seed'),
                'multiplicity': multiplicity,
                'timepoints': [tp[0] or '' for tp in timepoints],
                'timepoint_multiplicities': [tp[1] for tp in timepoints]}
        yield result_row


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
                timepoint_multiplicities[timepoint] += timepoint_multiplicity
        centroid_data.update({
            'cluster_duplicates': cluster_duplicates,
            'cluster_multiplicity': sum(seq['multiplicity'] for seq in sequences),
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
    parser.add_argument('--partis-seqmeta', type=csv_reader('unique_id'))
    parser.add_argument('--upstream-seqmeta', type=csv_reader('unique_id'))
    parser.add_argument('output', type=argparse.FileType('w'))
    args = parser.parse_args()
    args.upstream_seqmeta = args.upstream_seqmeta or {}
    return args

def main():
    args = get_args()
    fieldnames = ['sequence', 'orig_seqid', 'timepoint', 'mut_freq', 'is_seed', 'multiplicity', 'timepoints',
            'timepoint_multiplicities', 'duplicates']
            #'timepoint_multiplicities']
    if args.cluster_mapping:
        fieldnames += ['cluster_multiplicity', 'cluster_timepoints', 'cluster_timepoint_multiplicities', 'cluster_duplicates']
        #fieldnames += ['cluster_multiplicity', 'cluster_timepoints', 'cluster_timepoint_multiplicities']

    out_writer = csv.DictWriter(args.output, fieldnames=fieldnames)
    out_writer.writeheader()
    results = merge(args)
    if args.cluster_mapping:
        results = aggregate_clusters(results, args.cluster_mapping)
    out_writer.writerows(format_results(results))
    args.output.close()


if __name__ == '__main__':
    main()

