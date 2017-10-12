#!/usr/bin/env python

import argparse
import csv
import re
import collections
#import itertools


# The core merge/processing logic here

duplicate_seqid_regex = re.compile('\d+-(\d+)')

def upstream_seqmeta(args, seqid):
    default = {'original': seqid, 'timepoint': args.timepoint} if args.timepoint else {}
    return args.upstream_seqmeta.get(seqid, default)
    

def merge(args):
    """"Initial merge of upstream (pre-partis) metadata, including timepoint and multiplicity info coded in orig
    seqids, with the metadata output of process_partis (partis_seqmeta)."""
    # For each row of our partis sequence metadata (as output from process_partis)
    for seqid, row in args.partis_seqmeta.items():
        # We get the corresponding row of the upstream metadata (which contains our full lenght sequence duplicities...)
        upstream_row = upstream_seqmeta(args, seqid)
        duplicates = filter(lambda x: x, row['duplicates'].split(':'))
        seqids = [seqid] + duplicates
        timepoints_dict = collections.defaultdict(lambda: 0)
        for dup_seqid in seqids:
            dup_upstream_row = upstream_seqmeta(args, seqid)
            # The original seqid and orig multiplicity are the duplicites from vlad's pre-partis processing
            orig_seqid = dup_upstream_row.get('original')
            orig_seqid_match = duplicate_seqid_regex.match(orig_seqid) if isinstance(orig_seqid, str) else False
            if orig_seqid_match:
                orig_multiplicity = int(orig_seqid_match.groups()[0])
                #print "  multiplicity ", orig_multiplicity
            else:
                # TODO For seeds this is good but not for naive...
                #print "no multiplicity for", orig_seqid
                orig_multiplicity = 1
            timepoints_dict[dup_upstream_row.get('timepoint')] += orig_multiplicity
        timepoints = sorted(timepoints_dict.items())
        multiplicity = sum(t[1] for t in timepoints)
        result_row = {
                'sequence': seqid,
                'timepoint': upstream_row.get('timepoint'),
                'duplicates': seqids,
                'mut_freqs': row['mut_freqs'],
                'is_seed': not duplicate_seqid_regex.match(seqid),
                'orig_seqid': upstream_row.get('original'),
                'multiplicity': multiplicity,
                'timepoints': [tp[0] or '' for tp in timepoints],
                'timepoint_duplicities': [tp[1] for tp in timepoints]
                }
        yield result_row


def aggregate_clusters(merge_results, cluster_mapping):
    merge_results = {row['sequence']: row for row in merge_results}
    for centroid_id, sequence_ids in cluster_mapping.items():
        centroid_data = merge_results[centroid_id]
        sequences = [merge_results[sequence_id] for sequence_id in sequence_ids]
        cluster_duplicates = reduce(set.union, (seq['duplicates'] for seq in sequences), set())
        timepoint_duplicities = collections.defaultdict(lambda: 0)
        # Aggregate our timepoint duplicities
        for seq in sequences:
            for timepoint, timepoint_multiplicity in zip(seq['timepoints'], seq['timepoint_duplicities']):
                timepoint_duplicities[timepoint] += timepoint_multiplicity
        centroid_data.update({
            'cluster_duplicates': cluster_duplicates,
            'cluster_multiplicity': sum(seq['multiplicity'] for seq in sequences),
            'cluster_timepoints': timepoint_duplicities.keys(),
            'cluster_timepoint_duplicities': timepoint_duplicities.values(),
            })
        yield centroid_data


def format_list(row, key):
    if key in row:
        row[key] = ':'.join(map(str, row[key]))


def format_results(results):
    for row in results:
        format_list(row, 'timepoints')
        format_list(row, 'timepoint_duplicities')
        format_list(row, 'duplicates')
        format_list(row, 'cluster_duplicates')
        format_list(row, 'cluster_timepoints')
        format_list(row, 'cluster_timepoint_duplicities')
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


def timepoint_filter(args):
    mappings = []
    def f(x):
        mappings.append(x)
        return x['timepoint'] == args.timepoint
    if args.timepoint:
        return f


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--timepoint', help="optionally specify a specific timepoint on which to restrict")
    parser.add_argument('--cluster-mapping', type=cluster_reader)
    parser.add_argument('partis_seqmeta', type=csv_reader('unique_ids'))
    parser.add_argument('upstream_seqmeta')
    parser.add_argument('output', type=argparse.FileType('w'))
    args = parser.parse_args()
    index_by = 'original' if args.timepoint else 'new'
    args.upstream_seqmeta = \
            csv_reader(index=index_by,
                       filter_by=timepoint_filter(args))(
                               args.upstream_seqmeta
                               )
    return args

def main():
    args = get_args()
    fieldnames = ['sequence', 'orig_seqid', 'timepoint', 'mut_freqs', 'is_seed', 'multiplicity', 'timepoints',
            'timepoint_duplicities', 'duplicates']
            #'timepoint_duplicities']
    if args.cluster_mapping:
        fieldnames += ['cluster_multiplicity', 'cluster_timepoints', 'cluster_timepoint_duplicities', 'cluster_duplicates']
        #fieldnames += ['cluster_multiplicity', 'cluster_timepoints', 'cluster_timepoint_duplicities']

    out_writer = csv.DictWriter(args.output, fieldnames=fieldnames)
    out_writer.writeheader()
    results = merge(args)
    if args.cluster_mapping:
        results = aggregate_clusters(results, args.cluster_mapping)
    out_writer.writerows(format_results(results))
    args.output.close()


if __name__ == '__main__':
    main()

