#!/usr/bin/env python

import argparse
import csv
import re
import collections


# The core merge/processing logic here

duplicate_seqid_regex = re.compile('\d+-(\d+)')

def merge(args):
    for seqid, row in args.partis_seqmeta.items():
        upstream_row = args.upstream_seqmeta.get(seqid, {})
        duplicates = filter(lambda x: x, row['duplicates'].split(':'))
        seqids = [seqid] + duplicates
        result_row = {'sequence': seqid, 'timepoint': upstream_row.get('timepoint'), 'duplicates':
                row['duplicates'], 'mut_freqs': row['mut_freqs']}
        result_row['is_seed'] = not duplicate_seqid_regex.match(seqid)
        result_row['orig_seqid'] = upstream_row.get('original')
        timepoints = collections.defaultdict(lambda: 0)
        for dup_seqid in seqids:
            upstream_row = args.upstream_seqmeta.get(dup_seqid, {})
            orig_seqid = upstream_row.get('original')
            orig_seqid_match = duplicate_seqid_regex.match(orig_seqid) if isinstance(orig_seqid, str) else False
            if orig_seqid_match:
                orig_duplicity = int(orig_seqid_match.groups()[0])
            else:
                # TODO For seeds this is good but not for naive...
                orig_duplicity = 1
            timepoints[upstream_row.get('timepoint')] += orig_duplicity
        result_row['duplicity'] = sum(timepoints.values())
        timepoints = sorted(timepoints.items())
        result_row['timepoints'] = ':'.join(tp[0] or '' for tp in timepoints)
        result_row['timepoint_duplicities'] = ':'.join(str(tp[1]) for tp in timepoints)
        yield result_row


# Wrapping together the CLI

def csv_reader(index=None):
    def f(filename):
        with open(filename) as fh:
            if index:
                return {r[index]: r for r in csv.DictReader(fh)}
            else:
                return list(csv.DictReaer(fh))
    return f


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('partis_seqmeta', type=csv_reader('unique_ids')),
    parser.add_argument('upstream_seqmeta', type=csv_reader('new'))
    parser.add_argument('output', type=argparse.FileType('w'))
    return parser.parse_args()

def main():
    args = get_args()
    out_writer = csv.DictWriter(args.output,
            fieldnames=['sequence', 'orig_seqid', 'timepoint', 'mut_freqs', 'duplicity', 'is_seed', 'orig_seqid', 'duplicates', 'timepoints', 'timepoint_duplicities'])
    out_writer.writeheader()
    out_writer.writerows(merge(args))
    args.output.close()


if __name__ == '__main__':
    main()

