#!/usr/bin/env python

# WARNING: This file appears to be deprecated!

import collections
import argparse
from os import path
import csv
import re


# Processing code
# ---------------

# The rows are from a "translations" csv file with columns `sample,new,original,timepoint`

# String ending in digit
multiplicity_match = re.compile('.*-(\d*)\Z').match

def multiplicity(row):
    m = multiplicity_match(row['original'])
    return int(m.group(1)) if m else 1

def seqmeta_row(trans_row, merged=False):
    return {'sequence': trans_row['new' if merged else 'original'],
            'timepoint': trans_row['timepoint'],
            'multiplicity': multiplicity(trans_row)}


# File reading stuff
# ------------------

def translation_reader(filename):
    with open(filename) as fh:
        reader = csv.DictReader(fh)
        by_sample = collections.defaultdict(list)
        for row in reader:
            by_sample[row['sample']].append(row)
        return by_sample

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('translation_file', type=translation_reader)
    parser.add_argument('--outdir')
    args = parser.parse_args()
    return args

def write_out(filename, rows, merged=False):
    with open(filename, 'w') as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=['sequence', 'timepoint', 'multiplicity'])
        writer.writeheader()
        writer.writerows(seqmeta_row(row, merged=merged) for row in rows)


# Putting everything together
# ---------------------------

def main():
    args = get_args()
    for sample, sample_rows in args.translation_file.items():
        write_out(path.join(args.outdir, sample + '-seqmeta.csv'), sample_rows, merged=False)
    merged = (x for row in args.translation_file.values() for x in row)
    write_out(path.join(args.outdir, 'merged-seqmeta.csv'), merged, merged=True)


if __name__ == '__main__':
    main()

