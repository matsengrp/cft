#!/usr/bin/env python

import argparse
import collections
from os import path
import csv
import json
import re
import pprint
import yaml


def seqmeta_reader(filename):
    with open(filename) as fh:
        reader = csv.DictReader(fh)
        return {row["sequence"]: row for row in reader}


def yaml_reader(filename):
    with open(filename) as fh:
        try:
            result = yaml.load(fh)
        except yaml.YAMLError, e:
            raise
        return result if result else {}


def merge(d1, d2):
    d = d1.copy()
    d.update(d2)
    return d


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("tip_seqmeta", type=seqmeta_reader)
    parser.add_argument("sel_metric_info", type=yaml_reader)
    parser.add_argument("seqmeta_out", type=argparse.FileType("w"))
    parser.add_argument("metrics", nargs="?", default=['lbi', 'lbr', 'cons-dist-aa'])
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    if len(args.sel_metric_info) == 0:
        print 'merge_selection_metrics.py: no clusters in input file'
    elif len(args.sel_metric_info) > 1:
        print 'merge_selection_metrics.py warning: info for multiple clusters in input file, arbitrarily taking the first one'
    sinfo = args.sel_metric_info[0].get("lb", {})
    if len(set(args.metrics) - set(sinfo)) > 0:
        print 'merge_selection_metrics.py: requested metrics %s not found in input file' % ' '.join(set(args.metrics) - set(sinfo))
        args.metrics = [m for m in args.metrics if m in sinfo]

    out_fields = args.tip_seqmeta.values()[0].keys() + args.metrics
    writer = csv.DictWriter(
        args.seqmeta_out, fieldnames=out_fields, extrasaction="ignore"
    )
    writer.writeheader()

    seqmeta = collections.defaultdict(lambda: {"tip": False, "internal": True})
    seqmeta.update(
        {
            k: merge(v, {"tip": True, "internal": False})
            for k, v in args.tip_seqmeta.items()
        }
    )
    for metric in args.metrics:
        for seqid, val in sinfo[metric].items():
            seqmeta[seqid].update(
                {"sequence": seqid, "unique_id": seqid, metric: val}
            )
    for _, row in seqmeta.items():
        writer.writerow(row)


if __name__ == "__main__":
    main()
