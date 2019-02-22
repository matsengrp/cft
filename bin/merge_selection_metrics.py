#!/usr/bin/env python

import argparse
import collections
from os import path
import csv
import json
import re
import pprint
import process_partis
process_partis.insert_partis_path()
import utils

# File reading stuff
# ------------------

def seqmeta_reader(filename):
    with open(filename) as fh:
        reader = csv.DictReader(fh)
        return {row['sequence']: row for row in reader}

def seq_names_reader(filename):
    seqs_dict = utils.read_fastx(filename)
    return map(lambda x: x['name'], seqs_dict)

def find_partis_cluster_annotation(partis_output, seqmeta, seq_names):
    """
    Find the annotation that has the most overlapping sequence ids with the input sequences to "partis annotate --get-tree-metrics"
    This should be the same as the input sequences since we are passing it a cluster given by partis, but we can't assume this will 
    be the case. See https://github.com/matsengrp/cft/issues/264
    """
    seq_names = seqmeta.keys()
    _, annotation_list, _ = utils.read_output(partis_output)

    # this code comes mostly from https://github.com/psathyrella/partis/blob/dev/python/treeutils.py#L759
    max_in_common, annotation_idx_to_use = None, None
    # annotation_list should have length 1 here for the most part
    for i, annotation in enumerate(annotation_list):
        n_in_common = len(set(annotation['unique_ids']) & set(seq_names)) # take intersection of the two sets
        if max_in_common is None or n_in_common > max_in_common:
            annotation_idx_to_use = i
            max_in_common = n_in_common
        if max_in_common is None:
            raise Exception('cluster annotation with any of the following sequence ids: \'%s\' not found in annotations. Check what is being passed to \'partis annotate --get-tree-metrics\'' % seq_names)
        if max_in_common < len(seq_names):
            print( '    note: couldn\'t find annotated cluster that shared all sequences with input seqs (best was %d/%d)' % (max_in_common, len(seq_names)))
    
    return annotation_list[annotation_idx_to_use]


# Args
# ----

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tip_seqmeta', type=seqmeta_reader)
    parser.add_argument('selection_metrics_fname', type=str)
    parser.add_argument('asr_seqs', type=seq_names_reader)
    parser.add_argument('seqmeta_out', type=argparse.FileType('w'))
    args = parser.parse_args()
    return args


# Putting everything together
# ---------------------------


def merge(d1, d2):
    d = d1.copy()
    d.update(d2)
    return d


def main():
    args = get_args()
    out_fields = args.tip_seqmeta.values()[0].keys() + ['lbi', 'lbr']
    writer = csv.DictWriter(args.seqmeta_out, fieldnames=out_fields, extrasaction='ignore')
    writer.writeheader()
    
    
    seqmeta = collections.defaultdict(lambda: {'tip': False, 'internal': True})
    seqmeta.update({k: merge(v, {'tip': True, 'internal': False}) for k, v in args.tip_seqmeta.items()})
    
    cluster_annotation = find_partis_cluster_annotation(args.selection_metrics_fname, seqmeta, args.asr_seqs)
    selection_metrics = cluster_annotation.get('tree-info', {}).get('lb', {})  
    
    for metric in ['lbi', 'lbr']:
        if selection_metrics.get(metric):
            for seqid, val in selection_metrics[metric].items():
                seqmeta[seqid].update({'sequence': seqid, 'unique_id': seqid, metric: val})
    for _, row in seqmeta.items():
        writer.writerow(row)

    # for _, row in args.tip_seqmeta
    #     row = args.tip_seqmeta.get(seqid, {'sequence': seqid, 'unique_id': seqid})
    #     lb = args.selection_metrics
    #     #pprint.pprint(args.selection_metrics['lbi'])
    #     row['lbi'] = args.selection_metrics['lbi'][i]
    #     row['lbr'] = args.selection_metrics['lbr'][i]
    #     writer.writerow(row)

if __name__ == '__main__':
    main()

