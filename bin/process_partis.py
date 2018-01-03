#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import os.path
import csv
import json
import sys
import textwrap
import time
#import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Figure out where partis is so that partis utils and glutils ccan be loaded below
partis_path = os.environ.get('PARTIS')
if not partis_path or not os.path.exists(partis_path):
    msg = """\
    Must set environment variable PARTIS to point to an `partis` installation.
    You can clone the partis repo with,
        git clone --depth 1 git@github.com:psathyrella/partis.git
        export  PARTIS=$PWD/partis
    """
    print(textwrap.dedent(msg))
    sys.exit(1)

sys.path.insert(1, os.path.join(partis_path, 'python'))
import utils
import glutils
import clusterpath


# Make sure we can read the really big fields frequently found in partis output
csv.field_size_limit(sys.maxsize)


default_germline_sets = os.path.join(partis_path, '/data/germlines/human')



def indel_offset(indelfo):
    return sum(map(lambda indel: indel['len'] * (1 if indel['type'] == 'insertion' else -1),
                   indelfo['indels']))

def infer_frameshifts(line):
    "Infer frameshifts based on `in_frame` attr in partis output, or if more than one stop codon."
    attrs = ['in_frames', 'input_seqs']
    def infer_(args):
        in_frame, input_seq = args
        # Trim to multiple of 3 to avoid biopython warning
        input_seq = input_seq[0:((len(input_seq) / 3) * 3)]
        aa = Seq(input_seq).translate()
        stop_count = aa.count("*")
        return (not in_frame) or stop_count > 0
    return map(infer_, zip(*map(lambda x: line[x], attrs)))

def get_cluster_annotation(filename, unique_ids):
    unique_ids_string = ":".join(unique_ids)
    with open(filename) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row['unique_ids'] == unique_ids_string:
                return row
        raise ValueError("Cluster not found for unique_ids set {}".format(repr(unique_ids)))

def as_dict_rows(column_dict, columns=None):
    columns = columns or column_dict.keys()
    column_lengths = [len(column_dict[c]) for c in columns]
    assert min(column_lengths) == max(column_lengths), "columns can't be of different lengths"
    return [{c: column_dict[c][i] for c in columns}
            for i in range(max(column_lengths))]

def apply_filters(args, sequences):
    def filter_fn(sequence):
        if args.remove_frameshifts and sequence['frameshifted']:
            return False
        elif args.remove_stops and sequence['stops']:
            return False
        elif args.remove_mutated_invariants and sequence['mutated_invariants']:
            return False
        else:
            return True
    return filter(filter_fn, sequences)


def seqs(args, cluster_line):
    non_reversed_seqs = cluster_line['seqs']
    if args.indel_reversed_seqs:
        return [reversed_seq or non_reversed_seqs[i]
                for i, reversed_seq in enumerate(cluster_line['indel_reversed_seqs'])]
    else:
        return non_reversed_seqs


def process_cluster(args, cluster_line, seed_id):
    utils.process_input_line(cluster_line)
    utils.add_implicit_info(args.glfo, cluster_line)

    cluster_sequences = {
            'unique_id': cluster_line['unique_ids'] + ['naive'],
            'seq': seqs(args, cluster_line) + [cluster_line['naive_seq']],
            'is_seed': [(unique_id == seed_id) for unique_id in cluster_line['unique_ids']] + [False],
            'duplicates': [':'.join(x) for x in cluster_line['duplicates']] + [None],
            'frameshifted': [not x for x in cluster_line['in_frames']] + [False],
            'mutated_invariants': cluster_line['mutated_invariants'] + [False],
            'stops': cluster_line['stops'] + [False],
            'mut_freq': cluster_line['mut_freqs'] + [0.0]}

    for gene in 'vdj':
        for pos in ['start', 'end']:
            cluster_line[gene+'_'+pos] = cluster_line['regional_bounds'][gene][pos.startswith('e')]

    cluster_cols = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length']

    sequences = as_dict_rows(cluster_sequences)
    if args.remove_frameshifts or args.remove_stops or args.remove_mutated_invariants:
        sequences = apply_filters(args, sequences)
    cluster = {'sequences': sequences,
               'cdr3_start': cluster_line['codon_positions']['v'],
               'has_seed': seed_id in cluster_line['unique_ids'],
               'n_seqs' : len(sequences),
               'seed_id': seed_id}
    for k in cluster_cols:
        cluster[k] = cluster_line[k]
    for gene in 'vdj':
        for pos in ['start', 'end']:
            cluster[gene+'_'+pos] = cluster_line['regional_bounds'][gene][pos.startswith('e')]
    return cluster


def nth(sequence, n):
    for i, v in enumerate(sequence):
        if i == n:
            return  v

def nth_csv_row(filename, n):
    with open(filename) as handle:
        reader = csv.DictReader(handle)
        return nth(reader, n)


def processed_data(args):
    """Uses args to find the correct partition, cluster pair and all associated information. Cluster
    information is returned as by process_cluster."""

    cp = clusterpath.ClusterPath()
    cp.readfile(args.partition_file)

    # select partition, relative to best partition
    partition_index = cp.i_best + args.partition
    part = nth_csv_row(args.partition_file, partition_index)
    seed_id = cp.seed_unique_id

    # select cluster; unique_ids takes highest precedence
    if args.unique_ids:
        cluster_unique_ids = args.unique_ids
    # default to seed, when possibile
    elif cp.seed_unique_id and not args.cluster:
        cluster_unique_ids = next(cluster.split(':') for cluster in part['partition'].split(';')
                if seed_id in cluster.split(':'))
    # otherwise, assume we have args.cluster or default it to 0
    else:
        clusters = sorted((cluster.split(':') for cluster in part['partition'].split(';')), key=len, reverse=True)
        cluster_unique_ids = clusters[args.cluster or 0]

    # Get cluster annotation and put together into 
    cluster_annotation = get_cluster_annotation(args.cluster_annotation_file, cluster_unique_ids)
    data = {'n_clusters': part['n_clusters'],
            'logprob': part['logprob'],
            'partition_file': args.partition_file,
            'last_modified': time.ctime(os.path.getmtime(args.cluster_annotation_file)),
            'annotation_file': args.cluster_annotation_file}
    if args.seqs_out:
        data['seqs_file'] = os.path.relpath(args.seqs_out, args.paths_relative_to)
    # Process the annotation file specific details/data
    data.update(process_cluster(args, cluster_annotation, seed_id))
    return data


def write_cluster_meta(args, cluster_data):
    def attrs(base):
        return [base + '_' + k for k in ['gene', 'start', 'end']]
    to_keep = ['has_seed', 'seqs_file', 'n_seqs', 'last_modified', 'annotation_file', 'partition_file',
        'cdr3_start', 'cdr3_length'] + attrs('v') + attrs('d') + attrs('j')
    doc = {k: cluster_data[k] for k in to_keep}
    if args.namespace:
        doc = {args.namespace + ':' + k: v for k, v in doc.items()}
    with open(args.cluster_meta_out, 'w') as outfile:
        json.dump(doc, outfile, sort_keys=True, indent=4)

def write_seq_meta(args, cluster_data):
    to_keep = ['unique_id', 'is_seed', 'frameshifted', 'stops', 'mutated_invariants', 'duplicates', 'mut_freq']
    with open(args.seqmeta_out, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=to_keep, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(cluster_data['sequences'])

def write_seqs(args, cluster_data):
    with open(args.seqs_out, 'w') as outfile:
        SeqIO.write(
                (SeqRecord(Seq(sequence['seq']), id=sequence['unique_id'], description='')
                    for sequence in cluster_data['sequences']),
                outfile,
                'fasta')

def parse_args():
    def existing_file(fname):
        """Argparse type for an existing file"""
        if not os.path.isfile(fname):
            raise ValueError("Can't find file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_argument_group(title="Input files", description="(required)")
    inputs.add_argument(
        '--cluster-annotation-file',
        help='cluster annotations file as output by partis',
        type=existing_file, required=True)
    inputs.add_argument(
        '--partition-file',
        help='partitions file as output by partis',
        type=existing_file, required=True)

    outputs = parser.add_argument_group(title="Output files", description="(optional)")
    outputs.add_argument(
        '--seqmeta-out',
        help='per sequence metadata CSV file')
    outputs.add_argument(
        '--seqs-out',
        help='cluster sequences as a FASTA file')
    outputs.add_argument(
        '--cluster-meta-out',
        help='cluster sequences as a JSON file')
    # If we support a recursive option, we have to name I guess?
    #outputs.add_argument(
        #'--process-all-data-to',
        #help="writes all data for all partitions/clusters to the specified directory")

    partis_args = parser.add_argument_group(title="Partis args",
        description="""These arguments (as passed to partis) are required in order to process the data correctly.""")
    partis_args.add_argument(
        '--parameter-dir',
        help='parameter dir path, as passed to partis (if omitted, gls assumed to be ' + default_germline_sets + ')')
    partis_args.add_argument('--locus',
        help='again, as passed to partis',
        required=True)

    cluster_selection_args = parser.add_argument_group(title="Cluster selection args",
        description="""Given a partition file and associated cluster annotation file, there may be multiple
        clusters one might extract data for. These options allow you to specify a selection.""")
    cluster_selection_args.add_argument(
        '--partition',
        type=int, default=0,
        help='"best plus" index; defaults to 0 (best partition); 1 selects the next partition step, etc.')
    cluster_selection_args.add_argument(
        '--cluster',
        type=int,
        help="""index of cluster in partition-file after sorting by cluster size; defaults to seed cluster if
        seeded and 0 (the largest cluster) otherwise.""")
    # add a non sorted version?
    cluster_selection_args.add_argument(
        '--unique-ids',
        help='select a specific cluster using its unique_ids signature')
    cluster_selection_args.add_argument(
        '--unique-ids-file',
        help='select a specific cluster using its unique_ids signature in a single line in a file',
        type=lambda x: file(x).read().strip())


    other_args = parser.add_argument_group(title="Other options")
    other_args.add_argument(
        '--remove-frameshifts',
        help='if set, removes seqs with frameshifted indels from output',
        action="store_true")
    other_args.add_argument(
        '--remove-stops',
        help='if set, removes seqs with stop codons from output',
        action="store_true")
    other_args.add_argument(
        '--remove-mutated-invariants',
        help='if set, removes seqs with mutated "invariant" regions from output',
        action="store_true")
    other_args.add_argument(
        '--indel-reversed-seqs',
        help='if set, uses the "indel_reversed_seqs" output of partis instead of "seqs"',
        action="store_true")
    other_args.add_argument(
        '--paths-relative-to',
        default='.',
        help='files pointed to from metadata.json file will be specified relative to this path')
    other_args.add_argument(
        '--namespace',
        help='namespace to be applied to cluster meta attr names')
    # --indel-reversed-seqs
    # --remove-mutated-invariants

    # parse args and decorate with derived values
    args = parser.parse_args()
    # default paths_relative_to is just whatever the output dir is
    args.unique_ids = args.unique_ids or args.unique_ids_file
    # default germline set (discouraged)
    args.germline_sets = os.path.join(args.parameter_dir, 'hmm/germline-sets') if args.parameter_dir else default_germline_sets
    args.glfo = glutils.read_glfo(args.germline_sets, args.locus)

    return args


def main():
    """
    Run and save cluster file processing.
    """
    args = parse_args()

    cluster_data = processed_data(args)

    if args.seqmeta_out:
        write_seq_meta(args, cluster_data)
    if args.cluster_meta_out:
        write_cluster_meta(args, cluster_data)
    if args.seqs_out:
        write_seqs(args, cluster_data)


if __name__ == '__main__':
    main()


