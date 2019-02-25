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
import collections
import numpy
#import itertools as it
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


# Util for reading csv files from command line arguments
def csv_reader(index=None, filter_by=None):
    "Returns a function which reads csv data from a filename; constructor closes"
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



default_germline_sets = os.path.join(partis_path, 'data/germlines/human')

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


def get_upstream_row(upstream_seqmeta, seqid):
    upstream_seqmeta = upstream_seqmeta or {}
    row = upstream_seqmeta.get(seqid, {'timepoint': '', 'multiplicity': 1})
    row['sequence'] = seqid
    return row


def merge_upstream_seqmeta(partis_seqmeta, upstream_seqmeta):
    """"Merge upstream (pre-partis) metadata, indexed in a dict by unique_id, (potentially)
    including timepoint and multiplicity info, with the metadata output of process_partis (partis_seqmeta)."""
    n_seqs_multiplicity = 0
    # For each row of our partis sequence metadata (as output from process_partis)
    sequences = [] 
    for row in partis_seqmeta:
        seqid = row['unique_id']
        # We get the corresponding row of the upstream metadata (which contains our full lenght sequence multiplicities...)
        upstream_row = get_upstream_row(upstream_seqmeta, seqid)
        duplicates = filter(None, (row.get('duplicates') or '').split(':'))
        seqids = [seqid] + duplicates
        timepoints_dict = collections.defaultdict(lambda: 0)
        for dup_seqid in seqids:
            dup_upstream_row = get_upstream_row(upstream_seqmeta, seqid)
            # pre-partis filtering multiplicity
            dup_multiplicity = int(dup_upstream_row.get('multiplicity', 1))
            timepoints_dict[dup_upstream_row.get('timepoint')] += dup_multiplicity
        timepoints = sorted(timepoints_dict.items())
        multiplicity = sum(t[1] for t in timepoints)
        n_seqs_multiplicity += multiplicity
        result_row = {
                'unique_id': row['unique_id'],
                'sequence': seqid,
                'timepoint': upstream_row.get('timepoint'),
                'duplicates': seqids,
                'mut_freq': row['mut_freq'],
                'is_seed': row.get('is_seed'),
                'seq': row['seq'],
                'multiplicity': multiplicity,
                'timepoints': [tp[0] or '' for tp in timepoints],
                'timepoint_multiplicities': [tp[1] for tp in timepoints],
                'affinity': row.get('affinity')}
        # Currently arbitrary other upstream seqmeta isn't being merged in here, but could easily be
        sequences.append(result_row)
    return sequences, n_seqs_multiplicity 


def downsample_sequences(args, sequences):
    if args.max_sequences:
        always_include = set(args.always_include + [args.inferred_naive_name])
        always_include_seqs = filter(lambda x: x.get('unique_id') in always_include, sequences)
        rest_seqs = filter(lambda x: x.get('unique_id') not in always_include, sequences)
        # first take the always keep, then take as many as you can of the remaining seqs, in order of highest multiplicity
        return always_include_seqs + \
               sorted(rest_seqs,
                      # Sort by negative so we take the highest multiplicity (lowest neg value) first 
                      key=lambda seqmeta: - seqmeta['multiplicity'])[0:args.max_sequences - len(always_include_seqs)]
    else:
        return sequences


def subset_dict(d, keys):
    return {k: d[k] for k in keys if k in d}

def merge(d1, d2):
    d = d1.copy()
    d.update(d2)
    return d

def process_cluster(args, cluster_line, seed_id):
    #print("calling utils.add_implicit_info with args:", args.glfo, cluster_line)
    #utils.add_implicit_info(args.glfo, cluster_line)

    cluster_sequences = {
            'unique_id':             [args.inferred_naive_name] + cluster_line['unique_ids'],
            'seq': [cluster_line['naive_seq']] + seqs(args, cluster_line),
            'is_seed':               ['False'] + [(unique_id == seed_id) for unique_id in cluster_line['unique_ids']],
            'duplicates':               [None] + [':'.join(x) for x in cluster_line['duplicates']],
            'frameshifted':            [False] + [not x for x in cluster_line['in_frames']],
            'mutated_invariants':      [False] + cluster_line['mutated_invariants'],
            'stops':                   [False] + cluster_line['stops'],
            'mut_freq':                  [0.0] + cluster_line['mut_freqs'],
            'affinity':                 [None] + cluster_line['affinities']}

    for gene in 'vdj':
        for pos in ['start', 'end']:
            cluster_line[gene+'_'+pos] = cluster_line['regional_bounds'][gene][pos.startswith('e')]

    cluster_cols = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length']

    sequences = as_dict_rows(cluster_sequences)
    if args.remove_frameshifts or args.remove_stops or args.remove_mutated_invariants:
        sequences = apply_filters(args, sequences)

    n_seqs = len(sequences)

    # apply merging of multiplicity info here (or flesh out with default values otherwise)
    sequences, n_seqs_multiplicity = merge_upstream_seqmeta(sequences, args.upstream_seqmeta)
    # apply sequence downsampling here
    sequences = downsample_sequences(args, sequences)

    cluster = merge(
            subset_dict(cluster_line, ['naive_seq', 'v_per_gene_support', 'd_per_gene_support', 'j_per_gene_support']),
            {'sequences': sequences,
             'cdr3_start': cluster_line['codon_positions']['v'],
             'has_seed': seed_id in cluster_line['unique_ids'],
             # total in cluster output from partis
             'n_seqs' : n_seqs,
             # n_seqs_multiplicity represents the sum of all sequence multiplicities in the cluster
             'n_seqs_multiplicity' : n_seqs_multiplicity,
             # Should also add mean_mut_freq etc here
             'mean_mut_freq': numpy.mean(cluster_line['mut_freqs']),
             # Should be equal unless downsampled
             'n_sampled_seqs': len(sequences),
             'seed_id': seed_id})
    for k in cluster_cols:
        cluster[k] = cluster_line[k]
    for gene in 'vdj':
        for pos in ['start', 'end']:
            cluster[gene+'_'+pos] = cluster_line['regional_bounds'][gene][pos.startswith('e')]
    return cluster



def processed_data(args):
    """Uses args to find the correct partition, cluster pair and all associated information. Cluster
    information is returned as by process_cluster."""

    print("calling utils.read_output with args:", args.partition_file, args.glfo)
    file_glfo, annotation_list, cpath = utils.read_output(args.partition_file, glfo=args.glfo)
    if annotation_list is None:
        raise Exception('cluster annotation file not found')
    if file_glfo:  # will only be set if we're reading a yaml file
        args.glfo = file_glfo

    # select partition, relative to best partition
    ipart = cpath.i_best + args.partition

    # select cluster; unique_ids takes highest precedence
    if args.unique_ids:
        cluster_unique_ids = args.unique_ids
    # default to seed, when possibile
    elif cpath.seed_unique_id and not args.cluster:
        cluster_unique_ids = next(cluster for cluster in cpath.partitions[ipart] if cpath.seed_unique_id in cluster)
    # otherwise, assume we have args.cluster or default it to 0
    else:
        clusters = sorted(cpath.partitions[ipart], key=len, reverse=True)
        cluster_unique_ids = clusters[args.cluster or 0]

    # Get cluster annotation and put together into
    annotations = [l for l in annotation_list if l['unique_ids'] == cluster_unique_ids]
    if len(annotations) == 0:
        raise ValueError('requested uids %s not found in %s' % (cluster_unique_ids, args.partition_file))  # it was a value error before, so I'm leaving it at that
    elif len(annotations) > 1:
        print '%s more than one annotation with requested uids %s found in %s' % (utils.color('red', 'warning'), cluster_unique_ids, args.partition_file)  # shouldn't be possible
    cluster_annotation = annotations[0]
    data = {'n_clusters': len(cpath.partitions[ipart]),
            'logprob': cpath.logprobs[ipart],
            'partition_file': args.partition_file,
            'last_modified': time.ctime(os.path.getmtime(args.partition_file))}
    if args.seqs_out:
        data['seqs_file'] = os.path.relpath(args.seqs_out, args.paths_relative_to)
    # Process the annotation file specific details/data
    data.update(process_cluster(args, cluster_annotation, cpath.seed_unique_id))
    return data



def write_cluster_meta(args, cluster_data):
    def attrs(base):
        return [base + '_' + k for k in ['gene', 'start', 'end', 'per_gene_support']]
    to_keep = ['naive_seq', 'has_seed', 'seqs_file', 'n_seqs', 'n_seqs_multiplicity', 'last_modified', 'partition_file',
        'cdr3_start', 'cdr3_length', 'mean_mut_freq'] + attrs('v') + attrs('d') + attrs('j')
    doc = subset_dict(cluster_data, to_keep)
    for gene in 'vdj':
        attr = gene + '_per_gene_support'
        base = 'cft.gene_support:' if args.namespace else ''
        doc[attr] = [{base + 'gene': k,
                      base + 'prob': v}
                     for k, v in doc[attr].items()]
    if args.namespace:
        doc = {args.namespace + ':' + k: v for k, v in doc.items()}
    with open(args.cluster_meta_out, 'w') as outfile:
        json.dump(doc, outfile, sort_keys=True, indent=4)


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

def write_seq_meta(args, cluster_data):
    to_keep = ['unique_id', 'sequence', 'is_seed', 'frameshifted', 'stops', 'mutated_invariants',
            'mut_freq', 'timepoint', 'multiplicity', 'timepoints', 'timepoint_multiplicities', 'duplicates', 'affinity']
    with open(args.seqmeta_out, 'w') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=to_keep, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(format_results(cluster_data['sequences']))

def write_seqs(args, cluster_data):
    with open(args.seqs_out, 'w') as outfile:
        SeqIO.write(
                (SeqRecord(Seq(sequence.get('seq', '')), id=sequence['unique_id'], description='')
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
        '--partition-file',
        help='partitions file as output by partis',
        type=existing_file, required=True)
    inputs.add_argument(
        '--upstream-seqmeta',
        help='optionally, specify upstream seqmeta as a csv with cols: unique_id,timepoint,multiplicity',
        # Index rows by unique id
        type=csv_reader('unique_id'))

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
        '--max-sequences',
        help="""if set, downsamples semi-randomly, with preference towards sequences with higher multiplicity
        and order output by partis""",
        type=int)
    other_args.add_argument(
        '--always-include',
        type=lambda x: x.split(','), help='comma separated list of ids to keep if --max-sequences is set', default=[])
    other_args.add_argument(
        '--paths-relative-to',
        default='.',
        help='files pointed to from metadata.json file will be specified relative to this path')
    other_args.add_argument(
        '--namespace',
        help='namespace to be applied to cluster meta attr names')
    other_args.add_argument(
        '--inferred-naive-name',
        help='see scons option help')
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

    for seq in cluster_data['sequences']:
        if not seq.get('seq'):
            print seq

    if args.seqmeta_out:
        write_seq_meta(args, cluster_data)
    if args.cluster_meta_out:
        write_cluster_meta(args, cluster_data)
    if args.seqs_out:
        write_seqs(args, cluster_data)


if __name__ == '__main__':
    main()


