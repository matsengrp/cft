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
import indelutils
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



default_glfo_dir = os.path.join(partis_path, 'data/germlines/human')  # this should only be used as a last resort (e.g. you've completely lost the germline sets corresponding to your deprecated csv output files)

def subset_dict(d, keys):
    return {k: d[k] for k in keys if k in d}

def merge(d1, d2):
    d = d1.copy()
    d.update(d2)
    return d

def as_dict_rows(column_dict, columns=None):
    columns = columns or column_dict.keys()
    column_lengths = [len(column_dict[c]) for c in columns]
    assert min(column_lengths) == max(column_lengths), "columns can't be of different lengths"
    return [{c: column_dict[c][i] for c in columns}
            for i in range(max(column_lengths))]

def apply_filters(args, cluster_line):
    def filter_fn(iseq):
        if args.remove_frameshifts and not cluster_line['in_frames'][iseq]:
            return False
        elif args.remove_stops and cluster_line['stops'][iseq]:
            return False
        elif args.remove_mutated_invariants and cluster_line['mutated_invariants'][iseq]:
            return False
        else:
            return True
    return filter(filter_fn, range(len(cluster_line['seqs'])))

def add_regional_bounds(cluster_line):
    keys_to_add = []
    for gene in 'vdj':
        for pos in ['start', 'end']:
            key = gene+'_'+pos
            keys_to_add.append(key)
            cluster_line[key] = cluster_line['regional_bounds'][gene][pos.startswith('e')]
    return cluster_line, keys_to_add

def seqs(args, cluster_line):
    non_reversed_seqs = cluster_line['seqs']
    if args.indel_reversed_seqs:
        return [reversed_seq or non_reversed_seqs[i]
                for i, reversed_seq in enumerate(cluster_line['indel_reversed_seqs'])]
    else:
        return non_reversed_seqs

def get_cluster_seqs_dict(cluster_line, seed_id, args):
    """ Dict format from cluster line format """
    # add naive values to beginning of per-seq fields
    cluster_sequences = {
            'unique_id':             [args.inferred_naive_name] + cluster_line['unique_ids'],
            'sequence':             [args.inferred_naive_name] + cluster_line['unique_ids'],
            'seq': [cluster_line['naive_seq']] + seqs(args, cluster_line),
            'is_seed':               ['False'] + [(unique_id == seed_id) for unique_id in cluster_line['unique_ids']],
            'duplicates':               [None] + cluster_line['duplicates'],
            'frameshifted':            [False] + [not x for x in cluster_line['in_frames']],
            'mutated_invariants':      [False] + cluster_line['mutated_invariants'],
            'stops':                   [False] + cluster_line['stops'],
            'mut_freq':                  [0.0] + cluster_line['mut_freqs'],
            'affinity':                 [None] + cluster_line.get('affinities', [None for _ in cluster_line['unique_ids']]),
            'timepoint':                [None] + cluster_line['timepoints'],
            'duplicates':               [[]] + cluster_line['duplicates'],
            'multiplicity':           [None] + cluster_line['multiplicities'],
            'timepoints':               [[]] + cluster_line['duplicate_timepoints'],
            'timepoint_multiplicities': [[]] + cluster_line['duplicate_multiplicities']}
    res = as_dict_rows(cluster_sequences)
    #print(cluster_line['unique_ids'][314],cluster_line['duplicate_multiplicities'][314] )
    #print([i for i in res if i['unique_id'] == 'MIG474-349'])
    return res

def get_cluster_meta_dict(cluster_line, seed_id, args):
    return {'sequences': get_cluster_seqs_dict(cluster_line, seed_id, args),
            'cdr3_start': cluster_line['codon_positions']['v'],
            'has_seed': seed_id in cluster_line['unique_ids'],
            'mean_mut_freq': numpy.mean(cluster_line['mut_freqs']),
            'seed_id': seed_id,
            'match_indels_in_uid': args.match_indels_in_uid}

def add_additional_info(cluster_line, additional_per_seq_info, iseqs_to_keep):
    for key in additional_per_seq_info:
        cluster_line[key] = [additional_per_seq_info[key][iseq] for iseq in iseqs_to_keep]
    return cluster_line

def downsample_iseqs_by_multiplicity(cluster_line, multiplicity_seqmeta, max_sequences_count, always_include_ids):
    """ First take the always keep, then take as many as you can of the remaining seqs, in order of highest multiplicity """
    if len(multiplicity_seqmeta['multiplicities']) != len(cluster_line['seqs']):
        raise Exception('Something went wrong internally, mutiplicities are calculated for each seq in the cluster annotation but the number of seqs in the annotation does not match the number of multiplicities')
    always_include_iseqs = [iseq for iseq in range(len(cluster_line['seqs'])) if cluster_line['unique_ids'][iseq] in always_include_ids]
    rest_iseqs = [iseq for iseq in range(len(cluster_line['seqs'])) if cluster_line['unique_ids'][iseq] not in always_include_ids] 
    remaining_seqs_to_take_count = max_sequences_count - len(always_include_ids)
    downsampled_iseqs = always_include_iseqs + \
                        sorted(rest_iseqs,
                        key=lambda iseq: multiplicity_seqmeta['multiplicities'][iseq], # Sort by multiplicity
                        reverse=True)[:remaining_seqs_to_take_count]           # Descending order
    return downsampled_iseqs

def get_upstream_row(upstream_seqmeta, seqid):
    """ Get the corresponding row of the upstream metadata (which contains our full lenght sequence multiplicities...) """
    upstream_seqmeta = upstream_seqmeta or {}
    default_row = {'multiplicity': 1}
    row = merge(default_row, upstream_seqmeta.get(seqid, {}))
    timepoint = row.get('timepoint')
    # if we do not have timepoint info (including when timepoint == ''), then we assign a dummy timepoint so multiplicity is calculated even in the absence of timepoint data.
    row['timepoint'] = timepoint if timepoint and timepoint is not '' else 'no-timepoint'
    return row

def timepoint_multiplicity_mapping(seqid, duplicates, upstream_seqmeta):
    timepoints_dict = collections.defaultdict(lambda: 0)
    for dup_seqid in duplicates:
        dup_upstream_row = get_upstream_row(upstream_seqmeta, seqid)
        # pre-partis filtering multiplicity
        dup_multiplicity = int(dup_upstream_row['multiplicity'])
        #print(seqid, dup_multiplicity)
        # we have handled any case with missing timepoint info in get_upstream_row, so this should always work
        timepoints_dict[dup_upstream_row['timepoint']] += dup_multiplicity
    return timepoints_dict

def get_multiplicity_seqmeta(cluster_line, upstream_seqmeta):
    """"Merge upstream (pre-partis) metadata, indexed in a dict by unique_id, (potentially)
    including timepoint and multiplicity info, with the metadata output of process_partis (partis_seqmeta)."""
    multiplicity_seqmeta = collections.defaultdict(list) #initialize empty array for new per seq fields
    for iseq, seqid in enumerate(cluster_line['unique_ids']):
        upstream_row = get_upstream_row(upstream_seqmeta, seqid)
        duplicates = [seqid] + cluster_line['duplicates'][iseq] 
        timepoints_dict = timepoint_multiplicity_mapping(seqid, duplicates, upstream_seqmeta)
        timepoints = sorted(timepoints_dict.items())
        multiplicity = sum(t[1] for t in timepoints)
        multiplicity_seqmeta['timepoints'].append(upstream_row['timepoint']) #this represents the timepoint this exact sequence was sampled
        multiplicity_seqmeta['duplicates'].append(duplicates)
        multiplicity_seqmeta['multiplicities'].append(multiplicity)
        multiplicity_seqmeta['duplicate_timepoints'].append([tp[0] for tp in timepoints]) #this represents the timepoints duplicates (indentical seqs) were sampled
        multiplicity_seqmeta['duplicate_multiplicities'].append([tp[1] for tp in timepoints])
    return multiplicity_seqmeta

def match_indels_in_uid_seq(cluster_line, match_indels_in_uid):
    iseq_to_match = cluster_line['unique_ids'].index(match_indels_in_uid)
    ifos_to_match = cluster_line['indelfos'][iseq_to_match]['indels']
    if len(ifos_to_match) > 1:
        raise Exception('{} has more than 1 indel. We don\'t have a good way of matching more than one indel between two seqs right now. Use an id of a sequence with exactly one indel or make a suggestion for a good default in this case, or adapt get_iseqs_with_compatible_indels in partis/python/indelutils.py'.format(match_indels_in_uid))
    elif len(ifos_to_match) < 1:
        raise Exception('{} has no indel. Use an id of a sequence with exactly one indel.'.format(match_indels_in_uid))
    else:
        ifo_to_match = ifos_to_match[0]
    iseqs_to_keep = get_iseqs_with_compatible_indels(cluster_line, ifo_to_match)
    if len(iseqs_to_keep) < 1:
        raise Exception('No indels in cluster annotation for cluster containing {} matched {}'.format(match_indels_in_uid, ifo_to_match))
    return iseqs_to_keep 

def check_seed_for_indels(cluster_line, seed_id, partition_file):
    iseq_seed = cluster_line['unique_ids'].index(seed_id)
    ifos = cluster_line['indelfos'][iseq_seed]['indels']
    if len(ifos) > 0:
        print([indelutils.get_dbg_str(ifo) for ifo in ifos])
        raise Exception('indel in seed sequence {}. Options are 1. Look at the annotation for this cluster and find the indel in the seed. Rerun process_partis.py with --match-indels-in-uid <uid-of-seq-containing-indel-of-interest> to process only sequences containing that specific indel for further analysis of the indel 2. Run with --ignore-seed-indels. PS check out {}'.format(seed_id, partition_file))

def process_cluster(args, cluster_line, seed_id, glfo):
    if seed_id is not None and not args.match_indels_in_uid and not args.ignore_seed_indels:
        check_seed_for_indels(cluster_line, seed_id, args.partition_file)
    #assume we want all seqs in cluster
    iseqs_to_keep = set(range(len(cluster_line['seqs'])))
    # various cases where we downsample cluster sequences
    if args.match_indels_in_uid:
        #cluster_line['unique_ids'] = map(lambda i: '{}_indel_filtered'.format(i), cluster_line['unique_ids']) 
        iseqs_to_keep = iseqs_to_keep & set(match_indels_in_uid_seq(cluster_line, args.match_indels_in_uid, glfo))
    
    if args.largest_cluster_across_partitions:
        '''
        Deduplicate sequence records. When using largest_cluster_across_partitions for seeded clusters, we may end up with duplicate sequences in 
        these clusters because of how partis partitions seed clusters. If this option used, beware that this deduplication pays no respect to which 
        duplicate record is preserved of two with the same unique id.
        '''
        iseqs_to_keep = iseqs_to_keep & set({unique_id: iseq for iseq, unique_id in enumerate(cluster_line['unique_ids'])}.values())
    if args.remove_frameshifts or args.remove_stops or args.remove_mutated_invariants:
        iseqs_to_keep = iseqs_to_keep & set(apply_filters(args, cluster_line))

    # apply merging of multiplicity info here (or flesh out with default values otherwise)
    multiplicity_seqmeta = get_multiplicity_seqmeta(cluster_line, args.upstream_seqmeta)
    
    #print([cluster_line['duplicate_multiplicities'][iseq] for iseq, uid, in enumerate(cluster_line['unique_ids']) if uid == 'MIG474-349'])
    # apply sequence downsampling here
    cluster_line['n_unique_seqs'] = len(iseqs_to_keep) + 1 # total in cluster output from partis
    always_include = set(args.always_include + [args.inferred_naive_name])
    if args.max_sequences:
        iseqs_to_keep = downsample_iseqs_by_multiplicity(cluster_line, multiplicity_seqmeta, args.max_sequences, always_include)
    cluster_line['n_sampled_seqs'] = len(iseqs_to_keep)
    #print([cluster_line['duplicate_multiplicities'][iseq] for iseq, uid, in enumerate(cluster_line['unique_ids']) if uid == 'MIG474-349'])
    #filter cluster line to iseqs_to_keep
    cluster_line = utils.restrict_to_iseqs(cluster_line, iseqs_to_keep, glfo)
    
    # add the additional info computed in above for the iseqs we care about
    cluster_line = add_additional_info(cluster_line, multiplicity_seqmeta, iseqs_to_keep) 

    cluster_line['n_total_reads'] = sum(cluster_line['multiplicities']) + 1 #total reads accounting for multiplicity (must be calculated after subsetting cluster in restrict_to_iseqs if it should correspond to total reads represented by subset of cluster returned by restrict_to_iseqs)

    # this needs to happen after restrict_to_iseqs re-adds implicit partis linekeys including 'regional_bounds'
    cluster_line, regional_bounds_keys = add_regional_bounds(cluster_line)
    print([cluster_line['duplicate_multiplicities'][iseq] for iseq, uid, in enumerate(cluster_line['unique_ids']) if uid == 'MIG474-349'])
    return merge(
            subset_dict(cluster_line, regional_bounds_keys + ['n_total_reads', 'n_sampled_seqs', 'n_unique_seqs', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'naive_seq', 'v_per_gene_support', 'd_per_gene_support', 'j_per_gene_support']),
            get_cluster_meta_dict(cluster_line, seed_id, args))

def find_largest_cluster_across_partitions(cpath, annotation_list):
    '''
    Sometimes we'd like to choose the largest cluster across all partitions (not just within a given partition such as the most likely one). 
    This does that, and makes sure to restrict this to seed containing clusters if is a seed unique id.
    '''
    seed = cpath.seed_unique_id
    largest_cluster_len = 0
    for i, partition in enumerate(cpath.partitions):
        clusters_by_size = sorted(partition, key=lambda cluster: len(set(cluster)), reverse=True)
        if seed is not None:
            clusters_by_size = list(filter(lambda c: seed in c, clusters_by_size))
            if len(clusters_by_size) == 0:
                raise Exception(' --largest-cluster-across-partitions specified for a seeded partition and no clusters contain the seed. This should not happen, as both the seed info and the cluster ids are coming from partis here. Make sure the partition file specified is a valid partition that includes the seed sequence.')
        uids_largest_cluster_in_partition = clusters_by_size[0]
        unique_id_count = len(set(uids_largest_cluster_in_partition))
        if unique_id_count > largest_cluster_len:
            uids_largest_cluster = uids_largest_cluster_in_partition
            largest_cluster_len = unique_id_count
            ipart = i
    return uids_largest_cluster, ipart
 
def choose_cluster(partition_file, annotation_list, cpath, ipart=None, i_cluster=None, unique_ids=None):
    """Given a partition file and associated cluster annotation file, there may be multiple
    clusters one might extract data for. These options allow you to specify a selection."""
    
    # partition index is i_best unless specified
    if ipart is None:
        ipart = cpath.i_best 

    # select cluster; unique_ids takes highest precedence
    if unique_ids:
        cluster_unique_ids = unique_ids
    # default to seed, when possibile
    elif cpath.seed_unique_id and not i_cluster:
        cluster_unique_ids = next(cluster for cluster in cpath.partitions[ipart] if cpath.seed_unique_id in cluster)
    # otherwise, assume we have args.cluster or default it to 0
    else:
        clusters = sorted(cpath.partitions[ipart], key=len, reverse=True)
        cluster_unique_ids = clusters[i_cluster or 0]

    # Get cluster annotation and put together into
    annotations = [l for l in annotation_list if l['unique_ids'] == cluster_unique_ids]
    if len(annotations) == 0:
        raise ValueError('requested uids %s not found in %s' % (cluster_unique_ids, partition_file))  # it was a value error before, so I'm leaving it at that
    elif len(annotations) > 1:
        print '%s more than one annotation with requested uids %s found in %s' % (utils.color('red', 'warning'), cluster_unique_ids, partition_file)  # shouldn't be possible
    return annotations[0]

def processed_data(args):
    """Uses args to find the correct partition, cluster pair and all associated information. Cluster
    information is returned as by process_cluster."""

    glfo = None if utils.getsuffix(args.partition_file) == '.yaml' else glutils.read_glfo(args.glfo_dir if args.glfo_dir else default_glfo_dir, args.locus)
    #print("calling utils.read_output with args:", args.partition_file, glfo)
    glfo, annotation_list, cpath = utils.read_output(args.partition_file, glfo=glfo)  # returns glfo from the file if it's there, otherwise it returns the one we passed in
    if annotation_list is None:
        raise Exception('no annotations in %s (probably because cluster annotation file wasn\'t found)' % args.partition_file)

    ipart = args.partition if args.partition is not None else cpath.i_best
    unique_ids = args.unique_ids

    #find the largest cluster across all partitions
    if args.largest_cluster_across_partitions:
        if any([arg is not None for arg in (args.cluster, args.unique_ids, args.partition)]):
            raise Exception('Doesn\'t make sense to specify any of --partition, --cluster, --unique-ids when you use --largest-cluster-across-partitions as it will disregard those options and choose the largest cluster across all partitions.')
        unique_ids, ipart = find_largest_cluster_across_partitions(cpath, annotation_list)
    
    cluster_annotation = choose_cluster(args.partition_file, annotation_list, cpath, ipart, args.cluster, unique_ids)
    
    data = {'n_clusters': len(cpath.partitions[ipart]),
            'logprob': cpath.logprobs[ipart],
            'partition_file': args.partition_file,
            'last_modified': time.ctime(os.path.getmtime(args.partition_file))}    
    if args.seqs_out:
        data['seqs_file'] = os.path.relpath(args.seqs_out, args.paths_relative_to)
    # Process the annotation file specific details/data
    data.update(process_cluster(args, cluster_annotation, cpath.seed_unique_id, glfo))
    return data

def write_cluster_meta(args, cluster_data):
    def attrs(base):
        return [base + '_' + k for k in ['gene', 'start', 'end', 'per_gene_support']]
    to_keep = ['naive_seq', 'has_seed', 'seqs_file', 'n_unique_seqs', 'n_total_reads', 'last_modified', 'partition_file',
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
        yield row

def write_seq_meta(args, cluster_data):
    to_keep = ['unique_id', 'sequence', 'is_seed', 'frameshifted', 'stops', 'mutated_invariants',
            'mut_freq', 'timepoint', 'multiplicity', 'timepoints', 'timepoint_multiplicities', 'duplicates', 'affinity']
    #    to_keep = ['unique_id', 'is_seed', 'frameshifted', 'stops', 'mutated_invariants',
    #            'mut_freq', 'timepoint', 'multiplicity', 'duplicate_timepoints', 'duplicate_multiplicities', 'duplicates', 'affinity']
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
        help='DEPRECATED use --glfo-dir instead')
    partis_args.add_argument(
        '--glfo-dir',
        help='path to germline info, only necessary for deprecated .csv output files (if omitted, defaults to %s, which may or may not work)' %  default_glfo_dir)
    partis_args.add_argument('--locus',
        help='again, as passed to partis',
        required=True)

    cluster_selection_args = parser.add_argument_group(title="Cluster selection args",
        description="""Given a partition file and associated cluster annotation file, there may be multiple
        clusters one might extract data for. These options allow you to specify a selection.""")
    cluster_selection_args.add_argument(
        '--partition',
        type=int,
        help='partition step index.')
    cluster_selection_args.add_argument(
        '--cluster',
        type=int,
        help="""index of cluster in partition-file after sorting by cluster size; defaults to seed cluster if
        seeded and 0 (the largest cluster) otherwise.""")
    cluster_selection_args.add_argument(
        '--largest-cluster-across-partitions',
        help='select the largest cluster across all partitions. Must include seed if partis was run with a seed.',
        action="store_true")
    # add a non sorted version?
    cluster_selection_args.add_argument(
        '--unique-ids',
        help='select a specific cluster using its unique_ids signature')
    cluster_selection_args.add_argument(
        '--unique-ids-file',
        help='select a specific cluster using its unique_ids signature in a single line in a file',
        type=lambda x: file(x).read().strip())

    seqs_args = parser.add_argument_group(title="Options regarding methods for choosing sequences to be included or left out of the cluster of interest.")
    seqs_args.add_argument(
        '--match-indels-in-uid',
        help='process only sequences matching the one indel in the sequence corresponding to the uid passed here in the annotation chosen in choose_cluster()',
        type=str)
    seqs_args.add_argument(
        '--ignore-seed-indels',
        help='If --match-indels-in-uid has not been set, this allows processing of a seed cluster (without filtering) where there is an indel in the seed sequence.',
        action="store_true")
    seqs_args.add_argument(
        '--remove-frameshifts',
        help='if set, removes seqs with frameshifted indels from output',
        action="store_true")
    seqs_args.add_argument(
        '--remove-stops',
        help='if set, removes seqs with stop codons from output',
        action="store_true")
    seqs_args.add_argument(
        '--remove-mutated-invariants',
        help='if set, removes seqs with mutated "invariant" regions from output',
        action="store_true")
    seqs_args.add_argument(
        '--indel-reversed-seqs',
        help='if set, uses the "indel_reversed_seqs" output of partis instead of "seqs"',
        action="store_true")
    seqs_args.add_argument(
        '--max-sequences',
        help="""if set, downsamples semi-randomly, with preference towards sequences with higher multiplicity
        and order output by partis""",
        type=int)
    seqs_args.add_argument(
        '--always-include',
        type=lambda x: x.split(','), help='comma separated list of ids to keep if --max-sequences is set', default=[])
    
    other_args = parser.add_argument_group(title="Other options")
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

    # parse args and decorate with derived values
    args = parser.parse_args()
    # default paths_relative_to is just whatever the output dir is
    args.unique_ids = args.unique_ids or args.unique_ids_file
    if args.parameter_dir is not None:
        if args.glfo_dir is None:
            args.glfo_dir = os.path.join(args.parameter_dir, 'hmm', glutils.glfo_dir)  # thre reason you should be passing --glfo-dir instead is to avoid these path shenanigans, since partis writes the glfo-dir to the meta info yaml file, so you can just use that
        else:
            raise Exception('doesn\'t make sense to pass both --parameter-dir and --glfo-dir (just use the latter)')
        delattr(args, 'parameter_dir')

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


