#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python modules for processing partis output.

Input annotations csv file and output fasta files.
"""

import pandas as pd
import argparse
import os
import os.path
import itertools
import json
import time
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import textwrap


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


def parse_args():

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError('Invalid file: ' + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--annotations',
        help='input cluster annotations csv',
        type=existing_file, required=True)
    parser.add_argument(
        '--partition',
        help='input cluster partition csv',
        type=existing_file, required=True)
    parser.add_argument(
        '--cluster_base',
        help='basename for clusters',
        default='cluster')
    parser.add_argument(
        '--output_dir',
        default='.',
        help='directory for output files')
    log_or_param_dir = parser.add_mutually_exclusive_group(required=True)
    log_or_param_dir.add_argument(
        '--partis_log',
        help='log file containing relevant information about partis run (required if --param_dir not specified)',
        type=existing_file)
    log_or_param_dir.add_argument(
        '--param_dir',
        help='parameter directory passed to partis call',
        type=str)
    #parser.add_argument('--select_clustering', dest='select_clustering',
    #        help='choose a row from partition file for a different cluster',
    #        default=0, type=int)

    return parser.parse_args()


def calculate_bounds(line, glfo):
    """
    Calculate index of cdr3 region as well as start and end index
    for v, d and j genes.

    For light chains, partis outputs bounds for the d region as a zero-length
    set of bounds starting and ending at the same place (i.e., absent any indels
    this will mean the v region ends where the j region begins).
    """

    utils.process_input_line(line)
    utils.add_implicit_info(glfo, line)
    return line['codon_positions']['v'], line['regional_bounds']


def process_log_file(log_file):
    """
    Get function call from log file that will include information about
    location of inferred germlines, which chain we're using, etc.
    """

    # assumes function call is first line of log file
    # eventually i think duncan wants to append "CALL" or something like it
    # to the function call so it can be reliably searched for.
    # but for now we'll assume it's the first line, which it is.
    with open(log_file, 'r') as partis_log:
        call = partis_log.readline()

    call_args = call.split()
    if not any([list_item.startswith('--') for list_item in call_args]):
        # if there are command line arguments in the first line then
        # it is probably a command. otherwise this won't have any good
        # information for us.
        raise Exception('First line of provided log file not a valid partis command: {}'.format(call))

    if not '--chain' in call_args:
        # partis default is heavy chain
        chain = 'h'
    else:
        chain = call_args[1+call_args.index('--chain')]

    if not '--parameter-dir' in call_args:
        # currently we use IMGT germlines if no cached parameters provided.
        # we have no other way of getting this information since it's printed
        # to stdout if it's not provided. should we assume it's always
        # provided?
        inferred_gls = partis_path + '/data/germlines/human'
    else:
        inferred_gls = call_args[1+call_args.index('--parameter-dir')] + \
                '/hmm/germline-sets'

    # if the parameter file is not an absolute path then we don't know where
    # to look since we don't know where partis was run from!
    # so for now we'll spit an error
    if not os.path.isabs(inferred_gls):
        raise ValueError('Parameter directory must be an absolute path: ' \
                + str(inferred_gls))

    # even if the path is absolute it might have been deleted...
    if not os.path.isdir(inferred_gls):
        raise ValueError('Invalid parameter directory: ' + str(inferred_gls))

    return chain, inferred_gls


def process_data(annot_file, part_file, chain, glpath):
    """
    Melt data into dataframe from annotations and partition files
    """

    seed_ids = []
    part_df = pd.read_csv(part_file, converters={'seed_unique_id':str})
    if 'seed_unique_id' in part_df.columns:
        seed_ids = part_df.loc[0]['seed_unique_id'].split(':')

    #if args.select_clustering > 0:
    #    # we'll need to do a little poking at the partition file and to
    #    # reassign columns 'unique_ids', 'seqs' and 'naive_seq' to their
    #    # preferred partition values
    #    # basically create a dictionary with the annotations and then
    #    # repartition based on the selected partition
    #    annotations = select_different_cluster(args, annotations)

    output_df = pd.DataFrame()
    annotations = pd.read_csv(annot_file, dtype=object)
    glfo = glutils.read_glfo(glpath, chain=chain)
    for idx, cluster in annotations.fillna('').iterrows():
        current_df = pd.DataFrame()
        seq_ids = cluster['unique_ids'].split(':')
        unique_ids = [('seed_' if seq_id in seed_ids else '')+seq_id for seq_id in seq_ids]
        unique_ids.append('naive'+str(idx))
        seqs = cluster['seqs'].split(':')
        seqs.append(cluster['naive_seq'])

        # sometimes we'll have duplicate IDs, which is a no-no for
        # FastTree, so we'll group them together with set()
        idseqs = list(set(zip(unique_ids, seqs)))
        unique_ids, seqs = zip(*idseqs)
        current_df['unique_ids'] = unique_ids
        current_df['seqs'] = seqs
        current_df['cluster'] = str(idx)
        current_df['has_seed'] = any(seed_id in seq_ids for seed_id in \
                seed_ids)
        current_df['v_gene'] = cluster['v_gene']
        current_df['d_gene'] = cluster['d_gene']
        current_df['j_gene'] = cluster['j_gene']
        current_df['cdr3_length'] = str(cluster['cdr3_length'])
        current_df['seed_ids'] = ':'.join(seed_ids)
        cdr3, region_bounds = \
                calculate_bounds(cluster.to_dict(), glfo)
        current_df['cdr3_start'] = cdr3
        for gene in 'vdj':
            for pos in ['start', 'end']:
                current_df[gene+'_'+pos] = region_bounds[gene][pos.startswith('e')]
        output_df = pd.concat([output_df, current_df])

    return output_df


def write_json(df, fname, mod_date, cluster_base, annotations, partition, meta):
    """
    Write metatdata to json file from dataframe
    """

    # We want to know the individual (patient) identifier, the
    # timepoint (when sampled), the seed sequence name, and the
    # gene name.  Ideally this would be explicitly provided in the
    # json data, but instead we must extract it from the
    # paths provided in the json data.
    #
    # Given a partial path like, "QA255.016-Vh/Hs-LN2-5RACE-IgG-new"
    # "QA255"       - individual (patient) identification
    # "016"         - seed sequence name
    # "Vh"          - gene name.  Vh and Vk are the V heavy chain (IgG) vs. V kappa (there is also Vl which mean V lambda)
    # "LN2"         - timepoint.  LN1 and LN2 are early timepoints and LN3 and LN4 are late timepoints.
    #
    # "Hs-LN2-5RACE-IgG-new" is the name of the sequencing run,
    #
    # This currently will only work if the output files are spat into a directory
    # of the form"/path/to/output/QA255.016-Vh/Hs-LN2-5RACE-IgG-new/*.fa"
    # Otherwise no new fields will be added.

    def merge_two_dicts(dict1, dict2):
        """
        Merge two dictionaries into a copy
        """
        merged_dict = dict1.copy()
        merged_dict.update(dict2)
        return merged_dict


    def jsonify(df, cluster_id, mod_date, cluster_base, meta):
        data = df.iloc[0]
        return merge_two_dicts({
            'file': cluster_base+cluster_id+'.fa',
            'cluster_id': cluster_id,
            'v_gene': data['v_gene'],
            'v_start': data['v_start'],
            'v_end': data['v_end'],
            'd_gene': data['d_gene'],
            'd_start': data['d_start'],
            'd_end': data['d_end'],
            'j_gene': data['j_gene'],
            'j_start': data['j_start'],
            'j_end': data['j_end'],
            'cdr3_length': data['cdr3_length'],
            'cdr3_start': data['cdr3_start'],
            'seed': data['seed_ids'],
            'has_seed': str(data['has_seed']),
            'last_modified': time.ctime(mod_date),
            'annotation_file': annotations,
            'partition_file': partition
        }, meta)

    arr = [jsonify(g, k, mod_date, cluster_base, meta) \
            for k,g in df.groupby(['cluster'])]
    with open(fname, 'wb') as outfile:
        json.dump(arr,
                  outfile,
                  sort_keys=True,
                  indent=4,
                  separators=(',', ': '))


def iter_seqs(df):
    for row in df.itertuples():
        yield SeqRecord(Seq(row.seqs), id=row.unique_ids, description='')


def write_fasta(df, fname):
    print("writing {}".format(fname))
    with open(fname, 'w') as fh:
        SeqIO.write(iter_seqs(df), fh, "fasta")

    
def write_separate_fasta(df, output_dir, cluster_base):
    for k,g in df.groupby(['cluster']):
        fname = os.path.join(output_dir, cluster_base+'{}.fa'.format(k))
        write_fasta(g, fname)


def write_melted_partis(df, fname):
    """
    Cassie and other Overbaugh group members sometimes prefer output
    in a csv file with cluster assignments so they can be sorted by other
    metadata related to the read.

    Here we'll just scrape out the cluster information and read IDs and
    output them into a flattened csv file.
    """

    df.to_csv(
        fname,
        index=False,
        columns=[
            'unique_ids', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length',
            'cluster', 'has_seed'
        ])


def main():
    """
    Run and save cluster file processing.
    """

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # get metadata
    regex = re.compile(r'^(?P<pid>[^.]*).(?P<seedid>[0-9]*)-(?P<gene>[^/]*)/[^-]*-(?P<timepoint>[^-]*)')
    m = regex.match('/'.join(args.output_dir.split('/')[-2:]))
    meta = {}
    if m:
        meta = m.groupdict()

    if args.partis_log is not None:
        chain, inferred_gls = process_log_file(args.partis_log)
    else:
        chain = meta['gene'][1].lower()
        inferred_gls = args.param_dir + '/hmm/germline-sets'

    melted_annotations = process_data(args.annotations,
                                      args.partition,
                                      chain,
                                      inferred_gls)

    write_separate_fasta(melted_annotations,
                         args.output_dir,
                         args.cluster_base)

    write_melted_partis(melted_annotations,
                        os.path.join(args.output_dir, 'melted.csv'))

    write_json(melted_annotations,
               os.path.join(args.output_dir, 'metadata.json'),
               os.path.getmtime(args.annotations),
               args.cluster_base,
               args.annotations,
               args.partition,
               meta)


if __name__ == '__main__':
    main()
