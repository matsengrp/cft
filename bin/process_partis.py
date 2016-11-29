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
        '''
        Argparse type for an existing file
        '''
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
    parser.add_argument(
        '--baseline',
        action='store_true',
        help='output fasta files in baseline format',
        default=False)
    parser.add_argument(
        '--separate',
        action='store_true',
        help='output to per-cluster fasta files',
        default=False)
    parser.add_argument(
        '--chain',
        help='type of chain data used (h, k, l)',
        default='h')
    #parser.add_argument('--select_clustering', dest='select_clustering',
    #        help='choose a row from partition file for a different cluster',
    #        default=0, type=int)

    return parser.parse_args()


def calculate_cdr3_bounds(line, glfo):
    utils.process_input_line(line)
    line['v_gene'] = line['v_gene'].split('+')[0]
    utils.add_implicit_info(glfo, line)
    return line['codon_positions']['v']


def process_data(annot_file, part_file, chain):
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
    glfo = glutils.read_glfo(partis_path + '/data/germlines/human', chain=chain)
    for idx, cluster in annotations.fillna('').iterrows():
        current_df = pd.DataFrame()
        seq_ids = cluster['unique_ids'].split(':')
        unique_ids = [('seed_' if seq_id in seed_ids else '')+seq_id for seq_id in seq_ids]
        unique_ids.append('naive'+str(idx))
        current_df['unique_ids'] = unique_ids
        seqs = cluster['seqs'].split(':')
        seqs.append(cluster['naive_seq'])
        current_df['seqs'] = seqs
        current_df['cluster'] = str(idx)
        current_df['has_seed'] = any(seed_id in seq_ids for seed_id in \
                seed_ids)
        current_df['v_gene'] = cluster['v_gene']
        current_df['d_gene'] = cluster['d_gene']
        current_df['j_gene'] = cluster['j_gene']
        current_df['cdr3_length'] = str(cluster['cdr3_length'])
        current_df['seed_ids'] = ':'.join(seed_ids)
        current_df['cdr3_start'] = calculate_cdr3_bounds(cluster.to_dict(), glfo)
        output_df = pd.concat([output_df, current_df])

    return output_df


def write_json(df, fname, mod_date, cluster_base, annotations, partition):
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
    # Otherwise entries will be empty.

    regex = re.compile(r'^(?P<pid>[^.]*).(?P<seedid>[0-9]*)-(?P<gene>[^/]*)/[^-]*-(?P<timepoint>[^-]*)')
    m = regex.match('/'.join(fname.split('/')[-3:-1]))
    if m:
        meta = {
            'pid': m.group('pid'),
            'seedid': m.group('seedid'),
            'gene': m.group('gene'),
            'timepoint': m.group('timepoint')
            }
    else:
        meta = {
            'pid': None,
            'seedid': None,
            'gene': None,
            'timepoint': None
            }

    def jsonify(df, cluster_id, mod_date, cluster_base, meta):
        data = df.iloc[0]
        return {
            'file': cluster_base+cluster_id+'.fa',
            'cluster_id': cluster_id,
            'v_gene': data['v_gene'],
            'd_gene': data['d_gene'],
            'j_gene': data['j_gene'],
            'cdr3_length': data['cdr3_length'],
            'cdr3_start': data['cdr3_start'],
            'pid': meta['pid'],
            'seedid': meta['seedid'],
            'gene': meta['gene'],
            'timepoint': meta['timepoint'],
            'seed': data['seed_ids'],
            'has_seed': str(data['has_seed']),
            'last_modified': time.ctime(mod_date),
            'annotation_file': annotations,
            'partition_file': partition
        }

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

    melted_annotations = process_data(args.annotations,
                                      args.partition,
                                      args.chain)

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
               args.partition)



if __name__ == '__main__':
    main()
