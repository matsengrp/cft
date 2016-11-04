#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Python modules for processing partis output.

Input annotations csv file and output fasta files.
'''

import pandas as pd
import argparse
import os
import os.path
import itertools
import json


def parse_args():
    ''' parse input arguments '''

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--incsv', dest='incsv',
            help='input cluter annotations csv', type=str)
    parser.add_argument('--cluster_base', dest='cluster_base',
            help='basename for clusters', default='cluster', type=str)
    parser.add_argument('--input_dir', dest='input_dir',
            help='directory with data', type=str)
    parser.add_argument('--output_dir', dest='output_dir',
            help='directory for fasta files', type=str)
    parser.add_argument('--baseline', action='store_true', dest='baseline',
            help='output fasta files in baseline format', default=False)
    parser.add_argument('--separate', action='store_true', dest='separate',
            help='output to per-cluster fasta files', default=False)
    #parser.add_argument('--select_clustering', dest='select_clustering',
    #        help='choose a row from partition file for a different cluster',
    #        default=0, type=int)

    return parser.parse_args()


def create_dir(args):
    ''' create directories if they don't exist '''
    input_base = args.incsv[:args.incsv.index('-cluster-annotations.csv')]
    output_base = os.path.join(args.output_dir,
            input_base.replace(args.input_dir, ''))
    dirname = os.path.dirname(output_base)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    return output_base


def output_clusters(pandas_df, col, id1, id2):
    ''' split clusters and unnest them '''
    nested = [[id1[idx]] + [id2[idx]] + \
            clus.split(':') for idx, clus in enumerate(pandas_df[col])]
    return sum(nested, [])


def interleave_lists(pandas_df, cols, id1, id2, seed_id):
    ''' interleave two lists and prepend '>' to headers '''

    # we'll only have two lists (headers and sequences) so concatenate
    # headers and then dummy sequences
    out_list = []
    for item in range(2):
        out_list += [output_clusters(pandas_df, cols[item],
            id1[item], id2[item])]

    # we hope that the lists constructed are of the same length, but you
    # never know
    assert len(out_list[0]) == len(out_list[1])

    interleaved = [singleton for two_tuple in zip(out_list[0], out_list[1]) \
        for singleton in two_tuple]

    # add '>' to headers for fasta files
    interleaved[::2] = ['>'+('seed_'+item if item==seed_id else item) \
            for item in interleaved[::2]]

    return interleaved


def output_cluster_stats(args, fname, annotations, seed_id):
    ''' output to file various cluster statistics given a seed '''

    for idx, cluster in enumerate(annotations['unique_ids']):
        ids = cluster.split(':')
        if seed_id in ids:
            seed_cluster = idx
            with open(fname+'-cluster-stats.csv', 'wb') as stats:
                stats.write('seed_id,n_clones,cluster_id,clone_ids\n')
                stats.write('{},{},{},{}\n'.format(seed_id, len(ids),
                        args.cluster_base+str(seed_cluster), cluster))
    return seed_cluster


def process_data(args, output_base):
    ''' read data and interleave lists '''

    # Read data
    annotations = pd.read_csv(args.incsv)

    part_file = args.incsv.replace('-cluster-annotations.csv', '.csv')
    seed_id = pd.read_csv(part_file).loc[0]['seed_unique_id']
    seed_cluster = \
            output_cluster_stats(args, output_base, annotations, seed_id)

    #if args.select_clustering > 0:
    #    # we'll need to do a little poking at the partition file and to
    #    # reassign columns 'unique_ids', 'seqs' and 'naive_seq' to their
    #    # preferred partition values
    #    # basically create a dictionary with the annotations and then
    #    # repartition based on the selected partition
    #    annotations = select_different_cluster(args, annotations)

    # define headers for clusters
    nrow = annotations.shape[0]
    blanks = [''] * nrow
    cluster_ids = ['>>'+('seed_' if row==seed_cluster else '')+ \
            args.cluster_base+str(row) for row in range(nrow)]

    # get naive seqs
    naive_ids = ['>naive'+str(row) for row in range(nrow)]
    naive_seq = annotations['naive_seq'].tolist()

    # combine lists
    return interleave_lists(annotations,
            ['unique_ids', 'seqs'],
            [cluster_ids, blanks],
            [naive_ids, naive_seq],
            seed_id), seed_id

def get_json(args, fname, cluster, data, seed_id):
    ''' get metatdata for writing json file '''

    # this only works for heavy chain data, though currently partis
    # does not take light chains, right?
    return {'file': '-'.join([fname, cluster, 'seqs.fa']),
            'cluster_id': cluster,
            'v-gene': data['v_gene'],
            'd-gene': data['d_gene'],
            'j-gene': data['j_gene'],
            'cdr3-length': data['cdr3_length'],
            'seed': seed_id,
            'has_seed': cluster.startswith('seed_')
           }


def write_file(args, fname, in_list, seed_id):
    '''
    Output 1: Kleinstein lab--style fasta files. Group names are given as
    '>>>GROUPNAME', naive seqs are given as '>>NAIVE' and usual sequence
    headers are given as '>SEQUENCENAME'.
    This file can be read directly into BASELINe.

    Output 2: Separate fasta files for each group for tree building and all.
    '''

    annotations = pd.read_csv(args.incsv)
    cluster_id = -1

    list_of_dicts = []
    if args.separate:
        # separate fastas
        for header, seq in itertools.izip(in_list[::2], in_list[1::2]):
            if header.startswith('>>>'):
                # '>>>' denotes start of a separate cluster
                cluster_id += 1
                current_file = '-'.join([fname, header[3:], 'seqs.fa'])
                open(current_file, 'wb').close()
                cluster_dict = get_json(args, fname, header[3:],
                        annotations.iloc[cluster_id], seed_id)
                list_of_dicts.append(cluster_dict)
            else:
                with open(current_file, 'a') as seqs:
                    seqs.write('%s\n%s\n' % (header, seq))
        with open('-'.join([fname, 'data.json']), 'wb') as outfile:
            json.dump(list_of_dicts, outfile)
    else:
        # baseline or all in one
        write_naive = False
        with open('-'.join([fname, 'all-seqs.'+ \
                ('bl' if args.baseline else 'fa')]), 'wb') as seqs:
            for header, seq in itertools.izip(in_list[::2], in_list[1::2]):
                flag = (not header.startswith('>>>') and not \
                        header.startswith('>>') or write_naive)
                if flag:
                    seqs.write('%s\n%s\n' % (header, seq))
                write_naive = header.startswith('>>>seed')


def main():
    ''' run and save cluster file processing '''

    args = parse_args()
    if not os.path.isfile(args.incsv):
        raise ValueError('Invalid file: ' + str(args.incsv))

    output_base = create_dir(args)

    interleaved, seed_id = process_data(args, output_base)

    write_file(args, output_base, interleaved, seed_id)


if __name__ == '__main__':
    main()

