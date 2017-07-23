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
import warnings
import json
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tripl import tripl

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
        '--cluster-base',
        help='basename for clusters',
        default='cluster')
    parser.add_argument(
        '--output-dir',
        default='.',
        help='directory for output files')
    parser.add_argument(
        '--paths-relative-to',
        help='files pointed to from metadata.json file will be speicfied relative to this path; defaults to --output-dir')
    parser.add_argument(
        '--melted-base',
        help='basename for melted data output',
        default='melted')
    parser.add_argument(
        '-F', '--remove-frameshifts',
        help='if set, tries to remove seqs with frameshift indels from the output',
        action="store_true")
    parser.add_argument(
        '--partis-log',
        help='log file containing relevant information about partis run (required if --param_dir not specified)',
        required=True,
        type=existing_file)
    #parser.add_argument('--select_clustering', dest='select_clustering',
    #        help='choose a row from partition file for a different cluster',
    #        default=0, type=int)

    # default paths_relative_to is just whatever the output dir is
    args = parser.parse_args()
    args.paths_relative_to = args.paths_relative_to or args.output_dir

    return parser.parse_args()


def process_log_file(log_file):
    """
    Get function call from log file that will include information about
    location of inferred germlines, which locus we're using, etc.
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

    if not '--locus' in call_args:
        # partis default is heavy chain
        warnings.warn("Couldn't infer locus from log file; assuming igh")
        locus = 'igh'
    else:
        locus = call_args[1+call_args.index('--locus')]

    if not '--parameter-dir' in call_args:
        # currently we use IMGT germlines if no cached parameters provided.
        # we have no other way of getting this information since it's printed
        # to stdout if it's not provided. should we assume it's always
        # provided?
        inferred_gls = partis_path + '/data/germlines/human'
    else:
        inferred_gls = call_args[1+call_args.index('--parameter-dir')] + \
                '/hmm/germline-sets'
    # Hacky temp fix for the dataset rename 2017/04/03
    inferred_gls = inferred_gls.replace('kate-qrs-2016-09-09', 'kate-qrs')
    inferred_gls = inferred_gls.replace('laura-mb-2016-12-22', 'laura-mb')

    # if the parameter file is not an absolute path then we don't know where
    # to look since we don't know where partis was run from!
    # so for now we'll spit an error
    # This may not make as much sense now that we've moved into datascripts, and is partof why this shouold be
    # here
    if not os.path.isabs(inferred_gls):
        raise ValueError('Parameter directory must be an absolute path: ' \
                + str(inferred_gls))

    # even if the path is absolute it might have been deleted...
    if not os.path.isdir(inferred_gls):
        raise ValueError('Invalid parameter directory: ' + str(inferred_gls))

    return locus, inferred_gls


def indel_offset(indelfo):
    return sum(map(lambda indel: indel['len'] * (1 if indel['type'] == 'insertion' else -1),
                   indelfo['indels']))

def infer_frameshifts(line):
    "Infer frameshifts based on `in_frame` attr in partis output, or if more than one stop codon."
    attrs = ['in_frames', 'input_seqs']
    def infer_(args):
        in_frame, input_seq = args
        aa = Seq(input_seq).translate()
        stop_count = aa.count("*")
        return (not in_frame) or stop_count > 1
    return map(infer_,
               zip(*map(lambda x: line[x], attrs)))


def process_data(annot_file, part_file, locus, glpath):
    """
    Melt data into dataframe from annotations and partition files
    """

    seed_ids = []
    part_df = pd.read_csv(part_file, converters={'seed_unique_id': str})
    # I _think_ this is the only place we use the the part_df, at least here
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
    glfo = glutils.read_glfo(glpath, locus)
    to_keep = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length']
    for idx, cluster in annotations.fillna('').iterrows():
        current_df = pd.DataFrame()
        line = cluster.to_dict()
        utils.process_input_line(line)
        utils.add_implicit_info(glfo, line)
        current_df['unique_ids'] = line['unique_ids'] + ['naive{}'.format(idx)]
        current_df['seqs'] = line['input_seqs'] + [line['naive_seq']]
        current_df['duplicates'] = [':'.join(x) for x in line['duplicates']] + [None]
        current_df['frameshifts'] = infer_frameshifts(line) + [False] # mocking last entry
        current_df['mut_freqs'] = line['mut_freqs'] + [0.0] # mocking last entry
        for col in to_keep:
            current_df[col] = line[col]
        current_df['cluster'] = idx
        current_df['has_seed'] = any(seed_id in line['unique_ids'] for seed_id in \
                seed_ids)
        current_df['seed_ids'] = ':'.join(seed_ids)
        #current_df['seed_ids'] = map(seed_ids)
        current_df['cdr3_start'] = line['codon_positions']['v']
        for gene in 'vdj':
            for pos in ['start', 'end']:
                current_df[gene+'_'+pos] = line['regional_bounds'][gene][pos.startswith('e')]
        output_df = pd.concat([output_df, current_df])

    return output_df


def write_json(df, fname, mod_date, cluster_base, annotations, partition, outdir, paths_relative_to):
    """
    Write metatdata to json file from dataframe
    """

    def merge_two_dicts(dict1, dict2):
        """
        Merge two dictionaries into a copy
        """
        merged_dict = dict1.copy()
        merged_dict.update(dict2)
        return merged_dict


    def jsonify(df, cluster_id, mod_date, cluster_base):
        data = df.iloc[0]
        def attrs(base):
            return [base + '_' + k for k in ['gene', 'start', 'end']]
        data_ = {k: data[k] for k in attrs('v') + attrs('d') + attrs('j') + ['cdr3_length', 'cdr3_start']}
        data_['has_seed'] = bool(data['has_seed'])
        return tripl.namespaced('cft.cluster',
            id = cluster_id,
            seqs_file = os.path.relpath(os.path.join(outdir, cluster_base+str(cluster_id)+'.fa'), paths_relative_to),
            n_seqs = len(df), # Note... this count naive; good idea?
            last_modified = time.ctime(mod_date),
            annotation_file = annotations,
            partition_file = partition,
            **data_)

    arr = [jsonify(g, k, mod_date, cluster_base) \
            for k,g in df.groupby(['cluster'])]
    with open(fname, 'wb') as outfile:
        # Should fix so we can have multiple clusters more easily
        json.dump(arr[0],
                  outfile,
                  sort_keys=True,
                  indent=4)


def iter_seqs(df):
    for row in df.itertuples():
        yield SeqRecord(Seq(row.seqs), id=row.unique_ids, description='')


def write_fasta(df, fname):
    print("writing {}".format(fname))
    with open(fname, 'w') as fh:
        SeqIO.write(iter_seqs(df), fh, "fasta")

    
def write_separate_fasta(df, output_dir, cluster_base):
    for k, g in df.groupby(['cluster']):
        fname = os.path.join(output_dir, cluster_base+'{}.fa'.format(k))
        write_fasta(g, fname)


def handle_frameshifts(df, remove_frameshifts=False):
    """Remove frameshited sequences from the data frame, if remove_frameshifts is set; otherwise, just warn.
    Note that a sequence could be considered frameshifted merely by having more than 1 stop codon, but this
    might not actually be due to a frameshift."""
    frameshifted_seqs = list(df[df.frameshifts].unique_ids)
    if frameshifted_seqs:
        if remove_frameshifts:
            df_ = df[- df.frameshifts]
            # This check is to make sure we don't end up with any clusters with fewer than 3 sequences, which
            # would stimy our attempts at making trees
            if len(df_) > 2:
                print "Removing frameshifted seqs:", frameshifted_seqs
                df = df_
            else:
                warnings.warn("Not removing frameshifted seqs from cluster to avoid < 3 issue: " +
                        str(frameshifted_seqs))
        else:
            warnings.warn("Found possible frameshifted input_seqs: " + str(frameshifted_seqs))
    return df

def seqmeta_path(output_dir, cluster_base, k):
    return os.path.join(output_dir, cluster_base+'{}.seqmeta.csv'.format(k))


def write_melted_partis(df, fname):
    """
    Spits out a csv with rows corresponding to sequences, and columns corresonding to data about seqs.
    Cassie and other Overbaugh group members sometimes prefer this read-level format. It's also useful in
    various parts of the cft pipeline.
    """
    df.to_csv(
        fname,
        index=False,
        columns=[
            'unique_ids', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length',
            'cluster', 'has_seed', 'frameshifts', 'duplicates', 'mut_freqs'])


def melted_path(output_dir, melted_base):
    return os.path.join(output_dir, melted_base + '.csv')


def main():
    """
    Run and save cluster file processing.
    """

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print("Inferring chain and gls from partis log")
    locus, inferred_gls = process_log_file(args.partis_log)

    melted_annotations = process_data(args.annotations,
                                      args.partition,
                                      locus,
                                      inferred_gls)

    write_melted_partis(melted_annotations,
                        melted_path(args.output_dir, args.melted_base))

    melted_annotations = handle_frameshifts(melted_annotations, args.remove_frameshifts)

    write_separate_fasta(melted_annotations,
                         args.output_dir,
                         args.cluster_base)

    write_json(melted_annotations,
               os.path.join(args.output_dir, 'partis_metadata.json'),
               os.path.getmtime(args.annotations),
               args.cluster_base,
               args.annotations,
               args.partition,
               args.output_dir,
               args.paths_relative_to)


if __name__ == '__main__':
    main()

