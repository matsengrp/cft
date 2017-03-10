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
    parser.add_argument(
        '--paths-relative-to',
        help='files pointed to from metadata.json file will be speicfied relative to this path; defaults to --output-dir')
    parser.add_argument(
        '--melted_base',
        help='basename for melted data output',
        default='melted')
    parser.add_argument(
        '-F', '--remove-frameshifts',
        help='if set, tries to remove seqs with frameshift indels from the output',
        action="store_true")
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

    # default paths_relative_to is just whatever the output dir is
    args = parser.parse_args()
    args.paths_relative_to = args.paths_relative_to or args.output_dir

    return parser.parse_args()


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


def indel_offset(indelfo):
    return sum(map(lambda indel: indel['len'] * (1 if indel['type'] == 'insertion' else -1),
                   indelfo['indels']))

def infer_frameshifts(line):
    attrs = ['stops', 'indelfos', 'input_seqs']
    def infer_(args):
        stop, indelfo, input_seq = args
        aa = Seq(input_seq).translate()
        stop_count = aa.count("*")
        # We say it's a frameshift if the indell offsets don't leave us with a multiple of three, and if there
        # are stop codons. Can tweak this down the road, but for now...
        return bool(stop_count > 0 and indel_offset(indelfo) % 3)
    return map(infer_,
               zip(*map(lambda x: line[x], attrs)))


def process_data(annot_file, part_file, locus, glpath):
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
    glfo = glutils.read_glfo(glpath, locus)
    to_keep = ['v_gene', 'd_gene', 'j_gene', 'cdr3_length']
    for idx, cluster in annotations.fillna('').iterrows():
        current_df = pd.DataFrame()
        line = cluster.to_dict()
        utils.process_input_line(line)
        utils.add_implicit_info(glfo, line)
        current_df['unique_ids'] = line['unique_ids'] + ['naive{}'.format(idx)]
        current_df['seqs'] = line['input_seqs'] + [line['naive_seq']]
        current_df['frameshifts'] = infer_frameshifts(line) + [False] # mocking last entry
        current_df['mut_freqs'] = line['mut_freqs'] + [0.0] # mocking last entry
        for col in to_keep:
            current_df[col] = line[col]
        current_df['cluster'] = str(idx)
        current_df['has_seed'] = any(seed_id in line['unique_ids'] for seed_id in \
                seed_ids)
        current_df['seed_ids'] = ':'.join(seed_ids)
        current_df['cdr3_start'] = line['codon_positions']['v']
        for gene in 'vdj':
            for pos in ['start', 'end']:
                current_df[gene+'_'+pos] = line['regional_bounds'][gene][pos.startswith('e')]
        output_df = pd.concat([output_df, current_df])

    return output_df


def write_json(df, fname, mod_date, cluster_base, annotations, partition, meta, outdir, paths_relative_to):
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
    # of the form "/path/to/output/QA255.016-Vh/Hs-LN2-5RACE-IgG-new/*.fa"
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
            'file': os.path.relpath(os.path.join(outdir, cluster_base+cluster_id+'.fa'), paths_relative_to),
            'cluster_id': cluster_id,
            'n_seqs': len(df), # Note... this count naive; good idea?
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
    # Remove frameshift seqs if requested
    for k, g in df.groupby(['cluster']):
        fname = os.path.join(output_dir, cluster_base+'{}.fa'.format(k))
        write_fasta(g, fname)


def handle_frameshifts(df, remove_frameshifts=False):
    frameshifted_seqs = list(df[df.frameshifts].unique_ids)
    if frameshifted_seqs:
        if remove_frameshifts:
            print "Removing frameshifted seqs:", frameshifted_seqs
            df = df[- df.frameshifts]
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
            'cluster', 'has_seed', 'frameshifts', 'mut_freqs'])


def melted_path(output_dir, melted_base):
    return os.path.join(output_dir, melted_base + '.csv')


def main():
    """
    Run and save cluster file processing.
    """

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # get metadata
    regex = re.compile(r'^(?P<subject_id>[^.]*).(?P<seed_id>[0-9]*)-(?P<gene>[^/]*)/[^-]*-(?P<timepoint>[^-]*)')
    path = '/'.join(args.output_dir.split('/')[-3:])
    m = regex.match(path)
    meta = {}
    if m:
        meta = m.groupdict()
    else:
        raise Exception('--output_dir needs to be of the form "/path/to/output/QA255.016-Vh/Hs-LN2-5RACE-IgG"')

    if args.partis_log is not None:
        print("Inferring chain and gls from partis log")
        chain, inferred_gls = process_log_file(args.partis_log)
    else:
        print("Inferring chain and gls from path and param dir cl arg")
        chain = meta['gene'][1].lower()
        inferred_gls = args.param_dir + '/hmm/germline-sets'

    locus = dict(h='igh', k='igk', l='igl', a='tra', b='trb', d='trd', g='trg')[chain]

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
               os.path.join(args.output_dir, 'metadata.json'),
               os.path.getmtime(args.annotations),
               args.cluster_base,
               args.annotations,
               args.partition,
               meta,
               args.output_dir,
               args.paths_relative_to)


if __name__ == '__main__':
    main()
