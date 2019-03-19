#!/usr/bin/env python
from __future__ import division
import argparse
import json
import csv
from tripl import tripl
from ete3 import PhyloTree

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', nargs='+')
    parser.add_argument('-C', '--csv', action="store_true")
    parser.add_argument('-c', '--clonal-families-out')
    parser.add_argument('-d', '--datasets-out')
    parser.add_argument('-s', '--sequences-out')
    return parser.parse_args()

# Some generic data processing helpers helpers

def comp(f, g):
    def h(*args, **kw_args):
        return f(g(*args, **kw_args))
    return h

def strip_ns(a):
    return a.split(':')[-1]

inf = float("inf")
neginf = float("-inf")

def clean_record(d):
    if isinstance(d, list):
        return map(clean_record, d)
    elif isinstance(d, dict):
        return {strip_ns(k): clean_record(v)
                for k, v in d.items()}
    # can't have infinity in json
    elif d == inf or d == neginf:
        return None
    else:
        return d


# Pulling datasets information out

datasets_pull_pattern = [
    "db:ident",
    "tripl:type",
    "cft.dataset:id",
    "cft.cluster:_dataset",
    "cft.subject:_dataset",
    "cft.timepoint:_dataset",
    {"cft.dataset:build": ["cft.build:id", "cft.build:time", "cft.build:commit"]}]


def clean_dataset_record(d):
    d = d.copy()
    d['n_clonal_families'] = len(d['cft.cluster:_dataset'])
    d['n_subjects'] = len(d['cft.subject:_dataset'])
    d['n_timepoints'] = len(d['cft.timepoint:_dataset'])
    del d['cft.cluster:_dataset'], d['cft.subject:_dataset'], d['cft.timepoint:_dataset']
    return d

def pull_datasets(t):
    records = list(t.pull_many(datasets_pull_pattern, {'tripl:type': 'cft.dataset'}))
    return map(comp(clean_record, clean_dataset_record), records)


# Pulling clonal families information out

clonal_family_pull_pattern = [
   {    
     "cft.reconstruction:seqmeta": [{"tripl.csv:data": ["bio.seq:id", "cft.seq:cluster_multiplicity", "cft.seq:multiplicity"]}],
     "cft.reconstruction:cluster_aa": [{"bio.seq:set": ["*"]}],
     "cft.reconstruction:asr_tree": ["*"],
     "cft.reconstruction:asr_seqs": [{'bio.seq:set': ['bio.seq:id', 'bio.seq:seq']}],
     "cft.reconstruction:cluster": 
     [
      "db:ident",
      "tripl:type",
      "cft.cluster:id",
      "cft.cluster:naive_seq",
      "cft.cluster:has_seed",
      "cft.cluster:n_unique_seqs",
      "cft.cluster:n_sampled_seqs",
      "cft.cluster:size",
      "cft.cluster:sorted_index",
      "cft.cluster:v_end",
      "cft.cluster:v_start",
      "cft.cluster:v_gene",
      "cft.cluster:d_end",
      "cft.cluster:d_gene",
      "cft.cluster:d_start",
      "cft.cluster:j_end",
      "cft.cluster:j_gene",
      "cft.cluster:j_start",
      "cft.cluster:cdr3_length",
      "cft.cluster:cdr3_start",
      "cft.cluster:naive_seq",
      "cft.cluster:mean_mut_freq",
      {"cft.cluster:sample": ["cft.sample:id", "cft.sample:timepoint"],
       "cft.cluster:dataset": ["cft.dataset:id"],
       "cft.cluster:partition": ["cft.partition:id", "cft.partition:logprob", "cft.partition:step"],
       "cft.cluster:subject": ["cft.subject:id"],
       "cft.cluster:v_per_gene_support": ["cft.gene_support:gene", "cft.gene_support:prob"],
       "cft.cluster:d_per_gene_support": ["cft.gene_support:gene", "cft.gene_support:prob"],
       "cft.cluster:j_per_gene_support": ["cft.gene_support:gene", "cft.gene_support:prob"],
      }
     ]
   } 
]

def create_node_records(tree, nt_seqs_dict, aa_seqs_dict, seqmeta_dict):
    records = []
    leaves_counter = 1
    for n in tree.traverse('postorder'):
        n.label = n.id = n.name
        n.nt_seq = nt_seqs_dict[n.name]
        n.aa_seq = aa_seqs_dict[n.name]
        mult = None
        cluster_mult = None
        if n.name in seqmeta_dict.keys():
            mult = seqmeta_dict[n.name]["cft.seq:multiplicity"].pop()
            clust_mult = seqmeta_dict[n.name]["cft.seq:cluster_multiplicity"].pop()
        n.multiplicity = int(mult) if mult else mult
        n.cluster_multiplicity = int(clust_mult) if clust_mult else clust_mult
        n.type = "node"
        if n.is_leaf():
            # get height for leaves
            n.type = "leaf"
            n.height = leaves_counter
            leaves_counter +=1
        else:
            # get height for non leaves
            total_height = 0
            for child in n.children:
                total_height += child.height
            avg_height = total_height/len(n.children)
            n.height = avg_height
        if n.up:
            # get parent info, distance for non root
            n.parent = n.up.name
            n.length = n.get_distance(n.up)
            n.distance = n.get_distance("naive")
        else:
            # n is root
            n.type = "root"
            n.parent = None
            n.length = 0.0
            n.distance = 0.0
        records.append({'id': n.id, 'label': n.label, 'type': n.type, 'parent': n.parent, 'length': n.length, 'distance': n.distance, 'height': n.height, 'nt_seq': n.nt_seq, 'aa_seq': n.aa_seq, 'multiplicity': n.multiplicity, 'cluster_multiplicity': n.cluster_multiplicity})

    return records

def parse_tree_data(s, nt_seqs_dict, aa_seqs_dict, seqmeta_dict):
    t = PhyloTree(s, format=1)
    return create_node_records(t, nt_seqs_dict, aa_seqs_dict, seqmeta_dict)

def create_seqs_dict(seq_records):
    d = dict()
    for record in seq_records:
        d[record["bio.seq:id"].pop()] = record["bio.seq:seq"].pop()
    return d

def create_seqmeta_dict(seqmeta_records):
    d = dict()
    for record in seqmeta_records:
        seq_id = record["bio.seq:id"].pop()
        d[seq_id] = record
    return d

def clean_clonal_family_record(d):
    c = d.copy()
    c['cft.reconstruction:cluster'] = c['cft.reconstruction:cluster'][0]
    #c['cft.reconstruction:cluster']['cft.reconstruction:seqmeta'] = c['cft.reconstruction:seqmeta']
    aa_seqs_dict = create_seqs_dict(list(c['cft.reconstruction:cluster_aa'])[0]['bio.seq:set'])
    nt_seqs_dict = create_seqs_dict(list(c['cft.reconstruction:asr_seqs'])[0]['bio.seq:set'])
    seqmeta_dict = create_seqmeta_dict(list(c['cft.reconstruction:seqmeta'])[0]['tripl.csv:data'])
    if(c['cft.reconstruction:asr_tree'][0].get('tripl.file:contents')):
        c['cft.reconstruction:cluster']['cft.reconstruction:asr_tree'] = parse_tree_data(list(c['cft.reconstruction:asr_tree'][0]['tripl.file:contents'])[0], nt_seqs_dict, aa_seqs_dict, seqmeta_dict)
    c = c['cft.reconstruction:cluster']
    try:
        del c['cft.cluster:unique_ids']
    except Exception:
        pass    
    return c

def pull_clonal_families(t):
    result = map(comp(clean_record, clean_clonal_family_record),
            t.pull_many(clonal_family_pull_pattern, {'tripl:type': 'cft.reconstruction'}))
    #result[0]['cft.reconstruction:cluster']['cft.reconstruction:asr_tree'] = parse_tree_data(list(result['cft.reconstruction:cluster']['cft.reconstruction:asr_tree']['tripl.file:contents'])[0])
    return result

def write_out(data, filename, args):
    with open(filename, 'w') as fh:
        if args.csv:
            data = [{k: v for k, v in d.items()}
                    for d in data]
            writer = csv.DictWriter(fh, fieldnames=sorted(data[0].keys()))
            writer.writeheader()
            writer.writerows(data)
        else:
            # Then assume json
            json.dump(data, fh, default=list,
                indent=4,
                #allow_nan=False
                )

def main():
    args = get_args()
    t = tripl.TripleStore.loads(args.inputs)
    if args.datasets_out:
        write_out(pull_datasets(t), args.datasets_out, args)
    if args.clonal_families_out:
        write_out(pull_clonal_families(t), args.clonal_families_out, args)
    #if args.sequences_out:
        #write_out(pull_sequences(t), args.sequences_out, args)


if __name__ == '__main__':
    main()


