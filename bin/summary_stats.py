#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
summary stats of the seq data going into tree inference
'''

from __future__ import division, print_function
from Bio import AlignIO
import scipy, matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import argparse
import seaborn as sns
sns.set(style="white", color_codes=True)
from matplotlib.backends.backend_pdf import PdfPages

def hamming_distance(seq1, seq2):
    '''Hamming distance between two sequences of equal length'''
    return sum(x != y for x, y in zip(seq1, seq2))

parser = argparse.ArgumentParser(description='summary statistics of pre-tree data')
parser.add_argument('input', type=str, nargs='+', help='inseqs_trimmed.fa files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

for i, fname in enumerate(args.input):
    print(fname)
    seqs = {seq.id:str(seq.seq) for seq in AlignIO.read(fname, 'fasta')}
    nseqs = len(seqs)
    print(nseqs)
    if nseqs <= 2: continue

    distance_from_naive, degree = zip(*[(hamming_distance(seqs[seqid], seqs['naive0']),
                                         min(hamming_distance(seqs[seqid], seqs[seqid2]) for seqid2 in seqs if seqid2 != 'naive0' and seqid2 != seqid))
                                                   for seqid in seqs if seqid != 'naive0'])
    df = pd.DataFrame({'distance to naive sequence':distance_from_naive,
                       'nearest neighbor distance':degree})
    df['data set'] = i + 1
    if i == 0:
        aggdat = df
    else:
        aggdat = aggdat.append(df, ignore_index=True)

datasets = set(aggdat['data set'])
ndatasets = len(datasets)

# bw = .3
alpha = min([.9, 20/ndatasets])
bins = range(aggdat['distance to naive sequence'].max() + 2)

plt.figure(figsize=(6, 3))
plt.subplot(1, 2, 1)
for dataset, dataset_aggdat in aggdat.groupby('data set'):
    sns.distplot(dataset_aggdat['distance to naive sequence'],
                 bins=bins,
                 kde=False,
                 hist_kws={'histtype':'step', 'cumulative':True, 'alpha':alpha, 'lw':1})
plt.xlabel('distance to naive sequence')
plt.xlim([0, bins[-1]])
plt.ylabel('observed sequences')
plt.tight_layout()

bins = range(aggdat['nearest neighbor distance'].max() + 2)

plt.subplot(1, 2, 2)
for dataset, dataset_aggdat in aggdat.groupby('data set'):
    sns.distplot(dataset_aggdat['nearest neighbor distance'],
                 bins=bins,
                 kde=False,
                 hist_kws={'histtype':'step', 'cumulative':True, 'alpha':alpha, 'lw':1})
plt.xlabel('nearest neighbor distance')
plt.xlim([0, bins[-1]])
plt.ylabel('')
plt.tight_layout()

plt.savefig(args.outbase+'.pdf')
