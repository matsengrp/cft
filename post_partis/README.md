# Python modules for processing partis output

Run `bash scripts/demo.sh` for an example.
It will create the following directory structure:

```
_output/
└── seeds
    ├── QA255.006-Vh
    │   ├── Hs-LN2-5RACE-IgG-new-all-seqs.fa
    │   ├── Hs-LN2-5RACE-IgG-new-cluster0-seqs.fa
    │   ├── Hs-LN2-5RACE-IgG-new-cluster1-seqs.fa
    │   ├── Hs-LN2-5RACE-IgG-new-cluster3-seqs.fa
    │   ├── Hs-LN2-5RACE-IgG-new-cluster4-seqs.fa
    │   ├── Hs-LN2-5RACE-IgG-new-cluster5-seqs.fa
    │   ├── Hs-LN2-5RACE-IgG-new-cluster-stats.csv
    │   └── Hs-LN2-5RACE-IgG-new-seed_cluster2-seqs.fa
    ├── QB850.043-Vk
    │   ├── Hs-LN1-5RACE-IgK-100k-all-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster1-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster2-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster3-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster4-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster5-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster6-seqs.fa
    │   ├── Hs-LN1-5RACE-IgK-100k-cluster-stats.csv
    │   └── Hs-LN1-5RACE-IgK-100k-seed_cluster0-seqs.fa
    └── QB850.424-Vk
        ├── Hs-LN1-5RACE-IgK-100k-all-seqs.fa
        └── Hs-LN1-5RACE-IgK-100k-cluster-stats.csv
```

The output for seed `QA255.006-Vh` has 6 clusters, the seed appearing in `seed_cluster2` as well as all sequences in a single `.fa` file.

The output for seed `QB850.043-Vk` is the same as above with the seed cluster being cluster 0.
This seed also has some sequences of different length than others within separate clusters.

The output for seed `QB850.424-Vk` has all sequences in a single `.fa` file, also with some sequences of different length.
