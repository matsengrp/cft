
Python scripts to implement a simplistic Toil pipeline that converts
FASTA files into tree images.  Runs the fasta files through Fasttree
to get a Newick format tree which is then run throughFogtree to
produce an SVG.

Usage:
```
$ python pipeline.py file:JobStore sample_dedup.fa sample_dedup.fa sample_dedup.fa sample_dedup.fa sample_dedup.fa sample_dedup.fa
```
This will process multiple copies of `sample_dedup.fa` as if they were different FASTA files.
As of this moment we are only running fast tree on the input files.

TODO:
* Add working figtree job
* Gather output onto GlobalFielStore.
* Base output name on the FASTA name to avoid name clashes.
* Use SLURM to run jobs in parallel
