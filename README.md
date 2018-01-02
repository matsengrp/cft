# CFT: Clonal Family Tree

A pipeline for producing clonal family trees and ancestral state reconstructions using partis output.

Output data can be run through [cftweb](https://github.com/matsengrp/cftweb) for visualization and exploration.


## Input data

For each sample you'd like to process, cft needs to know:

* `partition-file` - the main partition output file from partis
* `cluster-annotation-file` - the cluster annotation file output from partis
* `locus` - as used in running partis
* `parameter-dir` - as used in running partis
* `per-sequence-meta-file` (optional) - as applicable, with noted columns `unique_id,timepoint,multiplicity`

CFT requires that you organize this information in a _dataset file_ as follows:

```yaml
# The dataset-id identifies the collection of samples in downstream organization
id: laura-mb-v14
samples:
  # each sample must be keyed by a per-dataset unique identifier
  Hs-LN-D-5RACE-IgG:
    locus: igh
    parameter-dir: /path/to/Hs-LN-D-5RACE-IgG/parameter-dir
    per-sequence-meta-file: /path/to/Hs-LN-D-5RACE-IgG/seqmeta.csv

    # Unseeded partitions go here
    partition-file: /path/to/Hs-LN-D-5RACE-IgG/partition.csv
    cluster-annotation-file: /path/to/Hs-LN-D-5RACE-IgG/partition-cluster-annotations.csv

    # seed partition runs should be organized under a `seeds` key as follows
    seeds:
      # seed sequence id
      BF520.1-igh:
        partition-file: /path/to/Hs-LN-D-5RACE-IgG/partition.csv
        cluster-annotation-file: /path/to/Hs-LN-D-5RACE-IgG/partition-cluster-annotations.csv
      # other seeds, as applicable...

  # another sample in our dataset...
  Hs-LN-D-5RACE-IgK:
    # etc.

  # etc.
```

Some notes about this:

* This file can either be in `.yaml` format (shown above), or `.json`.
* You may specify a `meta` attribute within a particular sample with keys from `[isotype, locus, shorthand, species, subject, timepoint]`.
* If a sample has multiple unseeded partis runs, these can be nested within an `other-partitions` key, much as `seeds` are specified

A more fleshed out example, as well as json examples, and python snippets can be seen [on the wiki](https://github.com/matsengrp/cft/wiki).


## Running the pipeline

_Note: Before you run the pipeline, you must follow the build environment setup section below._

CFT uses the `scons` build tool to execute the build process.
Running `scons` from within the `cft` checkout directory loads the `SConstruct` file, which specifies how data is to be processed:

* initial processing of partis results using `bin/process_partis.py` to produce (among other things) a sequence file for each cluster/clonal-family
    * a quality filter is applied removing sequences with frameshifts, mutations in "invariant" regions, or stop codons
* build trees out of the sequences for each of these clusters
* subset the sequences according to a couple of different strategies (seed lineage selection and overall diversity selection)
* ancestral state reconstructions using this sequence subset, producing:
    * a tree file
    * an svg representation of the tree, with tips colored by timepoints and a highlighted seed lineage
    * a fasta file with sequences corresponding to internal nodes on the tree (the ancestral state reconstructions)
* finally, produce a `output/<dataset-id>/metadata.json` file consumable by the `cftweb` web application summarizing this information

This particular `SConstruct` takes several command line parameters.
Below are the most frequently used options:

* `--infiles`: A `:` separated list of partis output directories relative to `--base-datapath`
  Defaults (presently) to `laura-mb/latest:kate-qrs/latest`.
* `--base-datapath`: The location of `--datapaths`, if not specified as absolute paths.
  Defaults to `/fh/fast/matsen_e/processed-data/partis/`.
* `--test`: Run on a small subset of all the seeds, as defined in the `SConstruct`, rather than the whole dataset; Useful for testing new code.

For the most complete and up to date reference on these, look at the tail `Local Options` section of `scons -h`.

You may also wish to take note of the following basic `scons` options:

* `-j`: Specify the number of jobs (parallelism) for the build
* `-k`: keep building as much as possible if there is an error
* `-n`: perform a "dry run" of the pipeline, only printing out the commands that would be run without actually running any

A separate "dataset" directory and corresponding `metadata.json` file will be created for each `infile` and placed within the `output` directory, organized by the `id` attribute of the dataset infile.

### Typical example usage

```
# If you're using conda, as above, first activate the environment
source activate cft

# Next build the 
scons --infiles=info1.yaml:info2.yaml -j 12

# Check out output
tree output
# or if you don't have tree
find output
```

Note that you can install tree with `sudo apt-get install tree` on Ubuntu for a nice ASCII-art file tree display of the output contents.


## Setting up the build environment

First, you'll need a number bioinformatics executables installed in order to run the pipeline:

* `FastTree`
* `phylip`
* `seqmagick`
* `muscle`
* `rppr` (from the `pplacer` suite)
* `xvfb-run`
* `R` (and `Rscript`)

If you are running on Fredhutch's compute nodes, you can use scicomp's module system to gain access to these dependencies:

```
module use ~matsengrp/modules
module load phylip seqmagick FastTree
```

Otherwise, you'll need to install each as you wish on whatever system you run.

### Python dependencies

We'll also need a number of python libraries.
The instructions below illustrate how one might set these dependencies up via `conda`.
One should however be able to complete the setup with `virtualenv` or whatever else as well.

First, you'll need to have `conda` installed.
If you are at the hutch you can load conda with `module load matsengrp-anaconda`.
Otherwise, Google is your friend.

Next, initialize and activate the conda environment:

```
conda create -n cft
source activate cft
```

Now install the python packages.
These are mostly installed via the `conda` command; `pip` is only used for packages not yet available via conda.

```
conda install biopython nestly pyqt
pip install ete3 scons
```

Your mileage may vary, but Chris Warth found it necessary to reactivate the environment and rehash the path, especially after installing `scons` and `nestly`.

```
source deactivate
source activate cft
hash -r
```

### Git submodules

Finally, there is some python code needed for the build script to execute which can be found in a number of git submodules.
In particular, this repository has a partis submodule which should be kept in sync to avoid build issues.
To check out these submodules, execute `git submodule init` then `git submodule update`.

If there are updates to the submodules, you can have those reflected in your checkouts by executing `git submodule update`.
If you want to update the version/commit pointed to by a submodule, you can `cd <submodule> && git update-repo-as-desired && cd .. && git add <submodule> && git commit -m "Your commit message here"`.

Note that if you'd like to use a different partis installation you can do so using the `PARTIS` env variable.

### Running partis

At the moment, this part of the pipeline doesn't require running partis at all if you only need to operate on data output by partis.
However, it might be worth compiling so you can use (at the very least) partis' `view-annotations` and/or `view-partitions` subcommands for inspecting partis' output files.
See `$PARTIS/README.md` for instructions on this.

### Using Slurm

The build pipeline is set up to use `slurm` for job submission on a number of the more compute heavy, long-running tasks.
If you have a slurm environment set up to submit to a cluster, and are able to write from slurm nodes to a shared filesystem, you can potentially run with significantly higher parallelism than you would be able to on a single computer.

If you are running on Fredhutch's servers, this should all be set up for you, and you should be able to submit upwards of 50-70 jobs using the `-j` flag, as specified below.
If you're not at the Hutch, setting up such a cluster is _way_ out of scope for this document, but if you're inspired, good luck figuring it out!


## Executing CFTWeb

Once the data is built, you can consume the fruits of this labor by passing the data off to the cftweb application.
For instructions on doing this please see the [cftweb](https://github.com/matsengrp/cftweb) repository.


