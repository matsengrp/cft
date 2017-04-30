# cft: clonal family tree

A pipeline for digesting partis clusters and inferring clonal family trees.

[![wercker status](https://app.wercker.com/status/73265f18b3a63457ecbc79018da52162/s/master "wercker status")](https://app.wercker.com/project/byKey/73265f18b3a63457ecbc79018da52162)

* Input: FASTA file, partis clustering file
* Intermediate step: tree building, which gets a per-clonal-family FASTA file, and an inferred naive sequence
* Output: Trees, mutation maps, and a summary `metadata.json`

Additionally, the `cftweb` directory has a web application for serving
up the results of these analyses in your browser.


## Setting up the build environment

This section is specific to setting up the data build environment.
For running the CFTWeb server, see the [CFTWeb README]('/cftweb/README.md').

First, the pipeline will need to have access to various commands for building the data (bioinformatics utilities and such).
You can use the Hutch's module system for these:

```
module use ~matsengrp/modules
module load phylip seqmagick FastTree
```

Next we'll create a virtual environment using conda for all of this projects python dependencies.
Use `module load matsengrp-anaconda` if the `conda` command is not available in your path. 
You can also mimic these instructions using `virtualenv` and`pip` if you would prefer.

```
conda create -n cft
source activate cft
```

Now load the python packages.
Again these are mostly specified using `conda` commands; `pip` is used if the package is not yet available via conda.

```
conda install pandas biopython nestly pyqt
pip install ete3 scons
```

Your mileage may vary, but Chris Warth found it necessary to reactivate the environment and rehash my path, especially after installing `scons` and `nestly`.

```
source deactivate
source activate cft
hash -r
```

Finally, you'll have to have partis installed somewhere in order to execute `bin/process_partis.py`.
This repository has a partis submodule that should be kept in sync so that bugs as a result of mismatching assumptions should be avoided.
There's also a `datascripts` submodule which contains code that Duncan uses to organize data for execution of partis, some of which is useful to this application (cft).
To check out these submodules, execute `git submodule init`.
If there are updates to the submodules, you can have those reflected in your checkouts by executing `git submodule update`.
If you want to update the version/commit pointed to by a submodule, you can `cd <submodule> && git update-repo-as-desired && cd .. && git add <submodule> && git commit -m "Your commit message here"`.

Note that if you'd like to use a different partis installation you can do so using the `PARTIS` env variable.

### Running partis

At the moment, this part of the pipeline doesn't require running partis at all.
However, it might be worth compiling so you can use (at the very least) partis' `view-annotations` and/or `view-partitions` subcommands for inspecting partis' output files.
See `$PARTIS/README.md` for instructions on this.


## Running

The build process can be initiated by executing `scons`.
This:

* does some initial processing of partis results using `process_partis.py` to produce (among other things) a sequence file for each cluster/clonal-family
* builds trees out of the sequences for each of these clusters
* subsets the sequences to sequences along the seed sequence lineage (this may change in the future)
* does ancestral state reconstructions using this sequence subset, producing:
    * a tree file
    * an svg representation of the tree, with tips colored by timepoints and a highlighted seed lineage
    * a fasta file with sequences corresponding to internal nodes on the tree (the ancestral state reconstructions)
* finally, produces a `metadata.json` file consumable by the `cftweb` web application summarizing this information

This build script takes several parameters.
For the most complete and up to date reference on these, look at the tail `Local Options` section of `scons -h`.
Below are the most frequently used options:

* `--datapaths`: A `:` separated list of partis output directories relative to `--base-datapath`
  Defaults (presently) to `laura-mb/latest:kate-qrs/latest`.
* `--base-datapath`: The location of `--datapaths`, if not specified as absolute paths.
  Defaults to `/fh/fast/matsen_e/processed-data/partis/`.
* `--asr-progs`: Should be set to either `dnaml` or `dnapars` or `dnaml:dnapars`, depending on whether you want to use parsimony or ML ancestral state reconstruction.
* `--test`: Run on a small subset of all the seeds, as defined in the `SConstruct`, rather than the whole dataset; Useful for testing new code.

A separate "dataset" directory and corresponding `metadata.json` file will be created for each combination of `datapaths` and `asr-progs`, and placed within the `output` directory.
Later, additional dataset level parameters may be added that will further expand this space.

After the data has been built, stripped down extracts of the data are placed in `output/bulids/`.
These are 1/10 the size of the full data directories, and are recommended for passing off data to CFTWeb (see below).

### Typical example usage

```
scons --datapaths=kate-qrs/v10 --asr-progs=dnaml -j 50
```

This will create a metadata file and build data at `output/kate-qrs-v10-dnaml/metadata.json`

Note that the `-j` flag can set the parallelism of the build.
Additionally, the `-n` flag can be used to execute a "dry run" where commands are printed but not executed.


### Executing CFTWeb

Once the data is built, you can consume the fruits of this labor by passing the data off the cftweb application.
For instructions `cd cftweb && cat README.md`.


