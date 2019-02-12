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

    # seed partition runs should be organized under a `seeds` key as follows
    seeds:
      # seed sequence id
      BF520.1-igh:
        partition-file: /path/to/Hs-LN-D-5RACE-IgG/partition.csv
      # other seeds, as applicable...

  # another sample in our dataset...
  Hs-LN-D-5RACE-IgK:
    # etc.

  # etc.
```

### Some notes about this:

* This file can either be in `.yaml` format (shown above), or `.json`.
* You may specify a `meta` attribute within a particular sample with keys from `[isotype, locus, shorthand, species, subject, timepoint]`.
* If a sample has multiple unseeded partis runs, these can be nested within an `other-partitions` key, much as `seeds` are specified
* A functional example is found in `test.yaml`, which is run by default by scons if --infiles is not set

A more fleshed out example, as well as json examples, and python snippets can be seen [on the wiki](https://github.com/matsengrp/cft/wiki).

You may also wish to take a look at `bin/dataset_utils.py`, a little utility script for filtering and merging dataset files.
You can get a comprehensive help menu by running `bin/dataset_utils.py -h`.
You may also wish to directly use the script which does the initial extraction of data from partis, `bin/process_partis.py`.

Note that in order for the data to process correctly, the following must be true of the naming scheme for sequences:

* must not include any of the characters: `:`, `;`, `,`

One final note: Partis now by default outputs `yaml` files which contain all the data necessary to run this pipeline.
However, it was previously the case that partis output separate `csv` files for partition and annotation data, and that both of these had to be passed into the `yaml` dataset files sketched out above.
If you are running cft on this old `csv` data, things should work fine as long as you don't rename or move these partition/annotation files in relation to each other.
As long as you keep the naming path relationships the same, cft will infer the location of the corresponding annotation csv files based on the specified partition file paths.


## Running the pipeline

_Note: Before you run the pipeline, you must follow the build environment setup section below._

CFT uses the `scons` build tool to execute the build process.
Running `scons` from within the `cft` checkout directory loads the `SConstruct` file, which specifies how data is to be processed:

* initial processing of partis results using `bin/process_partis.py` to produce (among other things) a sequence file for each cluster/clonal-family
    * a quality filter is applied removing sequences with stop codons, with out of frame CDR3 regions, and with mutations in the invariant codons bounding the CDR3
* build trees out of the sequences for each of these clusters
* subset the sequences according to a couple of different strategies (seed lineage selection and overall diversity selection)
* ancestral state reconstructions using this sequence subset, producing:
    * a tree file
    * an svg representation of the tree, with tips colored by timepoints and a highlighted seed lineage
    * a fasta file with sequences corresponding to internal nodes on the tree (the ancestral state reconstructions)
* finally, produce a `output/<dataset-id>/metadata.json` file consumable by the `cftweb` web application summarizing this information

Running `scons` without modifying the `SConstruct` will run default tests on the partis output in `tests/`.
To check that the output thereby produced matches the expected test output, run `diff -ubr --exclude='*metadata.json' tests/test-output output`

This particular `SConstruct` takes several command line parameters.
Below are the most frequently used options, which _must_ include `=` in the format `--option=value`:

* `--infiles`: A `:` separated list of partis output directories relative to `--base-datapath`
* `--base-datapath`: The location of `--datapaths`, if not specified as absolute paths.
* `--test`: Run on a small subset of all the seeds, as defined in the `SConstruct`, rather than the whole dataset; Useful for testing new code.

A separate "dataset" directory and corresponding `metadata.json` file will be created for each `infile` and placed within the `output` directory, organized by the `id` attribute of the dataset infile.
For the most complete and up to date reference on these, look at the tail `Local Options` section of `scons -h`.

You may also wish to take note of the following basic `scons` build options options:

* `-j`: specify the number of jobs (parallelism) for the build
* `-k`: if there is an error, stop building targets downstream of failure, but continue to build all targets not downstream of such errors
* `-i`: if there is an error, continue running all jobs, including those downstream of failure
* `-n`: perform a "dry run" of the pipeline, only printing out the commands that would be run without actually running any
* `--debug explain`: scons print out why it's building each target (e.g. hasn't been built yet, updated command, changed upstream target, updated executable/script, etc.), useful to have in the logs for debugging 

In general, it's good to run with `-k` so that on a first pass, you end up building as much of the data as you can properly build.
If there are errors, try rerunning to make sure the problem isn't just an errant memory issue on your cluster, then look back at the logs and see if you can't debug the issue.
If it's just a few clusters failing to build properly and you don't want to hold out on getting the rest of the built data into `cftweb`, you can rerun the build with `-i`, which will take a little longer to run through all of the failed build branches with missing files etc, but which should successfully compile the final output `metadata.json` files necessary for passing along to `cftweb`.


### Typical example usage

```
# If you're using conda, as below, first activate the environment
source activate cft

# Build the data, running 12 jobs at a time (parallelism) and appending all stdout/stderr to a log file
scons --infiles=info1.yaml:info2.yaml -k -j 12 --debug explain &>> 2018-05-24.info1-build.log

# You can watch a live tail of the log file from another terminal window or tmux pane with
tail -f 2018-05-24.info1-build.log

# Once its done running, you can take a look at the output
tree output
# or if you don't have tree
find output
```

Note that you can install tree with `sudo apt-get install tree` on Ubuntu for a nice ASCII-art file tree display of the output contents.


## Setting up the build environment

First, you'll need a number bioinformatics executables installed in order to run the pipeline:

* `phylip`
* `seqmagick`
* `muscle`
* `rppr` (from the `pplacer` suite)
* `xvfb-run`
* `R` (and `Rscript`)
* `prank`

If you are running on Fredhutch's compute nodes, you can use scicomp's module system to gain access to these dependencies. Start by loading these modules, which should allow you most of what you might need. When in doubt, check if you can install it with conda before resorting to loading the module.

```
module load phylip seqmagick 
ml pplacer
ml unload python2/2.7.8 
ml unload GCCcore/5.4.0
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
conda install python=2.7
conda install biopython pyqt scons
conda install -c r r-rcolorbrewer 
conda install -c etetoolkit ete3
conda install -c conda-forge nestly
conda install -c conda-forge yaml
conda install --force scipy=0.17.0
conda install -c bioconda fasttree
conda install -c anaconda gcc
```

Some have found it necessary to add the anaconda-installed GNU compiler (GCC) to their paths similar to this:
```
export PATH=/miniconda2/envs/cft_min/lib/gcc:$PATH
```

Your mileage may vary, but some have found it necessary to reactivate the environment and rehash the path, especially after installing `scons` and `nestly`.

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


