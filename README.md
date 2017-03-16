# cft: clonal family tree

A pipeline for digesting partis clusters and inferring clonal family trees.

[![wercker status](https://app.wercker.com/status/73265f18b3a63457ecbc79018da52162/s/master "wercker status")](https://app.wercker.com/project/byKey/73265f18b3a63457ecbc79018da52162)

* Input: FASTA file, partis clustering file
* Intermediate step: tree building, which gets a per-clonal-family FASTA file, and an inferred naive sequence
* Output: Trees, mutation maps, and a summary `metadata.json`

Additionally, the `cftweb` directory has a web application for serving
up the results of these analyses in your browser.


## Setting up the environment

The pipeline is specified as python, but it relies upon some tools
that are available as modules.  Load these modules before you create
your virtual environment.

```
	$ module use ~matsengrp/modules
	$ module load phylip seqmagick FastTree
```

Now create a virtual environment.  These instructions assume you will
operating in a virtual environment using `conda`.  Use `module load
matsengrp-anaconda` if the `conda` command is not available in your
path.  You can also mimic these instructions using `virtualenv` and
`pip` if you would prefer.

```
	$ conda create -n cft
	$ source activate cft
```

Now load the python packages.  Again these are mostly specified using
`conda` commands; `pip` is used if the package is not yet available
via conda.

```
	$ conda install pandas biopython flask nestly pyqt slackclient
	$ pip install ete3 scons flask-breadcrumbs slacker-log-handler cottonmouth colorbrewer
 
```

Your mileage may vary, but I have found it necessary to 
reactivate the environment and rehash my path, especially after
installing `scons` and `nestly`.

```
	$ source deactivate
	$ source activate cft
	$ hash -r
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

Running is a two step process:

* Run the analyses and produce a summary `metadata.json` file
* Run the `cftweb` application based on this file


### Running scons

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

* `--datapath`: The output directory of the partis run, defaulting (currently) to `/fh/fast/matsen_e/processed-data/partis/kate-qrs-2016-09-09/latest`.
* `--test`: Run on a small subset of all the seeds, as defined in the `SConstruct`, rather than the whole dataset; Useful for testing new code.
* `--asr-prog`: Should be set to either `dnaml` or `dnapars`, depending on whether you want to use parsimony or ML ancestral state reconstruction.

Typical example usage:

```
scons --datapath=/some/partis-run-2017-09-18/build7 --asr-prog=dnaml
```

This will create a metadata file and build data at `output/partis-run-build7-dnaml/metadata.json`

### Running the CFT web server

To actually consume the fruits of our labor, we need to fire up the CFT web server, which will provide us with a web interface for exploring the data.

Before exiting, the build will prompt you with the directions for how to execute the cftweb server given the output `metadata.json` file.
This should look something like:

```
cd cftweb && python -m cftweb /path/to/output/metadata.json
```

You may specify multiple json files in this fashion as long as each has a unique `dataset_id` attribute.

The default port is `5000`.
If someone else is running the web server on the same machine (or something else using that port), you can set a different one using the `-P` flag.

### Running CFT web in production

The default dev server shipped with Flask is rather stupid and can't handle multiple web requests at the same time, and also gets thrown for a loop if a socket connection closes before all the data has been sent, periodically crashing the app.
If you install `gevent` via pip or conda, the `gevent.wsgi` module is used instead, which should resolve these issues.
Eventually we'll have a more clever Docker based server setup, but for now this will keep us limping along.

Also note that for production, you may want to specify either the `--email` or `--slack` flags to the invocation above so that folks can be notified of errors.
For slack notifications (recommended), you will also need to obtain an API token (see <https://api.slack.com/web#authentication>), and set the `SLACK_TOKEN` environment variable accordingly.

Finally, the data output directories can be cut down to less than a 1/10 of their full size by executing the `./bin/publish_output.py` script.
Execute that script with the `-h` flag for further details on this.
This is recommended for moving data around in production/deployment.

