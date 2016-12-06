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
	$ conda install pandas
	$ conda install biopython
	$ conda install flask
	$ conda install nestly
	$ conda install pyqt
	$ pip install ete3
	$ pip install scons
	$ pip install flask-breadcrumbs
```

Your mileage may vary, but I have found it necessary to 
reactivate the environment and rehash my path, especially after
installing `scons` and `nestly`.

```
	$ source deactivate
	$ source activate cft
	$ hash -r
```

Finally, you'll have to have partis installed somewhere in order to execute `bin/process_partis.py` (via the `bin/demo.sh` script; see below).
If you have partis installed somewhere on your `PATH` already, you should be good to go.
If not, you can simply clone the partis repo and set the `PARTIS` environment variable so `process_partis.py` knows where to find things.

```
git clone --depth 1 git@github.com:psathryella/partis.git
export PARTIS=$PWD/partis
```

## Running

Running comes in three phases.

### Process partis output

First we execute the `process_partis.py` sciprt via `bin/demo.sh`.
This takes the output from partis and does some pre-processing on the data, splitting clusters up into separate repositories.

You'll need to specify an output directory for this intermediate data using the `outdir` environment variable.
Additionally, you'll likely need to specify the input `datapath` (there is a default, but there's a good chance it will be stale...).

```
datapath=path/to/yer/data/ outdir=processpartis-output ./bin/demo.sh
```

### Running scons

The rest of the build process will be executed using `scons`.
This builds trees out of each of the clusters, and does ancestral state reconstructions using these trees.

You can specify the `datapath` used as the outdir in the last step by using the `--datapath` flag.

```
scons --datapath=processpartis-output
```

Currently, all data is placed in an `output` subdirectory of your current working directory.

### Running the CFT web server

To actually consume the fruits of our labor, we need to fire up the CFT web server, which will provide us with a web interface for exploring the data.

Before exiting, the build will prompt you with the directions for how to execute the cftweb server given the output `metadata.json` file.
This should look something like:

```
cd cftweb && python -m cftweb --file /path/to/output/metadata.json
```

The default port is `5000`.
If someone else is running the web server on the same machine (or something else using that port), you can set a different one using the `-P` flag.


