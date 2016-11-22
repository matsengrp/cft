# cft: clonal family tree

A pipeline for digesting partis clusters and inferring clonal family trees.

* Input: FASTA file, partis clustering file
* Intermediate step: tree building, which gets a per-clonal-family FASTA file, and an inferred naive sequence
* Output: Trees, mutation maps, and a summary `metadata.json`

Additionally, the `cftweb` directory has a web application for serving up the results of these analyses as a web app.


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
	$ conda install -c etetoolkit ete3 ete3_external_apps
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

## Running

Build using `scons`.
You can specify the datapath to be built using the `--datapath` flag.

```
scons --datapath=path/to/yer/data/
```

Currently, all data is placed in an `output` subdirectory of your current working directory.

Before exiting, the build will prompt you with the directions for how to execute the cftweb server given the output `metadata.json` file.
This should look something like

```
cd cftweb && python -m cftweb -d output/metadata.json
```

