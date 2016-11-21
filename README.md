# cft: clonal family tree

A pipeline for digesting partis clusters and inferring clonal family trees.

* Input: FASTA file, partis clustering file
* Intermediate step: tree building, which gets a per-clonal-family FASTA file, and an inferred naive sequence
* Output: Trees and mutation maps


## Setting up the environment

The pipeline is specified as python, but it relies upon some tools
that are available as modules.  Load these modules before you create
your virtual environment.

```
	$ module use ~matsengrp/modules
	$ module load phylip
	$ module load seqmagick
	$ module load FastTree
```

Now create a virtual environment.  These instructios assume you will
operating in a virtual environment using `conda`.  Use `module load
matsengrp-anaconda` if the `conda` command is not available in your
path.  You can also mimic these instructions using `virtualenv` and
`pip` if you would prefer.

```
	$ conda create -n cft
	$ source activate cft
```

Now load the python packages.  Again these are mostly speciofic using
`conda` commands, but you can do something similar with `pip`.

```
	$ conda install pandas
	$ conda install biopython
	$ conda install flask
	$ conda install -c etetoolkit ete3 ete3_external_apps
	$ pip install scons
	$ pip install nestly
```

Your mileage may vary, but I have found it necessary to 
reactivate the environment and rehash my path, especially after
installing `scons` and `nestly`.

```
	$ source deactivate
	$ source activate cft
	$ hash -r
```
