
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

## SLURM support

You can run the same pipeline on slurm.
```
$ python pipeline.py --batchSystem slurm --disableCaching --logDebug ./jobStore ./sample_dedup.fa
```

## Installing Toil support for CWL

Support for CWL V1.0 is only available in code from the Toil
repository - not in any released pip or conda packages.

Furthermore, support for runnig CWL jobs on `slurm` is not stable.  I recommend downloading
Toil from a forked repo that I maintain with patches for slurm.

Your experience may vary, but I found I also needed to upgrade
`setuptools` in order to complete installation of the Toil dependencies.

To install from [the `dev/slurm` branch of] my forked repo,
```
$ conda create toil
$ source activate toil
$ module load nodejs
$ module load FastTree
$ pip install -e git+https://github.com/cswarth/toil@dev/slurm#egg=toil[cwl]                                                                                                                                                                                                                                                                                          
```

To install from the official Toil repo,
```
$ pip install git+https://github.com/BD2KGenomics/toil#egg=toil[cwl]
```

## Testing the Toil support

Copy the [first example][ex1] from the CWL User Guide into
`hello.cwl`.  Copy the YAML input specification into `hello.yml`.
Run the example with,
```
cwl_runner hello.cwl hello.yml
```

You may also use,
```
cwltoil hello.cwl hello.yml
```

## Running the CWL workflow on slurm


```
cwltoil --logDebug --batchSystem slurm --disableCaching --jobStore ./jobStore  --preserve-environment LD_LIBRARY_PATH PATH --no-container --workDir $PWD/workdir  --cleanWorkDir never workflow.cwl workflow.yml 
```


---

[ToilCWL]: http://toil.readthedocs.io/en/latest/running.html#running-cwl-workflows

[cwlref]: https://github.com/common-workflow-language

[cwltool]: https://github.com/common-workflow-language/cwltool

[ex1]: http://www.commonwl.org/v1.0/UserGuide.html#First_example


