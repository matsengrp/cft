
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

## SLURM support

You can run the same pipeline on slurm.
```
$ python pipeline.py --batchSystem slurm --disableCaching --logDebug ./jobStore ./sample_dedup.fa
```
## CWL Support

The Toil documentation [indicates support][ToilCWL] for the Common
Workflow Language.

* There are a [group of github repos][cwlref] containing the CWL
  specification and reference implementations of CWL parsers and
  tools.
* The Toil implementation of CWL depends upon [cwltool][cwltool] from the CWL
refrence project.
* CWL went through a revision process this Summer and is now at V1.0
* Virtually all the CWL workflow files I have discovered so far target
  CWL V1.0.  Parsers for earlier draft versions of CWL will not work on
  v1.0 files.
* There may be a `cwltool` commmand already available on your system, but it
  is almost certainly incompatible with Toil.
* The released Toil package on CRAN, available via `pip install toil`,
  does not have support for CWL V1.0
* To get Toil support for the CWL V1.0 specification, you must
  install the Toil package directly from sources. 
* CWL in Toil also requires an extra package that must be installed
  seperately from the main Toil support.
* I have not found any way to install Toil with `conda` .

## Installing Toil support for CWL

```
$ pip install --upgrade setuptools                                                                                                                                                                                                                                                                                                                                    
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

```
$ cwltoil hello.cwl hello.yml 
stoat 2016-10-25 09:40:47,363 MainThread INFO toil.lib.bioio: Root logger is at level 'INFO', 'toil' logger at level 'INFO'.
stoat 2016-10-25 09:40:48,699 MainThread INFO toil.lib.bioio: Root logger is at level 'INFO', 'toil' logger at level 'INFO'.
stoat 2016-10-25 09:40:48,721 MainThread INFO toil.jobStores.fileJobStore: Path to job store directory is '/tmp/tmpXIeNTX'.
stoat 2016-10-25 09:40:48,722 MainThread INFO toil.jobStores.abstractJobStore: The workflow ID is: '8286f2fd-fd46-4f8f-85c5-60d83ad86c61'
stoat 2016-10-25 09:40:48,754 MainThread INFO toil.common: Using the single machine batch system
stoat 2016-10-25 09:40:48,755 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxCores to CPU count of system (16).
stoat 2016-10-25 09:40:48,755 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxMemory to physically available memory (67378266112).
stoat 2016-10-25 09:40:48,755 MainThread INFO toil.common: Created the workflow directory at /tmp/toil-8286f2fd-fd46-4f8f-85c5-60d83ad86c61
stoat 2016-10-25 09:40:48,756 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxDisk to physically available disk (4628402176).
stoat 2016-10-25 09:40:48,756 MainThread INFO toil.batchSystems.singleMachine: Setting up the thread pool with 160 workers, given a minimum CPU fraction of 0.100000 and a maximum CPU value of 16.
stoat 2016-10-25 09:40:48,832 MainThread INFO toil.common: Written the environment for the jobs to the environment file
stoat 2016-10-25 09:40:48,832 MainThread INFO toil.common: Caching all jobs in job store
stoat 2016-10-25 09:40:48,833 MainThread INFO toil.common: 0 jobs downloaded.
stoat 2016-10-25 09:40:48,836 MainThread INFO toil.realtimeLogger: Real-time logging disabled
stoat 2016-10-25 09:40:48,858 MainThread INFO toil.leader: (Re)building internal scheduler state
stoat 2016-10-25 09:40:48,859 MainThread INFO toil.leader: Checked batch system has no running jobs and no updated jobs
stoat 2016-10-25 09:40:48,859 MainThread INFO toil.leader: Found 1 jobs to start and 0 jobs with successors to run
stoat 2016-10-25 09:40:48,859 MainThread INFO toil.leader: Starting the main loop
stoat 2016-10-25 09:40:48,861 Thread-1 INFO toil.batchSystems.singleMachine: Executing command: '/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/workflows/workflows/csw/venv/bin/_toil_worker file:/tmp/tmpXIeNTX V/9/jobXo1Kv7'.
stoat 2016-10-25 09:40:49,878 MainThread INFO toil.leader: No jobs left to run so exiting.
stoat 2016-10-25 09:40:49,879 MainThread INFO toil.leader: Finished the main loop
stoat 2016-10-25 09:40:49,879 MainThread INFO toil.leader: Waiting for stats and logging collator thread to finish ...
stoat 2016-10-25 09:40:50,364 MainThread INFO toil.leader: ... finished collating stats and logs. Took 0.484457969666 seconds
stoat 2016-10-25 09:40:50,364 MainThread INFO toil.leader: Waiting for service manager thread to finish ...
stoat 2016-10-25 09:40:50,860 MainThread INFO toil.leader: ... finished shutting down the service manager. Took 0.495107889175 seconds
stoat 2016-10-25 09:40:50,860 MainThread INFO toil.leader: Finished toil run successfully
stoat 2016-10-25 09:40:50,904 MainThread INFO toil.common: Attempting to delete the job store
stoat 2016-10-25 09:40:50,906 MainThread INFO toil.common: Successfully deleted the job store
```


[ToilCWL]: http://toil.readthedocs.io/en/latest/running.html#running-cwl-workflows

[cwlref]: https://github.com/common-workflow-language

[cwltool]: https://github.com/common-workflow-language/cwltool

[ex1]: http://www.commonwl.org/v1.0/UserGuide.html#First_example


