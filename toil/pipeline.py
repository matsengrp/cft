#!/usr/bin/env python
"""
Demonstration of a simple Toil pipeline
to convert FASTA files into SVG tree file.

Usage:
  $ python pipeline.py file:JobStore sample_dedup.fa sample_dedup.fa sample_dedup.fa sample_dedup.fa sample_dedup.fa sample_dedup.fa

Please see the README.md in the same directory
Tree Structure of pipeline (per sample)

fasttree -> figtree -> cleanup
Dependencies
FastTree:       module load FastTree
Figtree:        module load matsengrp-local
Toil:           module load matsengrp-anaconda
				source activate toil

"""
from __future__ import print_function

import argparse
import base64
import errno
import glob
import hashlib
import multiprocessing
import os
import shutil
import subprocess
import logging

from toil.job import Job
from toil.common import Toil, jobStoreLocatorHelp, Config

log = logging.getLogger(__name__)


class Fasttree(Job):

    def __init__(self, inputid):
        Job.__init__(self,  memory="2G", cores=1, disk="100M")
        self.inputid = inputid

    def run(self, fileStore):
        fileStore.logToMaster("inside Fasttree run()")
        input_file = fileStore.readGlobalFile(
            self.inputid, cache=False)

        with fileStore.writeGlobalFileStream() as (fileHandle, output_id):
            cmd = ["FastTree", input_file]
            fileStore.logToMaster(" ".join(cmd))
            subprocess.check_call(cmd, stdout=fileHandle)

        return output_id


def which(executable):
    log.info("PATH = {}".format( os.environ['PATH']))
    log.info("ENV = {}".format( os.environ))
    paths = (os.path.join(p, executable)
             for p in os.environ['PATH'].split(':'))
    values = (p for p in paths if os.path.isfile(p) and os.access(p, os.X_OK))
    try:
        return next(values)
    except StopIteration:
        raise OSError("{0} not found in PATH".format(executable))


def guess_figtree_jar_path():
    """
    Try to guess the path to the figtree jar.

    Should be in $(which figtree)/../lib/figtree.jar

    raises OSError if path cannot be found.
    """
    figtree_path = which('figtree')

    jar_path = os.path.abspath(os.path.join(os.path.dirname(figtree_path),
                                            '..', 'lib', 'figtree.jar'))
    if os.path.exists(jar_path):
        return jar_path
    raise OSError("FigTree jar could not be found.")


class Figtree(Job):

    def __init__(self, inputid, width=800, height=600):
        Job.__init__(self,  memory="2G", cores=1, disk="100M")
        self.figtree_path = guess_figtree_jar_path()
        self.inputid = inputid
        self.width = width
        self.height = height

    def run(self, fileStore):
        fileStore.logToMaster("inside Fasttree run()")
        tree_file = fileStore.readGlobalFile(
            self.inputid, cache=False)

        work_dir = fileStore.getLocalTempDir()
        svg_file = os.path.join(work_dir, 'output.svg')
        cmd = ['java', '-client', '-Djava.awt.headless=true', '-Xms64m', '-Xmx512m',
               '-jar', self.figtree_path,
               '-graphic', 'SVG', '-width', str(
                   self.width), '-height', str(self.height),
               tree_file, svg_file]
        fileStore.logToMaster(" ".join(cmd))
        subprocess.check_call(cmd)

        return fileStore.writeGlobalFile(svg_file)


def start_batch(job, datafiles, input_args):
    """
    Dispatch a fasttree job for each input file.
    """
    job.fileStore.logToMaster("start_batch")
    for input_file in datafiles:
        input_id = job.fileStore.writeGlobalFile(input_file, True)
        job.fileStore.logToMaster(" Starting the run fasttree")

        job.fileStore.logToMaster(
            "creating FastTree object for {}".format(input_file))
        fasttree = Fasttree(input_id)
        treefile = job.addChild(fasttree).rv()

        figtree = Figtree(treefile)
        figfile = fasttree.addChild(figtree).rv()

        root, _ = os.path.splitext(input_file)
        svgfile = root + ".svg"
        job.addFollowOnJobFn(cleanup,
                             figfile,
                             svgfile)


def cleanup(job, temp_output_id, output_file):
    """Copies back the temporary file to input once we've successfully sorted the temporary file.
    """
    job.fileStore.logToMaster("inside cleanup job")
    tempFile = job.fileStore.readGlobalFile(temp_output_id)
    job.fileStore.logToMaster("Copying {} to {}".format(tempFile, output_file))
    shutil.copyfile(tempFile, output_file)
    job.fileStore.logToMaster(
        "Finished copying sorted file to output: %s" % output_file)


def main():
    """
    This is a Toil pipeline to transfer TCGA data into an S3 Bucket

    Data is pulled down with Genetorrent and transferred to S3 via S3AM.
    """
    # Define Parser object and add to toil
    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(
        description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--sudo', dest='sudo', default=None, action='store_true',
                        help='Docker usually needs sudo to execute locally, but not when running Mesos or when '
                             'the user is a member of a Docker group.')
    Job.Runner.addToilOptions(parser)
    parser.add_argument('datafiles', nargs='+',
                        help='FASTA input', type=existing_file)

    args = parser.parse_args()

    assert args.jobStore is not None
    config = Config()
    config.setOptions(args)

    # Store inputs from argparse
    inputs = {'sudo': args.sudo}
    datafiles = [os.path.abspath(d) for d in args.datafiles]
    # Start Pipeline
    options = Job.Runner.getDefaultOptions("./toilWorkflow")

    Job.Runner.startToil(
        Job.wrapJobFn(start_batch, datafiles, inputs), options)


if __name__ == '__main__':
    main()
