#!/usr/bin/env python
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

log = logging.getLogger(__name__)


class Fasttree(Job):

    def __init__(self, inputid):
        log.info("inside Fasttree init()")
        Job.__init__(self,  memory="2G", cores=1, disk="100M")
        self.inputid = inputid

    def run(self, fileStore):
        log.info("inside Fasttree run()")
        input_file = fileStore.readGlobalFile(
            self.inputid, cache=False)

        work_dir = fileStore.getLocalTempDir()
        with open(os.path.join(work_dir, 'output.tree'), "wb") as stdout:
            cmd = ["fasttree", input_file]
            subprocess.check_call(cmd, stdout=stdout)
            log.info(cmd)

        return fileStore.writeGlobalFile(os.path.join(work_dir, 'output.tree'))


class Figtree(Job):

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

    def __init__(self, inputid, width=800, height=600):
        Job.__init__(self,  memory="2G", cores=1, disk="100M")
        self.figtree_path = guess_figtree_jar_path()
        self.inputid = inputid
        self.width = width
        self.height = height

    def run(self, fileStore):
        log.info("inside Fasttree run()")
        input_file = fileStore.readGlobalFile(
            self.inputid, cache=False)

        work_dir = fileStore.getLocalTempDir()
        tree_path = os.path.join(work_dir, 'output.tree')
        svg_path = os.path.join(work_dir, 'output.svg')
        with open(os.path.join(work_dir, 'output.tree'), "wb") as stdout:
            cmd = ['java', '-client', '-Djava.awt.headless=true', '-Xms64m', '-Xmx512m',
                   '-jar', self.figtree_path,
                   '-graphic', 'SVG', '-width', str(self.width), '-height', str(self.height),
                   tree_path, svg_path]

            subprocess.check_call(cmd, stdout=stdout)
            log.info(cmd)

        return fileStore.writeGlobalFile(os.path.join(work_dir, 'output.svg'))


def start_batch(job, datafiles, input_args):
    """
    Dispatch a fasttree job for each input file.
    """
    for input_file in datafiles:
        input_filestore_id = job.fileStore.writeGlobalFile(input_file, True)
        job.fileStore.logToMaster(" Starting the run fasttree")

        log.info("creating FastTree object for {}".format(input_file))
        job.addChild(Fasttree(input_filestore_id))


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
    # Store inputs from argparse
    inputs = {'sudo': args.sudo}
    datafiles = [os.path.abspath(d) for d in args.datafiles]
    # Start Pipeline
    Job.Runner.startToil(
        Job.wrapJobFn(start_batch, datafiles, inputs), args)


if __name__ == '__main__':
    main()
