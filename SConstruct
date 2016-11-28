#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Carry partially processed partis output through steps needed 
to propare files for display via cftweb application.

This SCons file expects to work on FASTA files produced by  
`post_partis/scripts/process_partis.py`

Those fasta files are run through `dnaml` and a custom python script
to generate a SVG tree file and a fasta file containing sequences at
the internal nodes of the tree.

Certain fasta files cannot be processed by this pipeline as they have
too few sequences for dnaml to handle.  Those short fasta files are
detected and skipped automatically.

You need to have `seqmagick` and `dnaml` available in your path to use this script.
Typical usage would be,

    $ module use ~matsengrp/modules
    $ module load phylip
    $ module load seqmagick
    $ scons -j 20


'''


from __future__ import print_function
import os
import re
import sys
import subprocess
import copy
import glob

from os import path
from warnings import warn
from nestly import Nest
from nestly.scons import SConsWrap
from SCons.Script import Environment
from sconsutils import Wait

from Bio import SeqIO  # for validating fasta files.
import json


# Directory where metadata files live.
AddOption('--datapath',
              dest='datapath',
              type='string',
              nargs=1,
              action='store',
              metavar='DIR',
              default="/home/cwarth/src/matsen/cft/post_partis/output",
              help='directory where JSON metadata can be found.')

environ = os.environ.copy()

env = Environment(ENV=environ, )
env.PrependENVPath('PATH', 'bin')
env.PrependENVPath('PATH', 'post_partis/scripts')
env.PrependENVPath('PATH', 'tree')

# simpler...
outdir_base = 'output'
#env['build'] = 'output'



def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env['ENV']) == 0

# check if  tools are available, but only if we aren't in a dry-run.
if not GetOption('no_exec') and not cmd_exists('process_partis.py'):
    msg = '''
    `process_partis.py` not found on PATH
    '''
    print(msg)
    sys.exit(1)


# 
def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# check if dependencies commands are available, but only if we aren't in a dry-run.
if not GetOption('no_exec'):
    if not cmd_exists('dnaml'):
        msg = '''
           Required dependency command, 
           `dnaml` not found on PATH
           Consider using,
                $ module use ~matsengrp/modules
                $ module load phylip
            '''
        warn(msg)
        sys.exit(1)

    if not cmd_exists('seqmagick'):
        msg = '''
           Required dependency command, 
           `seqmagick` not found on PATH
           Consider using,
                $ module load seqmagick
            '''
        warn(msg)
        sys.exit(1)
        
    if not cmd_exists('FastTree'):
        msg = '''
           Required dependency command, 
           `FastTree` not found on PATH
           Consider using,
                $ module load FastTree
            '''
        warn(msg)
        sys.exit(1)


def valid_fasta(f):
    ok = False
    for i,_ in enumerate(SeqIO.parse(str(f), "fasta")):
        ok = i >= 2
        if ok:
            break
    return ok



# Iterate through metadata entries for each cluster.
#
# Metadata is contained in json files below `datapath`, one per seed.
# Each json file holds an array of metadata entries, one per cluster.
#
# Some metadata entries are filtered out because they have too few
# sequences (which will cause dnaml to crash!)


datapath = GetOption('datapath')
print("datapath = {}".format(datapath))


def iter_meta(datapath):
    filter = re.compile('^.*\.json$').search
    for meta_path, dirs, files in os.walk(datapath):
        jsonfiles = [j for j in files if filter(j)]
        for j in jsonfiles:
            with open(path.join(meta_path, j), 'r') as fh:
                meta = json.load(fh)
                for m in meta:
                    # fix-up the path to the fasta file
                    # Assume fasta files are in the same directory as the metadata
                    m['fasta'] = path.join(meta_path, path.basename(m['file']))

                    # assume the last component of the path is the
                    # name of the run that this partis info came from.
                    # It would be better to have this explicitly
                    # included in the json info.
                    #
                    # This is used to distinguish multiple partis runs
                    # on the same seed.
                    m['run'] = path.basename(path.normpath(meta_path))
                    if valid_fasta(m['fasta']):
                        yield m


# Set up nestly
nest = Nest()
w = SConsWrap(nest, outdir_base, alias_environment=env)

# TODO These should now be aggregates now
#metadata = []
#svgfiles = []

w.add_aggregate('metadata', list)
w.add_aggregate('svgfiles', list)


# Our first nest; By metadata file
w.add('cluster',
        iter_meta(datapath),
        # Label the nodes as run/seed/fasta-filename/
        label_func=lambda m: path.join(m["run"], m["seed"], path.basename(path.splitext(m['fasta'])[0])))

@w.add_target()
def padn(outdir, c):
    return env.Command(
        path.join(outdir, "padN.fa"),
        c['cluster']['fasta'],
        "python bin/padseq.py  $SOURCE >$TARGET")

@w.add_target()
def newick(outdir, c):
    return env.Command(
        path.join(outdir, "fasttree.nwk"),
        c['padn'],
        "FastTree -nt -quiet $SOURCE > $TARGET")

@w.add_target()
def pruned_ids(outdir, c):
    return env.Command(
        path.join(outdir, "pruned_ids.txt"),
        c['newick'],
        "python bin/prune.py $SOURCE > $TARGET")

@w.add_target()
def fasta(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.fa"),
        [c['pruned_ids'], c['padn']],
        "seqmagick convert --include-from-file $SOURCES $TARGET")

@w.add_target()
def phy(outdir, c):
    return env.Command(
        path.join(outdir, "pruned.phy"),
        c['fasta'],
        "seqmagick convert $SOURCE $TARGET")

@w.add_target()
def dnaml_config(outdir, c):
    return env.Command(
        path.join(outdir, "dnaml.cfg"),
        c['phy'],
        "python bin/mkconfig.py $SOURCE > $TARGET")

@w.add_target()
def dnaml(outdir, c):
    return env.Command(
        map(lambda x: path.join(outdir, x), ["outtree", "outfile", "dnaml.log"]),
        c['dnaml_config'],
        ['srun --time=30 --chdir=' + outdir + ' --output=${TARGETS[2]} dnaml < $SOURCE',
            Wait("${TARGETS[1]}")])

def extended_metadata(metadata, dnaml_tree_tgt):
    m = copy.copy(metadata)
    m['svg'] = path.relpath(str(dnaml_tree_tgt[0]), outdir_base)
    m['fasta'] = path.relpath(str(dnaml_tree_tgt[1]), outdir_base)
    m['seedlineage'] = path.relpath(str(dnaml_tree_tgt[2]), outdir_base)
    del m['file']
    return m

@w.add_target()
def dnaml_tree(outdir, c):
    tgt = env.Command(
            [path.join(outdir, "dnaml.svg"),
                path.join(outdir, "dnaml.fa"),
                path.join(outdir, "dnaml.seedLineage.fa")],
            c['dnaml'],
            "xvfb-run -a bin/dnaml2tree.py --dnaml ${SOURCES[1]} --outdir ${TARGETS[0].dir} --basename dnaml")
    # Do aggregate work
    c['svgfiles'].append(tgt[0])
    m = extended_metadata(c['cluster'], tgt)
    c['metadata'].append(m)
    return tgt

def write_metadata(metadata):
    def build_fn(target, source, env):
        with open(str(target[0]), "w") as fh:
            json.dump(metadata, fh, sort_keys=True,
                           indent=4, separators=(',', ': '))
    return build_fn

w.pop('cluster')

@w.add_target()
def metadata(outdir, c):
    tgt = env.Command(
        [path.join(outdir, "metadata.json")],
        c['svgfiles'],
        write_metadata(c['metadata']))
    AlwaysBuild(tgt)
    return tgt

import textwrap

def print_hints(target, source, env):
    msg = """
        hint: to run the cft web interface,
            $ cd cftweb && python -m cftweb -d {}
    """.format(path.abspath(str(source[0])))
    print(textwrap.dedent(msg))

@w.add_target()
def hints(outdir, c):
    hints = env.Command(
        None,
        c['metadata'],
        print_hints)
    AlwaysBuild(hints)
    return hints


