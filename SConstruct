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
import os.path
import subprocess

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
env['build'] = 'output'

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
    # Exit if any of the prerequisites are missing.
    msg = ""
    if not cmd_exists('dnaml'):
        msg += '''
           Required dependency command, 
           `dnaml` not found on PATH
           Consider using,
                $ module use ~matsengrp/modules
                $ module load phylip
            '''

    if not cmd_exists('seqmagick'):
        msg += '''
           Required dependency command, 
           `seqmagick` not found on PATH
           Consider using,
                $ module load seqmagick
            '''
        
    if not cmd_exists('FastTree'):
        msg += '''
           Required dependency command, 
           `FastTree` not found on PATH
           Consider using,
                $ module load FastTree
            '''
        
    if len(msg):
        warn(msg)
        sys.exit(1)


# Return True if `f` (scons File object) is a "valid" fasta file.
# "valid" in this case means the file exists, is readable, has the expected fasta format, and
# includes 2 or more sequences.
# This test is necessary to filter out short fasta files that will cause `dnaml` to crash.
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
    for path, dirs, files in os.walk(datapath):
        jsonfiles = [j for j in files if filter(j)]
        for j in jsonfiles:
            with open(os.path.join(path, j), 'r') as fh:
                meta = json.load(fh)
                for m in meta:
                    # fix-up the path to the fasta file
                    # Assume fasta files are in the same directory as the metadata
                    m['fasta'] = os.path.join(path, os.path.basename(m['file']))

                    # assume the last component of the path is the
                    # name of the run that this partis info came from.
                    # It would be better to have this explicitly
                    # included in the json info.
                    #
                    # This is used to distinguish multiple partis runs
                    # on the same seed.
                    m['run'] = os.path.basename(os.path.normpath(path))
                    if valid_fasta(m['fasta']):
                        yield m
metadata = []
svgfiles = []

for m in iter_meta(datapath):
    # create a new environment populated with metadata from the json file.
    senv = env.Clone(**m)
    senv['outdir'] = os.path.join("$build", '${seed}', '${run}', os.path.splitext(os.path.basename(senv['fasta']))[0])
    # pad the sequences so all the same length.
    padn = senv.Command(
        [ os.path.join("$outdir", "padN.fa") ],
        [ '${fasta}' ],
        "python bin/padseq.py  $SOURCE >$TARGET")

    # use fasttree to make newick tree from sequences
    newick = senv.Command(
        [ os.path.join("$outdir", "fasttree.nwk") ],
        [ padn ],
        "FastTree -nt -quiet $SOURCE >$TARGET")

    # calculate list of sequences to be pruned
    pruned_ids = senv.Command(
        [ os.path.join("$outdir", "pruned_ids.txt") ],
        [ newick ],
        "python bin/prune.py $SOURCE >$TARGET")

    # prune out sequences to reduce taxa
    fasta = senv.Command(
        [ os.path.join("$outdir", "pruned.fa") ],
        [ pruned_ids, padn ],
        "seqmagick convert --include-from-file ${SOURCES[0]} ${SOURCES[1]} ${TARGET}")
    
    phy = senv.Command(
        [ os.path.join("$outdir", "pruned.phy") ],
        [fasta],
        "seqmagick convert ${SOURCE} ${TARGET}")

    config = senv.Command(
        [ os.path.join("$outdir", "dnaml.cfg") ],
        [phy],
        "python bin/mkconfig.py ${SOURCES[0]} >${TARGETS[0]}")


    # run dnaml (from phylip package) to create tree with inferred sequences at internal nodes
    dnaml = senv.Command(
        [ # targets
            os.path.join("$outdir", "outtree"),
            os.path.join("$outdir", "outfile"),
            os.path.join("$outdir", "dnaml.log")
        ],
        [ # sources
            config
        ],
        [ # commands
        'srun --time=30 --chdir=${outdir} --output=${TARGETS[2]} dnaml <${SOURCES[0]}',
        Wait("${TARGETS[1]}")
        # "cd ${TARGETS[0].dir} && dnaml <${SOURCES[0].file} >${TARGETS[2].file}"
        ])

    # parse dnaml output into fasta and newick files, and make SVG format tree with ETE package.
    # xvfb-run is needed because of issue https://github.com/etetoolkit/ete/issues/101
    svg = senv.Command(
        [ # targets
            os.path.join("$outdir", "dnaml.svg"),
            os.path.join("$outdir", "dnaml.fa"),
            os.path.join("$outdir", "dnaml.seedLineage.fa"),
            os.path.join("$outdir", "dnaml.newick")
        ],
        [dnaml],
        "xvfb-run -a bin/dnaml2tree.py --dnaml ${SOURCES[1]} --outdir ${TARGETS[0].dir} --basename dnaml")

    m['svg'] = os.path.relpath(str(svg[0]), senv.subst('${build}'))
    m['fasta'] = os.path.relpath(str(svg[1]), senv.subst('${build}'))
    m['seedlineage'] = os.path.relpath(str(svg[2]), senv.subst('${build}'))
    m['newick'] = os.path.relpath(str(svg[3]), senv.subst('${build}'))
    del m['file']
    metadata.append(m)
    svgfiles.append(svg[0])

def write_metadata(target, source, env):
    with open(str(target[0]), "w") as fh:
        json.dump(metadata, fh, sort_keys=True,
                       indent=4, separators=(',', ': '))


metafile = env.Command(
    [ # targets
          os.path.join("$build", "metadata.json")
    ],
    svgfiles,
    write_metadata)
AlwaysBuild(metafile)

import textwrap

def print_hints(target, source, env):
    msg = """\
		hint: to run the cft web interface,
			$ cd cftweb && python -m cftweb --file {}
	""".format(os.path.abspath(str(source[0])))
    print(textwrap.dedent(msg))


hints = env.Command(
    None,
    metafile,
    print_hints)
AlwaysBuild(metafile)

