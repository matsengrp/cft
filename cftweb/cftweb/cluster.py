#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Description
   A cluster instance holds al the information associated with a cluster of 
   b-cell receiptors.  This would include the sequences and the tree built from those
   sequences.

'''
from __future__ import print_function

import os
import uuid
import sys
import json
from jinja2 import Environment, FileSystemLoader

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def load_template(name):
    # Capture parent directory above where this script lives.  
    parent = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    print("__file__ = {}".format(__file__))
    print("os.path.dirname(__file__) = {}".format(os.path.dirname(__file__)))
    print("parent = {}".format(parent))
    templatedir = os.path.join(parent, 'templates')
    env = Environment(loader=FileSystemLoader(templatedir),
                          trim_blocks=True)
    # Alias str.format to strformat in template
    env.filters['strformat'] = str.format
    template = env.get_template(name)
    return template

# Cluster is a container for the json definition of a cluster.
# each json attribute will become an attribute of the Cluster instance.
# potentially this will cause problems if a JSON attribute conflicts with an existing
# instance attribute.
class Cluster(object):
    @classmethod
    def fromfile(cls, filename):
        # create a new cluster from json definition
        c = []
        array = []
        print("fromfile attempting to open {}".format(filename))
        try:
            dir = os.path.dirname(filename)
            with open(filename, 'rb') as fp:
                print("fromfile reading {}".format(filename))
                array = json.load(fp)

            c = [cls(data, dir) for data in array if type(data) is dict]
        except ValueError as e:
            print("Missing or improperly formatted JSON file \"{}\" - ignored.".format(filename))
            

        return c
    
    def __init__(self, data, dir):
        print("creating cluster")
        super(Cluster, self).__init__()
        self.__dict__.update(data)
        self.svg = os.path.join(dir, self.svg)
        self.tree = os.path.join(dir, self.tree)
        self.fasta = os.path.join(dir, self.fasta)
        self.id = str(uuid.uuid4())

    def fasta(self):
        # return a string containing the fasta sequences
        fasta = ""
        with open(self.fasta, 'rb') as fh:
            fasta = fh.read()
        return fasta

    def sequences(self):
        # return the sequences associated with this cluster
        records = []
        with open(self.fasta, "rU") as fh:
            records = list(SeqIO.parse(fh, "fasta"))
        return records
    
    def svgstr(self):
        # return the svg tree associated with this cluster
        svgstr = ""
        with open(self.svg, 'rb') as fh:
            svgstr = fh.read()
        return svgstr

    def seeds(self):
        # return the set of seed sequences (if any ) associated with this cluster.
        pass

