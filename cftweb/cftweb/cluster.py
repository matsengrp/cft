#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description
   A cluster instance holds all the information associated with a cluster of
   b-cell receptors.  This would include the sequences and the tree built from those
   sequences.

"""
from __future__ import print_function

import re
import os
import uuid
import sys
import json
import pprint
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node
from jinja2 import Environment, FileSystemLoader
from utils import iter_names_inorder, find_node, fake_seq

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

pp = pprint.PrettyPrinter(indent=4)


def load_template(name):
    # Capture parent directory above where this script lives.
    parent = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    templatedir = os.path.join(parent, 'templates')
    env = Environment(loader=FileSystemLoader(templatedir), trim_blocks=True)
    # Alias str.format to strformat in template
    env.filters['strformat'] = str.format
    template = env.get_template(name)
    return template


# Cluster is a container for the json definition of a cluster.
# Each json attribute will become an attribute of the Cluster instance.
# Potentially this will cause problems if a JSON attribute conflicts with an existing
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
            print(
                "Missing or improperly formatted JSON file \"{}\" - ignored.".
                format(filename))

        return c

    def __init__(self, data, dir):
        super(Cluster, self).__init__()
        self.__dict__.update(data)


        # Paths that come in via json are relative to the directory in
        # which the json file is found.  Make relative paths
        # absolute so we can access these paths later.
        #
        # As the json format has evolved we have added more
        # attributes.  Test for the presence of the attribute before
        # accessing in case we are working with an older version of
        # the json file.
        
        self.svg = os.path.join(dir, self.svg) if hasattr(self, 'svg') else "'svg' json attribute missing"
        self.fasta = os.path.join(dir, self.fasta) if hasattr(self,
                                                            'fasta') else "'fasta' json attribute missing"
        self.newick = os.path.join(dir, self.newick) if hasattr(self,
                                                            'newick') else "'newick' json attribute missing"

        # TEMPORARY HACK!
        #
        # We want to know the individual (patient) identifier, the
        # timepoint (when sampled), the seed sequence name, and the
        # gene name.  Ideally this would be explicitly provided in the
        # json data, but instead we must extract it from the
        # paths provided in the json data.
        #
        # Given a partial path like, "QA255.016-Vh/Hs-LN2-5RACE-IgG-new"
        # "QA255" 	- individual (patient) identification
        # "016" 	- seed sequence name
        # "Vh" 		- gene name.  Vh and Vk are the V heavy chain (IgG) vs. V kappa (there is also Vl which mean V lambda)
        # "LN2" 	- timepoint.  LN1 and LN2 are early timepoints and LN3 and LN4 are late timepoints.
        #
        # "Hs-LN2-5RACE-IgG-new" is the name of the sequencing run,
        # 

        path = self.seedlineage
        regex = re.compile(r'^(?P<pid>[^.]*).(?P<seedid>[0-9]*)-(?P<gene>[^/]*)/[^-]*-(?P<timepoint>[^-]*)')
        m = regex.match(path)
        if m:
            self.pid = m.group('pid')
            self.seedid = m.group('seedid')
            self.gene = m.group('gene')
            self.timepoint = m.group('timepoint')

        
        # Generate a unique identifier.  This id is used to index into
        # a dict of clusters and it will be visible in URLs.
        # Generating ids in this way means that the resulting URLs are
        # not persistent; they will change every time this server is
        # restarted.
        #
        # Whatever id is used here, it must be unique across all the
        # clusters being managed by the web interface.  If the id is
        # persistent across time, all the better!
        #
        # One alternative would be generate an id by
        # concatenating parts of the cluster data, something like
        # "<individual> + <timepoint> + <seed>", but it is hard to
        # guarantee that would actually be unique.  Researchers might
        # rerun the exact same run again because of some sequencing
        # error (happens all the time!).
        #
        # Another alternative would be to use random unique ids like uuid.uuid4().hex,
        # but generate them earlier in the pipeline and stick them in the json file.
        
        self.id = uuid.uuid4().hex

        self.nseq = 0
        with open(self.fasta, "rU") as fh:
            self.nseq = len(list(SeqIO.parse(fh, 'fasta')))


        
    def fasta(self):
        # return a string containing the fasta sequences
        fasta = ""
        with open(self.fasta, 'rb') as fh:
            fasta = fh.read()
        return fasta

    def sequences(self):
        # return the sequences records associated with this cluster
        records = []
        with open(self.fasta, "rU") as fh: 
            tree = self.tree()
            record_dict = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
            records = [record_dict[id] if id in record_dict else fake_seq() for id in iter_names_inorder(tree)]
        return records

    def tree(self):
        # return the sequences records associated with this cluster
        tree = Tree(self.newick, format=1)
        naive_node = find_node(tree, '.*naive.*')
        if naive_node:
            tree.set_outgroup(naive_node)
        return tree

    def svgstr(self):
        # return the svg tree associated with this cluster
        svgstr = ""
        try:
            with open(self.svg, 'rb') as fh:
                svgstr = fh.read()
        except Exception as e:
            svgstr = "%s".format(str(e))
        return svgstr

    def seeds(self):
        # return the set of seed sequences (if any ) associated with this cluster.
        pass
