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
from utils import iter_names_inorder, find_node, fake_seq, sort_tree

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
            dir = os.path.abspath(os.path.dirname(filename))
            with open(filename, 'rb') as fp:
                print("fromfile reading {}".format(filename))
                array = json.load(fp)["clusters"]
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
        #
        # Should really trigger a warning if any of these things happens...

        self.svg = os.path.join(dir, self.svg) if hasattr(self, 'svg') else "'svg' json attribute missing"
        self.fasta = os.path.join(dir, self.fasta) if hasattr(self, 'fasta') else "'fasta' json attribute missing"
        self.newick = os.path.join(dir, self.newick) if hasattr(self, 'newick') else "'newick' json attribute missing"
        self.seedlineage = os.path.join(dir, self.seedlineage) if hasattr(self, 'seedlineage') else "'seedlineage' json attribute missing"
        self.cluster_aa = os.path.join(dir, self.cluster_aa) if hasattr(self, 'cluster_aa') else "'cluster_aa' json attribute missing"
        self.seedlineage_aa = os.path.join(dir, self.seedlineage_aa) if hasattr(self, 'seedlineage_aa') else "'seedlineage_aa' json attribute missing"

        # Make sure that has_seed is a boolean attribute
        self.has_seed = self.has_seed == 'True'
        self.clustering_step = int(self.clustering_step)

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
        # "Hs-LN4-5RACE-IgK-100k/QB850.043-Vk/cluster4/dnaml.seedLineage.fa"
        # Note; using the data seedlineage here since we've made the self.seedlineage absolute in path
        path = data['seedlineage']
        regex = re.compile(r'^(?P<pid>[^.]*).(?P<seedid>[0-9]*)-(?P<gene>[^/]*)/[^-]*-(?P<timepoint>[^-]*)')
        m = regex.match(path)
        if m:
            self.pid = m.group('pid')
            self.seedid = m.group('seedid')
            self.gene = m.group('gene')
            self.timepoint = m.group('timepoint')

        
        # Generate a unique identifier.  This id is used to index into a dict of clusters and it will be visible in URLs.
        #
        # Whatever id is used here, it must be unique across all the clusters being managed by the web interface.
        # Our particular approach here is to use a hash of "<individual> + <timepoint> + <seed>". In theory, this could 
        # clash if (e.g.) researchers reran after some sequencing error, but because we only accept a single
        # json file these days, this scenario seems unlikely. The benefit of doing this versus just creating a
        # random uuid (or some such) is that we get persistent cluster urls across reboots.
        #
        # These identifiers may have to change if addition information is needed to maintain uniqueness, and
        # when possible care should be taken to retain existing ids. But for now this is only a soft guarantee.
        
        self.id = "c{}".format(hash((self.seed, self.timepoint, self.clustering_step)))

        # Compute and cache number of seqs by reading through the seq file
        self.nseq = 0
        with open(self.fasta, "rU") as fh:
            self.nseq = len(list(SeqIO.parse(fh, 'fasta')))


    # Public methods
        
    def fasta(self):
        # return a string containing the fasta sequences
        fasta = ""
        with open(self.fasta, 'rb') as fh:
            fasta = fh.read()
        return fasta

    def sequences(self, seq_mode="dna", as_dict=False):
        "return the sequences records associated with this cluster"
        assert(seq_mode in set(["dna", "aa"]))
        records = []
        with open(self.fasta if seq_mode == "dna" else self.cluster_aa, "rU") as fh: 
            records = (SeqIO.to_dict if as_dict else list)(SeqIO.parse(fh, "fasta"))
        return records

    def lineage_seqs(self, focus_node_name, seq_mode="dna"):
        tree = self.tree()
        focus_node = tree.search_nodes(name=focus_node_name)[0]
        lineage = [focus_node] + focus_node.get_ancestors()
        cluster_seqs = self.sequences(seq_mode=seq_mode, as_dict=True)
        seqs = [cluster_seqs[n.name] for n in lineage if n.name in cluster_seqs]
        # Add naive to the very bottom (we should really just rewrite things to ensure naive is the root)
        seqs.append(cluster_seqs["naive0"])
        return seqs

    def seed_seq(self, seq_mode="dna"):
        return self.sequences(seq_mode=seq_mode, as_dict=True)[self.seed]

    def tree(self):
        "return ETE tree structure arranged to seed node to the right of all other nodes"
        tree = Tree(self.newick, format=1)
        #naive_node = find_node(tree, '.*naive.*')
        #if naive_node:
            #tree.set_outgroup(naive_node)
        #sort_tree(tree, direction=1, predicate=lambda n: 'seed' in n.name)
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


