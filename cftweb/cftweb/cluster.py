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
import json
import pprint
import copy
from ete3 import Tree
from jinja2 import Environment, FileSystemLoader
from utils import find_node, sort_tree

from Bio import SeqIO

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
        #self.seedlineage = os.path.join(dir, self.seedlineage) if hasattr(self, 'seedlineage') else "'seedlineage' json attribute missing"
        self.cluster_aa = os.path.join(dir, self.cluster_aa) if hasattr(self, 'cluster_aa') else "'cluster_aa' json attribute missing"
        #self.seedlineage_aa = os.path.join(dir, self.seedlineage_aa) if hasattr(self, 'seedlineage_aa') else "'seedlineage_aa' json attribute missing"

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
        path = data['newick']
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



    # Public methods

    def get(self, attr):
        "Genneral purpose getter that takes an attr string and returns the corresponding attribute"
        return self.__dict__.get(attr)
        
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

    def multi_lineage_seqs(self, focus_node_names, seq_mode="dna"):
        "Note; order of focus_node_names is ignored; leaf order of self.tree is used."
        # Note also that here we are requiring that these then be leaf node names, which we didn't used to.
        # Consider relaxing this or submitting an issue.
        #
        # firt , sort the focus_nodes, as above
        tree = self.tree()
        #print("all tips:", tree.get_leaves())
        sorted_lineage_tips = [n for n in tree.get_leaves() if n.name in focus_node_names]
        #print('focus node 2:', find_node(tree, focus_node_names[1]))
        #print("focus node names", focus_node_names)
        #print("sorted_lineage_tips", sorted_lineage_tips)
        # get naive_node for reference, (and later appending, as of this writing)
        naive_node = find_node(tree, '.*naive.*')
        # Compute lineages in a few different shapes for different needs
        lineages = map(lambda n: list(reversed([n] + n.get_ancestors() + [naive_node])), sorted_lineage_tips)
        lineage_names = map(lambda ns: map(lambda n: n.name, ns), lineages)
        print("lineages", lineage_names)
        lineage_name_sets = map(set, lineage_names)
        # Get a set of all node names in subtree
        all_nodes = set.union(*lineage_name_sets)
        #print("all nodes:", all_nodes)
        # Now we have to sort the nodes (sequences), first by lineage order (when in same lineage), then by
        # sort_fn when in separate lineages.
        def shared_lineage_i(node1, node2):
            "Takes node names node1 and node2"
            print("Finding shared_lineage for:", node1, node2,)
            shared_lineages = [i for i, name_set in enumerate(lineage_name_sets) if set([node1, node2]).issubset(name_set)]
            try:
                ans = shared_lineages[0]
                print("val=", ans)
                return ans
                #return shared_lineages[0]
            except IndexError:
                print("val=", 'None')
                return None
        def node_comparator(node1, node2, sort_fn=lambda n: find_node(tree, n).get_distance(naive_node)):
            if node1 == "naive0":
                return -1
            if node2 == "naive0":
                return 1
            lineage_i = shared_lineage_i(node1, node2)
            if lineage_i == None:
                return cmp(sort_fn(node1), sort_fn(node2))
            else:
                lineage = lineage_names[lineage_i]
                return cmp(lineage.index(node1), lineage.index(node2))
        # And finally, here we have our sorted_nodes
        sorted_nodes = sorted(all_nodes, cmp=node_comparator)
        print("sorted nodes:", [x for x in sorted_nodes])
        # load up our sequences
        cluster_seqs = self.sequences(seq_mode=seq_mode, as_dict=True)
        # Now we process our nodes (in sort order) returning seqrecords with annotations about what other
        # lineages have been seen, where we have branches, which lineage the node is on, etc. Everything
        # needed to build a git tree.
        def lineage_stack_index(node_name):
            #print("lineage stack called on:", node_name)
            lineage = lineage_names[shared_lineage_i(node_name, node_name)]
            return lineage.index(node_name)
        def get(xs, n, missing=None):
            try:
                return xs[n]
            except IndexError:
                return missing
        def process_lineage_node(processed_nodes, node_name):
            try:
                sequence = cluster_seqs[node_name]
            except KeyError:
                # If we can't find the node in seqs, we leave it out; only root?
                #print("Can't find seq for:", node_name)
                return processed_nodes
            lineage_stack_i = lineage_stack_index(node_name) or 0
            current_lineage_i = shared_lineage_i(node_name, node_name)
            last_node = get(processed_nodes, -1)
            #print("testing node name", last_node)
            #print("testing result", last_node == None)
            last_seen = [0] if not last_node else last_node.lineage_annotations['lineages_seen']
            last_dead_lineages = [] if not last_node else last_node.lineage_annotations['dead_lineages']
            seen = copy.copy(last_seen)
            dead_lineages = copy.copy(last_dead_lineages)
            branch_to = []
            #lineage_stack = [get(names, lineage_stack_i) for names in lineage_names]
            next_lineage_stack = [get(names, lineage_stack_i + 1) for names in lineage_names]
            for i, _node_name in enumerate(next_lineage_stack):
                if i not in seen:
                    last_node_name = get(next_lineage_stack, i - 1)
                    if last_node_name != _node_name and i -1 == current_lineage_i:
                        seen.append(i)
                        branch_to.append(i)
            if node_name in focus_node_names:
                dead_lineages.append(current_lineage_i)
            annotations = {'node_name': node_name,
                           'lineage_index': current_lineage_i,
                           'lineage_stack_index': lineage_stack_i,
                           'branch_to': branch_to,
                           'lineages_seen': seen,
                           'dead_lineages': dead_lineages}
            sequence.lineage_annotations = annotations
            return processed_nodes + [sequence]
        # Apply process_lineage_nodes as a reduction, so we have access to the last_processed value
        processed_nodes = list(reversed(reduce(process_lineage_node, sorted_nodes, [])))
        print("processed_nodes:", [x.lineage_annotations for x in processed_nodes])
        return processed_nodes

    def lineage_seqs(self, focus_node_name, seq_mode="dna"):
        return self.multi_lineage_seqs([focus_node_name], seq_mode=seq_mode)

    def seed_seq(self, seq_mode="dna"):
        return self.sequences(seq_mode=seq_mode, as_dict=True)[self.seed]

    def tree(self):
        "return ETE tree structure arranged to seed node to the right of all other nodes"
        tree = Tree(self.newick, format=1)
        # I'm not sure we actually need this, but there are cases where we might. Better safe than sorry.
        naive_node = find_node(tree, '.*naive.*')
        if naive_node:
            tree.set_outgroup(naive_node)
        sort_tree(tree, direction=1, predicate=lambda n: 'seed' in n.name)
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


# This should really be coupled to the mode of id hash creation for clusters
def default_sort_key(c):
    return (c.pid, c.timepoint, c.seed, c.clustering_step)


class ClusterDB(object):
    "This class acts as a DB or Index of cluster data, which we can query based on indexed parameters"
    @classmethod
    def fromfile(cls, filename):
        "Load up clusters from a json file"
        dir = os.path.abspath(os.path.dirname(filename))
        with open(filename, 'rb') as fp:
            print("fromfile reading {}".format(filename))
            clusters = map(lambda x: Cluster(x, dir), json.load(fp)["clusters"])
            return cls(clusters)

    def __build_index__(self, attr):
        """Builds an index for a cluster attribute, so that clusters can be queried for more efficiently by said
        attributes."""
        attr_index  = {}
        for id, c in self.clusters.iteritems():
            c_attr_val = c.get(attr)
            try:
                attr_index[c_attr_val].add(c.id)
            except KeyError:
                attr_index[c_attr_val] = set([c.id])
        self.indices[attr] = attr_index

    def __init__(self, clusters):
        # We index clusters primarily by their primary id (should they/the-clusters be deciding ids or should
        # we (the cluster-sets)? Seems like it's our consistency concern...)
        self.indices = {}
        self.clusters = dict((c.id, c) for c in clusters if c)
        # Build indices;
        # Right now We don't really _need_ to index timepoints and patients, but the tool may need to
        # eventually
        for attr in ["seed", "seedid", "patient", "timepoints"]:
            self.__build_index__(attr)

    def __get_index_ids__(self, attr, vals):
        if not(type(vals) == list or type(vals) == set):
            vals = [vals]
        return set.intersection(*(self.indices[attr].get(v) for v in vals if self.indices[attr].get(v)))

    def get_by_id(self, id):
        return self.clusters[id]

    def get_by_ids(self, ids):
        return map(self.get_by_id, ids)

    def __query_param_ids__(self, attr, vals):
        if not(type(vals) == list or type(vals) == set):
            vals = [vals]
        # If we're getting by our cluster ids, just select by our primary key
        if attr == "id":
            return self.get_by_ids(vals)
        # If we have an index, use that
        if self.indices.get(attr):
            return self.__get_index_ids__(attr, vals)
        # Otherwise, manually filter through
        else:
            return set(id for id, c in self.clusters.iteritems() if c.get(attr) in vals)

    def query(self, params, sort_by=default_sort_key, page_n=None, per_page=10, first=False):
        """Return cluster records matching the query `params` argument (cluster_attr -> specific values), and
        optionally apply sorting and/or pagination to the results. Sort by can either be a fn or an attr str.
        First returns the first matching value."""
        # This returns the cluster ids which match the query params dictionary
        if params:
            ids = set.intersection(*(self.__query_param_ids__(k, v) for k, v in params.iteritems()))
        else:
            ids = self.clusters.keys()
        # Use these to get the clusters
        clusters = self.get_by_ids(ids)
        # Rewire as a fn if we are passed a sort_by attr str
        if type(sort_by) == str:
            sort_by = lambda c: c.get(sort_by)
        # If requested, sort by sort_by fn
        if sort_by:
            clusters = sorted(clusters, key=sort_by)
        if page_n:
            start = page_n * per_page
            clusters = clusters[start : start + per_page]
        if first:
            return clusters[0]
        return clusters

    def all(self, **kw_args):
        return self.query({}, **kw_args)

    def clustering_step_siblings(self, cluster_id):
        cluster = self.get_by_id(cluster_id)
        # Again, as mentioned above, this set of things necessary for identification should be factored out
        return self.query({'pid': cluster.pid, 'timepoint': cluster.timepoint, 'seed': cluster.seed})



