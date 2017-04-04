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
import warnings
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
        for attr in ['svg', 'fasta', 'newick', 'cluster_aa']:
            val = self.__dict__.get(attr)
            if val:
                self.__dict__[attr] = os.path.join(dir, val)
            else:
                warnings.warn("Attr '{}' is missing for cluster:\n{}".format(attr, self.__dict__))

        # Make sure that has_seed is a boolean attribute
        self.has_seed = self.has_seed == 'True'
        self.clustering_step = int(self.clustering_step)

        # TEMPORARY HACK!
        #
        # We want to know the individual (subject) identifier, the
        # the seed sequence name, and the
        # gene name.  Ideally this would be explicitly provided in the
        # json data, but instead we must extract it from the
        # paths provided in the json data.
        #
        # Given a partial path like, "QA255.016-Vh/Hs-LN2-5RACE-IgG-new"
        # "QA255" 	- individual (subject) identification
        # "016" 	- seed sequence name
        # "Vh" 		- gene name.  Vh and Vk are the V heavy chain (IgG) vs. V kappa (there is also Vl which mean V lambda)
        #
        # "Hs-LN2-5RACE-IgG-new" is the name of the sequencing run,
        # 
        # "Hs-LN4-5RACE-IgK-100k/QB850.043-Vk/cluster4/dnaml.seedLineage.fa"
        # Note; using the data seedlineage here since we've made the self.seedlineage absolute in path
        path = data['newick']
        regex = re.compile(r'^(?P<subject_id>[^.]*).(?P<seedid>[0-9]*)-(?P<gene>[^/]*)')
        m = regex.match(path)
        if m:
            self.subject_id = m.group('subject_id')
            self.seedid = m.group('seedid')
            self.gene = m.group('gene')

        
        # Generate a unique identifier.  This id is used to index into a dict of clusters and it will be visible in URLs.
        #
        # Whatever id is used here, it must be unique across all the clusters being managed by the web interface.
        # Our particular approach here is to use a hash of "<individual> + <seed>". In theory, this could 
        # clash if (e.g.) researchers reran after some sequencing error, but because we only accept a single
        # json file these days, this scenario seems unlikely. The benefit of doing this versus just creating a
        # random uuid (or some such) is that we get persistent cluster urls across reboots.
        #
        # These identifiers may have to change if addition information is needed to maintain uniqueness, and
        # when possible care should be taken to retain existing ids. But for now this is only a soft guarantee.
        
        self.id = "c{}".format(hash((self.dataset_id, self.seed, self.clustering_step)))



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
        sorted_lineage_tips = [n for n in tree.get_leaves() if n.name in focus_node_names]
        # get naive_node for reference, (and later appending, as of this writing)
        naive_node = find_node(tree, '.*naive.*')
        # Compute lineages in a few different shapes for different needs
        lineages = map(lambda n: [naive_node] + list(reversed(n.get_ancestors())) + [n], sorted_lineage_tips)
        sequences = self.sequences(seq_mode=seq_mode, as_dict=True)
        lineages = map(lambda ns: filter(lambda n: sequences.get(n.name), ns), lineages)
        lineage_names = map(lambda ns: map(lambda n: n.name, ns), lineages)
        lineage_name_sets = map(set, lineage_names)
        # Get a set of all node names in subtree
        all_nodes = set.union(*lineage_name_sets)
        # Now we have to sort the nodes (sequences), first by lineage order (when in same lineage), then by
        # sort_fn when in separate lineages.
        dists = dict((node.name, node.get_distance(naive_node)) for nodes in lineages for node in nodes)
        # We're going to imperatively iterate through, like a savage animal.
        # We'll keep track of sorted_nodes and _stack_cursor, the former being our end result, the latter
        # being our position in the stack (keeping track of what positions in each lineage have been visited).
        def get(xs, n, missing=None):
            try:
                return xs[n]
            except IndexError:
                return missing
        # Need to be able to grab the node pointed to by a cursor easily
        def node_name_by_cursor(cursor):
            i = cursor['lineage']
            j = cursor['lineage_slice'][i]
            return lineage_names[i][j]
        # This iterator gives us the next possible positions
        def next_cursors(cursor):
            for i, j in enumerate(cursor['lineage_slice']):
                if j + 1 < len(lineage_names[i]):
                    previous_cursor_node_name = None if i == 0 else get(lineage_names[i-1], j + 1)
                    this_cursor_node_name = lineage_names[i][j + 1]
                    # We only yield a new cursor if 0th lineage, or if we share with the last lineage (already
                    # yeilded)
                    if previous_cursor_node_name != this_cursor_node_name:
                        c = copy.deepcopy(cursor)
                        c['lineage'] = i
                        c['node_name'] = this_cursor_node_name
                        # We increment every value in the cursor lineage_slice till the names change
                        for k in range(i, len(cursor['lineage_slice'])):
                            if get(lineage_names[k], j + 1) == this_cursor_node_name:
                                c['lineage_slice'][k] += 1
                        yield c
        # Our stack cursor structure
        sorted_nodes = []
        _stack_cursor = {'lineage': 0, 'lineage_slice': [0 for _ in lineage_names]}
        # precompute distances for sorting between lineages
        while len(sorted_nodes) < len(all_nodes):
            cursors_by_dist = [(dists[node_name_by_cursor(cursor)], cursor) for cursor in next_cursors(_stack_cursor)]
            if not cursors_by_dist:
                break
            dist, cursor = min(cursors_by_dist)
            node_name = node_name_by_cursor(cursor)
            _stack_cursor = cursor
            sequence = sequences[node_name]
            lineage_index = cursor['lineage'] # this is which lineage in the stack we have
            lineage_stack_index = cursor['lineage_slice'][lineage_index] # this is the index in the lineage

            # computing branch points, lineages seen, and dead lineages
            last_node = get(sorted_nodes, -1)
            last_seen = [0] if not last_node else last_node.lineage_annotations['lineages_seen']
            last_dead_lineages = [] if not last_node else last_node.lineage_annotations['dead_lineages']
            seen = copy.copy(last_seen)
            dead_lineages = copy.copy(last_dead_lineages)
            branch_to = []
            next_lineage_stack = [get(names, lineage_stack_index + 1) for names in lineage_names]
            # This is the hacky part of this; the way we iterate through things and figure out when we're
            # branching; Wish we could use a pruned ete tree for this.
            for i, _node_name in enumerate(next_lineage_stack):
                if i not in seen:
                    last_node_name = get(next_lineage_stack, i - 1)
                    # Also need to make sure this is the branch from where this actually splits... splitting
                    # is happening at wrong time todo
                    next_this_stack = get(next_lineage_stack, lineage_index)
                    if last_node_name != _node_name and last_node_name == next_this_stack and next_this_stack:
                        seen.append(i)
                        branch_to.append(i)
            if node_name in focus_node_names:
                dead_lineages.append(lineage_index)

            # target
            annotations = {'node_name': node_name,
                           'lineage_index': lineage_index,
                           'lineage_stack_index': lineage_stack_index,
                           'branch_to': branch_to,
                           'lineages_seen': seen,
                           'dead_lineages': dead_lineages}
            sequence.lineage_annotations = annotations
            sorted_nodes.append(sequence)

        processed_nodes = list(reversed(sorted_nodes))
        return processed_nodes


    def lineage_seqs(self, focus_node_name, seq_mode="dna"):
        return self.multi_lineage_seqs([focus_node_name], seq_mode=seq_mode)

    def seed_seq(self, seq_mode="dna"):
        return self.sequences(seq_mode=seq_mode, as_dict=True)[self.seed]

    def tree(self):
        "return ETE tree structure arranged to seed node to the right of all other nodes"
        try:
            tree = Tree(self.newick, format=1)
            sort_tree(tree, direction=1, predicate=lambda n: 'seed' in n.name)
            return tree
        except:
            print("Error loading tree for self:", self.__dict__)
            return Tree()

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
    return (c.subject_id, c.seed, c.clustering_step)


class ClusterDB(object):
    "This class acts as a DB or Index of cluster data, which we can query based on indexed parameters"
    @classmethod
    def fromfiles(cls, filenames):
        def fromfile(filename):
            dir = os.path.abspath(os.path.dirname(os.path.dirname(filename)))
            with open(filename, 'rb') as fp:
                print("fromfiles reading {}".format(filename))
                dataset = json.load(fp)
                dataset['clusters'] = map(lambda x: Cluster(x, dir), dataset["clusters"])
                return dataset
        return cls(map(fromfile, filenames))

    def __build_index__(self, attr):
        """Builds an index for a cluster attribute, so that clusters can be queried for more efficiently by said
        attributes."""
        attr_index  = {}
        for id, c in self._clusters.iteritems():
            c_attr_val = c.get(attr)
            try:
                attr_index[c_attr_val].add(c.id)
            except KeyError:
                attr_index[c_attr_val] = set([c.id])
        self._indices[attr] = attr_index

    def __init__(self, datasets):
        # We index clusters primarily by their primary id (should they/the-clusters be deciding ids or should
        # we (the cluster-sets)? Seems like it's our consistency concern...)
        self._indices = {}
        self._clusters = dict((c.id, c) for dataset in datasets for c in dataset['clusters'] if c)
        self._build_info = dict((d['dataset_id'], d['build_info']) for d in datasets)
        # Build indices;
        # Right now We don't really _need_ to index patients, and datasets but the tool may need to
        # eventually support larger numbers of these things, and whereas with a real file based db it's not
        # worth indexing if your number of unique values is on the order of the page block size, since
        # everything is in memory in our implementation here, we still gain some efficiency, at the cost of a
        # little extra memory consumption, so what the hey.
        for attr in ["seed", "seedid", "subject_id", "dataset_id"]:
            self.__build_index__(attr)

    def dataset_ids(self):
        return sorted(self._build_info.keys())

    def build_info(self, dataset_id):
        return copy.copy(self._build_info.get(dataset_id, {}))

    def __get_index_ids__(self, attr, vals):
        if not(type(vals) == list or type(vals) == set):
            vals = [vals]
        result_sets = [self._indices[attr].get(v) for v in vals if self._indices[attr].get(v)]
        if result_sets:
            return set.intersection(*result_sets)
        else:
            warnings.warn("No matches for sub-query:\n{}".format({'attr': attr, 'vals': vals}))
            return set()

    def get_by_id(self, id):
        return self._clusters[id]

    def get_by_ids(self, ids):
        return map(self.get_by_id, ids)

    def __query_param_ids__(self, attr, vals):
        if not(type(vals) == list or type(vals) == set):
            vals = [vals]
        # If we're getting by our cluster ids, just select by our primary key
        if attr == "id":
            return self.get_by_ids(vals)
        # If we have an index, use that
        if self._indices.get(attr):
            return self.__get_index_ids__(attr, vals)
        # Otherwise, manually filter through
        else:
            return set(id for id, c in self._clusters.iteritems() if c.get(attr) in vals)

    def query(self, params, sort_by=default_sort_key, page_n=None, per_page=10, first=False):
        """Return cluster records matching the query `params` argument (cluster_attr -> specific values), and
        optionally apply sorting and/or pagination to the results. Sort by can either be a fn or an attr str.
        First returns the first matching value."""
        # This returns the cluster ids which match the query params dictionary
        if params:
            ids = set.intersection(*(self.__query_param_ids__(k, v) for k, v in params.iteritems()))
        else:
            ids = self._clusters.keys()
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
        return self.query({'dataset_id': cluster.dataset_id, 'subject_id': cluster.subject_id, 'seed': cluster.seed})



