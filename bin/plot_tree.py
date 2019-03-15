#!/usr/bin/env python
import argparse
import ete3
import math
import os
import csv
import re
import warnings
import colorbrewer

def find_node(tree, pattern):
    regex = re.compile(pattern).search
    nodes =  [ node for node in tree.traverse() for m in [regex(node.name)] if m]
    if not nodes:
        # warnings.warn("Cannot find matching node; looking for name matching '{}'".format(pattern))
        return
    else:
        if len(nodes) > 1:
            warnings.warn("multiple nodes found; using first one.\nfound: {}".format([n.name for n in nodes]))
        return nodes[0]


# reroot the tree on node matching pattern.
# Usually this is used to root on the naive germline sequence
# NOTE duplicates fcn in process_asr.py
def reroot_tree(tree, pattern):
    # find all nodes matching pattern
    node = find_node(tree, pattern)
    if tree != node:
        try:
            tree.remove_child(node)
            node.add_child(tree)
            tree.dist = node.dist
            node.dist = 0
            tree = node
        except:
            # fuckit
            pass
    return tree


# iterate up a lineage toward the root.
# lineage starts at node whose name matches pattern.
def iter_lineage(tree, pattern):
    node = find_node(tree, pattern)
    while node:
        yield node
        node = node.up


def timepoint_colors(annotations):
    # Should really improve this sort so that we're fully chronological
    def annotation_timepoints(annotation):
        return set(filter(lambda x: x, annotation.get('cluster_timepoints',
            annotation['timepoints']).split(':')))
    timepoints = sorted(reduce(set.union, (annotation_timepoints(x) for x in annotations.values()), set()))
    palette = ['#%02x%02x%02x' % c for c in colorbrewer.RdYlBu[min(max(len(timepoints), 3), 11)]] 
    # note:      ^^^ this bit    converts into a hex
    colors = [palette[i%11] for i in range(len(timepoints))]
    return dict(zip(timepoints, colors))


def timepoint_legend(ts, tp_colors):
    ts.legend.add_face(ete3.faces.TextFace(""), 0)
    ts.legend.add_face(ete3.faces.TextFace("Timepoints"), 1)
    for timepoint, color in tp_colors.iteritems():
        ts.legend.add_face(ete3.faces.CircleFace(12, color, "circle"), 0)
        ts.legend.add_face(ete3.faces.TextFace(timepoint), 1)
    ts.legend.add_face(ete3.faces.CircleFace(8, 'brown', "circle"), 0)
    ts.legend.add_face(ete3.faces.TextFace("seedlineage"), 1)


def scale_node(size):
    return 6 + (math.log(size) * 2)


def multiplicity_legend(ts):
    ts.legend.add_face(ete3.faces.TextFace(""), 0)
    ts.legend.add_face(ete3.faces.TextFace("Multiplicity"), 1)
    for multiplicity in [1, 10, 100]:
        size = scale_node(multiplicity)
        pie_face = ete3.PieChartFace([100], width=size, height=size, colors=['grey'])
        ts.legend.add_face(pie_face, 0)
        ts.legend.add_face(ete3.faces.TextFace(str(multiplicity)), 1)


def leaf_style(node, seqmeta, tp_colors, highlight_node=None):
    name = node.name + " (mf={}) ".format(round(float(seqmeta['mut_freq']), 3))
    F = ete3.faces.TextFace(name, fsize=9)
    ete3.add_face_to_node(F, node, column=1, position='branch-right')
    ## Style the node with color corresponding to timepoint
    nstyle = ete3.NodeStyle()
    #if node.name == highlight_node:
        #nstyle['fgcolor'] = 'brown'
    #else:
        #nstyle['fgcolor'] = tp_colors[seqmeta['timepoint']]
    nstyle['size'] = 0
    node.set_style(nstyle)
    # Style the node with color corresponding to timepoint
    if False:
        nstyle = ete3.NodeStyle()
        if node.name == highlight_node:
            nstyle['fgcolor'] = 'brown'
        else:
            nstyle['fgcolor'] = tp_colors[seqmeta['timepoint']]
        nstyle['size'] = 14
        node.set_style(nstyle)
    timepoints = seqmeta.get('cluster_timepoints', seqmeta['timepoints']).split(':')
    duplicities = [int(n) if n else 0 for n in seqmeta.get('cluster_timepoint_multiplicities', seqmeta['timepoint_multiplicities']).split(':')]
    multiplicity = int(seqmeta.get('cluster_multiplicity', seqmeta['multiplicity']))
    if len(timepoints) != len(duplicities):
        raise Exception('number of timepoints: {} doesn\'t match number of multiplicities: {}'.format(len(timepoints), len(duplicities)))
    percents = [d * 100 / multiplicity for d in duplicities]
    colors = [tp_colors[t] for t in timepoints]
    pie_node = ete3.PieChartFace(percents, width=scale_node(multiplicity), height=scale_node(multiplicity),
            colors=colors, line_color='black')
    ete3.add_face_to_node(pie_node, node, column=0)


def render_tree(fname, tree, annotations, highlight_node, inferred_naive_name, supports=False, support_cuttof=0.8):
    "render tree SVG"
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    tp_colors = timepoint_colors(annotations)

    seed_lineage = [n for n in iter_lineage(tree, highlight_node)]

    def my_layout(node):
        seqmeta = annotations.get(node.name)
        # handle leaves
        if seqmeta and not re.compile("naive").match(node.name):
            leaf_style(node, seqmeta, tp_colors, highlight_node)
        # Deal with naive and true internal nodes
        else:
            # Have to set all node sizes to 0 or we end up distorting the x-axis
            nstyle = ete3.NodeStyle()
            nstyle['size'] = 0
            node.set_style(nstyle)
            # Highlight just those nodes in the seedlineage which have nonzero branch lengths, or bad things
            if node in seed_lineage and node.dist > 0:
                position = "float"
                # Note; we use circleface instead of node style to avoid borking the x-axis
                circle_face = ete3.CircleFace(radius=5, color='brown', style="circle")
                ete3.add_face_to_node(circle_face, node, column=0, position=position)

    ts.layout_fn = my_layout
    if supports:
        for node in tree.get_descendants():
            if node.support and node.support < 0.8:
                node.support = None
        ts.show_branch_support = True
    timepoint_legend(ts, tp_colors)
    multiplicity_legend(ts)
    # whether or not we had rerooted on naive before, we want to do so for the SVG tree
    if inferred_naive_name not in tree.name:
        tree = reroot_tree(tree, inferred_naive_name)
    ts.scale = 2300
    tree.render(fname, tree_style=ts)



def get_args():
    def seqmeta_input(fname):
        with open(fname, 'r') as fhandle:
            return dict((row['sequence'], row) for row in csv.DictReader(fhandle))

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'tree', type=existing_file,
        help='tree file, in newick format')
    parser.add_argument(
        '--supports', action='store_true',
        help='tree file, in newick format')
    parser.add_argument(
        'seqmeta', type=seqmeta_input,
        help="Sequence metadata file for annotating mut_freqs")
    parser.add_argument('svg_out', help="output file")
    parser.add_argument(
        '--inferred-naive-name', type=str, required=True)
    parser.add_argument(
        '--seed', type=str, help="id of leaf [default 'seed']", default='seed')
    return parser.parse_args()



def main():
    args = get_args()

    # get the tree
    with open(args.tree) as fh:
        tree = ete3.Tree(fh.read(), format=(0 if args.supports else 1))

    tree = reroot_tree(tree, args.inferred_naive_name)

    # write newick file
    render_tree(args.svg_out, tree, args.seqmeta, args.seed, args.inferred_naive_name, supports=args.supports)


if __name__ == '__main__':
    main()


