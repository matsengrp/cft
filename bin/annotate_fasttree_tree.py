#!/usr/bin/env python

import argparse
import os
from ete3 import Tree, TextFace, TreeStyle


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate FastTree tree with prune ids."
    )
    parser.add_argument("tree_path", type=str, help="Path to FastTree tree file.")
    parser.add_argument("ids_path", type=str, help="Path to prune id file.")
    parser.add_argument(
        "--naive", type=str, required=True, help="The name of the naive sequence."
    )
    parser.add_argument("--seed", type=str, help="The name of the seed sequence.")
    parser.add_argument(
        "--output-path", type=str, required=True, help="The PNG output file path."
    )
    parser.add_argument("--size", type=int, default=900, help="size in pixels of png")

    args = parser.parse_args()

    tree = Tree(args.tree_path, format=1)
    tree.set_outgroup(tree & args.naive)

    with open(args.ids_path) as f:
        ids = f.readlines()
    ids = [id.rstrip("\n") for id in ids]

    for leaf_node in tree.get_leaves():
        if leaf_node.name in ([args.naive, args.seed] if args.seed else [args.naive]):
            color = "blue"
        elif leaf_node.name in ids:
            color = "red"
        else:
            color = "black"

        node_face = TextFace(leaf_node.name, fsize=6, fgcolor=color)
        leaf_node.add_face(node_face, column=0, position="float")

    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.show_scale = False

    tree.render(args.output_path, tree_style=ts, w=args.size)
