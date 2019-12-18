#!/usr/bin/env python

import argparse
import csv
import dendropy as dendro


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "tree",
        type=lambda x: dendro.Tree.get(
            path=x, schema="newick", preserve_underscores=True
        ),
    )
    parser.add_argument("centroid_ids", type=argparse.FileType("r"))
    parser.add_argument("cluster_mapping", type=argparse.FileType("w"))
    return parser.parse_args()


def dendro_mapping(tree, centroid_ids):
    dist_matrix = tree.phylogenetic_distance_matrix()
    terminals = [n.taxon for n in tree.leaf_node_iter()]
    centroids = [tree.find_node_with_taxon_label(cid).taxon for cid in centroid_ids]
    for node in terminals:
        distance, centroid = min(
            (dist_matrix.distance(node, centr), centr) for centr in centroids
        )
        yield (node.label, centroid.label, distance)


def main():
    args = get_args()
    writer = csv.writer(args.cluster_mapping)
    writer.writerow(["sequence", "centroid", "distance"])
    centroids = [seqid.strip() for seqid in args.centroid_ids.readlines()]
    for row in dendro_mapping(args.tree, centroids):
        writer.writerow(row)
    args.centroid_ids.close()
    args.cluster_mapping.close()


if __name__ == "__main__":
    main()
