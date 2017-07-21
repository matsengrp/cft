#!/usr/bin/env python

import argparse
import csv
import ete3


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tree', type=lambda x: ete3.Tree(x, format=0))
    parser.add_argument('centroid_ids', type=argparse.FileType('r'))
    parser.add_argument('cluster_mapping', type=argparse.FileType('w'))
    return parser.parse_args()


def mapping(tree, centroid_ids):
    for node in tree.get_leaves():
        distance, centroid = min((node.get_distance(cid), cid) for cid in centroid_ids)
        yield (node.name, centroid, distance)


def main():
    args = get_args()
    writer = csv.writer(args.cluster_mapping)
    writer.writerow(['sequence', 'centroid', 'distance'])
    centroids = [seqid.strip() for seqid in args.centroid_ids.readlines()]
    for row in mapping(args.tree, centroids):
        writer.writerow(row)
    args.centroid_ids.close()
    args.cluster_mapping.close()



if __name__ == '__main__':
    main()


