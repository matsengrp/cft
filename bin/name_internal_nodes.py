#!/usr/bin/env python

import argparse
import ete3

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('intree')
    parser.add_argument('outtree')
    return parser.parse_args()


def main():
    args = get_args()
    tree = ete3.Tree(args.intree, format=0)
    i = 0
    for node in tree.traverse():
        if not node.name:
            node.name = 'in-' + str(i)
            i += 1
    tree.write(outfile=args.outtree, format=1)


if __name__ == '__main__':
    main()

