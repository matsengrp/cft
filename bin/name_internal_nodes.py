#!/usr/bin/env python
# this will probably get removed when raxml-ng is implemented unless we need it
import argparse
import ete3

def get_args():
    parser = argparse.ArgumentParser(description="gives blank internal nodes names; if --splits is specified")
    parser.add_argument('intree')
    parser.add_argument('outtree')
    parser.add_argument('--prefix', default='in-')
    parser.add_argument('--splits', action='store_true')
    return parser.parse_args()

def splits_name(tree_node):
    new_name = '-'.join(sorted(n.name for n in tree_node.get_leaves()))
    return new_name

def main():
    args = get_args()
    # should have optiion for always rebuilding all internal nodes (rename?)
    # should be able to do this with the format opton in load below?
    tree = ete3.Tree(args.intree, format=1)
    i = 0
    for node in tree.traverse():
        #if not node.name:
        if not node.is_leaf():
            if args.splits:
                node.name = splits_name(node)
            else:
                node.name = args.prefix + str(i)
                i += 1
    tree.write(outfile=args.outtree, format=1)

if __name__ == '__main__':
    main()
