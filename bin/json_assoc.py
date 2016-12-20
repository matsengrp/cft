#!/usr/bin/env python

import argparse
import json


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('base_document', type=argparse.FileType('r'))
    parser.add_argument('extended_document', type=argparse.FileType('w'))
    parser.add_argument('kw_pairs', nargs='+')
    args = parser.parse_args()
    # Make sure we have an even number of kw pairs
    assert((len(args.kw_pairs) % 2) == 0)
    return args


def main():
    args = get_args()
    document = json.load(args.base_document)[0]
    kw_args_iter = iter(args.kw_pairs)
    for k in kw_args_iter:
        v = next(kw_args_iter)
        document[k] = v
    json.dump(document, args.extended_document)


if __name__ == '__main__':
    main()

