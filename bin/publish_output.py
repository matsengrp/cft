#!/usr/bin/env python

import argparse
import os
import subprocess


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('-p', '--prefix', help="add a prefix")
    parser.add_argument('inputs', nargs='+')
    return parser.parse_args()

def main():
    args = get_args()
    for inpath in args.inputs:
        name_parts = inpath.split('/')[-2:]
        if args.prefix:
            name_parts = [args.prefix] + name_parts
        outpath = os.path.join(args.outdir, '-'.join(name_parts))
        subprocess.check_call(['cp', inpath, outpath])


if __name__ == '__main__':
    main()

