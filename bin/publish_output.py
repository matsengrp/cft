#!/usr/bin/env python

import argparse
#import json
import os
import shutil


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('inputs', nargs='+')
    return parser.parse_args()

def main():
    args = get_args()
    for inpath in args.inputs:
        outpath = os.path.join(args.outdir, '-'.join(inpath.split('/')[-2:]))
        shutil.copy(inpath, outpath)


if __name__ == '__main__':
    main()

