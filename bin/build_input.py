#!/usr/bin/env python

import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_metadata')
    parser.add_argument('generated_infile')
    return parser.parse_args()

def main():
    args = get_args()
