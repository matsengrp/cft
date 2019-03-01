#!/usr/bin/env python

import argparse
import subprocess
import weblogolib as w

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create logo to represent alternative-annotations naive probabilities rom partition-cluster-annotations.csv.")
    parser.add_argument(
        'input_fasta_path', type=str,
        help="Path to file with all potential partis naive seqs; their frequency in this file should represent their probability.")
    parser.add_argument(
        '--output-base', type=str, required=True,
        help="The output basename.")

    args = parser.parse_args()
    
    with open(args.output_base + ".fasta", "rU") as f:
        seqs = w.read_seq_data(f)
    data = w.LogoData.from_seqs(seqs)
    subprocess.check_call(("rm " + args.output_base + ".fasta").split())

    options = w.LogoOptions()
    options.unit_name = "probability"
    options.yaxis_label = "Probability"
    options.xaxis_label = "Site Position"
    options.show_fineprint = False
    options.stacks_per_line = 500
    options.tic_length = 10

    format = w.LogoFormat(data, options)
    with open(args.output_base + ".png", 'w') as f:
        f.write(w.png_print_formatter(data, format))
