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
    '''
    SymbolColor wont import
    chemistry = w.colorscheme.ColorScheme(
        [
            w.colorscheme.SymbolColor("GSTYC", "green", "polar"),
            w.colorscheme.SymbolColor("NQ", "purple", "neutral"),
            w.colorscheme.SymbolColor("KRH", "blue", "basic"),
            w.colorscheme.SymbolColor("DE", "red", "acidic"),
            w.colorscheme.SymbolColor("PAWFLIMV", "black", "hydrophobic")
        ],
        alphabet=w.seq.unambiguous_protein_alphabet
    )
    bloom = ColorScheme(
        [
            SymbolColor("GSTYC", "green", "polar"),
            SymbolColor("NQ", "purple", "neutral"),
            SymbolColor("KRH", "blue", "basic"),
            SymbolColor("DE", "red", "acidic"),
            SymbolColor("PAWFLIMV", "black", "hydrophobic")
        ],
        alphabet=seq.unambiguous_protein_alphabet
    )
    '''
    options.color_scheme = w.colorscheme.chemistry
    options.unit_name = "probability"
    options.yaxis_label = "Probability"
    options.xaxis_label = "Site Position"
    options.show_fineprint = False
    options.stacks_per_line = 500
    options.tic_length = 10

    format = w.LogoFormat(data, options)
    with open(args.output_base + ".png", 'w') as f:
        f.write(w.png_print_formatter(data, format))
