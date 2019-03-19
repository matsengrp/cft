#!/usr/bin/env python

import argparse
import subprocess
import weblogolib as w
import os
import sys

default_partis_path = os.path.join(os.getcwd(), 'partis')
partis_path = os.environ.get('PARTIS', default_partis_path)
sys.path.append(os.path.join(partis_path, 'python'))
import utils as partisutils

def write_seq_according_to_probability(seq, probability, f):
    '''
    Write the sequences into fastas according to their probability; this happens because the logo package
    draws the height of the character according to its multiplicity at that site among the sequences in the file you give it
    '''
    sig_figs = 3
    for _ in range(int(probability*10**sig_figs)):
        f.write('>%s\n%s\n' % ('naive_w_probability_{}'.format(probability), seq))

def alternative_naives_with_probabilities(f):
    '''
    Create seq, probability tuples by reading the ranked naive probabilities fasta
    '''
    seqfos = partisutils.read_fastx(f)
    return [(sfo['seq'], float(sfo['name'].split('_probability_')[1])) for sfo in seqfos]

def write_logo_input(alternative_naives, logo_input_fname, cdr3_logo_input_fname, aa_cdr3_start, aa_cdr3_end, aa_naive_len):
    '''
    Create logo input fastas for full seqs and CDR3 seqs. See write_seq_according_to_probability() for details on how 
    weblogo wants these files to look.
    '''
    with open(logo_input_fname, 'w') as logo_input, open(cdr3_logo_input_fname, 'w') as cdr3_logo_input:
        for aa_seq, probability in alternative_naives:
            write_seq_according_to_probability(aa_seq, probability, logo_input)
            # We are using the v and j codon position annotations from the most probable naive, so warn if one of the
            # alternatives is not the same length since this could mess up the cdr3 start and end.
            if len(aa_seq) != aa_naive_len:
                warn('Skipping: alternative naive sequence with length %d; differs in length from most probable naive sequence with length %d' % (len(aa_seq), aa_naive_len))
                continue
            aa_cdr3 = aa_seq[aa_cdr3_start : aa_cdr3_end]
            write_seq_according_to_probability(aa_cdr3, probability, cdr3_logo_input)

def create_logo(input_seqs_fname, logo_fname, options):
    '''
    Create a logo plot png using weblogo from a fasta created with write_logo_input()
    '''
    with open(input_seqs_fname, "rU") as f:
        seqs = w.read_seq_data(f)
    data = w.LogoData.from_seqs(seqs)

    subprocess.check_call(("rm " + input_seqs_fname).split())
    
    format = w.LogoFormat(data, options)
    with open(logo_fname, 'w') as f:
        f.write(w.png_print_formatter(data, format))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create logo plot to represent uncertainty associated with alternative naive sequences.")
    parser.add_argument(
        'input_fasta_path', type=str,
        help="Path to file with all potential partis naive seqs; their frequency in this file should represent their probability.")
    parser.add_argument(
        '--aa-cdr3-start', type=int, required=True,
        help="Amino acid position of CDR3 start codon.")
    parser.add_argument(
        '--aa-cdr3-end', type=int, required=True,
        help="Amino acid position of CDR3 end codon.")
    parser.add_argument(
        '--aa-naive-len', type=int, required=True,
        help="Length of reference (most probable) naive sequence (number of amino acids).")
    parser.add_argument(
        '--logo-fname', type=str, required=True,
        help="The output file (must be .png) name for the (complete-sequence) logo plot.")
    parser.add_argument(
        '--cdr3-logo-fname', type=str, required=True,
        help="The output file (must be .png) name for the CDR3 logo plot.")
    parser.add_argument(
        '--outdir', type=str, required=True,
        help="The output directory to write the logo plots to.")

    args = parser.parse_args()
    
    alternative_naives = alternative_naives_with_probabilities(args.input_fasta_path)
    logo_input_fname = os.path.join(args.outdir, 'naive_logo.fasta')
    cdr3_logo_input_fname = os.path.join(args.outdir, 'naive_logo_cdr3.fasta')

    write_logo_input(alternative_naives, logo_input_fname, cdr3_logo_input_fname, args.aa_cdr3_start, args.aa_cdr3_end, args.aa_naive_len)

    options = w.LogoOptions()
    options.color_scheme = w.colorscheme.chemistry
    options.unit_name = "probability"
    options.yaxis_label = "Probability"
    options.xaxis_label = "Site Position"
    options.show_fineprint = False
    options.stacks_per_line = 500
    options.tic_length = 10

    create_logo(logo_input_fname, os.path.join(args.outdir, args.logo_fname), options)
    create_logo(cdr3_logo_input_fname, os.path.join(args.outdir, args.cdr3_logo_fname), options)

