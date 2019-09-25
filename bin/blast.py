#!/usr/bin/env python
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import os
import argparse
import subprocess

def write_blast_tsv(cline, outfname):
    stdout, stderr = cline()
    with open(outfname, 'r') as fh:
        lines = list(fh)
        fieldline = [l for l in lines if 'Fields:' in l].pop()
        fields = fieldline.split('Fields: ')[1].split(', ')
        datalines = [l for l in lines if '#' not in l]
    # overwrite the file with headers and data in a traditional tsv format
    with open(outfname, 'w') as fh:
        fh.write('\t'.join(fields))
        for l in datalines:
            fh.write(l)
        fh.close()

def blast(blast_constructor, query, db, evalue, outfile):
    blast_cline = blast_constructor(query=query, db=db, evalue=evalue, outfmt='7', out=outfile)
    write_blast_tsv(blast_cline, outfile)

def make_blast_db(infile, outfile, dbtype='nucl'):
    subprocess.check_call(['makeblastdb', '-in', infile, '-dbtype', dbtype, '-out', outfile])

def parse_args():
    parser = argparse.ArgumentParser(description='Find the closest N sampled sequences to each inferred ancestor.')
    parser.add_argument(
        'sampled_seqs', type=str,
        help='Fasta containing sampled sequences to match ancestors with.')
    parser.add_argument(
        'asr_seqs', type=str,
        help='Fasta containing the set of inferred ancestral sequences.')
    parser.add_argument(
        '--outdir', type=str,
        required=True,
        help='Name of JSON file to write the results to.')
    parser.add_argument(
        '--write-blast-alignments', action='store_true',
        default=False,
        help='If set, writes a fasta for each query sequence containing the aligned blast results for that query sequence.')
    parser.add_argument(
        '--evalue', type=float,
        default=0.001,
        help='''number of hits one can "expect" to see by chance when searching a database of a particular size. It decreases exponentially as the Score (S) of the match increases. Essentially, the E value describes the random background noise. For example, an E value of 1 assigned to a hit can be interpreted as meaning that in a database of the current size one might expect to see 1 match with a similar score simply by chance.''')
    return parser.parse_args()

def main():
    args = parse_args()
    outfname =  'sampled_ancestors.tsv'
    dbfname = os.path.join(args.outdir, 'sampled_seqs_blast_db')
    make_blast_db(args.sampled_seqs, dbfname)
    # nucleotide blast
    blast(NcbiblastnCommandline, args.asr_seqs, dbfname, args.evalue, os.path.join(args.outdir, 'blastn_' + outfname))
    # translated nucleotide (both db sequences and query sequences get translated using tblastx strategy) blast
    blast(NcbitblastxCommandline, args.asr_seqs, dbfname, args.evalue, os.path.join(args.outdir, 'tblastx_' + outfname))
    #TODO implement --write-blast-alignments

if __name__ == '__main__':
    main()
