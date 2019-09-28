#!/usr/bin/env python
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse
import subprocess
import warnings
import csv

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

def write_blast_alignment(blast_results_tsv, query_seqs, db_seqs):
    #TODO write fcn description, validate, and clean up
    query_dict = {record.id: record for record in SeqIO.parse(query_seqs, 'fasta')}
    db_seqs_dict = {record.id: record for record in SeqIO.parse(db_seqs, 'fasta')}
    with open(blast_results_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        #TODO probably break into a real for loop since reading this is hard
        blast_match_rows = [dict(r.items() + [('match_seq_record', db_seqs_dict[r['subject acc.ver']])]) for r in reader]
    for query_id, query_record in query_dict.items():
        query_match_rows = [r for r in blast_match_rows if r['query acc.ver'] == query_id]
        query_record.id += ' query'
        match_rows_to_write = sorted(query_match_rows, key=lambda r: float(r['% identity']), reverse=True)
        match_records_to_write = [SeqRecord(query_record.seq, id=query_record.id, description='', name='')]
        for match_row in match_rows_to_write:
            match_record = SeqRecord(match_row['match_seq_record'].seq, id='{} {}% identity'.format(match_row['subject acc.ver'], match_row['% identity']), description='', name='')
            match_records_to_write.append(match_record)
        fname = blast_results_tsv.split('.tsv')[0] + '.{}.'.format(query_id) + '.fasta'
        SeqIO.write(match_records_to_write, fname, 'fasta')

def blast(blast_constructor, query_seqs_fname, db_seqs_fname, db, evalue, outfile, write_blast_alignments):
    blast_cline = blast_constructor(query=query_seqs_fname, db=db, evalue=evalue, outfmt='7', out=outfile)
    write_blast_tsv(blast_cline, outfile)
    if write_blast_alignments:
        write_blast_alignment(outfile, query_seqs_fname, db_seqs_fname)

def make_blast_db(infile, outfile, dbtype='nucl'):
    try:
        subprocess.check_call(['makeblastdb', '-in', infile, '-dbtype', dbtype, '-out', outfile])
    except OSError, e:
        warnings.warn('blast command line tools not installed!')
        raise

def parse_args():
    parser = argparse.ArgumentParser(description='Create a BLAST database from <db_seqs>, query for <query_seqs>, and write results to a TSV (actually two TSVs: one for nucleotide and one for protein identitiy).')
    parser.add_argument(
        'db_seqs', type=str,
        help='Fasta containing sequences to create a BLAST database from.')
    parser.add_argument(
        'query_seqs', type=str,
        help='Fasta containing the set of query sequences.')
    parser.add_argument(
        '--outdir', type=str,
        required=True,
        help='Name of directory to write blast results to.')
    parser.add_argument(
        '--results-basename', type=str,
        default='blast_results',
        help='Base for the name of TSV file to write blast results to. Output TSVs will look like <results-basename>.{blastn, tblastx}.tsv')
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
    dbfname = os.path.join(args.outdir, 'blast_db')
    make_blast_db(args.db_seqs, dbfname)
    # nucleotide blast
    blast(NcbiblastnCommandline, args.query_seqs, args.db_seqs, dbfname, args.evalue, os.path.join(args.outdir, args.results_basename + '.blastn.tsv'), args.write_blast_alignments)
    # translated nucleotide (both db sequences and query sequences get translated using tblastx strategy) blast
    blast(NcbitblastxCommandline, args.query_seqs, args.db_seqs, dbfname, args.evalue, os.path.join(args.outdir, args.results_basename + '.tblastx.tsv'), args.write_blast_alignments)

if __name__ == '__main__':
    main()
