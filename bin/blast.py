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

def hit_sort_criteria(hit):
    '''
    sort by query, then by % identity, then by alignment length.
    '''
    return (hit['query acc.ver'], float(hit['% identity']), int(hit['alignment length']))

def write_blast_tsv(cline, outfname):
    '''
    parse BLAST output format 7 (TSV with weird headers) into a normal TSV
    '''
    stdout, stderr = cline()
    with open(outfname, 'r') as outfile:
        lines = list(outfile)
        headerline = [l for l in lines if 'Fields:' in l].pop()
        headers = headerline.split('\n')[0].split('Fields: ')[1].split(', ')
        datalines = [l for l in lines if '#' not in l]
    # overwrite the file with headers and data in a traditional tsv format
    with open(outfname, 'w') as outfile:
        dict_lines = [d for d  in csv.DictReader(['\t'.join(headers)] + datalines, delimiter='\t')]
        writer = csv.DictWriter(outfile, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        for l in sorted(dict_lines, key=hit_sort_criteria, reverse=True):
            writer.writerow(l)

def write_query_alignment(query_id, query_record, blast_matches, fname):
    '''
    write an alignment for a a query_sequence with with the query sequence
    followed by all of it's BLAST matches with their percent identities.
    '''
    query_matches = [d for d in blast_matches if d['query acc.ver'] == query_id]
    query_record.id += ' query'
    match_records_to_write = [SeqRecord(query_record.seq, id=query_record.id, description='', name='')]
    for match in sorted(query_matches, key=hit_sort_criteria, reverse=True):
        match_record = SeqRecord(match['match_seqrecord'].seq,
                                 id=match['subject acc.ver'],
                                 description=' nucleotide identity: {}% , aligment length: {}'.format(match['% identity'], match['alignment length']),
                                 name='')
        match_records_to_write.append(match_record)
    SeqIO.write(match_records_to_write, fname, 'fasta')

def write_all_query_alignments(blast_results_tsv, query_seqs, blast_matches):
    '''
    write an alignment for each query sequence. See write_query_alignment()
    '''
    query_dict = {record.id: record for record in SeqIO.parse(query_seqs, 'fasta')}
    for query_id, query_record in query_dict.items():
        write_query_alignment(query_id, query_record, blast_matches, blast_results_tsv.split('.tsv')[0] + '.{}.fasta'.format(query_id))
    
def blast_match_dicts(blast_results_tsv, db_seqs_fname):
    '''
    get a list of dicts corresponding to blast matches, including their Bio.SeqRecord
    '''
    db_seqs_dict = {record.id: record for record in SeqIO.parse(db_seqs_fname, 'fasta')}
    def blast_match_dict(tsv_row):
        return dict(tsv_row.items() + [('match_seqrecord', db_seqs_dict[tsv_row['subject acc.ver']])])
    with open(blast_results_tsv) as tsvfile:
        return [blast_match_dict(row) for row in csv.DictReader(tsvfile, delimiter='\t')]

def blast(blast_constructor, query_seqs_fname, db_seqs_fname, db, evalue, outfile, write_query_alignments=False):
    '''
    blast for sequences in <query_seqs_fname> among <db_seqs_fname> using <blast_constructor>, etc.
    optionally write query_alignments (see write_query_alignment()).
    '''
    blast_cline = blast_constructor(query=query_seqs_fname, db=db, evalue=evalue, outfmt='7', out=outfile)
    write_blast_tsv(blast_cline, outfile)
    if write_query_alignments:
        blast_matches = blast_match_dicts(outfile, db_seqs_fname)
        write_all_query_alignments(outfile, query_seqs_fname, blast_matches)

def make_blast_db(infile, outfile, dbtype='nucl'):
    '''
    make a blast database from <infile>.
    '''
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
        '--write-query-alignments', action='store_true',
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
    # nucleotide blast using NcbiblastnCommandline
    blast(NcbiblastnCommandline, args.query_seqs, args.db_seqs, dbfname, args.evalue, os.path.join(args.outdir, args.results_basename + '.blastn.tsv'), args.write_query_alignments)

if __name__ == '__main__':
    main()
