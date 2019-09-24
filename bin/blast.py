#!/usr/bin/env python
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os

def write_blast_tsv(cline, outfname):
    stdout, stderr = cline()
    with open(outfname, 'r') as fh:
        lines = list(fh)
        fieldline = [l for l in lines if "Fields:" in l].pop()
        fields = fieldline.split("Fields: ")[1].split(", ")
        print fields
        datalines = [l for l in lines if "#" not in l]
    with open(outfname, 'w') as fh:
        fh.write("\t".join(fields))
        for l in datalines:
            fh.write(l)
        fh.close()

def blast_results(fname):
    with open(fname) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        return list(blast_records)
    
def print_blast_record(record):
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.01:
                print 'score', hsp.score

if __name__ == '__main__':
    outdir = "/home/matsengrp/working/eharkins/cft/output/test_issue_186-qa255-synth-v18/QA255/QA255-g-merged/QA255.067-Vh/seed-part-2/seed-cluster/seed_lineage-dnaml"
    #outdir = "/home/matsengrp/working/eharkins/cft/output/test_issue_186-qa255-synth-v18/QA255/QA255-l-merged/QA255.067-VL/seed-part-1/seed-cluster/seed_lineage-dnaml"
    fname = os.path.join(outdir, "asr.ancestors_naive_and_seed.fa")
    outfname =  "blasttest.tsv"
    dbname = os.path.join(outdir, "blasttest")
    blastn_cline = NcbiblastnCommandline(query=fname, db=dbname, evalue=0.001, outfmt="7", out=os.path.join(outdir, "blastn_" + outfname))
    tblastx_cline = NcbitblastxCommandline(query=fname, db=dbname, evalue=0.001, outfmt=7, out=os.path.join(outdir, "tblastx_" + outfname))
    write_blast_tsv(blastn_cline, os.path.join(outdir, "blastn_" + outfname))
    write_blast_tsv(tblastx_cline, os.path.join(outdir, "tblastx_" + outfname))
    '''
    result = blast_results(outfname)
    for r in result:
        print_blast_record(r)
    '''
