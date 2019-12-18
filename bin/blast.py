#!/usr/bin/env python
import Bio
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
    """
    returns tuple used to sort by query, then by % identity, then by alignment length.
    """
    return (
        hit["query acc.ver"],
        float(hit["% identity"]),
        int(hit["alignment length"]),
    )


def write_blast_tsv(cline, outfname):
    """
    parse BLAST output format 7 (TSV with weird headers) into a normal TSV
    """

    def check_hitless_queries(lines):
        hitless_queries_count = len([l for l in lines if "# 0 hits found" in l])
        if hitless_queries_count > 0:
            print "BLAST found no hits for {} of the query sequences.".format(
                hitless_queries_count
            )

    def get_headers(lines):
        headerlines = [l for l in lines if "Fields:" in l]
        if len(headerlines) > 0:
            return headerlines.pop().split("\n")[0].split("Fields: ")[1].split(", ")
        return []

    try:
        stdout, stderr = cline()  # Run the actual blast CLI
    except Bio.Application.ApplicationError, e:
        raise
    with open(outfname, "r") as outfile:  # Read the blast CLI ouput
        lines = list(outfile)
        check_hitless_queries(lines)
        headers = get_headers(lines)
        datalines = [l for l in lines if "#" not in l]
    with open(
        outfname, "w"
    ) as outfile:  # Overwrite blast CLI output with headers in a traditional tsv format
        dict_lines = [
            d for d in csv.DictReader(["\t".join(headers)] + datalines, delimiter="\t")
        ]
        writer = csv.DictWriter(outfile, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for l in sorted(dict_lines, key=hit_sort_criteria, reverse=True):
            writer.writerow(l)


def write_query_alignment(query_id, query_record, blast_matches, fname):
    """
    write an alignment for a a query_sequence with with the query sequence
    followed by all of its BLAST matches with their percent identities.
    """
    query_matches = [d for d in blast_matches if d["query acc.ver"] == query_id]
    query_record.id += " query"
    match_records_to_write = [
        SeqRecord(query_record.seq, id=query_record.id, description="", name="")
    ]
    for match in sorted(query_matches, key=hit_sort_criteria, reverse=True):
        match_record = SeqRecord(
            match["match_seqrecord"].seq,
            id=match["subject acc.ver"],
            description=" nucleotide identity: {}% , aligment length: {}".format(
                match["% identity"], match["alignment length"]
            ),
            name="",
        )
        match_records_to_write.append(match_record)
    SeqIO.write(match_records_to_write, fname, "fasta")


def write_all_query_alignments(blast_results_tsv, query_seqs_fasta, blast_matches):
    """
    write an alignment for each query sequence. See write_query_alignment()
    """
    query_dict = {
        record.id: record for record in SeqIO.parse(query_seqs_fasta, "fasta")
    }
    for query_id, query_record in query_dict.items():
        write_query_alignment(
            query_id,
            query_record,
            blast_matches,
            blast_results_tsv.split(".tsv")[0] + ".{}.fasta".format(query_id),
        )


def top_hit_rows(query_matches, top_n_hits, query_record=None, include_seqs=False):
    """
    get the top <top_n_hits> matches for a given query. optionally include sequences in returned dicts
    """
    top_n_matches = query_matches[:top_n_hits]
    for match in top_n_matches:
        if include_seqs and query_record is not None:
            match["query seq"] = str(query_record.seq)
            match["subject (match) seq"] = str(match["match_seqrecord"].seq)
        del match["match_seqrecord"]
    return top_n_matches


def write_top_hits_summary_tsv(
    blast_results_tsv, query_seqs_fasta, blast_matches, top_n_hits
):
    """
    write out a tsv with info on the top <top_n_hits> matches among <blast_matches> for each query in <query_seqs_fasta> including sequences. duplicates some code from write_all_query_alignments and write_query_alignment but simpler to keep them separate
    """
    query_dict = {
        record.id: record for record in SeqIO.parse(query_seqs_fasta, "fasta")
    }
    hit_summaries = []
    for query_id, query_record in query_dict.items():
        query_matches = [d for d in blast_matches if d["query acc.ver"] == query_id]
        hit_summaries += top_hit_rows(
            query_matches, top_n_hits, query_record, include_seqs=True
        )
    with open(
        blast_results_tsv.split(".tsv")[0] + ".top_{}_hits.tsv".format(top_n_hits), "w"
    ) as outfile:
        headers = [
            "query acc.ver",
            "query seq",
            "subject acc.ver",
            "subject (match) seq",
            "% identity",
            "alignment length",
            "mismatches",
            "gap opens",
            "q. start",
            "q. end",
            "s. start",
            "s. end",
            "evalue",
            "bit score",
        ]
        writer = csv.DictWriter(outfile, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for r in hit_summaries:
            writer.writerow(r)


def blast_match_dicts(blast_results_tsv, db_seqs_fname):
    """
    get a list of dicts corresponding to blast matches, including their Bio.SeqRecord
    """
    db_seqs_dict = {record.id: record for record in SeqIO.parse(db_seqs_fname, "fasta")}

    def blast_match_dict(tsv_row):
        return dict(
            tsv_row.items()
            + [("match_seqrecord", db_seqs_dict[tsv_row["subject acc.ver"]])]
        )

    with open(blast_results_tsv) as tsvfile:
        return [
            blast_match_dict(row) for row in csv.DictReader(tsvfile, delimiter="\t")
        ]


def blast(
    blast_constructor,
    query_seqs_fname,
    db_seqs_fname,
    db,
    evalue,
    outfile,
    write_query_alignments=False,
    top_n_hits=None,
):
    """
    blast for sequences in <query_seqs_fname> among <db_seqs_fname> using <blast_constructor>, etc.
    optionally write query_alignments (see write_query_alignment()).
    """
    blast_cline = blast_constructor(
        query=query_seqs_fname, db=db, evalue=evalue, outfmt="7", out=outfile
    )
    write_blast_tsv(blast_cline, outfile)
    if write_query_alignments:
        blast_matches = blast_match_dicts(outfile, db_seqs_fname)
        write_all_query_alignments(outfile, query_seqs_fname, blast_matches)
    if top_n_hits:
        blast_matches = blast_match_dicts(outfile, db_seqs_fname)
        write_top_hits_summary_tsv(outfile, query_seqs_fname, blast_matches, top_n_hits)


def make_blast_db(infile, outfile, dbtype="nucl"):
    """
    make a blast database from <infile>.
    """
    try:
        subprocess.check_call(
            ["makeblastdb", "-in", infile, "-dbtype", dbtype, "-out", outfile]
        )
    except OSError, e:
        warnings.warn("blast command line tools not installed!")
        raise


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a BLAST database from <db_seqs>, query for <query_seqs>, and write results to a TSV (actually two TSVs: one for nucleotide and one for protein identitiy)."
    )
    parser.add_argument(
        "db_seqs",
        type=str,
        help="Fasta containing sequences to create a BLAST database from.",
    )
    parser.add_argument(
        "query_seqs", type=str, help="Fasta containing the set of query sequences."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Name of directory to write blast results to.",
    )
    parser.add_argument(
        "--results-basename",
        type=str,
        default="blast_results",
        help="Base for the name of TSV file to write blast results to. Output TSVs will look like <results-basename>.{blastn, tblastx}.tsv",
    )
    parser.add_argument(
        "--write-query-alignments",
        action="store_true",
        default=False,
        help="If set, writes a fasta for each query sequence containing the aligned blast results for that query sequence.",
    )
    parser.add_argument(
        "--top-n-hits",
        default=None,
        type=int,
        help="If set, writes a separate TSV file (including columns containing the query and match sequences) with just this many hits for each query.",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=0.001,
        help="""number of hits one can "expect" to see by chance when searching a database of a particular size. It decreases exponentially as the Score (S) of the match increases. Essentially, the E value describes the random background noise. For example, an E value of 1 assigned to a hit can be interpreted as meaning that in a database of the current size one might expect to see 1 match with a similar score simply by chance.""",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    dbfname = os.path.join(args.outdir, "blast_db")
    make_blast_db(args.db_seqs, dbfname)
    # nucleotide blast using NcbiblastnCommandline
    blast(
        NcbiblastnCommandline,
        args.query_seqs,
        args.db_seqs,
        dbfname,
        args.evalue,
        os.path.join(args.outdir, args.results_basename + ".blastn.tsv"),
        args.write_query_alignments,
        args.top_n_hits,
    )


if __name__ == "__main__":
    main()
