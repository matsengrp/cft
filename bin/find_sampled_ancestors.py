#!/usr/bin/env python
import argparse
import collections
import json
from Bio import SeqIO
from Bio import pairwise2

def percent_similarity(seq1, seq2):
    '''
    return percent similarity between two sequences according to the global alignment
    function pairwise2.align.globalxx, which returns a list of different possible
    alignments all with the same score. Score is number of matching positions. See http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc87
    '''
    alignments = pairwise2.align.globalxx(seq1, seq2)
    _, _, score, _, _ = alignments.pop()
    return score/float(len(seq1))

def match_summaries(query_record, ranked_match_records):
    def match_summary(record):
        return collections.OrderedDict([("id", record.id),
                                        ("dna_sequence", str(record.seq)),
                                        ("dna_similarity", record.dna_similarity),
                                        ("dna_similarity_over_aa_similarity", record.dna_similarity/float(record.aa_similarity)),
                                        ("aa_sequence", str(record.aa)),
                                        ("aa_similarity", record.aa_similarity)])
    return collections.OrderedDict([("inferred_ancestor", query_record.id),
                                    ("dna_sequence", str(query_record.seq)),
                                    ("aa_sequence", str(query_record.aa)),
                                    ("matches", [match_summary(record) for record in ranked_match_records])])

# Dont need to break out into another function quite yet:
#def find_seqrecord_closest_matches(query_record, records_to_match, match_count, rank_by='dna'

def get_sampled_ancestors(sampled_seqs, inferred_ancestors, match_count, rank_by='dna', aa=False):
    '''
    find the closest <match_count> sampled sequences in <sampled_seqs> to each inferred ancestor in <inferred_ancestors>.
    '''
    # TODO address time complexity? Majority of time is spent doing pairwise align function. At least make sure we aren't doing this more than we need to
    # TODO should we be able to take in user specified aa seqs or should we always translate here?
    sampled_ancestors = []
    for ancestor_record in inferred_ancestors:
        ancestor_record.aa = ancestor_record.seq.translate()
        for sampled_record in sampled_seqs:
            sampled_record.dna_similarity = percent_similarity(ancestor_record.seq, sampled_record.seq)
            sampled_record.aa = sampled_record.seq.translate()
            sampled_record.aa_similarity = percent_similarity(ancestor_record.aa, sampled_record.aa)
        sampled_seqs.sort(key=lambda record: record.dna_similarity)
        sampled_ancestors.append(match_summaries(ancestor_record, sampled_seqs[:match_count]))
    return sampled_ancestors

def parse_args():
    def get_seqrecords(fname):
        try:
            return [record for record in SeqIO.parse(fname, "fasta")]
        except IOError, e:
            print 'Error reading input fasta {}'.format(fname)
            raise 

    parser = argparse.ArgumentParser(description="Find the closest N sampled sequences to each inferred ancestor.")
    parser.add_argument(
        'sampled_seqs', type=get_seqrecords,
        help='Fasta containing sampled sequences to match ancestors with.')
    parser.add_argument(
        'asr_seqs', type=get_seqrecords,
        help='Fasta containing the set of inferred ancestral sequences.')
    parser.add_argument(
        'outfile', type=str,
        help='Name of JSON file to write the results to.')
    parser.add_argument(
        '--number-of-matches', type=int,
        default=3,
        help="Report this many of the closest sampled sequences for each inferred ancestor.")
    return parser.parse_args()

def main():
    args = parse_args()
    sampled_ancestors = get_sampled_ancestors(args.sampled_seqs, args.asr_seqs, args.number_of_matches)
    with open(args.outfile, 'w') as fh:
        json.dump(sampled_ancestors, fh)

if __name__ == "__main__":
    main()

