import os
import sys
import optparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import json

def parse_protospacers(protospacers_filename):
    protospacers = []
    with open(protospacers_filename) as f:
        for line in f:
            protospacer = line.strip()
            protospacer_len = len(protospacer)
            protospacers.append((protospacer, protospacer_len))
    return protospacers

def hamming_distance(s1, s2):
    distance = 0
    if len(s1) != len(s2):
        return 999
    for i in range(len(s1)):
        try:
            if s1[i] != s2[i]:
                distance += 1
        except IndexError:
           return 999 
    return distance

def find_genome_hits(protospacers, genome_filename, depth = 3):
    hits = {}
    with open(genome_filename) as f:
        for record in SeqIO.parse(f, 'fasta'):
            print(record)
            for i in range(len(record.seq)):
                substring = str(record.seq[i:i+23]).upper()
                substring_reverse_complement = str(record.seq[i:i+23].reverse_complement()).upper()
                if 'N' in substring:
                    continue
                for protospacer, plen in protospacers:
                    sp_spacer = protospacer[3::]
                    sp_substring = substring[3::]
                    sp_substring_reverse_complement = substring_reverse_complement[3::]
                    hd = hamming_distance(protospacer, substring)
                    sp_hd = hamming_distance(sp_spacer, sp_substring)
                    hd_rev_comp = hamming_distance(protospacer, substring_reverse_complement)
                    sp_hd_rev_comp = hamming_distance(sp_spacer, sp_substring_reverse_complement)
                    if hd <= depth:
                        if substring not in hits:
                            hits[substring] = {'protospacer_matches': [], 'mismatches': [], 'count': 0, 'sp_count': 0, 'nme_count':0}
                        hits[substring]['protospacer_matches'].append(protospacer)
                        hits[substring]['count'] += 1
                        hits[substring]['nme_count'] += 1
                        hits[substring]['mismatches'].append(hd)
                        print('Positive strand hit')
                        print(hits[substring])
                    if sp_hd <= depth:
                        if substring not in hits:
                            hits[substring] = {'protospacer_matches':[], 'mismatches': [], 'count':0, 'sp_count':0, 'nme_count':0}
                        hits[substring]['protospacer_matches'].append(sp_spacer)
                        hits[substring]['count'] += 1
                        hits[substring]['sp_count'] += 1
                        hits[substring]['mismatches'].append(sp_hd)
                        print('Positive strand hit')
                        print(hits[substring])
                    if hd_rev_comp <= depth:
                        if substring_reverse_complement not in hits:
                            hits[substring_reverse_complement] = {'protospacer_matches': [], 'mismatches': [], 'count': 0, 'sp_count': 0, 'nme_count':0}
                        hits[substring_reverse_complement]['protospacer_matches'].append(protospacer)
                        hits[substring_reverse_complement]['count'] += 1
                        hits[substring_reverse_complement]['nme_count'] += 1
                        hits[substring_reverse_complement]['mismatches'].append(hd_rev_comp)
                        print('Negative strand hit')
                        print(hits[substring_reverse_complement])
                    if sp_hd_rev_comp <= depth:
                        if substring_reverse_complement not in hits:
                            hits[substring_reverse_complement] = {'protospacer_matches': [], 'mismatches': [], 'count': 0, 'sp_count': 0, 'nme_count':0}
                        hits[substring_reverse_complement]['protospacer_matches'].append(sp_spacer)
                        hits[substring_reverse_complement]['count'] += 1
                        hits[substring_reverse_complement]['sp_count'] += 1
                        hits[substring_reverse_complement]['mismatches'].append(sp_hd_rev_comp)
                        print('Negative strand hit')
                        print(hits[substring_reverse_complement])

    return hits

def main():
    parser = optparse.OptionParser()
    parser.add_option('-g', '--genome-filename', action = 'store', dest = 'genome_filename', help = "Filename of genome")
    parser.add_option('-p', '--protospacers', action = 'store', dest = 'protospacers_filename', help = "Filename of protospacers")
    parser.add_option('-o', '--output-filename', action = 'store', dest = 'output_filename', help = "Filename for output")

    (options, args) = parser.parse_args()

    if not options.genome_filename or not options.protospacers_filename or not options.output_filename:
        print("Please specify the genome, protospacers, and output filenames.")
        parser.print_help()
        sys.exit(2)

    genome_filename = options.genome_filename
    protospacers_filename = options.protospacers_filename
    output_filename = options.output_filename

    protospacers = parse_protospacers(protospacers_filename)
    hits = find_genome_hits(protospacers, genome_filename, depth = 3)

    with open(output_filename, 'w', encoding='utf-8') as f:
        json.dump(hits, f, ensure_ascii=False, indent = 4)

if __name__ == '__main__':
    main()
