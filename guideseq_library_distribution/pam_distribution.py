import os
import sys
import optparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import csv
import re

prog = re.compile(r"[^N]([N]{15})[^N][A-Z-]*[^N](N{7})[^N]")

AMPLICONS = ["AAAGNNNNNNNNNNNNNNNACGCAATGGCCCAGACTGAGCACGTGANNNNNNNTTAAGCCAGCCCCGAC"]

# Scoring matrix for alignment
MATCH_DIC = {}
MATCH_DIC[('A', 'A')] = 2
MATCH_DIC[('A', 'N')] = 2
MATCH_DIC[('A', 'G')] = -1
MATCH_DIC[('A', 'C')] = -1
MATCH_DIC[('A', 'T')] = -1
MATCH_DIC[('G', 'A')] = -1
MATCH_DIC[('G', 'N')] = 2
MATCH_DIC[('G', 'G')] = 2
MATCH_DIC[('G', 'C')] = -1
MATCH_DIC[('G', 'T')] = -1
MATCH_DIC[('C', 'A')] = -1
MATCH_DIC[('C', 'N')] = 2
MATCH_DIC[('C', 'G')] = -1
MATCH_DIC[('C', 'C')] = 2
MATCH_DIC[('C', 'T')] = -1
MATCH_DIC[('T', 'A')] = -1
MATCH_DIC[('T', 'N')] = 2
MATCH_DIC[('T', 'G')] = -1
MATCH_DIC[('T', 'C')] = -1
MATCH_DIC[('T', 'T')] = 2

MATCH_DIC[('R','A')] = 2
MATCH_DIC[('R','C')] = -1
MATCH_DIC[('R','G')] = 2
MATCH_DIC[('R','T')] = -1

MATCH_DIC[('Y','A')] = -1
MATCH_DIC[('Y','C')] = 2
MATCH_DIC[('Y','G')] = -1
MATCH_DIC[('Y','T')] = 2

MATCH_DIC[('K','A')] = -1
MATCH_DIC[('K','C')] = -1
MATCH_DIC[('K','G')] = 2
MATCH_DIC[('K','T')] = 2

MATCH_DIC[('M','A')] = 2
MATCH_DIC[('M','C')] = 2
MATCH_DIC[('M','G')] = -1
MATCH_DIC[('M','T')] = -1

MATCH_DIC[('S','A')] = -1
MATCH_DIC[('S','C')] = 2
MATCH_DIC[('S','G')] = 2
MATCH_DIC[('S','T')] = -1

MATCH_DIC[('S','A')] = -1
MATCH_DIC[('S','C')] = 2
MATCH_DIC[('S','G')] = 2
MATCH_DIC[('S','T')] = -1

MATCH_DIC[('W','A')] = 2
MATCH_DIC[('W','C')] = -1
MATCH_DIC[('W','G')] = -1
MATCH_DIC[('W','T')] = 2

MATCH_DIC[('B','A')] = -1
MATCH_DIC[('B','C')] = 2
MATCH_DIC[('B','G')] = 2
MATCH_DIC[('B','T')] = 2

MATCH_DIC[('D','A')] = 2
MATCH_DIC[('D','C')] = -1
MATCH_DIC[('D','G')] = 2
MATCH_DIC[('D','T')] = 2

MATCH_DIC[('H','A')] = 2
MATCH_DIC[('H','C')] = 2
MATCH_DIC[('H','G')] = -1
MATCH_DIC[('H','T')] = 2

MATCH_DIC[('V','A')] = 2
MATCH_DIC[('V','C')] = 2
MATCH_DIC[('V','G')] = 2
MATCH_DIC[('V','T')] = -1

MATCH_DIC[('N','N')] = 2

def get_pam_distribution(starting_library_filename, pams_data):
    pam_distribution = {}
    with open(starting_library_filename) as f:
        for record in SeqIO.parse(f, 'fastq'):
            max_score = 0
            best_alignment = None
            for amplicon in AMPLICONS:
                alignments = pairwise2.align.localds(record.seq, amplicon, MATCH_DIC, -0.5, -0.1)
                for alignment in alignments:
                    if alignment.score > max_score:
                         best_alignment = alignment
                         max_score = alignment.score
            PAM_index = prog.search(str(best_alignment.seqB)).span(2)[0]
            PAM = best_alignment.seqA[PAM_index:PAM_index+7]
            if max_score < 125:
                continue
            if (not PAM.startswith('ACG') and not PAM.startswith('CAT')) or '-' in PAM:
                continue
            if PAM not in pam_distribution.keys():
                pam_distribution[PAM] = {'count':1}
            else:
                pam_distribution[PAM]['count'] = pam_distribution[PAM]['count'] + 1
    return pam_distribution

def save_pam_distribution(output_filename, pam_distribution):
    w = csv.writer(open(output_filename, "w"))
    for key, value in pam_distribution.items():
        w.writerow([key, value['count']])

def parse_pamsfile(pams_filename):
    pams_data = []
    with open(pams_filename) as f:
        reader = csv.DictReader(f, delimiter = '\t')
        for line in reader:
            pams_data.append(line)
    return pams_data

def main():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--starting-library', action = 'store', dest = 'starting_library', help = "Filename of starting library fastq")
    parser.add_option('-o', '--output', action = 'store', dest = 'output_filename', help = "Desired output filename")
    parser.add_option('-p', '--pams', action = 'store', dest = 'pams_filename', help = "Filename of pams file")

    (options, args) = parser.parse_args()

    if not options.starting_library or not options.output_filename:
        print('Please specify the starting library fastq filename and desired output filename.')
        parser.print_help()
        sys.exit(2)

    starting_library_fastq_filename = options.starting_library
    output_filename = options.output_filename

    pams_filename = options.pams_filename
    pams_data = parse_pamsfile(pams_filename)
    pam_distribution = get_pam_distribution(starting_library_fastq_filename, pams_data)
    save_pam_distribution(output_filename, pam_distribution)

if __name__ == '__main__':
    main()
