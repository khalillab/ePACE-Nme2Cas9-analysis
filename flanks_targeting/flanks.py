import sys
import requests
import csv
import optparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

VEP_URL = "https://rest.ensembl.org/sequence/region/human/{0}"
FLANK_SIZE = "250"

def parse_snp_report(filename, output_filename):
    with open(filename) as f:
        with open(output_filename, 'w') as output_handle:
            data_reader = csv.DictReader(f, dialect='excel-tab')
            count = 0
            for line in data_reader:
                chromosome = line['#chr']
                pos = line['pos']
                region = chromosome + ":" + pos + ".." + pos + ":1"
                headers= {"Content-type":"text/plain", "Accept": "text/plain"}
                query_params = {"format": "fasta", "expand_3prime": FLANK_SIZE, "expand_5prime": FLANK_SIZE}
                response = requests.get(VEP_URL.format(region),params=query_params,headers=headers)
                record = SeqRecord(Seq(response.text), id=line['snp_id'], description='dbSNP', annotations={'variation':line['variation']})
                SeqIO.write(record, output_handle, "fasta")
                count += 1
                if count % 100 == 0:
                    print(str(count) + " Records Fetched")

def main():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--snp-filename', action = 'store', dest = 'filename', help = "Filename of file containing dbSNP entries downloaded as tsv.")
    parser.add_option('-o', '--output-filename', action = 'store', dest = 'output_filename', help = "Desired output filename for flanking sequence fasta file.")

    (options, args) = parser.parse_args()

    if not options.filename or not options.output_filename:
        print('Please specify the input and output filenames.')
        parser.print_help()
        sys.exit(2)

    filename = options.filename
    output_filename = options.output_filename

    parse_snp_report(filename, output_filename)

if __name__ == '__main__':
    main()
