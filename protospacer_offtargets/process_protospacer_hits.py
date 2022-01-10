import os
import sys
import json
import optparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv

def process_protospacer_data(protospacer_data, output_filename):
    overlap_counts = 0
    sp_counts = {0:0,1:0,2:0,3:0}
    nme_counts = {0:0,1:0,2:0,3:0}
    spacer_hit_data = {}

    columns = ['Site', 'Sp-MM0', 'Sp-MM1', 'Sp-MM2', 'Sp-MM3', 'Nme-MM0', 'Nme-MM1', 'Nme-MM2', 'Nme-MM3']
    spacer_hit_data['CGCAAAGCTGCATCCACCCCCCG'] = [0,0,0,0]

    with open(output_filename, 'w') as f:
        w = csv.writer(f, delimiter=',')
        w.writerow(columns)
        for site, hit_info in protospacer_data.items():
            for i, pspacer in enumerate(hit_info['protospacer_matches']):
                if not pspacer in spacer_hit_data.keys():
                    spacer_hit_data[pspacer] = [0,0,0,0]
                mismatches = hit_info['mismatches'][i]
                spacer_hit_data[pspacer][mismatches] = spacer_hit_data[pspacer][mismatches] + 1

        print(spacer_hit_data)
        for spacer,hits in spacer_hit_data.items():
            if len(spacer) == 20:
                continue
            for second_spacer, second_hits in spacer_hit_data.items():
                if second_spacer == spacer[3:] and spacer != second_spacer:
                    all_hits = [spacer] + second_hits + hits
                    w.writerow(all_hits)
        
def main():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--data-filename', action = 'store', dest = 'data_filename', help = "Filename of protospacer hits")
    parser.add_option('-o', '--output-filename', action = 'store', dest = 'output_filename', help = "Filename for outptut")

    (options, args) = parser.parse_args()

    if not options.data_filename or not options.output_filename:
        print("Please specify the protospacer hits and output filenames.")
        parser.print_help()
        sys.exit(2)

    data_filename = options.data_filename
    output_filename = options.output_filename
    protspacer_data = {}

    with open(data_filename) as f:
        protospacer_data = json.load(f)

    process_protospacer_data(protospacer_data, output_filename)
    

if __name__ == "__main__":
    main()
