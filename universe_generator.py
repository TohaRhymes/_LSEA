import argparse
import sys
import csv
import json
import subprocess
from utils import get_snp_locations, get_overlapping_features, count_intervals, get_features_from_dir, read_features, read_gmt

def create_universe(dict, interval):
    with open("./tmp.bed", 'w', newline='') as bed_file:  # Here we write to new file
        my_writer = csv.writer(bed_file, delimiter='\t')
        id = 1
        for snp in dict:
            info = dict[snp]
            # Chromosome, start and end coordinates, id of an interval
            start = max(0, int(info[1]) - interval)
            bed_row = [info[0], start,
                       int(info[1]) + interval, id]
            id += 1
            my_writer.writerow(bed_row)
        subprocess.call("sort -k1,1 -k2,2n tmp.bed > universe.bed", shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Universe generator for LSEA')
    parser.add_argument('-variants', help='Path to TSV file with variant coordinates that will be used for universe generation. The file whould contain at least three columns: chromosome, position, and variant ID (must be unique). All other columns ar ignored.',
                        metavar='path', type=str, required=True)
    parser.add_argument('-interval', help='Size of the window around each target variant (Default: 500000)',
                        metavar='int', type=int, required=False, default=500000)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-features', help='Path to files with all feature annotations (one feature per line) in BED format and feature set description in GMT format. Two file paths separated by space should be provided in the following order: [BED] [GMT]. Cannot be used with -feature_files_dir',
                        metavar='path', nargs=2)
    group.add_argument('-feature_files_dir', help='A directory with feature files, one file per feature set, in BED format. File name will be used a the name of each feature set. Cannot be used with -features',
                        metavar='path', type=str)
    parser.add_argument('-o', help='Output path for json',
                        metavar='path', type=str, required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    variants = args.variants
    interval = args.interval
    if args.features is not None:
        bed, gmt = args.features
        gene_set_dict = read_gmt(gmt)
    else:
        feature_dir = args.feature_files_dir
        gene_set_dict = get_features_from_dir(feature_dir)
        bed = 'features.bed'
    out_path = args.o
    
    input_dict = get_snp_locations(variants)
    create_universe(input_dict, interval)
    features_in_universe = get_overlapping_features("./universe.bed", bed, "inter2.tsv")
    interval_counts_for_universe = count_intervals(gene_set_dict, features_in_universe, emit_raw=False)

    feature_dict = read_features(bed)
    out_dict = {"interval": interval,
            "universe_intervals_number" : len(input_dict), 
            "interval_counts" : interval_counts_for_universe, 
            "gene_set_dict" : gene_set_dict,
            "features" : feature_dict}

    base = open(out_path, "w")
    json.dump(out_dict, base)
    base.close()
