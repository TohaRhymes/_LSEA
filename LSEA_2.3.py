#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import csv
import json
import shutil
from collections import defaultdict
from typing import Counter
from scipy.stats import hypergeom
from rpy2 import robjects
from utils import get_overlapping_features, count_intervals


def tsv_len(tsv_file):
    with open(tsv_file) as f:
        for i, l in enumerate(f):
            pass
    return i


def run_plink(plink_path, bfile_path, tsv_file, input_dict, p, out_name):
#    print('Calculating independent loci with PLINK')
    tsv_plink = os.getcwd() + f"/{out_name}/" + \
        tsv_file.split('/')[-1].split('.')[0] + '_for_plink.tsv'
    out_plink = os.getcwd() + f"/{out_name}/" + \
        tsv_file.split('/')[-1].split('.')[0]
    with open(tsv_plink, 'w', newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter='\t')
        init_row = ["SNP", "Chr", "Pos", "P"]
        my_writer.writerow(init_row)
        for snp in input_dict:
            # Only extract necessary info from dictionary corresponding to csv
            row = [snp] + input_dict[snp][0:3]
            my_writer.writerow(row)
    subprocess.call('{0}/plink --bfile {1} --clump {2} --clump-field P --clump-p1 {3} --clump-p2 0.01 --clump-r2 0.1 --clump-snp-field SNP --clump-kb 500 --out {4} --allow-no-sex --allow-extra-chr \
      2> {5}'.format(plink_path, bfile_path, tsv_plink, p, out_plink, f'./{out_name}/PLINK_clumping.log'), 
      shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#    subprocess.call(f'md5sum {out_plink}.clumped', shell=True)
    return out_plink + ".clumped"


def normalize_path(path):
    if path[-1] == "/":
        path = path[:-1]
    return path


def get_snp_info(tsv_file, names):
    input_dict = defaultdict(list)
    names = list(map(lambda x: x.lower(), names))
    with open(tsv_file, 'r', newline='') as csvfile:  # Convert tsv to dictionary
        my_reader = csv.reader(csvfile, delimiter='\t')
        csv_headings = list(map(lambda x: x.lower(), next(my_reader)))
        try:  # We find necessary indices from header
            chr_index = csv_headings.index(names[0])
            pos_index = csv_headings.index(names[1])
            id_index = csv_headings.index(names[2])
            p_index = csv_headings.index(names[3])
        except ValueError:
            print("Check that your tsv file has headers and they are correct!")
            exit(1)
        for row in my_reader:  # Start from second row
            input_dict[row[id_index]] = [row[chr_index]] + \
                [row[pos_index]] + [row[p_index]]  # Name -> Chromosome, pos, pval
    return input_dict


def make_bed_file(clumped_file, interval, out_name):
#    print(f'Reading clumped file {clumped_file}')
    with open(f"./{out_name}/clumps.bed", 'w', newline='') as bed_file:  # Here we write to new file
        my_writer = csv.writer(bed_file, delimiter='\t')
        with open(clumped_file, 'r') as cl_file:  # Our result of clumping (SNPs sets)
            my_reader = csv.reader(cl_file, delimiter='\t')
            for row in my_reader:
#                print('Reading line')
                if len(row) != 0:
                    # What to do with NONE?
                    row = list(filter(lambda x: len(x) !=
                                      0 and x != " ", row[0].split(" ")))
                    if row[0] == "CHR":  # Skip first row
                        continue
#                    print('Got here')
                    bed_row = [row[0], max(int(row[3]) - interval, 0), int(row[3]) + interval]
                    my_writer.writerow(bed_row)
    # Sort file
    subprocess.call(
        f"bedtools sort -i ./{out_name}/clumps.bed > ./{out_name}/clumps_sorted.bed", shell=True)
    # Merge regions
    subprocess.call(
        f"bedtools merge -i ./{out_name}/clumps_sorted.bed > ./{out_name}/merged.bed", shell=True)
    with open(f"./{out_name}/merged_fixed_size.bed", 'w', newline='') as bed_file:  # Here we write to new file
        my_writer = csv.writer(bed_file, delimiter='\t', lineterminator='\n')
        with open(f'./{out_name}/merged.bed', 'r', newline='') as inter:
            my_reader = csv.reader(inter, delimiter='\t')
            for row in my_reader:
                middle_point = int(row[1]) + (int(row[2]) - int(row[1])) // 2
                new_row = [row[0], max(middle_point - interval, 0), middle_point + interval]
                my_writer.writerow(new_row)
    # Add numeration to intervals
    subprocess.call(
        "awk {'print $0\"\t\"FNR'}" + f" ./{out_name}/merged_fixed_size.bed > ./{out_name}/merged_with_line_numbers.bed", shell=True)


def p_val_for_gene_set(n_big, k_big, n, k):
    return hypergeom.sf(k, n_big, k_big, n)


def calcuate_qvals(pvals):
    x = robjects.r('''
            f <- function(pvals) {
                p.adjust(pvals, method = "fdr")
            }
            ''')
    r_f = robjects.globalenv['f']
    v = robjects.FloatVector(pvals)
    return list(r_f(v))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LSEA')
    parser.add_argument('-input', help='Input file in TSV format. The file should contain at least four columns: chromosome name, variant position, variant ID, and GWAS p-value. The file MUST have a header row.', metavar='file',
                        type=str, required=True)
    parser.add_argument('-universe', help='JSON file with universe. Several universe files may be specified separated by space',
                        metavar='path', nargs='+', type=str, required=True)
    parser.add_argument('-out', help='Relative path to output directory (Default: lsea_result). Will be created.',
                        metavar='name', type=str, required=False, default='lsea_result')
    parser.add_argument('-p', help='p-value cutoff to be used when identifying associated loci. If not specified, optimal cutoff will be estimated using a regression model (at least -n and -m options should be given)',
                        metavar='float', nargs='+', required=False, default=['1e-5', '5e-8'])
    parser.add_argument('-plink_dir', help='Path to a directory with PLINK executable', metavar='dir',
                        type=str, required=False)
    parser.add_argument('-bfile', help='Genotypes in PLINK .bed format that will be used for LD-based clumping. Only file prefix should be given.', metavar='prefix',
                        type=str, required=False)
    parser.add_argument('-qval_threshold', help='Q-value threshold for output (Default: 0.1)',
                        metavar='float', type=float, required=False, default="0.1")
    parser.add_argument('-column_names', help='Column names for input TSV. These names will be used if the column names do not match the default: chr, pos, id, p. Names should be given in the same order (chromosome, position, ID, p-value)',
                        metavar='name', nargs=4, type=str, required=False)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    
    print(f'INFO\tRunning LSEA with the following CMD: {" ".join(sys.argv[:])}')
    
    tsv_file = args.input
    path_to_plink_dir = args.plink_dir
    if path_to_plink_dir is not None:
        path_to_plink_dir = normalize_path(args.plink_dir)
    path_to_bfile = args.bfile
    qval_thresh = args.qval_threshold
    col_names = args.column_names
    json_files = args.universe
    p_cutoffs = [float(x) for x in args.p]
    out_name = args.out
    if col_names is None:
        col_names = ["chr", "pos", "id", "p"]
    print(f'INFO\tReading input file {tsv_file}...')
    input_dict = get_snp_info(tsv_file, col_names)

    if os.path.exists(f'./{out_name}'):
        print(f'WARN\tOutput diretory {out_name} exists, removing...')
        shutil.rmtree(f'./{out_name}')
    os.makedirs(f'./{out_name}')

    for universe_file in json_files:
        universe_name = os.path.basename(universe_file).replace('.json', '')
        print(f'INFO\tProcessing universe {universe_name}')
        universe= json.load(open(universe_file, "r"))
        interval = universe["interval"]
        with open(f'./{out_name}/features.bed', 'w') as feature_file:
            for feature in universe["features"]:
                bed_line = '\t'.join(universe["features"][feature])
                feature_file.write(f'{bed_line}\n')
        interval_counts_for_universe = universe["interval_counts"]
        print(f'INFO\tUniverse stats:')
        print(f'\tinterval size = {interval};\n\tinterval count = {universe["universe_intervals_number"]};')
        print(f'\ttotal number of features = {len(universe["features"])};')
        print(f'\tfeature sets in universe = {len(universe["gene_set_dict"])}')

        stats_file = open(f"./{out_name}/annotation_stats_{universe_name}.tsv", 'w')
        print('p_cutoff\tnum_loci\tannotated_loci\tunambiguoug_annotations\tsignificant_hits\tmin_qval', file=stats_file)

        for p_cutoff in p_cutoffs:
            print(f'INFO\tCalculating enrichment with p-value cutoff = {p_cutoff}')
            clumped_file = run_plink(path_to_plink_dir, path_to_bfile, tsv_file, input_dict, p_cutoff, out_name)
            make_bed_file(clumped_file, interval, out_name)
            n_intervals = sum(1 for line in open(f'./{out_name}/merged_with_line_numbers.bed'))
            target_features = get_overlapping_features(f"./{out_name}/merged_with_line_numbers.bed", 
                        f'./{out_name}/features.bed', f"./{out_name}/inter.tsv")    
            feature_set = universe["gene_set_dict"]
            interval_counts = count_intervals(feature_set, target_features)

            pvals = []
            for w in sorted(interval_counts, key=lambda item: len(interval_counts[item]), reverse=True):
                pvals.append(p_val_for_gene_set(universe["universe_intervals_number"], 
                    interval_counts_for_universe[w], n_intervals, len(interval_counts[w])))
            qvals = calcuate_qvals(pvals)

            with open(f"./{out_name}/{universe_name}_result_{p_cutoff}.tsv", 'w', newline='') as file:
                explained_loci = set()
                feature_names = defaultdict(set)
                hit_count = 0
                min_qval = 1
                my_writer = csv.writer(file, delimiter='\t')
                my_writer.writerow(["Gene_set", "Overlapping_loci", "Description", "p-value", "q-value"])
                for i, w in enumerate(sorted(interval_counts, key=lambda item: len(interval_counts[item]), reverse=True)):
                    if len(interval_counts[w]) >= 3 and qvals[i] <= qval_thresh:
                        min_qval = min(min_qval, qvals[i])
                        hit_count += 1
                        for locus in interval_counts[w]:
                            explained_loci.add(locus)
                            feature_names[locus].add(';'.join(interval_counts[w][locus]))
                        row = [w, len(interval_counts[w]), dict(interval_counts[w]), pvals[i], qvals[i]]
                        my_writer.writerow(row)
                # Calculating number of loci with ambiguous annotations
                ambiguous_loci = 0
                for _, feature_set in feature_names.items():
                    if len(feature_set) > 1:
                        expanded_set = sum([ftr.split(';') for ftr in feature_set], [])
                        feature_freq = dict(Counter(expanded_set))
                        if max(feature_freq.values()) < len(feature_set):
                            ambiguous_loci += 1
                unambiguous = len(explained_loci) - ambiguous_loci
                print(f'{p_cutoff}\t{n_intervals}\t{len(explained_loci)}\t{unambiguous}\t{hit_count}\t{min_qval}', file=stats_file)
        
        stats_file.close()
