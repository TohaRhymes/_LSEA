#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import csv
import json
import shutil
from collections import defaultdict, Counter
from typing import Dict, List

from scipy.stats import hypergeom
from utils import get_overlapping_features, count_intervals, log_message, check_and_create_dir, \
    get_filename_without_extension
from statsmodels.stats.multitest import fdrcorrection


def run_plink_clumping(plink_path: str,
                       bfile_path: str,
                       tsv_file: str,
                       input_dict: Dict,
                       p1: float,
                       p2: float,
                       r2: float,
                       kb: float,
                       out_name: str) -> str:
    """
    XXXXXXXXXXXXXXXXXX

    :param plink_path:
    :param bfile_path:
    :param tsv_file:
    :param input_dict:
    :param p:
    :param out_name:
    :return:
    """
    tsv_plink = os.path.join(out_name,
                             f'{get_filename_without_extension(tsv_file)}_for_plink.tsv')
    out_plink = os.path.join(out_name,
                             get_filename_without_extension(tsv_file))
    with open(tsv_plink, 'w', newline='') as csvfile:
        tsv_writer = csv.writer(csvfile, delimiter='\t')
        header = ["SNP", "Chr", "Pos", "P"]
        tsv_writer.writerow(header)
        for snp, values in input_dict.items():
            # todo modify and adjust
            # Only extract necessary info from dictionary corresponding to csv
            snp_info_row = [snp] + values[:3]
            tsv_writer.writerow(snp_info_row)
    # todo more params, instead of hardcode: p1, p2, r2, kb
    # The clumping procedure takes all SNPs that are significant at threshold p1 that have not already been clumped
    # (denoting these as index SNPs) and forms clumps of all other SNPs that are within a certain kb distance from the
    # index SNP (default 250kb) and that are in linkage disequilibrium with the index SNP,
    # based on an r-squared threshold (default 0.50). These SNPs are then subsetted based on the result for that SNP,
    # as illustrated below. This is a greedy algorithm and so each SNP will only appear in a single clump, if at all.
    subprocess.call(
        f'{os.path.join(plink_path, "plink")} '
        f'--bfile \"{bfile_path}\" '
        f'--clump \"{tsv_plink}\" '
        f'--clump-field P '
        f'--clump-p1 {p1} '  # Significance threshold for index SNPs
        f'--clump-p2 {p2} '  # Secondary significance threshold for clumped SNPs
        f'--clump-r2 {r2} '  # LD threshold for clumping
        f'--clump-kb {kb}'  # Physical distance threshold for clumping
        f'--clump-snp-field SNP '
        f'--out \"{out_plink}\" '
        f'--allow-no-sex '
        f'--allow-extra-chr '
        f'2> \"{out_name}/PLINK_clumping.log\"',
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out_plink + ".clumped"


def get_snp_info_from_tsv(tsv_file: str, names: List[str]) -> Dict[str, List[str]]:
    """
    XXXXXXXXXXXXXXXXXXXXXX

    :param tsv_file:
    :param names:
    :return:
    """
    input_dict = defaultdict(list)
    with open(tsv_file, 'r', newline='') as csvfile:  # Convert tsv to dictionary
        snp_reader = csv.reader(csvfile, delimiter='\t')
        csv_headings = next(snp_reader)
        try:  # We find necessary indices from header
            chr_index = csv_headings.index(names[0])
            pos_index = csv_headings.index(names[1])
            id_index = csv_headings.index(names[2])
            p_index = csv_headings.index(names[3])
        except ValueError:
            raise ValueError("Check that your tsv file has headers and they are correct!")
        for snp_info_row in snp_reader:  # Start from second row
            try:
                # Name -> Chr, pos, pval
                input_dict[snp_info_row[id_index]].extend(snp_info_row[chr_index])
                input_dict[snp_info_row[id_index]].extend(snp_info_row[pos_index])
                input_dict[snp_info_row[id_index]].extend(snp_info_row[p_index])
            except IndexError:
                raise IndexError(f"Error in string: {snp_info_row}. It s")
    return dict(input_dict)


def make_bed_file(clumped_file, interval, out_name, output_merged_file):
    """
    XXXXXXXXXXXXXXXXXXXXXXXXXXX

    :param output_merged_file:
    :param clumped_file:
    :param interval:
    :param out_name:
    :return:
    """
    clumps_file = os.path.join(out_name, "clumps.bed")
    sorted_file = os.path.join(out_name, "clumps_sorted.bed")
    merged_file = os.path.join(out_name, "merged.bed")
    merged_fixed_size_file = os.path.join(out_name, "merged_fixed_size.bed")
    with open(clumps_file, 'w', newline='') as bed_file:  # Here we write to new file
        clumps_writer = csv.writer(bed_file, delimiter='\t')
        with open(clumped_file, 'r') as cl_file:  # Our result of clumping (SNPs sets)
            clumped_reader = csv.reader(cl_file, delimiter='\t')
            for clump_info in clumped_reader:
                if len(clump_info) != 0:
                    # todo What to do with NONE?
                    clump_info = list(filter(lambda x: len(x) != 0 and x != " ", clump_info[0].split(" ")))
                    if clump_info[0].lower() == "chr":  # Skip first row
                        continue
                    bed_row = [clump_info[0], max(int(clump_info[3]) - interval, 0), int(clump_info[3]) + interval]
                    clumps_writer.writerow(bed_row)

    # Sort file
    subprocess.call(
        f"bedtools sort -i \"{clumps_file}\" > \"{sorted_file}\"", shell=True)
    # Merge regions
    subprocess.call(
        f"bedtools merge -i \"{sorted_file}\" > \"{merged_file}\"", shell=True)
    with open(merged_fixed_size_file, 'w', newline='') as bed_file:  # Here we write to new file
        fixed_merged_writer = csv.writer(bed_file, delimiter='\t', lineterminator='\n')
        with open(merged_file, 'r', newline='') as inter:
            merged_reader = csv.reader(inter, delimiter='\t')
            for row in merged_reader:
                middle_point = int(row[1]) + (int(row[2]) - int(row[1])) // 2
                new_row = [row[0], max(middle_point - interval, 0), middle_point + interval]
                fixed_merged_writer.writerow(new_row)
    # Add numeration to intervals
    subprocess.call(
        "awk {'print $0\"\t\"FNR'}" + f" \"{merged_fixed_size_file}\" > \"{output_merged_file}\"",
        shell=True)


def p_val_for_gene_set(n_big,
                       k_big,
                       n,
                       k):
    """

    :param n_big: features (genes)/intervals in total
    :param k_big: features (genes) from current set
    :param n: causal features (genes) in total
    :param k: causal features (genes) from current set
    :return:
    """
    return hypergeom.sf(k - 1, n_big, k_big, n)


def calculate_qvals(pvals):
    # todo other fdr corrections (?)
    return list(fdrcorrection(pvals)[1])


def count_lines(filename):
    """Count the number of lines in a given file."""
    with open(filename, 'r') as file_to_count:
        line_count = sum(1 for _ in file_to_count)
    return line_count


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LSEA')
    parser.add_argument('--input',
                        '-i',  # was: -i
                        help='Input file in TSV format. The file should contain at least four columns: '
                             'chromosome name, variant position, variant ID, and GWAS p-value. '
                             'NB: The file MUST have a header row.',
                        metavar='file',
                        type=str,
                        required=True)
    parser.add_argument('--column_names',
                        '-cols',  # was: -column_names
                        help='Column names for input TSV.'
                             'Names should be given in this order, separated by space: '
                             'chromosome, position, ID, p-value',
                        metavar='chr pos id p',
                        nargs=4,
                        type=str,
                        default=["chr", "pos", "id", "p"],
                        required=False)
    parser.add_argument('--universe',
                        '-u',  # was: -universe
                        help='JSON file(s) with universe. Several universe files may be specified separated by space',
                        metavar='path',
                        nargs='+',  # >=1 args
                        type=str,
                        required=True)
    parser.add_argument('--out',
                        '-o',  # was: -out
                        help='Relative path to output directory (default: lsea_result). '
                             'Will be created.',
                        metavar='name',
                        type=str,
                        required=False,
                        default='lsea_result')
    parser.add_argument('--clump_p1',
                        '-с_p',  # was: -p
                        help='p-value cutoff to be used when identifying associated loci '
                             '(corresponds to `--clump-p1` flag in PLINK).'
                             'If not specified, optimal cutoff will be estimated using a regression model '
                             '(from default: 1e-5 and 5e-8).',  # todo CHECKKK THIS IS TRUE?
                        metavar='float',
                        type=float,
                        nargs='+',
                        required=False,
                        default=['1e-5', '5e-8'])
    parser.add_argument('--clump_p2',
                        '-с_p2',
                        help='Secondary significance threshold for clumped SNPs'
                             '(corresponds to `--clump-p2` flag in PLINK).'
                             'Default: 0.01',
                        metavar='float',
                        type=float,
                        required=False,
                        default=0.01)
    parser.add_argument('--clump_r2',
                        '-с_r2',
                        help='LD threshold for clumping'
                             '(corresponds to `--clump-r2` flag in PLINK).'
                             'Default: 0.1',
                        metavar='float',
                        type=float,
                        required=False,
                        default=0.1)
    parser.add_argument('--clump_kb',
                        '-с_kb',
                        help='Physical distance threshold for clumping'
                             '(corresponds to `--clump-kb` flag in PLINK).'
                             'Default: 500',
                        metavar='int',
                        type=int,
                        required=False,
                        default=500)
    parser.add_argument('--plink_dir',
                        '-pd',  # was: -plink_dir
                        help='Path to a directory with PLINK executable',
                        metavar='dir',
                        type=str,
                        required=False)
    parser.add_argument('--bfile',
                        '-bf',  # was: -bfile
                        help='Genotypes in PLINK binary (bed+bim+fam) format that will be used for LD-based clumping. '
                             'Only file prefix should be given.',
                        metavar='prefix',
                        type=str,
                        required=False)
    parser.add_argument('--qval_threshold',
                        '-qt',  # was: -qval_threshold
                        help='Q-value threshold for output (default: 0.1)',
                        metavar='float',
                        type=float,
                        required=False,
                        default=0.05)
    parser.add_argument('--interval_count_threshold',
                        '-ict',
                        help='XXXXXXX',  # todo
                        metavar='float',
                        type=int,
                        required=False, default=3)
    parser.add_argument('--print_all',
                        '-a',
                        help='If set, save all pathways to resulting file, otherwise, only significant '
                             '(passing `qval_threshold` and `interval_count_threshold`)',
                        action='store_true',
                        required=False,
                        default=False)
    # todo add temp files deletion

    args = parser.parse_args()

    full_command = 'Running LSEA with the following CMD:\n'
    for arg in vars(args):
        value = getattr(args, arg)
        if value is not None:
            full_command += f"--{arg} {value}"
    log_message(full_command)

    tsv_file = args.input
    path_to_plink_dir = args.plink_dir
    if path_to_plink_dir is not None:
        path_to_plink_dir = os.path.normpath(path_to_plink_dir)
    path_to_bfile = args.bfile
    qval_thresh = args.qval_threshold
    interval_thresh = args.interval_count_threshold
    col_names = args.column_names
    if col_names is None:
        col_names = ["chr", "pos", "id", "p"]
    json_files = args.universe
    p_cutoffs = [float(x) for x in args.clump_p1]
    p2 = args.clump_p2
    r2 = args.clump_r2
    kb = args.clump_kb
    out_name = args.out
    check_and_create_dir(out_name)

    log_message(f'Reading input file {tsv_file}...')
    snp2chr_pos_pval = get_snp_info_from_tsv(tsv_file, col_names)

    for universe_file in json_files:
        universe_name = os.path.basename(universe_file).replace('.json', '')
        log_message(f'Processing universe: {universe_name}')
        universe = json.load(open(universe_file, "r"))
        interval = universe["interval"]
        with open(os.path.join(out_name, 'features.bed'), 'w') as feature_file:
            for feature in universe["features"]:
                bed_line = '\t'.join(universe["features"][feature])
                bed_line = bed_line.replace('chr', '')
                feature_file.write(f'{bed_line}\n')
        interval_counts_for_universe = universe["interval_counts"]
        log_message(
            f"""Universe stats:
            interval size = {interval};
            interval count = {universe["universe_intervals_number"]};
            total number of features = {len(universe["features"])};
            feature sets in universe = {len(universe["gene_set_dict"])}""")

        with open(os.path.join(out_name, f"annotation_stats_{universe_name}.tsv"), 'w', newline='') as stats_file:
            stats_writer = csv.writer(stats_file, delimiter='\t')
            stats_header = ['p_cutoff',
                            'num_loci',
                            'annotated_loci',
                            'unambiguous_annotations',
                            'significant_hits',
                            'min_qval']
            stats_writer.writerow(stats_header)

            for p_cutoff in p_cutoffs:
                output_merged_file = os.path.join(out_name, "merged_with_line_numbers.bed")
                features_file = os.path.join(out_name, "features.bed")
                inter_file = os.path.join(out_name, "inter.tsv")

                log_message(f'Calculating enrichment with p-value cutoff = {p_cutoff}')
                clumped_file = run_plink_clumping(path_to_plink_dir,
                                                  path_to_bfile,
                                                  tsv_file,
                                                  snp2chr_pos_pval,
                                                  p_cutoff,
                                                  p2,
                                                  r2,
                                                  kb,
                                                  out_name)
                make_bed_file(clumped_file, interval, out_name, output_merged_file)
                n_intervals = count_lines(output_merged_file)
                target_features = get_overlapping_features(output_merged_file,
                                                           features_file,
                                                           inter_file)  # return gene->[chr:s-e, chr:s-e, ...]
                feature_set = universe["gene_set_dict"]
                interval_counts = count_intervals(feature_set, target_features)

                pvals = []
                # todo CHECKKKK тут в двуз местах interval_counts[key] ТОЧНО ЛИ ВОЗВРАЩАЕТ уже длину ???
                # поглядеть как раньше
                for w in sorted(interval_counts, key=lambda item: interval_counts[item], reverse=True):
                    pvals.append(p_val_for_gene_set(universe["universe_intervals_number"],  # intervals in universe
                                                    interval_counts_for_universe[w],  # intervals in this set
                                                    n_intervals,  # significant clumps
                                                    interval_counts[w]  # significant clumps from this set
                                                    )
                                 )
                # todo instead of fdrcorrection, betted fdr for correlated features (???)
                qvals = calculate_qvals(pvals)

                with open(os.path.join(out_name, f"{universe_name}_result_{p_cutoff}.tsv"), 'w', newline='') as file:
                    explained_loci = set()
                    feature_names = defaultdict(set)
                    hit_count = 0
                    min_qval = 1
                    result_writer = csv.writer(file, delimiter='\t')
                    result_writer.writerow(["gene_set", "overlapping_loci", "p_value", "q_value", "description"])
                    for i, w in enumerate(
                            sorted(interval_counts, key=lambda item: len(interval_counts[item]), reverse=True)):
                        # todo почему ТРИ?? что за трешхолд такой может тоже в параметры вынести?
                        if len(interval_counts[w]) >= interval_thresh and qvals[i] <= qval_thresh or args.print_all:
                            min_qval = min(min_qval, qvals[i])
                            hit_count += 1
                            for locus in interval_counts[w]:
                                explained_loci.add(locus)
                                feature_names[locus].add(';'.join(interval_counts[w][locus]))
                            row = [w, len(interval_counts[w]), pvals[i], qvals[i], dict(interval_counts[w])]
                            result_writer.writerow(row)
                    # Calculating number of loci with ambiguous annotations
                    # todo а это вообще нужно?
                    # это же типа саммари для текущего пвала
                    # может както более разумно переименовать?
                    ambiguous_loci = 0
                    for _, feature_set in feature_names.items():
                        if len(feature_set) > 1:
                            expanded_set = sum([feature.split(';') for feature in feature_set], [])
                            feature_count = dict(
                                Counter(expanded_set))
                            if max(feature_count.values()) < len(feature_set):
                                ambiguous_loci += 1
                    unambiguous = len(explained_loci) - ambiguous_loci
                    stats_row = [p_cutoff, n_intervals, len(explained_loci), unambiguous, hit_count, min_qval]
                    stats_writer.writerow(stats_row)
