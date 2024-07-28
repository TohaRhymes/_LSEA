#!/usr/bin/env python

import argparse
import os
import sys
import csv
import json
import subprocess
from typing import Tuple, Dict

from utils import get_snp_locations, get_overlapping_features, count_intervals, get_features_from_dir, \
    read_features_from_bed, read_gmt, log_message, check_and_create_dir

import time
import random


def current_milli_time():
    return round(time.time() * 1000)


def current_random():
    cur_time = current_milli_time()
    return f"{cur_time}_{random.randint(1, cur_time)}"


def create_universe(snp2chrom_pos: Dict[str, Tuple],
                    interval: int,
                    universe_out: str = 'universe.bed',
                    tmp_file: str = 'tmp.bed') -> None:
    """
    Creates a universe of genomic intervals based on SNP positions and a specified interval size, outputs to a BED file.
    This function writes SNP intervals to a temporary file, sorts them, and outputs to a final BED file.
    The temporary file is deleted after use.


    :param snp2chrom_pos: (Dict[str, Tuple[str, int]]): Dictionary mapping SNP identifiers to tuples of chromosome and position.
    :param interval: (int): The interval size to extend around each SNP position, both upstream and downstream.
    :param universe_out: (str): The output path for the final sorted BED file. Defaults to 'universe.bed'.
    :param tmp_file: (str): The path for the temporary BED file used for sorting. Defaults to 'tmp.bed'.

    """
    with open(tmp_file, 'w', newline='') as bed_file:  # Here we write to new file
        bed_writer = csv.writer(bed_file, delimiter='\t')
        for cur_id, (snp, (chrom, pos)) in enumerate(snp2chrom_pos.items()):
            start = max(0, int(pos) - interval)
            end = int(pos) + interval
            # Chromosome, start, end coordinates, id
            bed_row = [chrom,
                       start,
                       end,
                       cur_id]
            bed_writer.writerow(bed_row)
    try:
        subprocess.call(f"sort -k1,1 -k2,2n {tmp_file} > {universe_out}", shell=True)
    finally:
        # Clean up the temporary file
        os.remove(tmp_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Universe generator for LSEA.')
    parser.add_argument('-v',
                        '--variants',
                        help='Path to TSV file with variant coordinates that will be used for universe generation. '
                             'The file should contain at least three columns: chromosome, position, and variant ID '
                             '(must be unique). All other columns ar ignored.',
                        metavar='path',
                        type=str,
                        required=True)
    parser.add_argument('-vc',
                        '--variants_colnames',
                        help='List of column names in the order of chromosome, position, and variant ID. '
                             'Provide exactly three names. Defaults are ("chr", "pos", "id").',
                        metavar='chr pos id',
                        nargs=3,
                        required=False,
                        default=("chr", "pos", "id"))
    parser.add_argument('-i',
                        '--interval',
                        help='Size of the window around each target variant (default: 500000)',
                        metavar='int',
                        type=int,
                        required=False,
                        default=500000)
    group_features = parser.add_mutually_exclusive_group(required=True)
    group_features.add_argument('-ft',
                                '--features',
                                help='Path to files with all feature annotations (one feature per line) in BED format '
                                     'and feature set description in GMT format. Two file paths separated by space '
                                     'should be provided in the following order: [BED] [GMT]. Cannot be used with '
                                     '--feature_files_dir',
                                metavar='path_bed path_gmt',
                                nargs=2)
    group_features.add_argument('-ffdir',
                                '--feature_files_dir',
                                help='A directory with feature files, one file per feature set, in BED format. File '
                                     'name will be used a the name of each feature set. Cannot be used with --features',
                                metavar='path',
                                type=str)
    parser.add_argument('-o',
                        '--out_json',
                        help='Output path for universe in json format.',
                        metavar='path',
                        type=str,
                        required=True)
    parser.add_argument('-tmp',
                        '--keep_temp',
                        help='If set, keep temporary files (mostly for debugging) (Default: False)',
                        action='store_true',
                        required=False,
                        default=False)
    # todo use keep_tmp in reality

    log_message("Processing input")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    log_message("Processing output dir")
    out_path = args.o
    out_dir = os.path.dirname(out_path)
    check_and_create_dir(out_dir)
    # temp files
    cur_id = current_random()
    universe_file = os.path.join(out_dir, f"universe_{cur_id}.bed")
    tmp_file = os.path.join(out_dir, f"tmp_{cur_id}.bed")
    inter_file = os.path.join(out_dir, f"inter2_{cur_id}.tsv")
    features_bed = os.path.join(out_dir, f"features_{cur_id}.bed")

    log_message("Read GMT/features input")
    variants = args.variants
    column_names = args.variants_colnames
    interval = args.interval
    if args.features is not None:
        bed, gmt = args.features
        set2features = read_gmt(gmt)
    else:
        feature_dir = args.feature_files_dir
        bed = features_bed
        set2features = get_features_from_dir(feature_dir,
                                             features_file_name=bed)

    log_message("Getting SNPs locations")
    snp2chrom_pos = get_snp_locations(variants,
                                      column_names)
    log_message("Creating universe")
    create_universe(snp2chrom_pos,
                    interval,
                    universe_out=universe_file,
                    tmp_file=tmp_file)
    log_message("Get overlaps")
    feature2intervals = get_overlapping_features(universe_file,
                                                 bed,
                                                 inter_file)
    log_message("Count intervals")
    interval_counts_for_universe = count_intervals(set2features,
                                                   feature2intervals)

    log_message("Creating features file")
    feature2pos = read_features_from_bed(bed)
    out_dict = {"interval": interval,
                "universe_intervals_number": len(snp2chrom_pos),
                "interval_counts": interval_counts_for_universe,
                "gene_set_dict": set2features,  # todo rename key
                "features": feature2pos}  # todo rename key

    log_message("Writing to out_path")
    with open(out_path, "w") as base:
        json.dump(out_dict, base)

    if args.keep_temp:
        log_message("Temporary files have not been deleted.")
    else:
        log_message("Deleting temporary files.")
        os.remove(universe_file)
        os.remove(tmp_file)
        os.remove(inter_file)
        os.remove(features_bed)
    log_message("Finished!")
