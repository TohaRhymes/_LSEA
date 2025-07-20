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

from tqdm import tqdm
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
                    tmp_file: str = 'tmp.bed',
                    keep_temp=False) -> None:
    """
    Creates a universe of genomic intervals based on SNP positions and a specified interval size, outputs to a BED file.
    This function writes SNP intervals to a temporary file, sorts them, and outputs to a final BED file.
    The temporary file is deleted after use.


    :param snp2chrom_pos: (Dict[str, Tuple[str, int]]): Dictionary mapping SNP identifiers to tuples of chromosome and position.
    :param interval: (int): The interval size to extend around each SNP position, both upstream and downstream.
    :param universe_out: (str): The output path for the final sorted BED file. Defaults to 'universe.bed'.
    :param tmp_file: (str): The path for the temporary BED file used for sorting. Defaults to 'tmp.bed'.
    :raises: OSError if file writing or sorting fails.
    """
    try:
        with open(tmp_file, 'w', newline='') as bed_file:  # Here we write to new file
            bed_writer = csv.writer(bed_file, delimiter='\t')
            for cur_id, (snp, (chrom, pos)) in tqdm(enumerate(snp2chrom_pos.items())):
                start = max(0, int(pos) - interval)
                end = int(pos) + interval
                # Chromosome, start, end coordinates, id
                bed_row = [chrom, start, end, cur_id]
                bed_writer.writerow(bed_row)
        log_message("Sorting Universe...")
        ret = subprocess.call(f"sort -k1,1 -k2,2n {tmp_file} > {universe_out}", shell=True)
        if ret != 0:
            raise RuntimeError(f"Sorting failed with exit code {ret}")
    except Exception as e:
        log_message(f"Error during universe creation: {e}", msg_type="ERROR")
        raise
    finally:
        # Clean up the temporary file
        if os.path.exists(tmp_file) and not keep_temp:
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
    out_path = args.out_json
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
    # Validate input files and columns
    if not os.path.isfile(variants):
        log_message(f"Variants file {variants} does not exist!", msg_type="ERROR")
        sys.exit(1)
    # Validate column names for variants
    with open(variants, 'r') as f:
        header = f.readline().strip().split('\t')
        for col in column_names:
            if col not in header:
                log_message(f"Column '{col}' not found in variants file header!", msg_type="ERROR")
                sys.exit(1)
    # Validate features/gmt or feature_files_dir
    if args.features is not None:
        bed, gmt = args.features
        if not os.path.isfile(bed):
            log_message(f"BED file {bed} does not exist!", msg_type="ERROR")
            sys.exit(1)
        if not os.path.isfile(gmt):
            log_message(f"GMT file {gmt} does not exist!", msg_type="ERROR")
            sys.exit(1)
        # Validate GMT format (at least 3 columns per row)
        with open(gmt, 'r') as gmtf:
            for i, row in enumerate(gmtf):
                if len(row.strip().split('\t')) < 3:
                    log_message(f"Row {i+1} in GMT file {gmt} does not have at least 3 columns!", msg_type="ERROR")
                    sys.exit(1)
        set2features = read_gmt(gmt)
    else:
        feature_dir = args.feature_files_dir
        if not os.path.isdir(feature_dir):
            log_message(f"Feature directory {feature_dir} does not exist!", msg_type="ERROR")
            sys.exit(1)
        bed = features_bed
        set2features = get_features_from_dir(feature_dir, features_file_name=bed)

    log_message("Getting SNPs locations")
    snp2chrom_pos = get_snp_locations(variants,
                                      column_names)
    log_message("Creating universe")
    create_universe(snp2chrom_pos,
                    interval,
                    universe_out=universe_file,
                    tmp_file=tmp_file,
                    keep_temp=args.keep_temp)
    log_message("Get overlaps")
    feature2intervals = get_overlapping_features(universe_file,
                                                 bed,
                                                 inter_file)
    log_message("Count intervals")
    interval_counts_for_universe = count_intervals(set2features,
                                                   feature2intervals,
                                                   return_set=False)

    log_message("Creating features file")
    feature2pos = read_features_from_bed(bed)
    out_dict = {"interval": interval,
                "universe_intervals_number": len(snp2chrom_pos),
                "interval_counts": interval_counts_for_universe,
                "gene_set_dict": set2features,  # todo rename key
                "features": feature2pos}  # todo rename key

    log_message(f"Writing to {out_path}")
    # Validate that SNPs were found
    if not snp2chrom_pos:
        log_message("No SNPs found in the variants file after parsing!", msg_type="ERROR")
        sys.exit(1)
    # Validate that features were found
    if not feature2pos:
        log_message("No features found in the BED file after parsing!", msg_type="ERROR")
        sys.exit(1)
    try:
        with open(out_path, "w") as base:
            json.dump(out_dict, base)
    except Exception as e:
        log_message(f"Failed to write output JSON: {e}", msg_type="ERROR")
        sys.exit(1)

    # Clean up temporary files with error handling
    if args.keep_temp:
        log_message("Temporary files have not been deleted.")
    else:
        log_message("Deleting temporary files.")
        for f in [universe_file, inter_file, features_bed]:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except Exception as e:
                log_message(f"Failed to delete temporary file {f}: {e}", msg_type="WARN")
    log_message("Finished!")
