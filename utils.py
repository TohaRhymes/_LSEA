from collections import defaultdict, Counter
import csv
import os
import subprocess
from typing import Dict, Tuple

import numpy as np

from datetime import datetime
import yaml


def count_intervals(set2features: Dict, features: Dict, emit_raw=True):
    res = defaultdict(int)
    explained_loci = []  # todo delete?
    for name in set2features.keys():
        int_count = 0  # todo delete?
        feature_list = set2features[name]  # List of genes for a trait
        interval_dict = defaultdict(set)  # Count every interval once
        for feature in features:
            if feature in feature_list:  # If gene belongs to the trait
                # Add every corresponding interval to set
                intervals = features[feature]
                int_count += len(intervals)
                for target_interval in intervals:
                    interval_dict[target_interval].add(feature)
        if emit_raw:
            res[name] = interval_dict
        else:
            # todo что из этого должно быть верным? если len, то тогда int_count можно удалить вообще
            # todo мы не используем интервал дикт, нам нужна по сути только длина ведь? len(intervals)
            # res[name] = int_count
            res[name] = len(interval_dict)
    return res


# Extract gene names that are in intersection
def get_overlapping_features(path_to_bed, path_to_gene_file, intersect_file):
    # todo delete inter_file
    features = defaultdict(list)  # Gene -> interval id
    subprocess.call(
        f"bedtools intersect -a {path_to_bed} -b {path_to_gene_file} -wo | perl -p -e 's/\r//g' > {intersect_file}",
        shell=True)  # сделать пресорт sort -k1,1 -k2,2n и -sorted
    with open(intersect_file, 'r', newline='') as inter:  # Our result of clumping (SNPs sets)
        intersect_reader = csv.reader(inter, delimiter='\t')
        for row in intersect_reader:
            # todo ['1', '0', '510177', '1', '1', '11869', '14409', 'DDX11L1', '2540']
            # todo правильные ли берем интервалы 1-2, not 5-6?
            chrom, start, end, _, _, _, _, gene, _ = row
            interval_name = f'{chrom}:{start}-{end}'
            features[gene].append(interval_name)
    features = dict(features)
    return features


# todo with header, any order
def get_snp_locations(tsv_file) -> Dict[Tuple]:
    """
    :param tsv_file:
    :return: dict: {variant_id:(chrom, pos)}
    """
    snp2chrom_pos = defaultdict(tuple)
    with open(tsv_file, 'r', newline='') as csvfile:  # Convert tsv to dictionary
        snps_reader = csv.reader(csvfile, delimiter='\t')
        _ = list(next(snps_reader))  # skip header (todo: поправить это недоразумение)
        for cur_row in snps_reader:  # Start from second row
            try:
                chrom = cur_row[0]
                pos = int(cur_row[1])
                variant_id = cur_row[2]
            except ValueError:
                print(f"Problem with: {cur_row}\nCheck that your tsv file's format is correct!")
                exit(1)
            snp2chrom_pos[variant_id] = (chrom, pos)
    return snp2chrom_pos


def read_gmt(path: str) -> Dict:
    """
    Read GMT file to dict
    :param path: path to GMT file.
    :return: Dict {set_name:list_of_features(e.g. genes)}
    """
    set2features = defaultdict(list)
    with open(path, 'r', newline='') as db:
        gmt_reader = csv.reader(db, delimiter='\t')
        for row in gmt_reader:
            gene_set = row[0]
            for gene in row[2:]:
                set2features[gene_set].append(gene)
    return dict(set2features)


# todo можно же поменять на такое (не повторяются же ключи?)
# def read_gmt(path: str) -> Dict:
#     set2features = dict()
#     with open(path, 'r', newline='') as db:
#         gmt_reader = csv.reader(db, delimiter='\t')
#         for row in gmt_reader:
#             gene_set = row[0]
#             set2features[gene_set] = row[2:]
#     return set2features


def get_features_from_dir(path, features_file_name='features.bed'):
    file_set = [x for x in os.listdir(path) if x.endswith('.bed')]
    feature_set_dict = defaultdict(list)
    with open(features_file_name, 'w') as features_file:
        for bed in file_set:
            bed_file = open(os.path.join(path, bed), 'r')
            bed_contents = [x.strip() for x in bed_file.readlines()]
            bed_file.close()
            set_name = os.path.splitext(bed)[0]
            for i, bed_entry in enumerate(bed_contents):
                elems = bed_entry.split('\t')
                feature_name = f'{set_name}:{i}'
                features_file.write(f'{elems[0]}\t{elems[1]}\t{elems[2]}\t{feature_name}\n')
                feature_set_dict[set_name].append(feature_name)
    return feature_set_dict


def read_features_from_bed(path):
    feature_dict = defaultdict(list)
    features = []
    with open(path, 'r') as feature_file:
        features = [line.strip().split('\t') for line in feature_file]
    for feature in features:
        feature_dict[feature[3]] = feature
    return feature_dict


def log_message(msg, msg_type="INFO"):
    # Get current date and time
    now = datetime.now()
    # Format the current date and time
    formatted_now = now.strftime("%Y-%m-%d %H:%M:%S")
    # Log the message with date and time
    print(f"[{msg_type} {formatted_now}] {msg}")


def check_and_create_dir(out_name: str):
    if os.path.exists(out_name):
        log_message(f'Output diretory {out_name} exists, writing there (files can be rewrote)...', msg_type="WARN")
        # shutil.rmtree(out_name)
    else:
        log_message(f'Creating directory {out_name} and writing there...', msg_type="WARN")
        os.makedirs(out_name)
