from collections import defaultdict, Counter
import csv
import os
import subprocess
from typing import Dict, Tuple, List

from tqdm import tqdm
from datetime import datetime


def count_intervals(set2features: Dict, features: Dict, return_set: bool = True) -> Dict:
    """
    Counts the unique intervals associated with feature sets.

    :param set2features: (dict) A dictionary mapping set names to lists of features (genes).
    :param features: (dict) A dictionary mapping feature names (genes) to their respective intervals.
    :param return_set: (bool) If True, for each set returns dict: {all associated intervals (key): features, covered by this interval}; if False, returns the count only.
    :return: (dict) A dictionary where each key is a set name and the value is either the dict of unique intervals (if return_set) or the count.
    """
    res = {}
    for name, feature_list in tqdm(set2features.items()):
        if return_set:
            interval_dict = defaultdict(set)
            for feature in feature_list:
                if feature in features:
                    for target_interval in features[feature]:
                        interval_dict[target_interval].add(feature)
            res[name] = interval_dict
        else:
            interval_set = set()
            for feature in feature_list:
                if feature in features:
                    interval_set.update(features[feature])
            res[name] = len(interval_set)
    return res


def get_overlapping_features(path_to_bed: str,
                             path_to_gene_file: str,
                             intersect_file: str) -> Dict[str, List]:
    """
    Computes the intersection of genomic features from a BED file with genes from another gene file,
    extracts overlapping gene names along with their intervals, and returns a dictionary mapping each gene
    to the list of intervals it overlaps.

    This function requires the BEDTools software to be installed on the system and available in the path.

    :param path_to_bed: (str) The file path to the BED file containing genomic features.
    :param path_to_gene_file: (str) The file path to the gene file in BED format.
    :param intersect_file: (str) The file path where the intersection results will be temporarily stored.
    :return: (dict) A dictionary where keys are gene names and values are lists of interval strings in the
                 format "chromosome:start-end" for each overlap with the genomic features.
    :raises: RuntimeError if BEDTools fails to execute properly.
    :raises: IOError if there is an error reading from the intersection file.
    """
    feature2intervals = defaultdict(list)  # Gene -> interval id
    try:
        ret = subprocess.call(
            f"bedtools intersect -a \"{path_to_bed}\" -b \"{path_to_gene_file}\" -wo | perl -p -e 's/\r//g' > \"{intersect_file}\"",
            shell=True)   # todo (??) сделать пресорт sort -k1,1 -k2,2n и -sorted
        if ret != 0:
            raise RuntimeError(f"BEDTools intersect failed with exit code {ret}")
        with open(intersect_file, 'r', newline='') as inter:  # The results of clumping (SNPs sets)
            intersect_reader = csv.reader(inter, delimiter='\t')
            for row in tqdm(intersect_reader):
                if len(row) < 8:
                    log_message(f"Skipping malformed line: {row} (the structure is not [chrom, start, end, _, _, _, _, gene, _])", msg_type="WARN")
                    continue  # skip malformed lines
                chrom, start, end, _, _, _, _, gene, *_ = row
                interval_name = f'{chrom}:{start}-{end}'
                feature2intervals[gene].append(interval_name)
    except Exception as e:
        raise RuntimeError(f"Error during feature intersection: {e}")
    finally:
        # Remove the intermediate file to clean up
        if os.path.exists(intersect_file):
            os.remove(intersect_file)
    return dict(feature2intervals)


def get_snp_locations(tsv_file: str,
                      column_names: Tuple[str] = ('chr', 'pos', 'id'),
                      ) -> Dict[str, Tuple[str, int]]:
    """
    Extracts SNP locations from a TSV file, relying on a provided list of column names.
    
    :param tsv_file: (str) The path to the TSV file.
    :param column_names: (List[str]) List of column names in the order of 'chromosome', 'position', and 'variant_id'.
    :return: (dict) A dictionary mapping variant IDs to a tuple of chromosome and position.
    :raises:  ValueError: If there is a problem parsing the TSV file, including missing data or data conversion issues.
    :raises:  FileNotFoundError: If the TSV file is not found.
    """
    # the dictionary to store the SNP information
    snp2chrom_pos = defaultdict(tuple)
    
    if not os.path.isfile(tsv_file):
        raise FileNotFoundError(f"TSV file {tsv_file} not found.")
    with open(tsv_file, 'r', newline='') as csvfile:
        snps_reader = csv.reader(csvfile, delimiter='\t')

        # Read the header and determine column indices
        header = next(snps_reader)
        for col in column_names:
            if col not in header:
                raise ValueError(f"Column '{col}' not found in TSV header!")
        idx_chrom, idx_pos, idx_variant_id = [header.index(col) for col in column_names]

        # Process each row according to the identified column indices
        for row in tqdm(snps_reader):
            try:
                chrom = row[idx_chrom]
                pos = int(row[idx_pos])
                variant_id = row[idx_variant_id]
            except ValueError:
                raise ValueError(
                    f"Problem with: {row}\nCheck that your TSV file's format is correct and all necessary columns are present (from {column_names})!")
            except IndexError:
                raise ValueError(
                    f"Problem with: {row}\nOne of the specified columns was not found in the file: {tsv_file}")
            snp2chrom_pos[variant_id] = (chrom, pos)
    return dict(snp2chrom_pos)


def read_gmt(path: str) -> Dict[str, List]:
    """
    Reads a GMT file and returns a dictionary mapping set names to lists of features.
    
    :param path: Path to GMT file.
    :return: Dict {set_name: list_of_features (e.g. genes)}
    :raises: FileNotFoundError if the file does not exist.
    :raises: ValueError if any row does not have at least 3 columns.
    """
    set2features = dict()
    if not os.path.isfile(path):
        raise FileNotFoundError(f"GMT file {path} not found.")
    with open(path, 'r', newline='') as db:
        gmt_reader = csv.reader(db, delimiter='\t')
        for i, row in enumerate(tqdm(gmt_reader)):
            if len(row) < 3:
                raise ValueError(f"Row {i+1} in GMT file {path} does not have at least 3 columns (gene set, id/link, features (genes))!")
            gene_set = row[0]
            set2features[gene_set] = row[2:]
    return set2features


def get_features_from_dir(path: str,
                          features_file_name: str = 'features.bed') -> dict:
    """
    Processes all BED files in a given directory, extracts features from each file, and compiles them into a single
    BED file with unique feature identifiers. Additionally, returns a dictionary mapping each set name to the list of
    its features.

    :param path: (str) The directory containing BED files.
    :param features_file_name:  (str) Name of the output BED file to write compiled features. Default is 'features.bed'.

    :return: (dict) A dictionary where each key is the name of the original BED file (without the extension) and the value
          is a list of unique feature identifiers for features in that file.

    :raises: FileNotFoundError: If the specified directory does not exist.
    :raises: IOError: If there's an issue reading a BED file or writing to the output file.

    """
    if not os.path.isdir(path):
        raise FileNotFoundError(f"No directory found at {path}")
    # find all BED files in the specified directory
    file_set = [x for x in os.listdir(path) if x.endswith('.bed')]
    if not file_set:
        raise FileNotFoundError(f"No BED files found in directory {path}")
    feature_set_dict = defaultdict(list)
    try:
        with open(features_file_name, 'w') as features_file:
            for bed in tqdm(file_set):
                bed_path = os.path.join(path, bed)
                if not os.path.isfile(bed_path):
                    log_message(f"Skipping missing BED file: {bed_path}", msg_type="WARN")
                    continue
                with open(bed_path, 'r') as bed_file:
                    bed_contents = [x.strip() for x in bed_file.readlines()]
                set_name = os.path.splitext(bed)[0]
                for i, bed_entry in enumerate(bed_contents):
                    elems = bed_entry.split('\t')
                    if len(elems) < 3:
                        log_message(f"Skipping malformed line in BED file {bed_path} at line {i+1}", msg_type="WARN")
                        continue  # skip malformed lines
                    feature_name = f'{set_name}:{i}'
                    features_file.write(f'{elems[0]}\t{elems[1]}\t{elems[2]}\t{feature_name}\n')
                    feature_set_dict[set_name].append(feature_name)
    except Exception as e:
        raise IOError(f"Error processing BED files in {path}: {e}")
    return dict(feature_set_dict)


def read_features_from_bed(path: str) -> Dict[str, List]:
    """
    Reads a BED file and maps each feature's name (assumed to be in the fourth column) to its corresponding row as a list.

    :param path: (str) The path to the BED file.
    :return: (dict) A dictionary mapping feature names to lists representing their corresponding rows in the BED file.
    :raises FileNotFoundError: If the BED file is not found.
    """
    feature2pos = defaultdict(list)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"BED file {path} not found.")
    with open(path, 'r') as feature_file:
        features = [line.strip().split('\t') for line in feature_file]
    for feature in tqdm(features):
        if len(feature) < 4:
            log_message(f"Skipping malformed line in BED file {path} at line {features.index(feature) + 1}", msg_type="WARN")
            continue  # skip malformed lines
        feature2pos[feature[3]] = feature
    return dict(feature2pos)


def log_message(msg, msg_type="INFO"):
    """
    Logs a message with the current date and time, prefixed by a specified message type.

    :param msg: (str) The message to log.
    :param msg_type: (str) The type of message, e.g., "INFO", "WARN". Default is "INFO".
    """
    now = datetime.now()
    formatted_now = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{msg_type} {formatted_now}] {msg}")


def check_and_create_dir(out_name: str):
    """
    Checks if a directory exists at the specified path and creates it if it does not. Logs messages about the actions taken.

    :param out_name: (str) The directory path to check or create.
    """
    if os.path.exists(out_name):
        log_message(f'Output diretory {out_name} exists, writing there (files can be  overwritten)...', msg_type="WARN")
    else:
        log_message(f'Creating directory {out_name} and writing there...', msg_type="WARN")
        os.makedirs(out_name)


def get_filename_without_extension(file_path: str) -> str:
    """
    Cuts just path and extension of file.
    :param file_path: path to file in any format.
    :return: (str) just base name of the file, without path and extension.
    """
    # Get the base name of the file (e.g., 'example.txt' from '/path/to/example.txt')
    base_name = os.path.basename(file_path)
    # Split the base name and get the first part (e.g., 'example' from 'example.txt')
    file_name_without_extension = os.path.splitext(base_name)[0]
    return file_name_without_extension
