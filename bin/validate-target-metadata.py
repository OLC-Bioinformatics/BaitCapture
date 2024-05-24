#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse
import sys
import pkg_resources


__author__ = 'Liam Brown'
__email__ = 'liam.brown@inspection.gc.ca'


def get_deflines(fasta_file):
    with open(fasta_file, 'r') as fasta:
        deflines = [record.id for record in SeqIO.parse(fasta, 'fasta')]
    return deflines


def get_targets(csv_file):
    df = pd.read_csv(csv_file)
    targets = df['target'].tolist()
    return targets


def compare_lists(deflines, targets):
    missing_deflines = {defline: defline in targets for defline in deflines}
    if not all(missing_deflines.values()):
        missing = [defline for defline, present in missing_deflines.items() if not present]
        raise ValueError('\nThe following FASTA deflines are missing in the CSV file: ' + ', '.join(missing) +
                         "\n\nPlease ensure that all deflines in the `--targets` FASTA file are present in the 'target' column of the provided `--target_metadata` CSV file.")


def main():
    parser = argparse.ArgumentParser(description='Check FASTA and CSV files for missing information.')
    parser.add_argument('fasta_file', help='The FASTA file to check.')
    parser.add_argument('csv_file', help='The CSV file to check.')
    args = parser.parse_args()

    # Output package versions to file
    with open('package-versions.txt', 'w') as f:
        for package in pkg_resources.working_set:
            f.write(f"{package.project_name}=={package.version}\n")

    deflines = get_deflines(args.fasta_file)
    targets = get_targets(args.csv_file)

    try:
        compare_lists(deflines, targets)
        print('All deflines in the FASTA file are present in the CSV file.')
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())