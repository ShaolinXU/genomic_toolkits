#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 20:17:52 2018

@author: xushaolin
"""

import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description="Summarize gene names from GenBank files.")
parser.add_argument("-i", "--infile", required=True, help="Input GenBank file")
parser.add_argument("-o", "--outfile", required=True, help="Output file for gene names")
parser.add_argument(
    "-s", "--strange_record", required=True, help="Output file for strange IDs"
)
parser.add_argument(
    "-k", "--only_test_key", type=int, default=0, help="Only test key set (default: 0)"
)

args = parser.parse_args()

infile = args.infile
outfile = args.outfile
strange_record = args.strange_record
only_test_key = args.only_test_key

genname_set = set()
strange_id = set()
key_set = set()

# Parse GenBank file and extract gene names or feature keys
for i, record in enumerate(SeqIO.parse(infile, "genbank")):
    for feature in record.features:
        # Collect all qualifier keys for testing
        key_set.update(feature.qualifiers.keys())

        if only_test_key == 0:
            # Define the order of qualifier keys to check for gene names
            qualifier_priority = ["gene", "product", "note", "label"]
            # Only process certain feature types
            if feature.type in {"tRNA", "rRNA", "gene", "CDS"}:
                for key in qualifier_priority:
                    if key in feature.qualifiers:
                        # Add the first value (lowercased) to the gene name set
                        genname_set.add(feature.qualifiers[key][0].lower())
                        break  # Stop after the first found key
    print(i)  # Print progress

# Output results
if only_test_key == 0:
    # Save strange IDs (currently always empty)
    pd.DataFrame(sorted(strange_id), columns=["strange_id"]).to_csv(
        strange_record, sep="\t", header=False, index=False
    )
    # Save sorted gene names
    pd.DataFrame(sorted(genname_set), columns=["gene_name"]).to_csv(
        outfile, sep="\t", header=False, index=False
    )
else:
    # Save all unique qualifier keys
    pd.DataFrame(sorted(key_set), columns=["key_set"]).to_csv(
        "feature_qualifiers_key.txt", sep="\t", header=False, index=False
    )
