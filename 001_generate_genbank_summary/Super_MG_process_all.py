#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:13:49 2018
get all 13 PCGs of the each animals
@author: xushaolin
"""
# fix the character: : . ' " in name of sequence

import argparse
import os
from typing import Dict, List


import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Process GenBank files and extract gene information."
)
parser.add_argument("-i", "--infile", required=True, help="Input GenBank file")
parser.add_argument("-g", "--ingene", required=True, help="Gene name file")
parser.add_argument(
    "-t", "--out_table", required=True, help="Output table name (without extension)"
)
parser.add_argument("-o", "--out_dir", required=True, help="Output directory")
parser.add_argument("-c", "--fcol_names", required=True, help="File with column names")
parser.add_argument("-f", "--ffor_fasta", required=True, help="File with fasta names")

args = parser.parse_args()

infile = args.infile
ingene = args.ingene
out_table = args.out_table
out_dir = args.out_dir
fcol_names = args.fcol_names
ffor_fasta = args.ffor_fasta

# create the directory for output
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# read the col_name from fcol_names
col_names = []
with open(fcol_names, "r") as fh:
    for line in fh:
        if len(line.strip()) > 0:
            col_names.append(line.strip())


# the function to parse the gene_name file
def parse_gene_name(in_file: str) -> Dict[str, List[str]]:
    the_dict: Dict[str, List[str]] = {}
    the_gene_now: str = ""
    with open(in_file, "r") as fh_2:
        for line_2 in fh_2:
            if line_2.startswith(">"):
                the_gene_now = line_2[1:].strip()
                the_dict[the_gene_now] = []
            elif len(line_2.strip()) == 0:
                continue
            else:
                if the_gene_now:  # Ensure gene header was found
                    the_dict[the_gene_now].append(line_2.strip())
    return the_dict


# extract the gene name info to dictionary
the_gene_name_dict = parse_gene_name(ingene)

# extract the dictionary's keys
the_key = list(the_gene_name_dict.keys())

# create the table for the final data
the_table = pd.DataFrame(columns=col_names)

# read the genbank file into a generator
all_records = SeqIO.parse(infile, "genbank")
# Process each GenBank record and extract relevant information
for position, record in enumerate(all_records):
    # Collect general record info
    new_info = [
        record.name,
        record.annotations.get("organism", ""),
        ":".join(record.annotations.get("taxonomy", [])),
        record.description,
    ]

    # Extract reference info if available
    ref = record.annotations.get("references", [{}])
    if ref and hasattr(ref[0], "journal"):
        journal = ref[0].journal
        PMED = getattr(ref[0], "pubmed_id", "")
        title = getattr(ref[0], "title", "")
        authors = getattr(ref[0], "authors", "")
    else:
        journal = PMED = title = authors = ""
    new_info.extend([journal, PMED, title, authors])

    # Prepare gene info placeholders
    the_gene_info = ["__NaN__"] * (len(the_key) * 2)

    # Helper function to extract gene name from feature qualifiers
    def get_gene_name(_feature):
        for key in ["gene", "product", "note", "label"]:
            if key in _feature.qualifiers:
                return _feature.qualifiers[key][0].lower().strip()
        return None

    # Iterate over features and fill gene info if matched
    for feature in record.features:
        gene_name = get_gene_name(feature)
        if gene_name:
            for j, name in enumerate(the_key):
                if gene_name in the_gene_name_dict[name]:
                    the_gene_info[j * 2] = f"{name}:{feature.location.strand}"
                    the_gene_info[j * 2 + 1] = str(feature.extract(record.seq))
                    break

    # Combine all info and append to the table
    new_info.extend(the_gene_info)
    the_table = pd.concat(
        [the_table, pd.DataFrame([new_info], columns=col_names)], ignore_index=True
    )
    print(position)

final_out_table_txt = out_dir + "/" + out_table + ".txt"
final_out_table_csv = out_dir + "/" + out_table + ".csv"

the_table.to_csv(final_out_table_txt, sep="\t", index=False)
the_table.to_csv(final_out_table_csv, sep=",", index=False)


for_fasta = list()

with open(ffor_fasta, "r") as fh:
    for line in fh:
        for_fasta.append(line.strip())


# Write gene sequences to separate FASTA files for each gene
for i, gene in enumerate(the_key):
    fasta_filename = os.path.join(out_dir, f"{for_fasta[i + 2]}.fasta")
    with open(fasta_filename, "w") as fh:
        for j in range(the_table.shape[0]):
            seq = the_table.iloc[j, i * 2 + 9]
            if seq != "__NaN__":
                # Use species name if available, else use record name
                species = the_table.iloc[j, 1].replace(" ", "_")
                header = (
                    f">{species}_{the_table.iloc[j, 0]}"
                    if species
                    else f">{the_table.iloc[j, 0]}"
                )
                fh.write(f"{header}\n{seq}\n")
