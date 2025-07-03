
# bamdam extract

import sys 
import re
import csv
import pysam
import math
import argparse
import os
import subprocess

from bamdam import utils
from bamdam import alignment_utils

def run_extract(args):
    formatted_tax = parse_tax(args)
    lca_file_type = utils.find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("Error: It looks like you're trying to run bamdam extract with a metaDMG-style lca file. Please use an ngsLCA-style lca file.")
        sys.exit()
    extract_reads(args.in_lca, args.in_bam, args.out_bam, formatted_tax, args.only_top_ref)


def parse_tax(args):
    if hasattr(args, 'tax') or hasattr(args, 'tax_file'):
        if args.tax and args.tax_file:
            raise ValueError("Please provide exactly one of --tax or --tax_file, not both.")
        
        tax = args.tax if args.tax else []
        if not isinstance(tax, list):
            tax = [tax]
        
        if args.tax_file:
            if not os.path.exists(args.tax_file):
                raise FileNotFoundError(f"tax_file path does not exist: {args.tax_file}")
            with open(args.tax_file, 'r') as f:
                lines = [line.strip() for line in f if line.strip()]
            if not lines:
                raise ValueError(f"tax_file is empty: {args.tax_file}")
            tax.extend(lines)
        
        tax = [kw.lstrip("'").lstrip('"').rstrip("'").rstrip('"') for kw in tax] # strip quotes 
        
        for t in tax:
            if not t.isdigit():
                raise ValueError(f"Invalid tax ID: '{t}' is not numeric.")

        # formatted_keywords = []
        # good to surround the digit-only tax ids with a :, so we don't accidentally hit substring tax ids
        #for keyword in exclude_keywords:
        #    if keyword.isdigit():
        #        keyword = f"{keyword}:"
        #    formatted_keywords.append(keyword)

        return tax
    else:
        return []


def extract_reads(in_lca, in_bam, out_bam, tax, only_top_ref = False):
    # extracts all reads with a tax path containing a certain keyword.
    # also optionally shortens the header to only necessary ids.
    # subsetting the header is kinda slow because it requires running through the input twice.

    if len(tax) > 1 and only_top_ref:
        print("Warning: You have set the --only_top_ref flag and given multiple tax IDs. Please be aware that your output bam will contain only reads to the most-hit reference genome among all the tax IDs given.")
    # surround the tax ids so that you don't get e.g. 2001 when you wanted 200
    tax_pattern = "|".join([f"[[:space:]]{t}:" for t in tax])

    # get all the read names matching the requested tax id(s) out of the lca file
    grep_command = (
        f"grep -E '{tax_pattern}' {in_lca} | "
        f"awk -F'\\t' '{{print $1}}' | "
        f"awk -F':' '{{for(i=1; i<=NF-3; i++) printf $i (i<NF-3 ? \":\" : \"\\n\")}}'"
    )
    result = subprocess.run(grep_command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: No matches found for any of the requested tax, or grep/awk command failed.")
        return
    read_names = set(result.stdout.strip().splitlines())

    # count ref hits among those read names
    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        header = bam_in.header.to_dict() # get original header; we are only subsetting the sq lines
        reference_count = {}
        for read in bam_in:
            if read.query_name in read_names:
                ref_name = bam_in.get_reference_name(read.reference_id)
                if ref_name:
                    reference_count[ref_name] = reference_count.get(ref_name, 0) + 1

    if only_top_ref:
        # find the most common reference
        most_common_reference = max(reference_count, key=reference_count.get)
        print(f"The most common reference is {most_common_reference} with {reference_count[most_common_reference]} alignments.")
        print(f"Your output bam will contain all alignments to this reference, even if there is more than one per read.")
        header['SQ'] = [sq for sq in header.get('SQ', []) if sq['SN'] == most_common_reference] 
        reference_names = {most_common_reference}
    else: # get all the headers matching all of the refs we hit with our reads
        reference_names = set(reference_count.keys())
        header['SQ'] = [sq for sq in header.get('SQ', []) if sq['SN'] in reference_names]

    # important : we have to re-link the reference IDs in each read row to the new header because of how the bam compression works
    ref_name_to_id = {sq['SN']: idx for idx, sq in enumerate(header.get('SQ', []))}
    # write the filtered reads with the updated header and re-linked reference IDs
    with pysam.AlignmentFile(out_bam, "wb", header=header) as bam_writer:
        with pysam.AlignmentFile(in_bam, "rb") as bam_reader_again:
            for read in bam_reader_again:
                if read.query_name in read_names:
                    # relink the reference_id to match the new header
                    if read.reference_id >= 0:
                        ref_name = bam_reader_again.get_reference_name(read.reference_id)
                        if ref_name in ref_name_to_id:
                            read.reference_id = ref_name_to_id[ref_name]
                        else:
                            continue  # skip reads with references not in the new header
                    bam_writer.write(read)
