
# bamdam extract

import sys 
import re
import csv
import pysam
import math
import argparse
import os
import subprocess
import random

from bamdam import utils
from bamdam import alignment_utils

def run_extract(args):
    formatted_tax = parse_tax(args)
    lca_file_type = utils.find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("Error: It looks like you're trying to run bamdam extract with a metaDMG-style lca file. Please use an ngsLCA-style lca file.")
        sys.exit()
    extract_reads(args.in_lca, args.in_bam, args.out_bam, formatted_tax, args.only_top_ref, args.only_top_alignment)


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

        return tax
    else:
        return []


def extract_reads(in_lca, in_bam, out_bam, tax, only_top_ref = False, only_top_alignment = False):
    # extracts all reads with a tax path containing a certain keyword.
    # also optionally shortens the header to only necessary ids.
    # subsetting the header is kinda slow because it requires running through the input twice.

    if len(tax) > 1 and only_top_ref:
        print("Warning: You have set the --only_top_ref flag and given multiple tax IDs. Please be aware that your output bam will contain only reads to the most-hit reference genome among all the tax IDs given.")
        # ^ is this true? 
    
    # catch it if there is no tax id or tax id file
    if not tax:
        if only_top_ref:
            print("No taxonomic nodes specified. Will extract all alignments to the best reference across all taxa.")
        elif only_top_alignment:
            print("No taxonomic nodes specified. Will extract the best alignment for each read.")
        else:
            print("Error: No taxonomic nodes specified. Please specify at least one of: a list of taxonomic nodes, only_top_ref or only_top_alignment.")
            sys.exit(1)
    else:
        # surround the tax ids so that you don't get e.g. 2001 when you wanted 200
        tax_pattern = "|".join([f"[[:space:]]{t}:" for t in tax])

    # get all the read names matching the requested tax id(s) out of the lca file
    if tax:
        grep_command = (
            f"grep -E '{tax_pattern}' {in_lca} | "
            f"awk -F'\\t' '{{print $1}}' | "
            f"awk -F':' '{{for(i=1; i<=NF-3; i++) printf $i (i<NF-3 ? \":\" : \"\\n\")}}'"
        )
    else:
        grep_command = (
            f"grep -E '{tax_pattern}' {in_lca} | "
            f"awk -F'\\t' '{{print $1}}')" )
    result = subprocess.run(grep_command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: No matches found for any of the requested tax, or grep/awk command failed.")
        return
    read_names = set(result.stdout.strip().splitlines())
    if(len(read_names) == 0):
        print(f"Error: No reads are assigned to tax id {tax[0]}. Please use numeric tax ids as input, and make sure they appear in your lca or tsv file.")
        sys.exit(1)

    # count ref hits among those read names (have to do this either way to relink ref IDs at the end, but as a bonus we can find the most common ref with this too)
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
        print(f"The most common reference is {most_common_reference} with {reference_count[most_common_reference]} total alignments.")
        if only_top_alignment:
            print(f"Your output bam will contain only reads with alignments to this reference, and only the best alignment to the reference per read.")
        else:
            print(f"Your output bam will contain all of the alignments to this reference, even if there is more than one per read.")
        header['SQ'] = [sq for sq in header.get('SQ', []) if sq['SN'] == most_common_reference] 
        reference_names = {most_common_reference}
    else: # get all the headers matching all of the refs we hit with our reads
        reference_names = set(reference_count.keys())
        header['SQ'] = [sq for sq in header.get('SQ', []) if sq['SN'] in reference_names]
    
    # we have to re-link the reference IDs in each read row to the new header because of how the bam compression works 
    ref_name_to_id = {sq['SN']: idx for idx, sq in enumerate(header.get('SQ', []))}

    # now do the extraction. different logic if you want the top alignments or not
    if only_top_alignment:
        with pysam.AlignmentFile(out_bam, "wb", header=header) as bam_writer:
            with pysam.AlignmentFile(in_bam, "rb") as bam_reader_again:
                current_read_name = "start"; current_best_alignment_list = []; current_best_alignment_score = int(-1000); relevant_ref_name = ""
                for alignment in bam_reader_again:
                    if alignment.query_name in read_names:
                        ref_name = bam_reader_again.get_reference_name(alignment.reference_id)
                        if ref_name in reference_names:
                            relevant_ref_name = ref_name
                            # here we are storing the best alignments for each read until the read is done and we know to write them. 
                            if alignment.query_name != current_read_name and current_read_name != "start":
                                # but first write the old read's best alignment. choose the best alignment if there are more than one
                                if len(current_best_alignment_list) > 1:
                                    current_best_alignment = random.choice(current_best_alignment_list) # there are multiple best alignments for this read; pick one randomly
                                elif len(current_best_alignment_list) == 1:
                                    current_best_alignment = current_best_alignment_list[0] # there is only one best alignment for this read
                                elif len(current_best_alignment_list) == 0:
                                    print(f"Error: No best alignment found for read {current_read_name}")
                                    sys.exit(1) # this should not happen
                                if current_best_alignment.reference_id >= 0: 
                                    current_best_alignment.reference_id = ref_name_to_id[ref_name] 
                                    # relink the reference_id to match the new header
                                else:
                                    print(f"Error: An alignment doesn't have a reference, for read name: {current_read_name}")
                                bam_writer.write(current_best_alignment)
                                # then re initialize the best alignment tracker
                                current_best_alignment_list = [alignment]
                                current_best_alignment_score = int(alignment.get_tag("AS"))
                                current_read_name = alignment.query_name
                            else: # add it to the tracker
                                try:
                                    alignment_score = int(alignment.get_tag("AS"))
                                except KeyError:
                                    print("Error: Can't find AS tags for at least one of your alignments in the BAM file, so it is not possible to determine the best alignments.")
                                    sys.exit(1)
                                if alignment_score > current_best_alignment_score:
                                    current_best_alignment_list = [alignment]
                                    current_best_alignment_score = alignment_score
                                if alignment_score == current_best_alignment_score:
                                    current_best_alignment_list.append(alignment)
                                if current_read_name == "start":
                                    current_read_name = alignment.query_name
                # writing is triggered by seeing a new relevant alignment to a new read name, so the last read with a relevant alignment hasn't been written yet
                # catch and write that here
                if len(current_best_alignment_list) > 1:
                    current_best_alignment = random.choice(current_best_alignment_list) # there are multiple best alignments for this read; pick one randomly
                elif len(current_best_alignment_list) == 1:
                    current_best_alignment = current_best_alignment_list[0] # there is only one best alignment for this read
                elif len(current_best_alignment_list) == 0:
                    print(f"Error: No best alignment found for read {current_read_name}")
                    sys.exit(1) # this should not happen
                if current_best_alignment.reference_id >= 0: 
                    current_best_alignment.reference_id = ref_name_to_id[relevant_ref_name] 
                    # relink the reference_id to match the new header
                else:
                    print(f"Error: An alignment doesn't have a reference, for read name: {current_read_name}")
                bam_writer.write(current_best_alignment)
                                
    # if you don't need the top alignment, we don't need to track them as we go, we just need to check the refs and relink
    else:
        with pysam.AlignmentFile(out_bam, "wb", header=header) as bam_writer:
            with pysam.AlignmentFile(in_bam, "rb") as bam_reader_again:
                current_read_name = "start"; current_best_alignment_list = []; current_best_alignment_score = int(-1000)
                for alignment in bam_reader_again:
                    if alignment.query_name in read_names:
                    # not just top alignments
                        if alignment.reference_id >= 0:
                            ref_name = bam_reader_again.get_reference_name(alignment.reference_id)
                            if ref_name in ref_name_to_id:
                                alignment.reference_id = ref_name_to_id[ref_name]
                            else:
                                continue  # skip reads with references not in the new header
                        bam_writer.write(alignment)




