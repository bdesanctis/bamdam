#!/usr/bin/env python3

# bamdam by Bianca De Sanctis, bddesanctis@gmail.com  
# last updated july 2 2025

import sys 
import re
import csv
import pysam
import math
import argparse
import os
import hyperloglog
import subprocess
from functools import lru_cache
try: # optional library for progress bars 
    from tqdm import tqdm
    tqdm_imported = True
except ImportError:
    tqdm_imported = False

# local
from bamdam import (
    utils,
    alignment_utils,
    shrink,
    compute,
    combine,
    extract,
    plotdamage,
    plotbaminfo,
    krona,
)

def main():

    # Initialize
    parser = argparse.ArgumentParser(
        description="Bamdam processes ancient metagenomic bam and lca files. Type bamdam [command] -h for more detailed help regarding a specific command.")
    
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Shrink
    parser_shrink = subparsers.add_parser('shrink', help="Filter the BAM and LCA files.")
    parser_shrink.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of CPU threads to use for BAM compression/decompression (default: 1)",
    )
    parser_shrink.add_argument("--in_lca", type=str, required=True, help="Path to the input LCA file (required)")
    parser_shrink.add_argument("--in_bam", type=str, required=True, help="Path to the input (read-sorted) BAM file (required)")
    parser_shrink.add_argument("--out_lca", type=str, required=True, help="Path to the short output LCA file (required)")
    parser_shrink.add_argument("--out_bam", type=str, required=True, help="Path to the short output BAM file (required)")
    parser_shrink.add_argument("--stranded", type=str, required=True, help="Either ss for single stranded or ds for double stranded (required)")
    parser_shrink.add_argument("--mincount", type=int, default=5, help="Minimum read count to keep a node (default: 5)")
    parser_shrink.add_argument("--upto", type=str, default="family", help="Keep nodes up to and including this tax threshold (default: family)")
    parser_shrink.add_argument("--minsim", type=float, default=0.9, help="Minimum similarity to reference to keep an alignment (default: 0.9)")
    parser_shrink.add_argument("--exclude_tax", type=str, nargs='+', default=[], help="Numeric tax ID(s) to exclude when filtering (default: none)")
    parser_shrink.add_argument("--exclude_tax_file", type=str, default=None, help="File of numeric tax ID(s) to exclude when filtering, one per line (default: none)")
    parser_shrink.add_argument("--annotate_pmd", action='store_true', help="Annotate output bam file with PMD tags  (default: not set)")
    parser_shrink.set_defaults(func=shrink.run_shrink) 

    # Compute
    parser_compute = subparsers.add_parser('compute', help="Compute tsv and subs files.")
    parser_compute.add_argument("--in_bam", type=str, required=True, help="Path to the BAM file (required)")
    parser_compute.add_argument("--in_lca", type=str, required=True, help="Path to the LCA file (required)")
    parser_compute.add_argument("--out_tsv", type=str, required=True, help="Path to the output tsv file (required)")
    parser_compute.add_argument("--out_subs", type=str, required=True, help="Path to the output subs file (required)")
    parser_compute.add_argument("--stranded", type=str, required=True, help="Either ss for single stranded or ds for double stranded (required)")
    parser_compute.add_argument("--k", type=int, default=29, help="Value of k for per-node counts of unique k-mers and duplicity (default: 29)")
    parser_compute.add_argument("--upto", type=str, default="family", help="Keep nodes up to and including this tax threshold; use root to disable (default: family)")
    parser_compute.add_argument("--plotdupdust", type=str, help="Path to create a duplicity-dust plot for this sample (default: not set)")
    parser_compute.set_defaults(func=compute.run_compute)

    # Extract
    parser_extract = subparsers.add_parser('extract', help="Extract bam alignments of reads assigned to specific taxa.")
    parser_extract.add_argument("--in_bam", type=str, required=True, help="Path to the BAM file (required)")
    parser_extract.add_argument("--in_lca", type=str, required=True, help="Path to the LCA file (required)")
    parser_extract.add_argument("--out_bam", type=str, required=True, help="Path to the filtered BAM file (required)")
    parser_extract.add_argument("--tax", type=str, nargs='+', default=[], help="Numeric tax ID(s) to extract (default: none)")
    parser_extract.add_argument("--tax_file", type=str, default=None, help="File of numeric tax ID(s) to extract, one per line (default: none)")
    parser_extract.add_argument("--only_top_ref", action='store_true', help="Only keep alignments to the most-hit reference (default: not set)")
    parser_extract.set_defaults(func=extract.run_extract)

    # Plot damage
    parser_plotdamage = subparsers.add_parser('plotdamage', help="Produces a postmortem damage plot for a specified taxonomic node using the subs file.")
    group_input_plotdamage = parser_plotdamage.add_mutually_exclusive_group(required=True)
    group_input_plotdamage.add_argument("--in_subs", nargs='+', help="Input subs file(s)")
    group_input_plotdamage.add_argument("--in_subs_list", help="Path to a text file contaning input subs files, one per line")
    parser_plotdamage.add_argument("--tax", type=str, required=True, help="Taxonomic node ID (required)")
    parser_plotdamage.add_argument("--outplot", type=str, default="damage_plot.png", help="Filename for the output plot, ending in .png or .pdf (default: damage_plot.png)")
    parser_plotdamage.add_argument("--ymax", type=str, default="0", help="Maximum for y axis (optional)")
    parser_plotdamage.set_defaults(func=plotdamage.run_plotdamage)

    # Plot bam info
    parser_plotbaminfo = subparsers.add_parser('plotbaminfo', help="Produces a mismatch and read length distribution plot for an input bam.")
    group_input_plotbaminfo = parser_plotbaminfo.add_mutually_exclusive_group(required=True)
    group_input_plotbaminfo.add_argument("--in_bam", nargs='+', help="Input bam file(s)")
    group_input_plotbaminfo.add_argument("--in_bam_list", help="Path to a text file containing input bams, one per line")
    parser_plotbaminfo.add_argument("--outplot", type=str, default="baminfo_plot.png", help="Filename for the output plot, ending in .png or .pdf (default: baminfo_plot.png)")
    parser_plotbaminfo.set_defaults(func=plotbaminfo.run_plotbaminfo)

    # Combine
    parser_combine = subparsers.add_parser(
        'combine', help="Combine multiple bamdam tsv files to generate a multi-sample tsv.")
    group_input_combine = parser_combine.add_mutually_exclusive_group(required=True)
    group_input_combine.add_argument("--in_tsv", nargs='+', help="List of input tsv file(s)")
    group_input_combine.add_argument("--in_tsv_list", help="Path to a text file containing paths to input tsv files, one per line")
    parser_combine.add_argument("--out_tsv", type=str, default="combined.tsv", help="Path to output tsv file name (default: combined.tsv)")
    parser_combine.add_argument("--minreads", type=float, default=50, help="Minimum reads across samples to include taxa (default: 50)")
    parser_combine.add_argument(
        "--include",
        nargs="*",
        choices=["damage", "duplicity", "dust", "taxpath", "gc", "all", "none"],
        default=["all"],
        help="Additional metrics to include in output file. Specify any combination of the options, 'all', or 'none'. Supports: damage, duplicity, dust, taxpath, gc (default: all)",
    )
    parser_combine.set_defaults(func=combine)

    # krona
    parser_krona = subparsers.add_parser(
        'krona', help="Generate Krona XML file from one or more tsv files.")
    group_input_krona = parser_krona.add_mutually_exclusive_group(required=True)
    group_input_krona.add_argument("--in_tsv", type=str,nargs='+', help="Path to tsv file(s)")
    group_input_krona.add_argument("--in_tsv_list", help="Path to a text file containing paths to input tsv files, one per line")
    parser_krona.add_argument("--out_xml", type=str, default="out.xml", help="Path to output xml file name (default: out.xml)")
    parser_krona.add_argument("--minreads", type=int, default=100, help="Minimum reads across samples to include taxa (default: 100)")
    parser_krona.add_argument("--maxdamage",type=float, default=None, help="Force a maximum value for the 5' C-to-T damage color scale. If not provided, the maximum value is determined from the data, with a minimum threshold of 0.3. (not recommended by default)")
    parser_krona.set_defaults(func=krona)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if '--help' in sys.argv or '-h' in sys.argv or len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    # Validation checks
    if hasattr(args, 'stranded') and args.stranded not in ["ss", "ds"]:
        parser.error(f"Invalid value for stranded: {args.stranded}. Valid values are ss or ds.")
    if hasattr(args, 'mincount') and not isinstance(args.mincount, int):
        parser.error(f"Invalid integer value for mincount: {args.mincount}")
    if hasattr(args, 'k') and (not isinstance(args.k, int) or not isinstance(args.k, int) or args.k > 50):
        parser.error(f"Invalid integer value for k : {args.k} (max 49, and that is much higher than recommended in any case)")
    if hasattr(args, 'upto') and not re.match("^[a-z]+$", args.upto):
        parser.error(f"Invalid value for upto: {args.upto}. Must be a string of only lowercase letters.")
    if hasattr(args, 'minsim') and not isinstance(args.minsim, float):
        parser.error(f"Invalid float value for minsim: {args.minsim}")
    if hasattr(args, 'in_lca') and not os.path.exists(args.in_lca):
        parser.error(f"Input LCA path does not exist: {args.in_lca}")
    if hasattr(args, 'upto') and args.upto=="clade":
        parser.error(f"Clade is not a valid taxonomic level in bamdam because there can be multiple clades in one taxonomic path.")
    if hasattr(args, 'upto') and "sub" in args.upto:
        parser.error(f"The taxonomic level cannot start with 'sub' (eg subfamily, subphylum) because this is inconsistently defined in taxonomy (not all species belong to a subfamily, but they all belong to a family).") 
    if hasattr(args, 'upto') and args.upto != args.upto.lower():
        parser.warning(f"Warning: {args.upto} as provided is not in lowercase, but it should be. Converting to lowercase and moving on.")
        args.upto = args.upto.lower()
    
    if hasattr(args, 'in_bam') and args.command != "plotbaminfo":
        sortorder = utils.get_sorting_order(args.in_bam)
        if sortorder != "queryname":
            print("Error: Your bam file does not appear to be read-sorted. Please try again with it once it has been read-sorted (samtools sort -n), which should be the same order as your lca file. \n")
            exit(-1) 

    if hasattr(args, 'in_tsv') and hasattr(args, 'in_tsv_list') and args.in_tsv is not None and args.in_tsv_list is not None:
        raise ValueError("You cannot specify both --in_tsv and --in_tsv_list at the same time. Use one or the other.")
    if hasattr(args, 'in_subs') and hasattr(args, 'in_subs_list') and args.in_subs is not None and args.in_subs_list is not None:
        raise ValueError("You cannot specify both --in_subs and --in_subs_list at the same time. Use one or the other.")
    if hasattr(args, 'in_bams') and hasattr(args, 'in_bam_list') and args.in_bam is not None and args.in_bam_list is not None:
        raise ValueError("You cannot specify both --in_bam and --in_bam_list at the same time. Use one or the other.")


    if hasattr(args, 'minreads') and args.minreads < 0:
        raise ValueError("Min reads must be a non-negative integer.")
    if hasattr(args, 'include'):
        invalid_metrics = set(args.include) - {'damage', 'duplicity', 'dust', 'taxpath', 'all', 'none'}
        if invalid_metrics:
            raise ValueError(f"Invalid metrics in --include: {', '.join(invalid_metrics)}. Allowed values are: 'damage', 'duplicity', 'dust', 'taxpath', 'all', 'none', or combinations of the first four.")

    if args.command == 'shrink':
        print("Hello! You are running bamdam shrink with the following arguments:")
        print(f"in_lca: {args.in_lca}")
        print(f"in_bam: {args.in_bam}")
        print(f"out_lca: {args.out_lca}")
        print(f"out_bam: {args.out_bam}")
        print(f"stranded: {args.stranded}")
        print(f"mincount: {args.mincount}")
        print(f"upto: {args.upto}")
        print(f"minsim: {args.minsim}")
        if hasattr(args, 'exclude_tax_file') and args.exclude_tax_file:
            print(f"exclude_tax: loaded from {args.exclude_taxfile}")
        if hasattr(args, 'exclude_tax') and args.exclude_tax:
            print(f"exclude_tax: {args.exclude_tax}")
        if hasattr(args, 'annotate_pmd') and args.annotate_pmd:
            print(f"annotate_pmd: {args.annotate_pmd}")

    elif args.command == 'compute':
        print("Hello! You are running bamdam compute with the following arguments:")
        print(f"in_bam: {args.in_bam}")
        print(f"in_lca: {args.in_lca}")
        print(f"out_tsv: {args.out_tsv}")
        print(f"out_subs: {args.out_subs}")
        print(f"stranded: {args.stranded}")
        print(f"k: {args.k}")
        print(f"upto: {args.upto}")
        if hasattr(args, 'plotdupdust'):
            print(f"printdupdust: {args.plotdupdust}")

    elif args.command == 'extract':
        print("Hello! You are running bamdam extract with the following arguments:")
        print(f"in_bam: {args.in_bam}")
        print(f"in_lca: {args.in_lca}")
        print(f"out_bam: {args.out_bam}")
        if args.tax:
            print(f"tax: {args.tax}")
        else:
            print(f"tax_file: {args.tax_file}")
        print(f"only_top_ref: {args.only_top_ref}")

    elif args.command == 'combine':
        print("Hello! You are running bamdam combine with the following arguments:")
        if hasattr(args, 'input') and args.in_tsv:
            print(f"Input files: {', '.join(args.in_tsv)}")
        if hasattr(args, 'input_files') and args.in_tsv_files:
            print(f"Input file list: {args.in_tsv_files}")
        print(f"Output file: {args.out_tsv}")
        print(f"Min reads: {args.minreads}")
        if args.include:
            print(f"Included metrics: {', '.join(args.include)}")   

    elif args.command == 'krona':
        print("Hello! You are running bamdam krona with the following arguments:")
        if hasattr(args, 'input') and args.in_tsv:
            print(f"Input files: {' '.join(args.in_tsv)}")
        if hasattr(args, 'input_files') and args.in_tsv_files:
            print(f"Input file list: {args.in_tsv_files}")
        print(f"Output file: {args.out_xml}")    
        print(f"Min reads: {args.minreads}")
        if hasattr(args, 'maxdamage') and args.maxdamage is not None:
            print(f"Max damage value for colour scale: {args.maxdamage}")

    if not tqdm_imported:
        print("The library tqdm is not available, so progress bars will not be shown. This will not impact performance.")

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()



if __name__ == "__main__":
    main()

