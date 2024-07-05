
# Simple script to quickly extract reads from a bam file given an associated lca file

import os
import argparse
import sys
import pysam

def main(in_lca, in_bam, out_bam, keyword):

    lcaheaderlines = 0

    with pysam.AlignmentFile(in_bam, "rb", check_sq=False, require_index=False) as infile, \
         pysam.AlignmentFile(out_bam, "wb", header=infile.header) as outfile, \
         open(in_lca, 'r') as lcafile:

        # Skip the LCA header lines
        for _ in range(lcaheaderlines):
            lcaline = next(lcafile)

        # Read the first LCA line and extract the read name and check if it's got the keyword in it
        lcaline = next(lcafile)
        iskeywordinlca = keyword in lcaline
        lcareadname = ":".join(lcaline.strip().split('\t')[0].split(':')[0:7])
        currentlymatching = False
        notdone = True

        # Get the first bam read
        try:
            bamread = next(infile)
        except StopIteration:
            print("Looks like your bam file is empty?")

        # iterate the lca file until you hit a keyword match. then, iterate the bam file until you hit the lca readname, and write all of those. 
        # then iterate the lca file again. 
        while notdone:
            if iskeywordinlca:
                # iterate the bam file until you find a match
                lcareadname = ":".join(lcaline.strip().split('\t')[0].split(':')[0:7])
                currentlymatching = bamread.query_name == lcareadname
                while not currentlymatching: # iterate until you hit the first instance of the read in the bam file
                    bamread = next(infile)
                    currentlymatching = bamread.query_name == lcareadname
                # now assume you are matching
                outfile.write(bamread)
                while currentlymatching:
                    try:
                        bamread = next(infile)
                        currentlymatching = bamread.query_name == lcareadname
                        if currentlymatching:
                            outfile.write(bamread)
                        else:
                            # done matching with this read name. increment the lca file
                            try:
                                lcaline = next(lcafile)
                                iskeywordinlca = keyword in lcaline
                            except StopIteration: 
                                notdone = False # done the lca file
                                break 
                    except StopIteration:
                        notdone = False
                        break # done the bam file
            else:
                while not iskeywordinlca:
                    try:
                        lcaline = next(lcafile)
                        iskeywordinlca = keyword in lcaline
                    except StopIteration:
                        notdone = False
                        break # done the lca file 

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('true', 't'):
        return True
    elif v.lower() in ('false', 'f'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    
    # Initialize
    parser = argparse.ArgumentParser(
        description="extractreads.py is a script to quickly extract reads from a bam file whose tax paths contain a keyword",
        epilog="\
          Written by Bianca De Sanctis in 2024. For assistance please contact bddesanctis@gmail.com.")
    
    # Mandatory arguments
    parser.add_argument("--in_lca", type=str, required=True, help="Path to the (sorted) LCA file")
    parser.add_argument("--in_bam", type=str, required=True, help="Path to the (sorted) BAM file")
    parser.add_argument("--out_bam", type=str, required=True, help="Path to the filtered BAM file")
    parser.add_argument("--keyword", type=str, required=True, help="Keyword or phrase to filter for")

    args = parser.parse_args()

    if '--help' in sys.argv or '-h' in sys.argv:
        parser.print_help()
        sys.exit()

    if not os.path.exists(args.in_lca):
        parser.error(f"Input LCA path does not exist: {args.in_lca}")
    if not os.path.exists(args.in_bam):
        parser.error(f"Input BAM path does not exist: {args.in_bam}")

    main(args.in_lca,args.in_bam,args.out_bam,args.keyword)




