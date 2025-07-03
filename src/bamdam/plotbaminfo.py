
# bamdam plotbaminfo

import sys 
import re
import csv
import pysam
import math
import argparse
import os
import subprocess


try: # optional library only needed for plotting 
    import matplotlib.pyplot as plt
    matplotlib_imported = True
except:
    matplotlib_imported = False

def run_plotbaminfo(args):
    make_baminfo_plot(args.in_bam, args.in_bam_list, args.outplot)


def make_baminfo_plot(in_bam, in_bam_list, plotfile):

    if matplotlib_imported == False:
        print(f"Error: Cannot find matplotlib library for plotting. Try: pip install matplotlib")
        return
    
    if in_bam:
        bamfiles = in_bam if isinstance(in_bam, list) else [in_bam]
    elif in_bam_list:
        with open(in_bam_list, 'r') as file:
            bamfiles = [line.strip() for line in file if line.strip()]
    else:
        raise ValueError("Either --in_bam or --in_bam_list must be provided.")
    if not all(isinstance(b, str) for b in bamfiles):
        raise TypeError(f"bamfiles contains non-string elements: {bamfiles}")

    mismatch_counts_all = []
    read_length_counts_all = []

    for bam_file in bamfiles:
        bamfile = pysam.AlignmentFile(bam_file, "rb", require_index=False)

        mismatch_counts = {}  
        read_length_counts = {}  
        current_readname = None
        mismatch_bins = {}  
        alignment_count = 0
        current_read_length = 0
        total_reads = 0

        for read in bamfile:
            readname = read.query_name

            if readname != current_readname:
                if current_readname is not None:
                    for mismatch_bin, count in mismatch_bins.items():
                        fraction = count / alignment_count
                        if mismatch_bin not in mismatch_counts:
                            mismatch_counts[mismatch_bin] = fraction
                        else:
                            mismatch_counts[mismatch_bin] += fraction

                    if current_read_length not in read_length_counts:
                        read_length_counts[current_read_length] = 1
                    else:
                        read_length_counts[current_read_length] += 1

                    total_reads += 1

                current_readname = readname
                mismatch_bins = {}  
                alignment_count = 0
                current_read_length = read.query_length

            mismatches = read.get_tag("NM") if read.has_tag("NM") else 0
            if mismatches not in mismatch_bins:
                mismatch_bins[mismatches] = 1
            else:
                mismatch_bins[mismatches] += 1

            alignment_count += 1

        if current_readname is not None:
            for mismatch_bin, count in mismatch_bins.items():
                fraction = count / alignment_count
                if mismatch_bin not in mismatch_counts:
                    mismatch_counts[mismatch_bin] = fraction
                else:
                    mismatch_counts[mismatch_bin] += fraction

            if current_read_length not in read_length_counts:
                read_length_counts[current_read_length] = 1
            else:
                read_length_counts[current_read_length] += 1

            total_reads += 1

        bamfile.close()

        # Add the results to the overall dictionaries for all files
        mismatch_counts_all.append(mismatch_counts)
        read_length_counts_all.append(read_length_counts)

    # Plotting 
    plt.figure(figsize=(10, 5)) 

    # Mismatch Frequency Plot
    ax1 = plt.subplot(1, 2, 1)
    mismatch_color = '#009E73'  # Consistent color for the mismatch plot
    for i, bam_file in enumerate(bamfiles):
        mismatch_dict = mismatch_counts_all[i]
        x_vals = sorted(mismatch_dict.keys())
        y_vals = [mismatch_dict[x] for x in x_vals]
        ax1.plot(x_vals, y_vals, linestyle='-', label=os.path.basename(bam_file), color=mismatch_color, linewidth=2)
    ax1.set_xlabel("Number of mismatches", fontsize=12)
    ax1.set_ylabel("Frequency (per read)", fontsize=12)
    ax1.set_title("Mismatch Frequency Plot", fontsize=16)
    ax1.grid(False)

    # Read Length Distribution
    ax2 = plt.subplot(1, 2, 2)
    readlength_color = '#CC79A7'  # Consistent color for the read length plot
    for i, bam_file in enumerate(bamfiles):
        readlen_dict = read_length_counts_all[i]
        x_vals = sorted(readlen_dict.keys())
        y_vals = [readlen_dict[x] for x in x_vals]
        ax2.plot(x_vals, y_vals, linestyle='-', label=os.path.basename(bam_file), color=readlength_color, linewidth=2)
    ax2.set_xlabel("Read length", fontsize=12)
    ax2.set_ylabel("Number of reads", fontsize=12)
    ax2.set_title("Read Length Distribution", fontsize=16)
    ax2.grid(False)

    # Save the plot
    plt.tight_layout()
    file_extension = os.path.splitext(plotfile)[1].lower()
    if file_extension == '.pdf':
        plt.savefig(plotfile, format='pdf')
    else:
        if file_extension != '.png':
            print("Warning: Invalid plot file suffix. Your plot file is being saved in .png format, using the filename you specified.")
        plt.savefig(plotfile)

    plt.close()
