
# utils.py

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

level_order = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus", "subgenus", "species", "subspecies"]

# useful functions used across modules

# little string processing functions
@lru_cache(maxsize=1024)
def parse_cigar_cached(cigar):
    """Cache CIGAR parsing since many reads have identical CIGAR strings"""
    return re.findall(r"\d+\D", cigar)

@lru_cache(maxsize=1024)
def parse_md_cached(md):
    """Cache MD parsing since many reads have identical MD strings"""
    return re.compile(r"\d+|\^[A-Za-z]+|[A-Za-z]").findall(md)

@lru_cache(maxsize=128)
def is_reverse_strand(flagsum):
    """Cache-enabled function to check if a read is on reverse strand"""
    bin_flagsum = bin(flagsum)[::-1]
    bit_position = 4  # 2^4 = 16 this flag means it's aligned in reverse
    return len(bin_flagsum) > bit_position and bin_flagsum[bit_position] == "1"

@lru_cache(maxsize=128)
def rev_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
    return ''.join(complement[base] for base in reversed(seq))

@lru_cache(maxsize=128)
def get_rep_kmer(seq): # representative canonical kmer representation for counting, gets the lexicographical min of a kmer and its rev complement
    rep_kmer = min(seq,rev_complement(seq))
    return rep_kmer


# common file checks


def line_count(file_path): # just a wc -l wrapper
    result = subprocess.run(['wc', '-l', file_path], stdout=subprocess.PIPE, text=True)
    line_count = int(result.stdout.split()[0])
    return line_count

def get_sorting_order(file_path):
    # bamdam needs query / read sorted bams
    # this is the fastest way to check if HD is the first line of the header 
    # a bam is basically a gzipped sam; take advantage of this so we don't have to read in the full header just for this check (pysam can't stream it)

    command = f"gunzip -dc {file_path}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    # read the output in binary and manually decode it 
    line = process.stdout.readline()
    hd_index = line.find(b'@HD')
    if hd_index != -1:
        line = line[hd_index:]
        line = line.decode('ascii')
        fields = line.strip().split("\t")
        for field in fields:
            if field.startswith("SO:"):
                sorting_order = field.split(":")[1]
                process.stdout.close()
                process.terminate()
                return sorting_order
    process.stdout.close()
    process.terminate()
    # here is a slower version, try this if HD line is not at the top
    try:
        with pysam.AlignmentFile(file_path, "rb",require_index=False) as bamfile:
            hdline = bamfile.header.get('HD', {})
            sorting_order = hdline.get('SO', 'unknown')
            return sorting_order
    except:
        return "unknown"
    
def find_lca_type(original_lca_path):
    # is the lca file from ngslca output, or from metadmg output?
    # let's detect it then act appropriately 

    lcaheaderlines = 0
    with open(original_lca_path, 'r') as lcafile:
        for lcaline in lcafile:
            if "root" in lcaline and "#" not in lcaline:
                break
            lcaheaderlines += 1
    
    with open(original_lca_path,'r') as file:
        for _ in range(lcaheaderlines):
            next(file)
        try:
            firstline = next(file)
        except StopIteration:
            print("\nError: LCA file contains no data lines after header.")
            sys.exit()
        line = firstline.split('\t')
        # line[1] is the first tax id in an ngslca-style format, and the full read in metadmg-style format.
        if ":" in line[1]:
            filetype = "ngslca"
        else:
            filetype = "metadmg"
        if filetype!= "ngslca" and filetype!= "metadmg":
            print("\n Error: Your input lca file doesn't look like it's in lca file format")
            sys.exit()

    return filetype


def calculate_dust(seq): # input seq as a string
    # parameters as given by sga and original publication: Morgulis A. "A fast and symmetric DUST implementation to Mask Low-Complexity DNA Sequences". J Comp Bio.
    # ! but scaling here between 0 and 100 inclusive is as implemented in Prinseq, not in sga or Morgulis 
    # throws a -1 if the read has an N in it
    
    readlength = len(seq) 
    if readlength < 3:
        print(f"Warning: Cannot calculate dust score for a very short sequence (wait, why do you have reads this short?)")
        return 0
    
    w = 64
    k = 3
    firstwindowend = min(readlength,w)
    l = firstwindowend - 2
    maxpossibledust = l*(l-1)/2

    kmer_counts = {}
    for i in range(firstwindowend - k + 1):
        kmer = seq[i:i + k]
        if not all(base in {'A', 'C', 'T', 'G'} for base in kmer):
            # print(f"Warning: Skipping DUST calculations for a read with non-ACGT characters.")
            return -1
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    currentdust = sum((count * (count - 1)) / 2 for count in kmer_counts.values()) 

    if firstwindowend == readlength:
        #  read is less than window size
        return currentdust * (100 / maxpossibledust) 

    # otherwise perform sliding window :
    maxdust = currentdust
    for i in range(1, readlength - w +1):
        oldkmer = seq[(i-1) : (i+2)]
        newkmer = seq[(i+w-3): (i+w)]
        if not all(base in {'A', 'C', 'T', 'G'} for base in newkmer):
            # print(f"Warning: Skipping DUST calculations for a read with non-ACGT characters.")
            return -1
        kmer_counts[oldkmer] += -1
        if kmer_counts[oldkmer] == 0:
            del kmer_counts[oldkmer]
        if newkmer not in kmer_counts:
            kmer_counts[newkmer] = 1
        else:
            kmer_counts[newkmer] += 1
        currentdust = sum((count * (count - 1)) / 2 for count in kmer_counts.values())
        if currentdust > maxdust:
            maxdust = currentdust

    return maxdust * 100 / (maxpossibledust) #  standardize so it's max 100 

