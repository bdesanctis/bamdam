
# bamdam shrink module


import sys 
import re
import csv
import pysam
import math
import argparse
import os
from functools import lru_cache
try: # optional library for progress bars 
    from tqdm import tqdm
    tqdm_imported = True
except ImportError:
    tqdm_imported = False

# local
from bamdam import utils
from bamdam import alignment_utils

def run_shrink(args):
    formatted_exclude_tax = parse_exclude_tax(args)
    lca_file_type = utils.find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("You are running bamdam shrink with a metaDMG-style lca file. This is ok, but be aware the output lca file will be in ngsLCA lca file format, as all the other functions in bamdam require this format.")
    shortlcalines = write_shortened_lca(args.in_lca, args.out_lca, args.upto, args.mincount, formatted_exclude_tax, lca_file_type)
    write_shortened_bam(args.in_bam, args.out_lca, args.out_bam, args.stranded, args.minsim, args.annotate_pmd, shortlcalines, args.threads)


def parse_exclude_tax(args):
    if hasattr(args, 'exclude_tax') or hasattr(args, 'exclude_tax_file'):
        if args.exclude_tax and args.exclude_tax_file:
            raise ValueError("Please only provide one of --exclude_tax or --exclude_tax_file, not both.")
        
        exclude_tax = args.exclude_tax if args.exclude_tax else []
        if not isinstance(exclude_tax, list):
            exclude_tax = [exclude_tax]
        
        if args.exclude_tax_file:
            if not os.path.exists(args.exclude_tax_file):
                raise FileNotFoundError(f"exclude_tax_file path does not exist: {args.exclude_tax_file}")
            with open(args.exclude_tax_file, 'r') as f:
                exclude_tax.extend([line.strip() for line in f if line.strip()])
        
        exclude_tax = [kw.lstrip("'").lstrip('"').rstrip("'").rstrip('"') for kw in exclude_tax] # strip quotes 
        
        # formatted_keywords = []
        # good to surround the digit-only tax ids with a :, so we don't accidentally hit substring tax ids
        #for keyword in exclude_keywords:
        #    if keyword.isdigit():
        #        keyword = f"{keyword}:"
        #    formatted_keywords.append(keyword)

        return exclude_tax
    else:
        return []


def write_shortened_lca(original_lca_path,short_lca_path,upto,mincount,exclude_tax,lca_file_type): 

    print("\nWriting a filtered lca file...")

    lcaheaderlines = 0
    with open(original_lca_path, 'r') as lcafile:
        for lcaline in lcafile:
            if "root" in lcaline and "#" not in lcaline:
                break
            lcaheaderlines += 1
    total_short_lca_lines = 0

    # pass 1: make a dictionary with all the tax ids and their counts 
    number_counts = {}
    with open(original_lca_path, 'r') as file:
        for _ in range(lcaheaderlines):
            next(file) 
        for line in file:
            if upto in line: 
                entry = line.strip().split('\t')
                if len(entry) > 1:  
                    if lca_file_type=="ngslca":
                        fields = entry[1:] # can ditch the read id 
                    elif lca_file_type=="metadmg":
                        fields = entry[6].split(';')
                    if exclude_tax:
                        # go check if you need to skip this line
                        taxidd = fields[0].split(":")[0].strip("'").strip('"')
                        matching_tax = [keyword for keyword in exclude_tax if taxidd == keyword]
                        if len(matching_tax)>0:
                            continue # this only checks the node the read is actually assigned to 
                        # note! we are skipping keywords BEFORE aggregation, which happens later. so you might still end up with a family-level line if you
                        # specified that family in the "exclude tax" list, for example if there was a species in the sample which was not itself in your "exclude tax" list
                    # moving on
                    # now explicitly check if upto is in the line on its own (e.g. we see "family", not just "subfamily" - yes this can happen rarely and weirdly and we will not include them)
                    if any(item.split(":")[2].strip("'").strip('"') == upto for item in fields):
                        keepgoing = True; field = 0
                        while keepgoing:
                            taxid = fields[field].split(':')[0].strip("'").strip('"')
                            taxlevel = fields[field].split(':')[2].strip("'").strip('"')
                            if taxid in number_counts:
                                number_counts[taxid] += 1
                            else:
                                number_counts[taxid] = 1
                            if taxlevel == upto:
                                keepgoing = False
                            field +=1 

    goodnodes = [key for key, count in number_counts.items() if count >= mincount]
    # these are the nodes that have at least the min count of reads assigned to them (or below them), and which are at most upto

    # pass 2: rewrite lines into a new lca file that pass the filter
    oldreadname = ""
    with open(original_lca_path, 'r') as infile, open(short_lca_path, 'w') as outfile:
            for _ in range(lcaheaderlines):
                next(infile) 
            for line in infile:
                tab_split = line.find('\t')
                if lca_file_type == "ngslca":
                    newreadname = line[:tab_split].rsplit(':', 3)[0]
                elif lca_file_type == "metadmg":
                    newreadname = line[:tab_split]
                if newreadname == oldreadname:
                    print("Error: You have duplicate entries in your LCA file, for example " + newreadname + ". You should fix this and re-run bamdam. Here is a suggested fix: awk '!seen[$0]++' input_lca > deduplicated_lca")
                    exit(-1)
                if upto in line: 
                    entry = line.strip().split('\t')
                    if len(entry) > 1:  
                        if lca_file_type == "ngslca":
                            fields = entry[1:] 
                        elif lca_file_type == "metadmg":
                            fields = entry[6].split(';')
                        if exclude_tax:
                            # go check if you need to skip this line
                            matching_tax = [keyword for keyword in exclude_tax if fields[0].split(":")[0].strip("'").strip('"') == keyword]
                            if len(matching_tax)>0:
                                continue # this only checks the node the read is actually assigned to 
                        for field in fields:
                            number = field.split(':')[0].strip("'").strip('"') # metadmg lca files can have quotation marks everywhere
                            level = field.split(':')[2].strip("'").strip('"')
                            if level == upto:
                                if number in goodnodes: # you only need to check the upto counts, as they will be higher than anything underneath them 
                                    if lca_file_type == "ngslca":
                                        outfile.write(line)
                                    elif lca_file_type == "metadmg":
                                        # reformat the output lca as an ngslca file format no matter how it came in 
                                        firstentry = ":".join(entry[0:4])
                                        restentry = "\t".join(fields)
                                        fullentry = "\t".join([firstentry, restentry]).replace('"', '') + "\n"
                                        outfile.write(fullentry)
                                    total_short_lca_lines += 1
                                    break

    print("Wrote a filtered lca file. \n")

    return total_short_lca_lines


def write_shortened_bam(
    original_bam_path,
    short_lca_path,
    short_bam_path,
    stranded,
    minsimilarity,
    annotate_pmd,
    totallcalines,
    threads=1,
):
    # runs through the existing bam and the new short lca file at once, and writes only lines to the new bam which are represented in the short lca file
    # also optionally annotates with pmd scores as it goes
    # takes in minsimilarity as a percentage, and will keep reads w/ equal to or greater than 1 - NM flag divided by read length to this percentage 

    if tqdm_imported:
        print(f"Writing a filtered bam file (a progress bar will initiate once the bam header has been written)...")
    else:
        print(f"Writing a filtered bam file...")

    # go and get get header lines in the OUTPUT lca, not the input (it will be 0, but just in case I modify code in the future)
    lcaheaderlines = 0
    with open(short_lca_path, 'r') as lcafile:
        for lcaline in lcafile:
            if "root" in lcaline and "#" not in lcaline:
                break
            lcaheaderlines += 1
    currentlcaline = lcaheaderlines

    with (
        pysam.AlignmentFile(
            original_bam_path,
            "rb",
            check_sq=False,
            require_index=False,
            threads=threads,
        ) as infile,
        pysam.AlignmentFile(
            short_bam_path, "wb", template=infile, threads=threads
        ) as outfile,
        open(short_lca_path, "r") as shortlcafile,
    ):

        for _ in range(lcaheaderlines): 
            lcaline = next(shortlcafile)

        lcaline = next(shortlcafile)
        tab_split = lcaline.find('\t')
        colon_split = lcaline[:tab_split].rsplit(':', 3)[0]
        lcareadname = colon_split
    
        currentlymatching = False
        notdone = True
        bamreadnumber = 0

        try:
            bamread = next(infile)
        except StopIteration:
            notdone = False

        progress_bar = tqdm(total=totallcalines-lcaheaderlines, unit='lines', disable=not tqdm_imported) if tqdm_imported else None
        if progress_bar:
            update_interval = 100 # this is arbitrary. could also be e.g. (totallcalines // 1000)

        while notdone:
            if bamread.query_name == lcareadname:
                # copy this line and all the rest until you hit a nonmatching LCA line
                readlength = bamread.query_length # same read length for all the alignments
                similarity = 1 - bamread.get_tag('NM') / readlength # not the same NM for all the alignments
                if similarity >= minsimilarity:
                    if annotate_pmd:
                        pmd = alignment_utils.get_pmd(bamread, stranded)
                        bamread.set_tag('DS','%.3f' % pmd)  
                    outfile.write(bamread) # write the read!
                currentlymatching = True
                while currentlymatching:
                    try:
                        bamread = next(infile)
                        if bamread.query_name == lcareadname:
                            similarity = 1 - bamread.get_tag('NM') / readlength
                            if similarity >= minsimilarity:
                                if annotate_pmd:
                                    pmd = alignment_utils.get_pmd(bamread, stranded)
                                    bamread.set_tag('DS','%.3f' % pmd) # replace a tag if it's already there
                                outfile.write(bamread) # write the read! 
                        else:
                            currentlymatching = False
                    except StopIteration:
                        notdone = False
                        break
                try:
                    lcaline = next(shortlcafile)
                    tab_split = lcaline.find('\t')
                    lcareadname = lcaline[:tab_split].rsplit(':', 3)[0]
                    currentlcaline +=1 
                    if progress_bar and currentlcaline % update_interval == 0 :
                        progress_bar.update(update_interval) 
                except StopIteration:
                    notdone = False
            else:
                try:
                    bamread = next(infile)
                    bamreadnumber +=1 
                except StopIteration:
                    notdone = False
    if progress_bar:
        progress_bar.close()

    print("Wrote a filtered bam file. Done! \n") 

