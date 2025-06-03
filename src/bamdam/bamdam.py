#!/usr/bin/env python3

# bamdam by Bianca De Sanctis, bddesanctis@gmail.com  
# WORKING VERSION; last updated may 5 2025
# since last git commit: added GC content and more averaging to the combine function, fixed ^Ms at ends of tsv file lines, added an empty lca error catch

import sys 
import re
import csv
import pysam
import math
import argparse
import os
import hyperloglog
import subprocess
try: # optional library only needed for plotting 
    import matplotlib.pyplot as plt
    matplotlib_imported = True
except:
    matplotlib_imported = False
try: # optional library for progress bars 
    from tqdm import tqdm
    tqdm_imported = True
except ImportError:
    tqdm_imported = False

level_order = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus", "subgenus", "species", "subspecies"]

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

def write_shortened_lca(original_lca_path,short_lca_path,upto,mincount,exclude_keywords,lca_file_type): 

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
                    if exclude_keywords:
                        # go check if you need to skip this line
                        taxidd = fields[0].split(":")[0].strip("'").strip('"')
                        matching_keywords = [keyword for keyword in exclude_keywords if taxidd == keyword]
                        if len(matching_keywords)>0:
                            continue # this only checks the node the read is actually assigned to 
                        # note! we are skipping keywords BEFORE aggregation, which happens later. so you might still end up with a family-level line if you
                        # specified that family in the "exclude keywords" list, for example if there was a species in the sample which was not itself in your "exclude keywords" list
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
                        if exclude_keywords:
                            # go check if you need to skip this line
                            matching_keywords = [keyword for keyword in exclude_keywords if fields[0].split(":")[0].strip("'").strip('"') == keyword]
                            if len(matching_keywords)>0:
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


def write_shortened_bam(original_bam_path,short_lca_path,short_bam_path,stranded,minsimilarity,annotate_pmd,totallcalines): 
    # runs through the existing bam and the new short lca file at once, and writes only lines to the new bam which are represented in the short lca file
    # does two passes, the first of which makes a shortened header as well and adds the command str to the end of the bam header
    # also annotates with pmd scores as it goes
    # now takes in minsimilarity as a percentage, and will keep reads w/ equal to or greater than NM flag to this percentage 

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

    with pysam.AlignmentFile(original_bam_path, "rb", check_sq=False, require_index=False) as infile, \
        pysam.AlignmentFile(short_bam_path, "wb", header=infile.header) as outfile, \
        open(short_lca_path, 'r') as shortlcafile:

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
                        pmd = get_pmd(bamread, stranded)
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
                                    pmd = get_pmd(bamread, stranded)
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

def get_mismatches(seq, cigar, md):  
    # parses a read, cigar and md string to determine mismatches and positions. 
    # does not output info on insertions/deletions, but accounts for them.
    # thanks jonas oppenheimer who wrote half of this function :)

    # goes in two passes: once w/ cigar and once w/ md 

    # the strategy is to reconstruct the alignment through the cigar string first, so read_seq and ref_seq are the same length
    # but might have "-"s and "N"s floating around, and next to inform the substitutions from the md string 

    cig = re.findall(r'\d+\D', cigar)
    md_list = re.compile(r'\d+|\^[A-Za-z]+|[A-Za-z]').findall(md)

    ref_seq = ''
    read_seq = ''
    query_pos = 0 # indexes the ref reconstruction (which might have "-"s or "N"s added if there is soft clipping, indels etc)
    read_pos = 0 # indexes the read reconstruction (which is the input read but with potential added "-"s if the ref has insertions)
    for x in cig:
        cat = x[-1]
        if cat == 'H': # doesn't consume reference or query
            print(f"Warning: You have cigar strings have Hs in them (hard clipping). These specific reads may not be parsed correctly")
            continue
            # ! i've never actually seen these in a cigar string so i'm not 100% sure this works

        bases = int(x[:-1])

        if cat == 'S': # soft clip: bases are present in the read, but not the reference. soft clipping is enabled by default in bwa but not bowtie2.
            # pad the reference with "-"s and keep the whole read (to get the positions right).
            # add to the md_list to deal with this manually below (md tags do not record soft clipped regions)
            md_extra = "softclip"+str(bases)
            if x == cig[len(cig)-1]: # the end of the read is soft clipped
                md_list.append(md_extra)
            elif x == cig[0]: # the start of the read is soft clipped
                md_list.insert(0,md_extra)
            else:
                sys.exit("Error: It looks like one of your cigar strings " + cigar + " encodes a soft clip internally to a read. That shouldn't be happening.")

            # now pad the sequences as needed
            read_seq += 's' * bases # soft clipped. i don't care what you are. # seq[read_pos:read_pos + bases]
            # i don't think we should include the original bases in the read, because i don't want to count them downstream,
            # e.g. in kmer counts, gc count etc, because they didn't actually hit the reference!
            ref_seq += 's' * bases 
            query_pos += bases
            read_pos += bases
            continue

        elif cat in ['M', '=', 'X']: # match
            ref_seq += seq[query_pos:query_pos + bases]
            read_seq += seq[read_pos:read_pos + bases]
            query_pos += bases
            read_pos += bases

        elif cat in ['D', 'N']: # 'D' : the reference has something there and the read doesn't, but pad them both 
            ref_seq += 'N' * bases # N means it's an ACTG, we just don't know which one right now
            read_seq += '-' * bases

        elif cat == 'I': # I: the read has something there and the reference doesn't, but pad them both 
            read_seq += seq[read_pos:read_pos + bases] # 'N' * bases
            ref_seq += '-' * bases
            query_pos += bases
            read_pos += bases

        else:
            sys.exit("Error: You've got some strange cigar strings here.")
    
    rec_pos = 0 # reconstruction position. refers to position in read_seq and ref_seq. 
    mismatch_list = [] # of format [ref, read, pos in alignment] 
    for x in md_list:
        if x.startswith('softclip'):
            # skip ahead
            num = int(x[len("softclip"):])
            rec_pos += num
        elif x.startswith('^'): # this can be arbitrarily long (a ^ and then characters)
            num_chars_after_hat = len(x) -1 
            rec_pos += num_chars_after_hat
        else: 
            # you're a number or a letter
            if x.isdigit(): # can be multiple digits 
                # if you're a number, you're matching. skip ahead. don't need to add to the mismatch list.
                rec_pos += int(x)
            else: # it will only be one character at a time ; a mismatch
                refhere = x
                readhere = read_seq[rec_pos]
                char_list = list(ref_seq)
                char_list[rec_pos] = x
                ref_seq = "".join(char_list)
                # moving on. get the position in the actual read itself, not the reconstructed alignment
                read_up_to_here = read_seq[0:rec_pos]
                read_pos = len(read_up_to_here) - read_up_to_here.count("-") # we padded it with "-"s earlier
                mismatch_list.append([refhere,readhere,read_pos +1])
                # in genetics we are 1-based (a->g at position 1 means the actual first position, whereas python is 0-based) 
                if refhere is None or readhere is None:
                    print("Warning: There appears to be an inconsistency with seq " + seq + " and cigar " + cigar + " and md " + md)
                    break
                rec_pos += 1 

    # returns same-length read and ref, padded as needed with "-" and "s" for indels and soft clips respectively,
    # and a mismatch list which is on the coordinates of the read sequence alone. 
    # for example: read AGTTCTGAG, cigar 1S6M1D2M and md 2A3^G2 should yield (note md tags don't include soft clipped regions):
    # sGTTCTG-AG read_seq
    #  |||||| ||
    # sGTACTGNAG ref_seq
    # mismatch_list: [A T 4]. 

    return [mismatch_list, ref_seq, read_seq] 

def mismatch_table(read,cigar,md,flagsum):
    # wrapper for get_mismatches that also reverse complements if needed, and mirrors around the middle of the read so you shouldn't have to keep the length

    # first parse the mismatches
    mms, refseq, readseq = get_mismatches(read,cigar,md)

    readlength = len(readseq)
    # the mismatch list mms is on the coordinates of the read itself, NOT the reconstructed refseq or readseq.
    # so first, get the reconstructions on the coordinates of the read too:
    non_dash_indices = [i for i, char in enumerate(readseq) if char != "-"]
    readseq_no_read_dashes = "".join(readseq[i] for i in non_dash_indices)
    refseq_no_read_dashes = "".join(refseq[i] for i in non_dash_indices)
    matchs = []
    for i, (read_char, ref_char) in enumerate(zip(readseq_no_read_dashes, refseq_no_read_dashes)):
        if read_char == ref_char and read_char in "ACTG":  # only count matches if both are an A C T or G 
            pos = non_dash_indices[i] + 1  # 1 based
            matchs.append([ref_char, read_char, pos])

    # some processing: first, figure out if it's forward or backwards
    # second field of a bam is a flagsum and it parses like this (see eg bowtie2 manual)
    bin_flagsum = bin(flagsum)[::-1]
    bit_position = 4 #  2^4 = 16 this flag means it's aligned in reverse
    backwards = len(bin_flagsum) > bit_position and bin_flagsum[bit_position] == '1' 

    if backwards: # flip all the positions and reverse complement all the nucleotides. the read in the bam is reverse-complemented if aligned to the reverse strand. 
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        mmsc = []
        matchsc = []

        for entry in range(0,len(mms)): # mismatches
            new_entry = [
                complement.get(mms[entry][0],"N"), 
                complement.get(mms[entry][1],"N"), 
                readlength - mms[entry][2] +1
            ]
            mmsc.append(new_entry)

        for entry in range(0, len(matchs)): # matches
            new_entry = [
                complement.get(matchs[entry][0],"N"), 
                complement.get(matchs[entry][1],"N"), 
                readlength - matchs[entry][2] + 1
            ]
            matchsc.append(new_entry)
    else:
        mmsc = mms
        matchsc = matchs
    # now accounts for everything EXCEPT unmerged but retained reverse mate pairs of paired end reads (which i should still try to catch later; maybe to do), should be 5' -> 3' 

    # mirroring for accurate cumulative damage counting later (end of read becomes -1, -2...)
    for entry in range(0,len(mmsc)): # mismatches
        pos = mmsc[entry][2]
        if pos > readlength/2:
            mmsc[entry][2] = -(readlength - pos +1)
            
    for entry in range(0, len(matchsc)): # matches
        pos = matchsc[entry][2]
        if pos > readlength / 2:
            matchsc[entry][2] = -(readlength - pos + 1)

    return mmsc, matchsc, refseq


def rev_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
    return ''.join(complement[base] for base in reversed(seq))

def get_rep_kmer(seq): # representative canonical kmer representation for counting, gets the lexicographical min of a kmer and its rev complement
    rep_kmer = min(seq,rev_complement(seq))
    return rep_kmer
 
def get_pmd(read, stranded):
    ### important!!! the original PMDtools implementation has a bug:
    # in lines 868 and 900 of the main python script, it multiplies the likelihood by the double-stranded models, and then there is an if statement that multiplies by the single-stranded models
    # resulting in totally incorrect pmd scores for single-stranded mode. the reimplementation here should be correct.
    
    # input is a pysam read object
    seq = read.query_sequence
    cigar = read.cigarstring
    md = read.get_tag('MD')
    rawphred = read.query_qualities
    flagsum = read.flag

    # set pmd score parameters . these are their original parameters, and i need to do some testing, but i think they are sensible enough in general. 
    P = 0.3
    C = 0.01
    pi = 0.001 

    # important note! in the PMDtools manuscript, they say "DZ=0" in the null model.
    # however, in the PMDtools code, Dz=0.001 in the null model.
    # here i am making the latter choice because i think it makes more biological sense.
    Dn = 0.001

    # find out if you're backwards
    bin_flagsum = bin(flagsum)[::-1]
    bit_position = 4 #  2^4 = 16 this flag means it's aligned in reverse
    backwards = len(bin_flagsum) > bit_position and bin_flagsum[bit_position] == '1' # backwards is a boolean 
    # do something if you are:
    mmsc, refseq, readseq = get_mismatches(seq, cigar, md)
    # adjust phred index if there are things happening in the read
    phred = [0 if base in ['-', 'N'] else rawphred.pop(0) for base in readseq]

    # run through both sequences to add up pmd likelihoods
    refseqlist = list(refseq)
    readseqlist = list(readseq)
    readlength = len(readseqlist)
    if backwards:
        refseqlist = rev_complement(refseqlist)
        readseqlist = rev_complement(readseqlist)
        phred = phred[::-1]
    pmd_lik = 1
    null_lik = 1
    pos = 0 # need a separate tracker to cope with indels. if there's a "-" in the read reconstruction because of an insertion in the ref, it should not count as a "position" in the read

    if stranded == "ss":
        for b in range(0,readlength):
            if readseqlist[b] == "-":
                continue # no pos +=1 
            # looking for c-> anything anywhere
            if refseqlist[b] == "C" and (readseqlist[b] == "T" or readseqlist[b] == "C"):
                # everything is relevant to 5 prime and 3 prime ends, get both distances
                epsilon = 1/3 * 10**(-phred[b] / 10)
                z = pos + 1 # pos 1 has distance 1, from pmd manuscript
                y = readlength - pos
                Dz = ((1-P)**(z-1))*P + C
                Dy = ((1-P)**(y-1))*P + C
                if readseqlist[b] == "T": # ss m
                    pmd_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dz)*(1-Dy) + (1-pi)*epsilon*Dz*(1-Dy) + (1-pi)*epsilon*Dy*(1-Dz) + pi*epsilon*(1-Dz)*(1-Dy))
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn)*(1-Dn) + (1-pi)*epsilon*Dn*(1-Dn) + (1-pi)*epsilon*Dn*(1-Dn) + pi*epsilon*(1-Dn)*(1-Dn))                 
                if readseqlist[b] == "C": # ss m
                    pmd_lik *= (1-pi)*(1-epsilon)*(1-Dz)*(1-Dy) + (1-pi)*epsilon*Dz*(1-Dy) + (1-pi)*epsilon*Dy*(1-Dz) + pi*epsilon*(1-Dz)*(1-Dy)
                    null_lik *= (1-pi)*(1-epsilon)*(1-Dn)*(1-Dn) + (1-pi)*epsilon*Dn*(1-Dn) + (1-pi)*epsilon*Dn*(1-Dn) + pi*epsilon*(1-Dn)*(1-Dn)
                pos +=1 

    if stranded == "ds":
        for b in range(0,readlength):
            if readseqlist[b] == "-":
                continue # no pos +=1 
            if refseqlist[b] == "C" and (readseqlist[b] == "T" or readseqlist[b] == "C"):
                # get distance and stuff to 5 prime end
                epsilon = 1/3 * 10**(-phred[b] / 10)
                z = pos + 1   # 5 prime
                Dz = ((1-P)**(z-1))*P + C

                if readseqlist[b] == "T": # ds mm 
                    pmd_lik *=  1 - ((1-pi)*(1-epsilon)*(1-Dz) + (1-pi)*epsilon*Dz + pi*epsilon*(1-Dz)) 
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn) + (1-pi)*epsilon*Dn + pi*epsilon*(1-Dn)) 
                
                if readseqlist[b] == "C": # ds match
                    pmd_lik *= (1-pi)*(1-epsilon)*(1-Dz) + (1-pi)*epsilon*Dz + pi*epsilon*(1-Dz)
                    null_lik *= (1-pi)*(1-epsilon)*(1-Dn) + (1-pi)*epsilon*Dn + pi*epsilon*(1-Dn)

            if refseqlist[b] == "G" and (readseqlist[b] == "A" or readseqlist[b] == "G"):
                # get distance and stuff to 3 prime end
                epsilon = 1/3 * 10**(-phred[b] / 10) # phred score 30 gives an error rate of 0.001 (then * 1/3)
                z = readlength - pos  # 3 prime
                Dz = ((1-P)**(z-1))*P + C
                if readseqlist[b] == "A": # ds mm
                    pmd_lik *=  1 - ((1-pi)*(1-epsilon)*(1-Dz) + (1-pi)*epsilon*Dz + pi*epsilon*(1-Dz)) 
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn) + (1-pi)*epsilon*Dn + pi*epsilon*(1-Dn)) 
                if readseqlist[b] == "G": # ds m
                    pmd_lik *= (1-pi)*(1-epsilon)*(1-Dz) + (1-pi)*epsilon*Dz + pi*epsilon*(1-Dz)
                    null_lik *= (1-pi)*(1-epsilon)*(1-Dn) + (1-pi)*epsilon*Dn + pi*epsilon*(1-Dn)
            
            pos +=1 
    
    if pmd_lik == 0 or null_lik == 0:
        pmd_score = 0
    else:
        pmd_score = math.log(pmd_lik/null_lik)

    return pmd_score

def line_count(file_path): # just a wc -l wrapper
    result = subprocess.run(['wc', '-l', file_path], stdout=subprocess.PIPE, text=True)
    line_count = int(result.stdout.split()[0])
    return line_count

def calculate_dust(seq):
    # parameters as given by sga and original publication: Morgulis A. "A fast and symmetric DUST implementation to Mask Low-Complexity DNA Sequences". J Comp Bio.
    # between 0 and 100 inclusive; throws a -1 if the read has an N in it
    
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

def get_hll_info(seq,k):
    # output to dump into hll objects
    rep_kmers = []
    total_kmers = 0
    if len(seq) > k:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if not all(base in {'A', 'C', 'T', 'G'} for base in kmer):         
                #print(f"Warning: Skipping k-mer calculations for a read with non-ACGT characters.")       
                continue # skip this k-mer, non ACTG characters are not allowed
            else:
                rep_kmers.append(get_rep_kmer(kmer))
                total_kmers +=1 
    else:
        print(f"Warning: One of your reads is shorter than k.")
    return rep_kmers, total_kmers

def gather_subs_and_kmers(bamfile_path, lcafile_path, kn, upto,stranded):
    print("\nGathering substitution and kmer metrics per node...")
    # this function is organized in an unintuitive way. it uses a bunch of nested loops to pop between the bam and lca files line by line. 
    # it matches up bam read names and lca read names and aggregates some things per alignment, some per read, and some per node, the last of which are added into a large structure node_data.
    # altogether this uses very little ram

    lcaheaderlines = 0
    with open(lcafile_path, 'r') as lcafile:
        for lcaline in lcafile:
            if "root" in lcaline and "#" not in lcaline:
                break
            lcaheaderlines += 1

    # initialize 
    node_data = {}
    bamfile = pysam.AlignmentFile(bamfile_path, "rb",require_index=False) 
    lcafile = open(lcafile_path, 'r')
    oldreadname = ""
    oldmd = ""
    oldcigar = ""
    oldflagsum = ""
    nodestodumpinto = []
    num_alignments = 0
    currentsubdict = {}
    nms = 0 
    pmdsover2 = 0
    pmdsover4 = 0
    are_pmds_in_the_bam = True # just assume true, then set to false quickly below if you notice otherwise 
    lcalinesskipped = 0
    readswithNs = 0

    currentlcalinenum = 0
    if tqdm_imported:
        totallcalines = line_count(lcafile_path)
        # should be super fast compared to anything else; probably worth it to initiate a progress bar.
        progress_bar = tqdm(total=totallcalines-lcaheaderlines, unit='lines', disable=not tqdm_imported) if tqdm_imported else None
        if progress_bar:
            update_interval = 100 # this is arbitrary ofc. could also be e.g. (totallcalines // 1000) 

    for _ in range(lcaheaderlines +1):
        currentlcaline = next(lcafile) 

    for read in bamfile:

        # get the basic info for this read
        readname = read.query_name

        # find out if it's a new read. if so, you just finished the last read, so do a bunch of stuff for it. 
        # the first read will skip this if statement because of the second condition
        if readname != oldreadname and oldreadname != "":

            # do k-mer things for this read
            dust = calculate_dust(seq)
            # then get all the rep kmers to dump into the hyperloglog for each relevant node below
            rep_kmers, total_kmers = get_hll_info(seq,kn)

            # get the lca entry and nodes we wanna update
            lcaentry = currentlcaline.split('\t')
            colon_split = lcaentry[0].rsplit(':', 3)[0]
            lcareadname = colon_split

            while oldreadname != lcareadname:
                # skip em but don't forget your progress bar
                # this must be because there are reads in the lca which are not in the bam. presumably because we didn't write them because none of them met the similarity cutoff.
                lcalinesskipped += 1
                currentlcaline = next(lcafile)
                currentlcalinenum += 1 
                lcaentry = currentlcaline.split('\t')
                colon_split = lcaentry[0].rsplit(':', 3)[0]
                lcareadname = colon_split
                if progress_bar and (currentlcalinenum % update_interval == 0) :
                    progress_bar.update(update_interval)  # update the progress bar

            fields = lcaentry[1:]
            nodestodumpinto = []
            for i in range(len(fields)):
                # i recently changed fields to level here, it should be safer
                splitfields = fields[i].split(":")
                level = splitfields[2].strip("'").strip('"')
                nodename = splitfields[0].strip("'").strip('"')
                if level == upto: 
                    nodestodumpinto.append(nodename)
                    break
                nodestodumpinto.append(nodename)

            # now update everything to all the relevant nodes
            for node in nodestodumpinto:

                # you will skip this if statement if your node already exists; otherwise just initialize it then move on
                if node not in node_data:
                    if are_pmds_in_the_bam:
                        node_data[node] = {'total_reads': 0,'pmdsover2': 0, 'pmdsover4': 0, 'meanlength': 0, 'total_alignments': 0, 
                            'ani': 0, 'avgdust' : 0, 'avgreadgc': 0, 'tax_path' : "", 'subs': {}, 'damagerelevantreadsp1': 0, 'damagerelevantreadsm1': 0,
                            'dp1' : 0, 'dm1' : 0,  'hll': hyperloglog.HyperLogLog(0.01), 'totalkmers' : 0}
                    else:
                        node_data[node] = {'total_reads': 0,'meanlength': 0, 'total_alignments': 0, 
                            'ani': 0, 'avgdust' : 0, 'avgreadgc': 0, 'tax_path' : "", 'subs': {}, 'damagerelevantreadsp1': 0, 'damagerelevantreadsm1': 0,
                            'dp1' : 0, 'dm1' : 0, 'hll': hyperloglog.HyperLogLog(0.01), 'totalkmers' : 0}
                        
                # now populate/update it        
                node_data[node]['meanlength'] = ((node_data[node]['meanlength'] * node_data[node]['total_reads']) + readlength) / (node_data[node]['total_reads'] + 1)

                if dust != -1:
                    node_data[node]['avgdust'] = ( (node_data[node]['avgdust'] * node_data[node]['total_reads']) + dust) / (node_data[node]['total_reads'] + 1)
                else:
                    readswithNs +=1 

                ani_for_this_read = (readlength - nms/num_alignments)/readlength 
                node_data[node]['ani'] = (ani_for_this_read + node_data[node]['ani'] * node_data[node]['total_reads']) / (node_data[node]['total_reads'] + 1)
                gc_content_for_this_read = (seq.count('C') + seq.count('G')) / readlength
                node_data[node]['avgreadgc'] = ((node_data[node]['avgreadgc'] * node_data[node]['total_reads']) + gc_content_for_this_read) / (node_data[node]['total_reads'] + 1)
                node_data[node]['total_alignments'] += num_alignments

                if are_pmds_in_the_bam:
                    node_data[node]['pmdsover2'] += pmdsover2 / num_alignments
                    node_data[node]['pmdsover4'] += pmdsover4 / num_alignments

                # update hyperloglogs
                for kmer in rep_kmers:
                    node_data[node]['hll'].add(kmer)
                node_data[node]['totalkmers'] += total_kmers

                # updates substitution tables similarly
                other_sub_count = 0
                if currentsubdict:
                    for sub, count in currentsubdict.items():
                        if not ((sub[0] == 'C' and sub[1] == 'T') or (sub[0] == 'G' and sub[1] == 'A')):
                            other_sub_count += count # don't include c>t or g>a in any case, regardless of library 
                        if sub in node_data[node]['subs']: 
                            node_data[node]['subs'][sub] += count / num_alignments
                        else:
                            node_data[node]['subs'][sub] = count / num_alignments # so, this can be up to 1 per node. 
                # add the tax path if it's not already there
                if node_data[node]['tax_path'] == "":
                    try:
                        lca_index = next(i for i, entry in enumerate(lcaentry) if entry.split(":")[0].strip("'").strip('"') == node)
                        tax_path = ';'.join(lcaentry[lca_index:]).replace('\n','')
                        node_data[node]['tax_path'] = tax_path
                    except StopIteration: # this should not happen
                        print(f"Error: Something weird has gone wrong. Cannot find node '{node}' in its supposed lca entry. Are there weird characters in your lca entries?")
                        print(f"The problematic line is {currentlcaline}")
                        print(f"Will try to continue.")
                
                # only at the end should you update total reads 
                node_data[node]['total_reads'] += 1

            # move on to the next lca entry. re initialize a bunch of things here 
            oldreadname = readname
            oldmd = ""
            oldcigar = ""
            oldflagsum = ""
            currentlcaline = next(lcafile)
            currentlcalinenum += 1 
            if progress_bar and (currentlcalinenum % update_interval == 0) :
                progress_bar.update(update_interval)  # Update the progress bar
            currentsubdict = {}
            num_alignments = 0
            nms = 0
            if are_pmds_in_the_bam:
                pmdsover2 = 0
                pmdsover4 = 0

        # now for the current alignment.
        # the following might change for different alignments of the same read: 
        seq = read.query_sequence
        readlength = len(seq)
        cigar = read.cigarstring
        md = read.get_tag('MD')
        nms += read.get_tag('NM')
        if are_pmds_in_the_bam: 
            try:
                pmd = float(read.get_tag('DS'))
            except KeyError:
                pmd = 0
            if(pmd>2):
                pmdsover2 += 1 
            if(pmd>4):
                pmdsover4 += 1
        flagsum = read.flag
        num_alignments += 1 

        # go and get the mismatch table for this read if the name/md/cigar/flagsum is different to before (this is expensive, so there is a catch to avoid it when possible)
        if (readname != oldreadname) or (cigar != oldcigar) or (md != oldmd) or (flagsum != oldflagsum):
            subs, matches, refseq = mismatch_table(seq,cigar,md,flagsum) 
            oldcigar = cigar; oldmd = md; oldflagsum = flagsum

        allsubs = subs + matches
        for sub in allsubs:
            key = "".join(str(sub))
            if key in currentsubdict:
                currentsubdict[key] +=1
            else:
                currentsubdict[key] = 1

        # quick catch for the starting read; check if the first read (and then presumably the whole bam) has a pmd score 
        if oldreadname == "":
            oldreadname = readname
            try:
                pmd = float(read.get_tag('DS'))
            except KeyError:
                are_pmds_in_the_bam = False

    if progress_bar:
        progress_bar.close()

    bamfile.close() 
    lcafile.close()

    if lcalinesskipped > 0:
        print("\nWarning: " + str(lcalinesskipped) + " reads in the input LCA file did not appear in the input bam file and so were not used. This may have happened if the minimum similarity used in bamdam shrink did not match that used in ngsLCA. This will not affect output statistics, except that these reads will not be included.")
        
    if readswithNs >0 :
        print("\nSkipped k-mer counting and DUST score computation for " + str(readswithNs) + " reads with non-ACGT characters. \n" +
            "If all reads assigned to a taxonomic node have non-ACGT characters, both the mean DUST and duplicity will be 0 for that node.")

    print("\nGathered substitution and kmer data for " + str(len(node_data)) + " taxonomic nodes. Now sorting and writing output files... ")

    return node_data, are_pmds_in_the_bam

def format_subs(subs, nreads):
    formatted_subs = []
    
    for key, value in subs.items():
        # extract position and check if it is within the range -15 to 15 (more is unnecessary for damage)
        # easy to remove the condition if you want to write all of the subs though!
        parts = key.strip("[]").replace("'", "").split(", ")
        pos = int(parts[2])
        if (-15 <= pos <= 15) and (parts[0] in {'A', 'C', 'T', 'G'}) and (parts[1] in {'A', 'C', 'T', 'G'}):
            formatted_key = "".join(parts)
            formatted_value = round(value / nreads, 3)  
            formatted_subs.append((pos, f"{formatted_key}:{formatted_value}"))
            formatted_subs.sort(key=lambda x: (x[0] > 0, (x[0])))

    return " ".join(sub[1] for sub in formatted_subs)

def calculate_node_damage(subs, stranded):
    # also in here calculate the avg gc content of the reference

    ctp1 = 0  # C>T at 5' position 1
    ctm1 = 0  # C>T at 3' position -1 for ss
    gam1 = 0  # G>A at 3' position -1 for ds
    c_p1 = 0  # total C at 5' position 1
    c_m1 = 0  # total C at 3' position -1 for ss
    g_m1 = 0  # total G at 3' position -1 for ds
    total_gc = 0 # gc content in ref 
    total_bases = 0

    for key, count in subs.items():
        parts = key.strip("[]").replace("'", "").split(", ")
        
        from_base, to_base, pos = parts[0], parts[1], int(parts[2])

        # 5' end: check for C>T and all Cs at position +1
        if from_base == 'C' and pos == 1:
            if to_base == 'T':
                ctp1 += count  # C>T at position 1
            c_p1 += count  # all Cs at position 1

        # 3' end (ss): check for C>T and all Cs at position -1
        if stranded == "ss" and from_base == 'C' and pos == -1:
            if to_base == 'T':
                ctm1 += count  # C>T at position -1
            c_m1 += count  # all Cs at position -1

        # 3' end (ds): check for G>A and all Gs at position -1
        if stranded == "ds" and from_base == 'G' and pos == -1:
            if to_base == 'A':
                gam1 += count  # G>A at position -1
            g_m1 += count  # all Gs at position -1

        # calculate total GC content on the reference
        if from_base == 'C' or from_base == 'G':
            total_gc += count  

        total_bases += count  

    avgrefgc = total_gc / total_bases if total_bases > 0 else 0

    dp1 = ctp1 / c_p1 if c_p1 > 0 else 0
    dm1 = (ctm1 / c_m1 if c_m1 > 0 else 0) if stranded == "ss" else (gam1 / g_m1 if g_m1 > 0 else 0)

    return dp1, dm1, avgrefgc

def parse_and_write_node_data(nodedata, tsv_path, subs_path, stranded, pmds_in_bam):
    # parses a dictionary where keys are node tax ids, and entries are total_reads, meanlength, total_alignments, etc 

    statsfile = open(tsv_path, 'w', newline='')
    subsfile = open(subs_path, 'w', newline='')
    if pmds_in_bam:
        header = ['TaxNodeID', 'TaxName', 'TotalReads','Duplicity', 'MeanDust','Damage+1', 'Damage-1','MeanLength', 'ANI','AvgReadGC','AvgRefGC',
                   'UniqueKmers', 'RatioDupKmers', 'TotalAlignments', 'PMDsover2', 'PMDSover4','taxpath'] 
    else:
        header = ['TaxNodeID', 'TaxName', 'TotalReads','Duplicity', 'MeanDust','Damage+1', 'Damage-1','MeanLength', 'ANI','AvgReadGC','AvgRefGC',
                   'UniqueKmers', 'RatioDupKmers', 'TotalAlignments', 'taxpath']
    statsfile.write("\t".join(header) + "\n")
    writer = csv.writer(
        statsfile,
        delimiter="\t",
        quotechar='"',
        quoting=csv.QUOTE_NONNUMERIC,
        lineterminator="\n",
    )
    subswriter = csv.writer(
        subsfile,
        delimiter="\t",
        quotechar='"',
        quoting=csv.QUOTE_NONE,
        lineterminator="\n",
    )

    rows = []
    subsrows = {}

    for node in nodedata:
        tn = nodedata[node]
        
        # get formatted subs
        fsubs = format_subs(tn['subs'], tn['total_reads'])

        # number of unique k-mers approximated by the hyperloglog algorithm 
        numuniquekmers = len(tn['hll'])
        if numuniquekmers > 0:
            duplicity = tn['totalkmers'] / numuniquekmers
            ratiodup = (tn['totalkmers'] - numuniquekmers) / tn['totalkmers']
        else:
            duplicity = 0
            ratiodup = 0
        if ratiodup < 0 :
            # you can get tiny negative numbers from essentially zero duplicity and error 
            # (there is up to 1% error in the uniq kmer counting)
            ratiodup = 0

        taxname = tn['tax_path'].split(";")[0].split(":")[1]

        # only now calculate damage per node from the subs dict (see output subs file)
        # do not use formatted subs ; these should be raw numbers: how many READS for this taxa have these matches/mismatches?
        # # (avg'd over all the alignments per read, so maybe not an integer) 
        # also calculate the gc content for the average ref (so, unbiased by damage), weighted equally by read, not alignment
        dp1, dm1, avgrefgc  = calculate_node_damage(tn['subs'], stranded)

        # write 
        if pmds_in_bam:
            row = [int(node), taxname, tn['total_reads'], round(duplicity, 3), round(tn['avgdust'], 2), 
                round(dp1, 4), round(dm1, 4), 
                round(tn['meanlength'], 2),  round(tn['ani'], 4), round(tn['avgreadgc'], 3),round(avgrefgc, 3), numuniquekmers, round(ratiodup,3),tn['total_alignments'],
                round(tn['pmdsover2'] / tn['total_reads'], 3), round(tn['pmdsover4'] / tn['total_reads'], 3), tn['tax_path']]
        else:
            row = [int(node), taxname, tn['total_reads'], round(duplicity, 2), round(tn['avgdust'], 2), 
                round(dp1, 4), round(dm1, 4), 
                round(tn['meanlength'], 2),  round(tn['ani'], 4), round(tn['avgreadgc'], 3),round(avgrefgc, 3),
                numuniquekmers, round(ratiodup,3), tn['total_alignments'], tn['tax_path']]
        rows.append(row)

        subsrows[int(node)] = [int(node), taxname, fsubs]
    
    rows.sort(key=lambda x: x[2], reverse=True)

    for row in rows:
        writer.writerow(row)
    
    for row in rows:
        subswriter.writerow(subsrows[row[0]])

    statsfile.close()
    subsfile.close()

    print("Wrote final tsv and subs files. Done!")

def extract_reads(in_lca, in_bam, out_bam, tax, subset_header = False, only_top_ref = False):
    # extracts all reads with a tax path containing a certain keyword.
    # also optionally shortens the header to only necessary ids.
    # subsetting the header is kinda slow because it requires running through the input twice.

    if only_top_ref and not subset_header:
        print("Error: The --only_top_ref flag requires the --subset_header flag to be set. Please set both flags and re-run.")
        sys.exit()

    # parse keyword.
    # if you gave in a tax id as a digit, you probably are referring to the tax id and don't want to also get paths with the keyword as a substring of another tax id (eg you gave 200 and you get 2001)
    tax_pattern = f"[[:space:]]{tax}:" if tax.isdigit() else tax

    # if you don't need to subset the header and you're happy with keeping all the refs, this is easy.
    # in this case, the fastest way to do this is to pass samtools a file of read names. 
    if not subset_header and not only_top_ref:
        # essentially this is a wrapper for grep / awk / samtools view (with the samtools inside pysam).
        # writes (then later deletes) a temp file with readnames, because that's what samtools view wants. 
        # actually i did try to do this more cleverly, using the known ordering of bam and lca, but it was slower.
        hashtax = abs(hash(tax_pattern)) # ensures no accidental file overwriting; the abs is because sometimes the hash is negative
        print("Writing a temp file to the current directory. Will delete when done.")
        tmp_file = f"{out_bam}.{hashtax}.readnames.tmp"
        # get all the associated read names with the tax keyword 
        grep_command = (
            f"grep '{tax_pattern}' {in_lca} | awk -F'\\t' '{{print $1}}' "
            f"| awk -F':' '{{for(i=1; i<=NF-3; i++) printf $i (i<NF-3 ? \":\" : \"\\n\")}}' > {tmp_file}"
        )
        # writes all the read names to a temp file 
        result = subprocess.run(grep_command, shell=True)
        if result.returncode != 0:
            print(f"No matches found for keyword: {tax} or grep/awk command failed.")
            return  
        pysam.view("-N", tmp_file, "-b", in_bam, "-o", out_bam, catch_stdout=False)
        try:
            os.remove(tmp_file)
        except OSError as e:
            print(f"Error removing temp file: {tmp_file} : {e.strerror}")
        return

    # otherwise, we gotta deal with keeping the read names in memory. go get all the relevant references
    grep_command = (
        f"grep '{tax_pattern}' {in_lca} | awk -F'\\t' '{{print $1}}' "
        f"| awk -F':' '{{for(i=1; i<=NF-3; i++) printf $i (i<NF-3 ? \":\" : \"\\n\")}}'"
    )
    result = subprocess.run(grep_command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"No matches found for keyword: {tax} or grep/awk command failed.")
        return
    read_names = set(result.stdout.strip().splitlines())

    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        header = bam_in.header.to_dict()
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
    else: # get all the headers matching all of the refs
        reference_names = set(reference_count.keys())
        header['SQ'] = [sq for sq in header.get('SQ', []) if sq['SN'] in reference_names]

    # important step: we have to re-link the reference IDs in each read row to the new header because of how the bam compression works
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

def calculate_damage_for_plot(items):
    # specific to plotdamage

    ctp_5prime = {i: 0 for i in range(1, 16)}  # C>T at 5' positions 1 to 15
    gap_5prime = {i: 0 for i in range(1, 16)}  # G>A at 5' positions 1 to 15
    total_c_5prime = {i: 0 for i in range(1, 16)}  # total C at 5' positions 1 to 15
    total_g_5prime = {i: 0 for i in range(1, 16)}  # total G at 5' positions 1 to 15

    ctp_3prime = {i: 0 for i in range(1, 16)}  # C>T at 3' positions -1 to -15
    total_c_3prime = {i: 0 for i in range(1, 16)}  # total C at 3' positions -1 to -15
    gap_3prime = {i: 0 for i in range(1, 16)}  # G>A at 3' positions -1 to -15
    total_g_3prime = {i: 0 for i in range(1, 16)}  # total G at 3' positions -1 to -15

    other_5prime = {i: 0 for i in range(1, 16)}  # other mismatches at 5' positions (not including C>T or G>A)
    other_3prime = {i: 0 for i in range(1, 16)}  # other mismatches at 3' positions (not including C>T or G>A)
    total_nondeam_5prime = {i: 0 for i in range(1, 16)} #  denominator part : all non C>T or G>A matches or mismatches at each position
    total_nondeam_3prime = {i: 0 for i in range(1, 16)} 

    for item in items:
        # split each item into mutation and proportion parts
        mutation, proportion = item.split(":")
        from_base, to_base, pos = mutation[:1], mutation[1:2], int(mutation[2:])
        count = float(proportion)  

        # 5'
        if 1 <= pos <= 15:
            if from_base == 'C':
                if to_base == 'T':
                    ctp_5prime[pos] += count  # C>T at 5' positions 1 to 15
                total_c_5prime[pos] += count  
            elif from_base == 'G':
                if to_base == 'A':
                    gap_5prime[pos] += count  # G>A at 5' positions 1 to 15
                total_g_5prime[pos] += count  

            if from_base == to_base:
                total_nondeam_5prime[pos] += count
            if from_base != to_base:
                if (from_base == 'C' and to_base == 'T') or (from_base == 'G' and to_base == 'A'):
                    continue 
                else:
                    other_5prime[pos] += count
                    total_nondeam_5prime[pos] += count

        # 3' 
        elif -15 <= pos <= -1:
            abs_pos = abs(pos)
            if from_base == 'C':
                if to_base == 'T':
                    ctp_3prime[abs_pos] += count  # C>T at 3' positions -1 to -15
                total_c_3prime[abs_pos] += count  
            elif from_base == 'G':
                if to_base == 'A':
                    gap_3prime[abs_pos] += count  # G>A at 3' positions -1 to -15
                total_g_3prime[abs_pos] += count 
            
            if from_base == to_base:
                total_nondeam_3prime[abs_pos] += count
            if from_base != to_base:
                if (from_base == 'C' and to_base == 'T') or (from_base == 'G' and to_base == 'A'):
                    continue 
                else:
                    other_3prime[abs_pos] += count
                    total_nondeam_3prime[abs_pos] += count

    proportion_ct_5prime = {i: ctp_5prime[i] / total_c_5prime[i] if total_c_5prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_ga_5prime = {i: gap_5prime[i] / total_g_5prime[i] if total_g_5prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_ct_3prime = {i: ctp_3prime[i] / total_c_3prime[i] if total_c_3prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_ga_3prime = {i: gap_3prime[i] / total_g_3prime[i] if total_g_3prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_other_5prime = {i: other_5prime[i] / total_nondeam_5prime[i] if total_nondeam_5prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_other_3prime = {i: other_3prime[i] / total_nondeam_3prime[i] if total_nondeam_3prime[i] > 0 else 0 for i in range(1, 16)}

    return proportion_ct_5prime, proportion_ga_5prime, proportion_ct_3prime, proportion_ga_3prime, proportion_other_5prime, proportion_other_3prime
    

def make_damage_plot(in_subs_list, in_subs, tax, plotfile, ymax=0):
    # just damage from the subs file; should be super fast 

    if not matplotlib_imported:
        print("Error: Cannot find matplotlib library for plotting. Try: pip install matplotlib")
        return

    # did we get one or more files? check they exist then parse the input style 
    subs_files = []
    if in_subs:
        subs_files = in_subs
    elif in_subs_list:
        with open(in_subs_list, 'r') as file:
            subs_files = [line.strip() for line in file if line.strip()]   

    # validate input files: check existence
    for file in subs_files:
        if not os.path.exists(file):
            print(f"Error: File {file} does not exist.")
            continue

    positions = list(range(-15, 0)) + list(range(1, 16))
    values_all_files = {key: [] for key in ['Other', 'CT', 'GA']}  # store data for all files

    for file in subs_files:
        with open(file, 'r') as f:
            file_content = f.readlines()

        tax_pattern = f"^{tax}\\t" if tax.isdigit() else tax
        matched_line = [line for line in file_content if re.match(tax_pattern, line)]

        if len(matched_line) == 0:
            print(f"Warning: {file} does not contain an entry for {tax}. Skipping this file.")
            continue
        elif len(matched_line) > 1:
            raise ValueError(f"More than one matching line found in {file}. Please be more specific, e.g., by using a tax ID instead of a name.")

        split_line = matched_line[0].split('\t')
        tax_id, tax_name, data_part = split_line[0], split_line[1], split_line[2]
        data_items = data_part.split()

        ctp, gap_5prime, ctm, gap_3prime, other_5prime, other_3prime = calculate_damage_for_plot(data_items)

        values = {'Other': [0] * 30, 'CT': [0] * 30, 'GA': [0] * 30}
        for i in range(1, 16):
            # 5' 
            values['CT'][positions.index(i)] = ctp[i]
            values['GA'][positions.index(i)] = gap_5prime[i]
            values['Other'][positions.index(i)] = other_5prime[i]

            # 3' 
            values['CT'][positions.index(-i)] = ctm[i]
            values['GA'][positions.index(-i)] = gap_3prime[i]
            values['Other'][positions.index(-i)] = other_3prime[i]

        for key in values:
            values_all_files[key].append(values[key])

    if not any(values_all_files['Other']):
        print("Warning: No valid rows found for the given tax ID.")
        return

    if ymax == 0 or ymax == "0":
        max_y = min(1.0, max(max(item for sublist in values_all_files[key] for item in sublist) for key in ['Other', 'CT', 'GA']) * 1.2)
    else:
        max_y = float(ymax)

    # do the plotting

    plt.figure(figsize=(10, 5))  
    color_palette = {'Other': '#009E73', 'CT': '#F8766D', 'GA': '#56B4E9'}  # colorblind-friendly palette

    # 5' 
    ax1 = plt.subplot(1, 2, 1)
    handles1 = []  # store the line handles for ax1
    labels1 = ['Other', 'C to T', 'G to A']  # labels for legend, one entry for each damage type
    for key in values_all_files:
        color = color_palette[key] 
        for j, file_values in enumerate(values_all_files[key]):
            line, = ax1.plot([str(pos) for pos in positions if pos > 0],
                             [file_values[i] for i, pos in enumerate(positions) if pos > 0],
                             label=f'{key} ({subs_files[j]})', color=color, linewidth=2)
            if j == 0: 
                handles1.append(line)
    ax1.set_ylim(0, max_y)
    ax1.set_xlabel('Position', fontsize=14)
    ax1.set_ylabel('Frequency', fontsize=14)
    ax1.yaxis.set_ticks_position('left')
    ax1.tick_params(right=False, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)

    # 3' 
    ax2 = plt.subplot(1, 2, 2)
    handles2 = []  
    for key in values_all_files:
        color = color_palette[key] 
        for j, file_values in enumerate(values_all_files[key]):
            line, = ax2.plot([str(pos) for pos in positions if pos < 0],
                             [file_values[i] for i, pos in enumerate(positions) if pos < 0],
                             label=f'{key} ({subs_files[j]})', color=color, linewidth=2)
            if j == 0:  
                handles2.append(line)
    ax2.set_ylim(0, max_y)
    ax2.set_xlabel('Position', fontsize=14)
    ax2.yaxis.set_ticks_position('right')
    ax2.tick_params(left=False, labelsize=12)

    ax2.legend(handles=handles2, labels=labels1, loc='upper left', fontsize=12) # legend

    plt.suptitle(f'Damage Plot for {tax_name} (tax ID {tax_id})', fontsize=18)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    file_extension = os.path.splitext(plotfile)[1].lower()
    if file_extension == '.pdf':
        plt.savefig(plotfile, format='pdf')
    else:
        if file_extension != '.png':
            print("Warning: Invalid plot file suffix. Your plot file is being saved in png format with the filename you requested.")
        plt.savefig(plotfile)

    plt.close()


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



def parse_exclude_keywords(args):
    if hasattr(args, 'exclude_keywords') or hasattr(args, 'exclude_keyword_file'):
        if args.exclude_keywords and args.exclude_keyword_file:
            raise ValueError("Please only provide one of --exclude_keywords or --exclude_keyword_file, not both.")
        
        exclude_keywords = args.exclude_keywords if args.exclude_keywords else []
        if not isinstance(exclude_keywords, list):
            exclude_keywords = [exclude_keywords]
        
        if args.exclude_keyword_file:
            if not os.path.exists(args.exclude_keyword_file):
                raise FileNotFoundError(f"exclude_keyword_file path does not exist: {args.exclude_keyword_file}")
            with open(args.exclude_keyword_file, 'r') as f:
                exclude_keywords.extend([line.strip() for line in f if line.strip()])
        
        exclude_keywords = [kw.lstrip("'").lstrip('"').rstrip("'").rstrip('"') for kw in exclude_keywords] # strip quotes 
        
        # formatted_keywords = []
        # good to surround the digit-only tax ids with a :, so we don't accidentally hit substring tax ids
        #for keyword in exclude_keywords:
        #    if keyword.isdigit():
        #        keyword = f"{keyword}:"
        #    formatted_keywords.append(keyword)

        return exclude_keywords
    else:
        return []


def tsvs_to_matrix(parsed_data, output_file, include='all', minreads=50):
    # for combine. pretty straightforward, read them all in and construct the matrix by pulling values out as needed
    # "meandamage" output is actually weighted by number of reads

    include_damage = 'damage' in include or 'all' in include
    include_duplicity = 'duplicity' in include or 'all' in include
    include_dust = 'dust' in include or 'all' in include
    include_taxpath = 'taxpath' in include or 'all' in include
    include_gc = "gc" in include or "all" in include

    tax_data = {}
    for sample_name, records in parsed_data.items():  # parsed_data is a dict: {sample_name: records}
        if sample_name.endswith('.tsv'):
            sample_name = sample_name.replace('.tsv', '')

        print(f"Processing sample: {sample_name} with {len(records)} records.")
        for record in records:
            taxpath = record[-1]
            tax = taxpath.split(';')[0].strip('"')
            if tax not in tax_data:
                tax_data[tax] = {
                    "TaxPath": taxpath,
                    "samples": {},
                    "TotalReads": 0,
                    "WeightedDamage": 0 if include_damage else None,
                    "WeightedDust": 0 if include_dust else None,
                    "WeightedDup": 0 if include_duplicity else None,
                    "WeightedReadGC": 0 if include_gc else None,
                    "WeightedRefGC": 0 if include_gc else None,
                    "WeightSum": 0,
                }

            sample_data = {
                "reads": int(record[2]),
                "damage": float(record[5]) if include_damage else None,
                "duplicity": float(record[3]) if include_duplicity else None,
                "dust": float(record[4]) if include_dust else None,
                "avg_read_gc": float(record[9]) if include_gc else None,
                "avg_ref_gc": float(record[10]) if include_gc else None,
            }
            tax_data[tax]['samples'][sample_name] = sample_data
            tax_data[tax]['TotalReads'] += sample_data['reads']
            tax_data[tax]["WeightSum"] += sample_data["reads"]
            if include_damage and sample_data['damage'] is not None:
                tax_data[tax]['WeightedDamage'] += sample_data['damage'] * sample_data['reads']
            if include_dust and sample_data["dust"] is not None:
                tax_data[tax]["WeightedDust"] += (
                    sample_data["dust"] * sample_data["reads"]
                )
            if include_duplicity and sample_data["duplicity"] is not None:
                tax_data[tax]["WeightedDup"] += (
                    sample_data["duplicity"] * sample_data["reads"]
                )
            if include_gc and sample_data["avg_read_gc"] is not None:
                tax_data[tax]["WeightedReadGC"] += (
                    sample_data["avg_read_gc"] * sample_data["reads"]
                )
            if include_gc and sample_data["avg_ref_gc"] is not None:
                tax_data[tax]["WeightedRefGC"] += (
                    sample_data["avg_ref_gc"] * sample_data["reads"]
                )

    for tax, data in tax_data.items():
        if include_damage and data['WeightSum'] > 0:
            data['MeanDamage'] = data['WeightedDamage'] / data['WeightSum']
        else:
            data['MeanDamage'] = 'NA'
        if include_dust and data["WeightSum"] > 0:
            data["MeanDust"] = data["WeightedDust"] / data["WeightSum"]
        else:
            data["MeanDust"] = "NA"
        if include_duplicity and data["WeightSum"] > 0:
            data["MeanDup"] = data["WeightedDup"] / data["WeightSum"]
        else:
            data["MeanDup"] = "NA"
        if include_gc and data["WeightSum"] > 0:
            data["MeanReadGC"] = data["WeightedReadGC"] / data["WeightSum"]
            data["MeanRefGC"] = data["WeightedRefGC"] / data["WeightSum"]
        else:
            data["MeanReadGC"] = "NA"
            data["MeanRefGC"] = "NA"

    tax_data = {tax: data for tax, data in tax_data.items() if data['TotalReads'] >= minreads}
    sorted_tax_data = sorted(tax_data.items(), key=lambda x: x[1]['TotalReads'], reverse=True)

    with open(output_file, "w") as outfile:
        header = ['Tax', 'TotalReads']
        if include_damage:
            header.append('MeanDamage')
        if include_duplicity:
            header.append("MeanDup")
        if include_dust:
            header.append("MeanDust")
        if include_gc:
            header.append("MeanReadGC")
            header.append("MeanRefGC")
        for sample_name in sorted(parsed_data.keys()):
            if sample_name.endswith('.tsv'):
                sample_name = sample_name.replace(".tsv", "")
            header.append(f"{sample_name}_reads")
            if include_damage:
                header.append(f"{sample_name}_damage")
            if include_duplicity:
                header.append(f"{sample_name}_duplicity")
            if include_dust:
                header.append(f"{sample_name}_dust")
            if include_gc:
                header.append(f"{sample_name}_avgreadgc")
                header.append(f"{sample_name}_avgrefgc")
        if include_taxpath:
            header.append("TaxPath")
        outfile.write("\t".join(header) + "\n")

        for tax, data in sorted_tax_data:
            row = [tax, str(data['TotalReads'])]
            if include_damage:
                row.append(str(round(data["MeanDamage"], 3)))
            if include_duplicity:
                row.append(
                    str(round(data["MeanDup"], 3))
                    if isinstance(data["MeanDup"], float)
                    else "NA"
                )
            if include_dust:
                row.append(
                    str(round(data["MeanDust"], 3))
                    if isinstance(data["MeanDust"], float)
                    else "NA"
                )
            if include_gc:
                row.append(
                    str(round(data["MeanReadGC"], 3))
                    if isinstance(data["MeanReadGC"], float)
                    else "NA"
                )
                row.append(
                    str(round(data["MeanRefGC"], 3))
                    if isinstance(data["MeanRefGC"], float)
                    else "NA"
                )

            for sample_name in sorted(parsed_data.keys()):
                if sample_name.endswith('.tsv'):
                    sample_name = sample_name.replace(".tsv", "")
                sample_data = data['samples'].get(sample_name, {})
                row.append(str(sample_data.get("reads", 0)))
                if include_damage:
                    row.append(str(sample_data.get('damage')) if sample_data.get('damage') is not None else 'NA')
                if include_duplicity:
                    row.append(str(sample_data.get('duplicity')) if sample_data.get('duplicity') is not None else 'NA')
                if include_dust:
                    row.append(str(sample_data.get('dust')) if sample_data.get('dust') is not None else 'NA')
                if include_gc:
                    row.append(
                        str(sample_data.get("avg_read_gc"))
                        if sample_data.get("avg_read_gc") is not None
                        else "NA"
                    )
                    row.append(
                        str(sample_data.get("avg_ref_gc"))
                        if sample_data.get("avg_ref_gc") is not None
                        else "NA"
                    )
            if include_taxpath:
                row.append(data["TaxPath"])
            outfile.write('\t'.join(row) + '\n')


def make_krona_xml(in_tsv, in_tsv_files, out_xml, minreads, maxdamage):
    # create an xml text file to be loaded into kronatools for visualization. adds damage as a colour option.

    level_order = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"] # detect "upto"

    # did we get one or more files? check they exist then parse the input style 
    input_files = []
    if in_tsv:
        input_files = in_tsv
    elif in_tsv_files:
        with open(in_tsv_files, 'r') as file:
            input_files = [line.strip() for line in file if line.strip()]   

    # validate input files: check existence, header format, and minreads threshold
    valid_files = []
    for file in input_files:
        if not os.path.exists(file):
            print(f"Error: File {file} does not exist.")
            continue
        with open(file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                print(f"Error: File {file} does not have enough lines to check reads.")
                continue
            header = lines[0].strip().split('\t')
            if header[0] != "TaxNodeID":
                print("Error: It looks like your input tsv file was not generated from bamdam compute or bamdam combine, or perhaps you edited the header.")
                continue
            fields = lines[1].strip().split('\t')
            try:
                reads = int(fields[2])
                if reads < minreads:
                    print(f"Warning: File {file} only contains taxa with at most ({reads}), which is below minreads ({minreads}). This file will be excluded.")
                    continue
            except (ValueError, IndexError):
                print(f"Error: File {file} has an unexpected format and will be excluded.")
                continue
        valid_files.append(file)
    input_files = valid_files

    # If no valid files remain, exit early
    if not input_files:
        print("Error: No valid input files remain after filtering. Exiting.")
        exit(-1)

    # checks complete

    # put everything into a tree structure in memory (these tsv files are usually pretty small) and then afterwards interpret it into krona xml format
    tree = {}
    sample_names = []
    sample_max_reads = {} # for colour scale
    sample_max_damage = {}  # for colour scale
    sample_reads_at_root = {} # krona format needs this

    for file in input_files:
        sample_name = file.split('/')[-1].replace('.tsv', '')
        sample_names.append(sample_name)

        with open(file, 'r') as file:
            lines = file.readlines()
        header = lines[0].strip().split('\t')

        # the first thing to do is get the top level; the "upto" that was used
        taxpath_col = header.index("TaxPath") if "TaxPath" in header else -1

        if not 'toplevel' in locals():
            # toplevel should be the same across tsv files!
            levels_found = set()
            for line in lines[1:]:  
                fields = line.strip().split('\t')
                if int(fields[2]) < minreads:
                    break
                taxpath = fields[taxpath_col]
                first_node = taxpath.split(';')[0]
                level = first_node.split(':')[2].strip('"')
                levels_found.add(level)
            toplevel = None
            for level in level_order:
                if level in levels_found:
                    toplevel = level
                    break

        max_reads = int(lines[1].strip().split('\t')[2])  # tsvs are ordered; top line is max reads
        sample_max_reads[sample_name] = max_reads

        for line in lines[1:]:  # skip the header
            fields = line.strip().split('\t')
            reads = float(fields[2]);
            if reads < minreads:
                break # you're done with this file
            dup = float(fields[3]); dust = float(fields[4]); damage = float(fields[5]); lengths = float(fields[7])
            taxpath = fields[-1].strip('"').strip("'") # 16 used to be -1
            taxpathsplit = taxpath.split(";")
            fullnode = taxpathsplit[0].split(":")
            taxid = fullnode[0]
            taxname = fullnode[1]
            level = fullnode[2]
            if toplevel not in taxpath:
                continue # skip this line, i don't need to deal with this rare weird case.
            if reads >= minreads:
                if taxid not in tree: # add a new node to the tree structure if it's not already there
                    tree[taxid] = {
                        'taxname': taxname,
                        'samples': {sample_name: {'reads': reads, 'duplicity': dup, 'dust': dust, 'damage': damage, 'length': lengths}},
                        'taxpath': taxpath,
                        'children': set()
                    }
                else: # otherwise just add this sample info in
                    if sample_name not in tree[taxid]['samples']:
                        tree[taxid]['samples'][sample_name] = {'reads': reads, 'duplicity': dup, 'dust': dust, 'damage': damage, 'length': lengths}
                    # else:  # this should not ever happen anymore
                        # print(f"Warning: The file {sample_name} has two lines for the same tax id {taxid}, or this file is in your list more than once.")
                if level != toplevel:
                        # attach it to its parent ; make the parent to attach it to if needed
                        parent = taxpathsplit[1].split(":")
                        parent_taxid = parent[0]; parent_taxname = parent[1]; parent_level = parent[2]
                        if parent_taxid not in tree:
                            # initialize it. we'll fill it in later. it must be there. 
                            tree[parent_taxid] = {'taxname': parent_taxname, 'samples': {}, 'taxpath': ", ".join(taxpathsplit[1:]), 'children': set()}
                        tree[parent_taxid]['children'].add(taxid)


        # sum up the reads at the nodes without children in all the nodes which have info for this sample
        root_reads = sum(tree[taxid]['samples'][sample_name]['reads'] for taxid in tree if 
                        sample_name in tree[taxid]['samples'] and not any(taxid in tree[child_id]['children'] for child_id in tree))
        sample_reads_at_root[sample_name] = root_reads
        
        # also get max damage of all the taxa for this sample  
        sample_max_damage[sample_name] = max(
            node['samples'][sample_name]['damage']
            for node in tree.values()
            if sample_name in node['samples'] and node['samples'][sample_name]['damage'] is not None
        )

    # a very handy debug point, i will leave it here in case anyone is ever changing this code and needs it (needs pprint package):
#    from pprint import pprint
#    treefile = "tree.tmp"
#    def write_tree_to_file(tree, filename):
#        with open(filename, "w") as f:
#            pprint(tree, stream=f)  # print the tree in a readable format
#        print(f"Tree structure written to file: {filename}")
#    write_tree_to_file(tree, treefile)

    # now build the xml manually by iterating through the tree 

    # first build a summary of everything if there is more than one thing; 
    # just iterate through existing nodes to do this then treat the summary like a bonus sample 
    if len(input_files) > 1:
        all_children = {child for node in tree.values() for child in node['children']}
        tree_items = list(tree.items())

        for taxid, node in tree_items:
            all_reads = 0
            weighted_damage_sum = 0
            total_reads_for_damage = 0
            duplicity_values = []
            dust_values = []
            length_values = []

            for sample in node['samples']:
                sample_data = node['samples'][sample]
                all_reads += sample_data['reads']
                duplicity_values.append(sample_data['duplicity'])
                dust_values.append(sample_data['dust'])
                length_values.append(sample_data['length'])
                if sample_data['damage'] is not None:
                    weighted_damage_sum += sample_data['reads'] * sample_data['damage']
                    total_reads_for_damage += sample_data['reads']

            avg_duplicity = sum(duplicity_values) / len(duplicity_values) if duplicity_values else 0
            avg_dust = sum(dust_values) / len(dust_values) if dust_values else 0
            damage = weighted_damage_sum / total_reads_for_damage if total_reads_for_damage > 0 else 0
            avg_length = sum(length_values) / len(length_values) if length_values else 0 
 
            node['samples']['Summary'] = {
                'reads': all_reads,
                'damage': round(damage,3),
                'duplicity': round(avg_duplicity,3),
                'dust': round(avg_dust,3),
                'length': round(avg_length,3)
            }

        summary_root_reads = sum(
            node['samples']['Summary']['reads']
            for taxid, node in tree_items
            if 'Summary' in node['samples'] and taxid not in all_children
        )
        sample_reads_at_root['Summary'] = summary_root_reads

        sample_names.insert(0, "Summary")

    # right now i am standardizing the colour scale for both damage and number of reads across all krona plots
    # but this is set up so that can be changed easily enough
    overall_max_reads = max(
        node['samples'][sample]['reads']
        for node in tree.values()
        for sample in node['samples']
    )
    overall_max_damage = max(
        node['samples'][sample]['damage']
        for node in tree.values()
        for sample in node['samples']
        if node['samples'][sample]['damage'] is not None
    )

    # allow user to specify their own max value for the damage colour scale
    if maxdamage is not None:
        overall_max_damage = maxdamage  # use user-specified value
    else:
        overall_max_damage = max(overall_max_damage, 0.3)  # otherwise auto-detect but ensure at least 0.3

    def add_node_to_xml(node_id, tree, indent_level=2):

        node_data = tree[node_id]
        indent = "\t" * indent_level
        node_str = f'{indent}<node name="{node_data["taxname"]}">\n'
        
        for attr in ["reads", "damage", "duplicity", "dust", "length"]:
            node_str += f"{indent}\t<{attr}>"
            for sample_name in sample_names:
                # add a 0 if this sample name doesn't appear for this node
                value = node_data['samples'].get(sample_name, {}).get(attr, 0)
                node_str += f"<val>{value}</val>"
            node_str += f"</{attr}>\n"

        node_str += f"{indent}\t<taxid><val>{node_id}</val></taxid>\n" # add the tax id too 

        for child_id in node_data["children"]:
            node_str += add_node_to_xml(child_id, tree, indent_level + 1)
        
        node_str += f"{indent}</node>\n"
        return node_str

    with open(out_xml, "w") as file:
        file.write('<krona>\n')
        file.write("\t<attributes magnitude=\"reads\">\n")
        file.write("\t\t<attribute display=\"Reads\">reads</attribute>\n")
        file.write("\t\t<attribute display=\"Damage\">damage</attribute>\n")
        file.write("\t\t<attribute display=\"Duplicity\">duplicity</attribute>\n")
        file.write("\t\t<attribute display=\"Dust\">dust</attribute>\n")
        file.write("\t\t<attribute display=\"Mean Read Length\">length</attribute>\n")
        file.write("\t\t<attribute display=\"Tax ID\">taxid</attribute>\n")
        file.write("\t</attributes>\n")
        file.write(f'\t<color attribute="reads" valueStart="0" valueEnd="{overall_max_reads}" hueStart="120" hueEnd="240"></color>\n')
        file.write(f'\t<color attribute="damage" valueStart="0" valueEnd="{overall_max_damage}" hueStart="0" hueEnd="240"></color>\n')
        file.write("\t<datasets>\n")
        for sample_name in sample_names:
            file.write(f"\t\t<dataset>{sample_name}</dataset>\n")
        file.write("\t</datasets>\n")
        file.write('\t<node name="Root">\n')
        file.write('\t\t<reads>')
        for sample_name in sample_names:
            root_reads = sample_reads_at_root.get(sample_name, 0)  
            file.write(f'<val>{root_reads}</val>')
        file.write('</reads>\n')

        for root_taxid in [taxid for taxid in tree if not any(taxid in tree[child_id]["children"] for child_id in tree)]:
            try:
                file.write(add_node_to_xml(root_taxid, tree))
            except RecursionError:
                print(f"There was a recursion error for {root_taxid}. We will skip it. If your tree is truly ginormous, you can increase the default recursion limit.")
                continue

        file.write("</krona>\n")
    print(f"Krona XML written to {out_xml}")


def shrink(args):
    formatted_exclude_keywords = parse_exclude_keywords(args)
    lca_file_type = find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("You are running bamdam shrink with a metaDMG-style lca file. This is ok, but be aware the output lca file will be in ngsLCA lca file format, as all the other functions in bamdam require this format.")
    shortlcalines = write_shortened_lca(args.in_lca, args.out_lca, args.upto, args.mincount, formatted_exclude_keywords, lca_file_type)
    write_shortened_bam(args.in_bam, args.out_lca, args.out_bam, args.stranded, args.minsim, args.annotate_pmd, shortlcalines)

def compute(args):
    lca_file_type = find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("Error: It looks like you're trying to run bamdam compute with a metaDMG-style lca file. Please use an ngsLCA-style lca file.")
        sys.exit()
    nodedata, pmds_in_bam = gather_subs_and_kmers(args.in_bam, args.in_lca, kn=args.k, upto=args.upto, stranded=args.stranded)
    parse_and_write_node_data(nodedata, args.out_tsv, args.out_subs, args.stranded, pmds_in_bam)  

def extract(args):
    lca_file_type = find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("Error: It looks like you're trying to run bamdam extract with a metaDMG-style lca file. Please use an ngsLCA-style lca file.")
        sys.exit()
    extract_reads(args.in_lca, args.in_bam, args.out_bam, args.keyword, args.subset_header, args.only_top_ref)

def plotdamage(args):
    make_damage_plot(args.in_subs_list, args.in_subs, args.tax, args.outplot, args.ymax)

def plotbaminfo(args):
    make_baminfo_plot(args.in_bam, args.in_bam_list, args.outplot)

def combine(args):
    # parse real quick
    input_files = []
    if args.in_tsv:
        input_files = args.in_tsv
    elif args.in_tsv_list:
        with open(args.in_tsv_list, 'r') as file:
            input_files = [line.strip() for line in file if line.strip()]   
    parsed_data = {}
    for file_path in input_files:
        sample_name = file_path.split('/')[-1].replace('.tsv', '')  
        parsed_data[sample_name] = []
        with open(file_path, 'r') as file:
            for i, line in enumerate(file):
                fields = line.strip().split('\t')
                if i == 0: # skip the header line
                    continue
                parsed_data[sample_name].append(fields)
    tsvs_to_matrix(parsed_data, args.out_tsv, args.include, args.minreads)

def krona(args):
    make_krona_xml(args.in_tsv, args.in_tsv_list, args.out_xml, args.minreads, args.maxdamage)

def main():

    # Initialize
    parser = argparse.ArgumentParser(
        description="Bamdam processes ancient metagenomic bam and lca files. Type bamdam [command] -h for more detailed help regarding a specific command.")
    
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Shrink
    parser_shrink = subparsers.add_parser('shrink', help="Filter the BAM and LCA files.")
    parser_shrink.add_argument("--in_lca", type=str, required=True, help="Path to the input LCA file (required)")
    parser_shrink.add_argument("--in_bam", type=str, required=True, help="Path to the input (read-sorted) BAM file (required)")
    parser_shrink.add_argument("--out_lca", type=str, required=True, help="Path to the short output LCA file (required)")
    parser_shrink.add_argument("--out_bam", type=str, required=True, help="Path to the short output BAM file (required)")
    parser_shrink.add_argument("--stranded", type=str, required=True, help="Either ss for single stranded or ds for double stranded (required)")
    parser_shrink.add_argument("--mincount", type=int, default=5, help="Minimum read count to keep a node (default: 5)")
    parser_shrink.add_argument("--upto", type=str, default="family", help="Keep nodes up to and including this tax threshold (default: family)")
    parser_shrink.add_argument("--minsim", type=float, default=0.9, help="Minimum similarity to reference to keep an alignment (default: 0.9)")
    parser_shrink.add_argument("--exclude_keywords", type=str, nargs='+', default=[], help="Keyword(s) to exclude when filtering (default: none)")
    parser_shrink.add_argument("--exclude_keyword_file", type=str, default=None, help="File of keywords to exclude when filtering, one per line (default: none)")
    parser_shrink.add_argument("--annotate_pmd", action='store_true', help="Annotate output bam file with PMD tags  (default: not set)")
    parser_shrink.set_defaults(func=shrink) 

    # Compute
    parser_compute = subparsers.add_parser('compute', help="Compute tsv and subs files.")
    parser_compute.add_argument("--in_bam", type=str, required=True, help="Path to the BAM file (required)")
    parser_compute.add_argument("--in_lca", type=str, required=True, help="Path to the LCA file (required)")
    parser_compute.add_argument("--out_tsv", type=str, required=True, help="Path to the output tsv file (required)")
    parser_compute.add_argument("--out_subs", type=str, required=True, help="Path to the output subs file (required)")
    parser_compute.add_argument("--stranded", type=str, required=True, help="Either ss for single stranded or ds for double stranded (required)")
    parser_compute.add_argument("--k", type=int, default=29, help="Value of k for per-node counts of unique k-mers and duplicity (default: 29)")
    parser_compute.add_argument("--upto", type=str, default="family", help="Keep nodes up to and including this tax threshold; use root to disable (default: family)")
    parser_compute.set_defaults(func=compute)

    # Extract
    parser_extract = subparsers.add_parser('extract', help="Extract alignments of reads containing a keyword in an associated lca file.")
    parser_extract.add_argument("--in_bam", type=str, required=True, help="Path to the BAM file (required)")
    parser_extract.add_argument("--in_lca", type=str, required=True, help="Path to the LCA file (required)")
    parser_extract.add_argument("--out_bam", type=str, required=True, help="Path to the filtered BAM file (required)")
    parser_extract.add_argument("--keyword", type=str, required=True, help="Keyword or phrase to filter for, e.g. a taxonomic node ID (required)")
    parser_extract.add_argument("--subset_header", action='store_true', help="Subset the header to only relevant references (default: not set)")
    parser_extract.add_argument("--only_top_ref", action='store_true', help="Only keep alignments to the most-hit reference (default: not set)")
    parser_extract.set_defaults(func=extract)

    # Plot damage
    parser_plotdamage = subparsers.add_parser('plotdamage', help="Produces a postmortem damage plot for a specified taxonomic node using the subs file.")
    group_input_plotdamage = parser_plotdamage.add_mutually_exclusive_group(required=True)
    group_input_plotdamage.add_argument("--in_subs", nargs='+', help="Input subs file(s)")
    group_input_plotdamage.add_argument("--in_subs_list", help="Path to a text file contaning input subs files, one per line")
    parser_plotdamage.add_argument("--tax", type=str, required=True, help="Taxonomic node ID (required)")
    parser_plotdamage.add_argument("--outplot", type=str, default="damage_plot.png", help="Filename for the output plot, ending in .png or .pdf (default: damage_plot.png)")
    parser_plotdamage.add_argument("--ymax", type=str, default="0", help="Maximum for y axis (optional)")
    parser_plotdamage.set_defaults(func=plotdamage)

    # Plot bam info
    parser_plotbaminfo = subparsers.add_parser('plotbaminfo', help="Produces a mismatch and read length distribution plot for an input bam.")
    group_input_plotbaminfo = parser_plotbaminfo.add_mutually_exclusive_group(required=True)
    group_input_plotbaminfo.add_argument("--in_bam", nargs='+', help="Input bam file(s)")
    group_input_plotbaminfo.add_argument("--in_bam_list", help="Path to a text file containing input bams, one per line")
    parser_plotbaminfo.add_argument("--outplot", type=str, default="baminfo_plot.png", help="Filename for the output plot, ending in .png or .pdf (default: baminfo_plot.png)")
    parser_plotbaminfo.set_defaults(func=plotbaminfo)

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
        sortorder = get_sorting_order(args.in_bam)
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
        if hasattr(args, 'exclude_keyword_file') and args.exclude_keyword_file:
            print(f"exclude_keywords: loaded from {args.exclude_keyword_file}")
        if hasattr(args, 'exclude_keywords') and args.exclude_keywords:
            print(f"exclude_keywords: {args.exclude_keywords}")
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

    elif args.command == 'extract':
        print("Hello! You are running bamdam extract with the following arguments:")
        print(f"in_bam: {args.in_bam}")
        print(f"in_lca: {args.in_lca}")
        print(f"out_bam: {args.out_bam}")
        print(f"keyword: {args.keyword}")

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

