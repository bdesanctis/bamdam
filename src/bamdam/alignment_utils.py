
# alignment parsing utils

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

def get_mismatches(seq, cigar, md):  
    # parses a read, cigar and md string to determine mismatches and positions. 
    # does not output info on insertions/deletions, but accounts for them.
    # thanks jonas oppenheimer who wrote half of this function :)

    # goes in two passes: once w/ cigar and once w/ md 

    # the strategy is to reconstruct the alignment through the cigar string first, so read_seq and ref_seq are the same length
    # but might have "-"s and "N"s floating around, and next to inform the substitutions from the md string

    cig = utils.parse_cigar_cached(cigar)
    md_list = utils.parse_md_cached(md)

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

    COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

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
    backwards = utils.is_reverse_strand(flagsum)

    if backwards:  # flip all the positions and reverse complement all the nucleotides. the read in the bam is reverse-complemented if aligned to the reverse strand.
        mmsc = []
        matchsc = []

        for entry in range(0,len(mms)): # mismatches
            new_entry = [
                COMPLEMENT.get(mms[entry][0], "N"),
                COMPLEMENT.get(mms[entry][1], "N"),
                readlength - mms[entry][2] + 1,
            ]
            mmsc.append(new_entry)

        for entry in range(0, len(matchs)): # matches
            new_entry = [
                COMPLEMENT.get(matchs[entry][0], "N"),
                COMPLEMENT.get(matchs[entry][1], "N"),
                readlength - matchs[entry][2] + 1,
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
        refseqlist = utils.rev_complement(refseqlist)
        readseqlist = utils.rev_complement(readseqlist)
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
