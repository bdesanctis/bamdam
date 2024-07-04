# BamDam by Bianca De Sanctis, bddesanctis@gmail.com
# Last updated July 3 2024

import sys 
import re
import csv
import pysam
import math
import argparse
import os
try:
    from tqdm import tqdm
except ImportError:
    pass  # tqdm is not available . it's not necessary, just nice for progress bars 


def write_shortened_lca(original_lca_path,short_lca_path,upto,mincount,lcaheaderlines,other_lca_keywords):

    print("\nWriting a filtered LCA file...")

    numheaderlines = lcaheaderlines # in the lca file to ignore 

    # pass 1: make a dictionary with all the tax ids and their counts 
    number_counts = {}
    with open(original_lca_path, 'r') as file:
        for _ in range(numheaderlines):
            next(file) 
        for line in file:
            if upto in line: 
                if not other_lca_keywords or all(keyword in line for keyword in other_lca_keywords):                        
                    entry = line.strip().split('\t')
                    if len(entry) > 1:  
                        fields = entry[1:] # can ditch the read id 
                        keepgoing = True; field = 0
                        while keepgoing:
                            number = fields[field].split(':')[0]
                            if number in number_counts:
                                number_counts[number] += 1
                            else:
                                number_counts[number] = 1
                            if upto in fields[field]:
                                keepgoing = False
                            field +=1 

    goodnodes = [key for key, count in number_counts.items() if count >= mincount]
    # these are the nodes that have at least the min count of reads assigned to them (or below them), and which are at most upto

    # pass 2: rewrite lines into a new lca file that pass the filter
    oldreadname = ""
    with open(original_lca_path, 'r') as infile, open(short_lca_path, 'w') as outfile:
            for _ in range(numheaderlines):
                next(infile) # assumes it has the usual 2 comment lines at the top from ngslca
            for line in infile:
                entry = line.strip().split('\t')
                readnamesplit = entry[0].split(':')[0:7]
                newreadname = ":".join(readnamesplit)
                if newreadname == oldreadname:
                    print("Alert! You have duplicate entries in your LCA file, for example " + newreadname + ". That's a problem. You should fix this (e.g. using uniq) and re-run BamDam.")
                if upto in line:
                    # you can just go straight to the upto level and check if that node has high enough count 
                    if not other_lca_keywords or all(keyword in line for keyword in other_lca_keywords):                        
                        for field in entry[1:]:
                            if f":{upto}" in field:
                                number = field.split(':')[0]
                                if number in goodnodes:
                                    outfile.write(line)
                                    break

    print("Wrote shortened lca file. \n")

def write_shortened_bam(original_bam_path,short_lca_path,short_bam_path,stranded,lcaheaderlines,minsimilarity,other_lca_keywords): 
    # runs through the existing bam and the new short lca file at once, and writes only lines to the new bam which are represented in the short lca file
    # also annotates with pmd scores as it goes

    # now takes in minsimilarity as a percentage, and will keep reads w/ equal to or greater than NM flag to this percentage 

    print("Writing a filtered bam file, annotated with PMD scores... ")

    # Count the number of reads in the BAM file for the progress bar
    with pysam.AlignmentFile(original_bam_path, "rb", check_sq=False, require_index=False) as infile, \
         pysam.AlignmentFile(short_bam_path, "wb", header=infile.header) as outfile, \
         open(short_lca_path, 'r') as shortlcafile:

        # Skip the LCA header lines
        for _ in range(lcaheaderlines):
            lcaline = next(shortlcafile)

        # Read the first LCA line and extract the read name
        lcaline = next(shortlcafile)
        lcareadname = ":".join(lcaline.strip().split('\t')[0].split(':')[0:7])

        currentlymatching = False
        notdone = True

        try:
            bamread = next(infile)
        except StopIteration:
            notdone = False

        # pbar_bam = tqdm(total=bam_line_count, desc='Processing BAM file')

        while notdone:
            if bamread.query_name == lcareadname:
                # Copy this line and all the rest until you hit a nonmatching LCA line
                similarity = 1 - bamread.get_tag('NM') / bamread.query_length
                if similarity >= minsimilarity:
                    pmd = get_pmd(bamread, stranded)
                    bamread.tags += [('DS','%.3f' % pmd)]
                    outfile.write(bamread)
                currentlymatching = True
                while currentlymatching:
                    try:
                        bamread = next(infile)
                        # pbar_bam.update(1)
                        if bamread.query_name == lcareadname:
                            similarity = 1 - bamread.get_tag('NM') / bamread.query_length
                            if similarity >= minsimilarity:
                                # write the read! 
                                pmd = get_pmd(bamread, stranded)
                                bamread.tags += [('DS','%.3f' % pmd)]
                                outfile.write(bamread)
                        else:
                            currentlymatching = False
                    except StopIteration:
                        notdone = False
                        break
                try:
                    lcaline = next(shortlcafile)
                    lcareadname = ":".join(lcaline.strip().split('\t')[0].split(':')[0:7])
                except StopIteration:
                    notdone = False
            else:
                try:
                    bamread = next(infile)
                    # pbar_bam.update(1)
                except StopIteration:
                    notdone = False

    # pbar_bam.close()
    print("Wrote a shortened bam file. \n")

def get_mismatches(seq, cigar, md):  
    # parses a read, cigar and md string to determine nucleotide mismatches and positions. 

    # does not output info on insertions/deletions but accounts for them.
    # thanks jonas oppenheimer who wrote half of this function :)

    # goes in two passes: once w/ cigar and once w/ md 
    '''
    reconstructs a nucleotide mismatch + position table given the query, cigar and md strings
    '''

    cig = re.findall('\d+\D', cigar)
    md_pattern = re.compile(r'\d+|\^[A-Za-z]+|[A-Za-z]')
    md_list = md_pattern.findall(md)

    ref_seq = ''
    read_seq = ''
    query_pos = 0 # indexes the ref reconstruction
    read_pos = 0 # indexes the read reconstruction (which is the input read but with potential added "-"s if the ref has insertions)
    for x in cig:
        cat = x[-1]
        if cat in ['H', 'P']: # doesn't consume reference or query
            continue

        bases = int(x[:-1])

        if cat == 'S': # soft clip
            read_seq += seq[read_pos:read_pos + bases] # include the bases in the reconstructed read 
            continue

        elif cat in ['M', '=', 'X']: # match
            ref_seq += seq[query_pos:query_pos + bases]
            read_seq += seq[read_pos:read_pos + bases]
            query_pos += bases
            read_pos += bases

        elif cat in ['D', 'N']: # 'D' : the reference has something there and the read doesn't, but pad them both 
            ref_seq += 'N' * bases
            read_seq += '-' * bases

        elif cat == 'I': # I: the read has something there and the reference doesn't, but pad them both 
            read_seq += seq[read_pos:read_pos + bases] # 'N' * bases
            ref_seq += '-' * bases
            query_pos += bases
            read_pos += bases

        else:
            sys.exit("you've got some strange cigar strings here")
    
    ref_pos = 0
    read_pos = 0
    mismatch_list = [] # of format [ref, read, pos in alignment]
    for x in md_list:
        if x.startswith('^'): # remember this can be arbitrarily long (a ^ and then characters)
            num_chars_after_hat = len(x) -1 
            for i in range(0,num_chars_after_hat):
                currchar = x[i+1]
                refhere = currchar
                readhere = "-" 
                ref_pos += 1 # but not the read pos
                # don't need to append to mismatch list 
        else: 
            if x.isdigit(): # can be multiple digits 
                # you're a number or a letter
                # if you're a number, you're matching. skip ahead. don't need to add to the mismatch list. 
                ref_pos += int(x)
                read_pos += int(x)
            else: # it will only be one character at a time 
                refhere = x
                readhere = read_seq[ref_pos]
                # update ref_seq
                char_list = list(ref_seq)
                char_list[ref_pos] = x
                ref_seq = "".join(char_list)
                # moving on
                # ALERT, i just changed the below to read_pos instead of ref_pos! i hope this fixes my problem
                mismatch_list.append([refhere,readhere,read_pos +1]) # in genetics we are 1-based (a->g at position 1 means the actual first position, whereas python is 0-based) 
                if refhere is None or readhere is None:
                    print("a problem arises with seq " + seq + " and cigar " + cigar + " and md " + md)
                    break
                ref_pos += 1
                read_pos += 1
                # if you're a letter, you're mismatching. the letter is given in the ref. you can find the letter in the read in the seq[ref_pos]. 

    return [mismatch_list, ref_seq, read_seq]

def mismatch_table(read,cigar,md,flagsum,phred):
    # wrapper for get_mismatches that also reverse complements if needed, and mirrors around the middle of the read 
    # so you shouldn't have to keep the length

    # also, a convenient place to calculate a pmd score without rewriting large pieces of code

    # "read" here is seq, a character vector

    readlength = len(read)

    # first parse the mismatches
    mms, refseq, readseq = get_mismatches(read,cigar,md)

    # some processing: first, figure out if it's forward or backwards
    # second field of a bam is the 0, 16, ... it's called a flag sum and it parses like this (see bowtie2 manual)
    bin_flagsum = bin(flagsum)[::-1]
    bit_position = 4 #  2^4 = 16 this flag means it's aligned in reverse
    backwards = len(bin_flagsum) > bit_position and bin_flagsum[bit_position] == '1' # backwards is a boolean 
    if backwards: # flip all the positions and reverse complement all the nucleotides. the read in the bam is reverse-complemented if aligned to the reverse strand. 
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        mmsc = []
        for entry in range(0,len(mms)):
            new_entry = [
                complement.get(mms[entry][0]), 
                complement.get(mms[entry][1]), 
                readlength - mms[entry][2] +1
            ]
            mmsc.append(new_entry)
        phred = phred[::-1]
    else:
        mmsc = mms
    # now everything, EXCEPT unmerged but retained reverse mate pairs of paired end reads (which i should still try to catch later; to do), should be 5' -> 3' 

    for entry in range(0,len(mmsc)):
        pos = mmsc[entry][2]
        if pos > readlength/2:
            mmsc[entry][2] = -(readlength - pos +1)

    return mmsc

def rev_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
    return [complement[base] for base in reversed(seq)]
 
def get_pmd(read, stranded):
    # input is a pysam read object
    seq = read.query_sequence
    cigar = read.cigarstring
    md = read.get_tag('MD')
    rawphred = read.query_qualities
    flagsum = read.flag

    # set pmd score parameters . these are their original parameters, but i think they are sensible enough in general.
    P = 0.3
    C = 0.01
    pi = 0.001 

    # important note! in the PMDtools manuscript, they say "DZ=0" in the null model.
    # however, in the PMDtools code, Dz=0.001 in the null model.
    # here i am making the latter choice because i think it makes more biological sense.
    Dn = 0.001

    # find out if you're backwards
    bin_flagsum = bin(flagsum)[::-1]
    bit_position = 4 #  2^4 = 16 this flag means it's aligned in reverse
    backwards = len(bin_flagsum) > bit_position and bin_flagsum[bit_position] == '1' # backwards is a boolean 
    # do something if you are:
    mmsc, refseq, readseq = get_mismatches(seq, cigar, md)
    # adjust phred index if there are things happening in the read
    phred = [0 if base in ['-', 'N'] else rawphred.pop(0) for base in readseq]

    # run through both sequences to add up pmd likelihoods
    refseqlist = list(refseq)
    readseqlist = list(readseq)
    readlength = len(readseqlist)
    if backwards:
        refseqlist = rev_complement(refseqlist)
        readseqlist = rev_complement(readseqlist)
        phred = phred[::-1]
    actual_read_pos = 0 # position in the read / position in the ref (they are aligned now)
    pmd_lik = 1
    null_lik = 1
    pos = 0 # need a separate tracker to cope with indels. if there's a "-" in the read reconstruction because of an insertion in the ref, it should not count as a "position" in the read

    if stranded == "ss":
        for b in range(0,readlength):
            if readseqlist[b] == "-":
                continue # no pos +=1 
            # looking for c-> anything anywhere
            if refseqlist[b] == "C" and (readseqlist[b] == "T" or readseqlist[b] == "C"):
                # everything is relevant to 5 prime and 3 prime ends, get both distances
                epsilon = 1/3 * 10**(-phred[b] / 10)
                z = pos + 1 # pos 1 has distance 1, from pmd manuscript
                y = readlength - pos
                Dz = ((1-P)**(z-1))*P + C
                Dy = ((1-P)**(y-1))*P + C
                if readseqlist[b] == "T": # ss m
                    pmd_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dz)*(1-Dy) + (1-pi)*epsilon*Dz*(1-Dy) + (1-pi)*epsilon*Dy*(1-Dz) + pi*epsilon*(1-Dz)*(1-Dy))
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn)*(1-Dn) + (1-pi)*epsilon*Dn*(1-Dn) + (1-pi)*epsilon*Dn*(1-Dn) + pi*epsilon*(1-Dn)*(1-Dn))                 
                if readseqlist[b] == "C": # ss m
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
                z = pos + 1   # 5 prime
                Dz = ((1-P)**(z-1))*P + C

                if readseqlist[b] == "T": # ds mm 
                    pmd_lik *=  1 - ((1-pi)*(1-epsilon)*(1-Dz) + (1-pi)*epsilon*Dz + pi*epsilon*(1-Dz)) 
                    null_lik *= 1 - ((1-pi)*(1-epsilon)*(1-Dn) + (1-pi)*epsilon*Dn + pi*epsilon*(1-Dn)) 
                
                if readseqlist[b] == "C": # ds match
                    pmd_lik *= (1-pi)*(1-epsilon)*(1-Dz) + (1-pi)*epsilon*Dz + pi*epsilon*(1-Dz)
                    null_lik *= (1-pi)*(1-epsilon)*(1-Dn) + (1-pi)*epsilon*Dn + pi*epsilon*(1-Dn)

            if refseqlist[b] == "G" and (readseqlist[b] == "A" or readseqlist[b] == "G"):
                # get distance and stuff to 3 prime end
                epsilon = 1/3 * 10**(-phred[b] / 10) # phred score 30 gives an error rate of 0.001 (then * 1/3)
                z = readlength - pos  # 3 prime
                Dz = ((1-P)**(z-1))*P + C
                if readseqlist[b] == "A": # ds mm
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

# bunch of handy kmer functions
def generate_kmers_recursive(prefix, k, kmers):
    if k == 0:
        kmers.append(prefix)
        return
    for nucleotide in "ACGT":
        generate_kmers_recursive(prefix + nucleotide, k - 1, kmers)
def generate_kmers(k):
    kmers = []
    generate_kmers_recursive("", k, kmers)
    return kmers
def create_kmer_index(kmers):
    kmer_index = {kmer: index for index, kmer in enumerate(kmers)}
    return kmer_index
def generate_kmer_table(read, k):
    kmer_table = {}
    for i in range(len(read) - k + 1):
        kmer = read[i:i + k]
        if kmer in kmer_table:
            kmer_table[kmer] += 1
        else:
            kmer_table[kmer] = 1
    return kmer_table
def map_kmers_to_index(kmer_table, kmer_index):
    mapped_kmers = {kmer_index[kmer]: count for kmer, count in kmer_table.items() if kmer in kmer_index}
    return mapped_kmers
def calculate_kmer_complexity(kmers, k=5):
    # this takes in a dict of kmer counts, keyed by their indices, and returns an entropy
    num_possible_kmers = 4**k
    entropy = 0
    for key in kmers:
        entropy += - kmers[key]/num_possible_kmers * math.log(kmers[key]/num_possible_kmers)
    pairwise_differences = 0
    keys = list(kmers.keys())
    sumkeys = 0
    for i in range(len(keys)):
        sumkeys += kmers[keys[i]]
        for j in range(i + 1, len(keys)):
            pairwise_differences += abs(kmers[keys[i]] - kmers[keys[j]])
    gini = pairwise_differences / (2*len(kmers)*sumkeys)
    # "most equally distributed" gini is 0, "least equally distributed" gini is 1
    # entropy is less well defined, so let's use gini for now 
    return gini

def gather_subs_and_kmers(bamfile_path, lcafile_path, k, upto,stranded,lcaheaderlines):
    # takes in a list of node ids as determined from shorten_files, and outputs the damage distribution at each of these nodes. 
    # Run through the bam/lca once and assign kmer complexity and damage to each node you care about. 
    print("Gathering substitution and kmer metrics per node...")
    numheaderlines = lcaheaderlines 

    # initialize kmer list so that indexing is fast and easy
    kmers = generate_kmers(k)
    kmer_index = create_kmer_index(kmers)
    
    # initialize other stuff
    # node_data = {node: {'total_reads': 0, 'pmdsover2': 0, 'pmdsover4': 0, 'dp1': 0, 'dpm1' = 0, 'meanlength': 0, 'total_alignments': 0, 
    # 'ani' : 0, 'avgperreadgini' : 0, 'avggc': 0, 'tax_path': "", 'subs': {}, 'kmer_count': {}}
    node_data = {}
    bamfile = pysam.AlignmentFile(bamfile_path, "rb") 
    lcafile = open(lcafile_path, 'r')
    oldreadname = ""
    nodestodumpinto = []
    num_alignments = 0
    currentsubdict = {}
    nms = 0 
    pmdsover2 = 0
    pmdsover4 = 0
    stop = False # for debugging

    for _ in range(numheaderlines +1):
        currentlcaline = next(lcafile) # assumes it has the usual 2 comment lines at the top from ngslca

    for read in bamfile:

        # get the basic info for this read
        readname = read.query_name

        # immediately find out if it's a new read. if so, you JUST finished the last read, so do a bunch of stuff for it.
        if readname != oldreadname and oldreadname != "":

            # go get the k-mer table now
            kmer_table = generate_kmer_table(seq,k)
            kmer_indices = map_kmers_to_index(kmer_table,kmer_index)

            # get the lca entry and nodes we wanna update
            lcaentry = currentlcaline.split('\t')
            fields = lcaentry[1:]
            nodestodumpinto = []
            for i in range(len(fields)):
                if upto in fields[i]:
                    nodestodumpinto.append(fields[i].split(':')[0])
                    break
                nodestodumpinto.append(fields[i].split(':')[0])
            lcareadnamesplit = lcaentry[0].split(':')[0:7]
            lcareadname = ":".join(lcareadnamesplit)
            if oldreadname != lcareadname:
                print("hey there is a mismatch between your lca and bam files at read " + oldreadname + " in the bam and " + lcareadname + " in the lca. ")
                print("this mismatch could have some bad downstream implications, best to find out what's causing it. perhaps you have some reads in one file but not in the other? or you forgot to order the bam file? \n")
                break

            # now update everything to all the relevant nodes
            for node in nodestodumpinto:
                if node not in node_data:
                    # that's ok! add it. 
                    node_data[node] = {'total_reads': 0,'pmdsover2': 0, 'pmdsover4': 0, 'meanlength': 0, 'total_alignments': 0, 
                                       'ani': 0, 'avgperreadgini' : 0, 'avggc': 0, 'tax_path' : "", 'subs': {}, 'kmer_count': {},
                                       'dp1' : 0, 'dm1' : 0}
                node_data[node]['meanlength'] = ((node_data[node]['meanlength'] * node_data[node]['total_reads']) + readlength) / (node_data[node]['total_reads'] + 1)
                node_data[node]['avgperreadgini'] = ( (node_data[node]['avgperreadgini'] * node_data[node]['total_reads']) + calculate_kmer_complexity(kmer_table)) / (node_data[node]['total_reads'] + 1)
                ani_for_this_read = (readlength - nms/num_alignments)/readlength 
                node_data[node]['ani'] = (ani_for_this_read + node_data[node]['ani'] * node_data[node]['total_reads']) / (node_data[node]['total_reads'] + 1)
                gc_content_for_this_read = (seq.count('C') + seq.count('G')) / readlength
                node_data[node]['avggc'] = ((node_data[node]['avggc'] * node_data[node]['total_reads']) + gc_content_for_this_read) / (node_data[node]['total_reads'] + 1)
                node_data[node]['total_reads'] += 1
                node_data[node]['total_alignments'] += num_alignments
                node_data[node]['pmdsover2'] += pmdsover2 / num_alignments
                node_data[node]['pmdsover4'] += pmdsover4 / num_alignments
                # only consider a transition snp "damage" if it's in every alignment of a read! 
                ctp1 = currentsubdict.get("['C', 'T', 1]", 0) # c -> t on the pos 1 
                if ctp1 == num_alignments: 
                    node_data[node]['dp1'] += 1 
                if stranded == "ss":
                    ctm1 = currentsubdict.get("['C', 'T', -1]", 0) # c -> t on the minus 1 
                    if ctm1 == num_alignments:
                        node_data[node]['dm1'] += 1 
                if stranded == "ds":
                    gam1 = currentsubdict.get("['G', 'A', -1]", 0)
                    if gam1 == num_alignments:
                        node_data[node]['dm1'] += 1 

                # updates kmer counts
                for kmer, count in kmer_indices.items(): 
                    if kmer in node_data[node]['kmer_count']:
                        node_data[node]['kmer_count'][kmer] += count
                    else:
                        node_data[node]['kmer_count'][kmer] = count
                # updates substitution tables similarly
                if currentsubdict:
                    for sub, count in currentsubdict.items():
                        if sub in node_data[node]['subs']: 
                            node_data[node]['subs'][sub] += count / num_alignments
                        else:
                            node_data[node]['subs'][sub] = count / num_alignments # so, this can be up to 1 per node. 
                # add the tax path if it's not already there
                if node_data[node]['tax_path'] == "":
                    lca_index = next(i for i, entry in enumerate(lcaentry) if entry.startswith(node))
                    tax_path = ','.join(lcaentry[lca_index:]).replace('\n','')
                    node_data[node]['tax_path'] = tax_path

            # move on. re initialize a bunch of things here 
            oldreadname = readname
            currentlcaline = next(lcafile)
            currentsubdict = {}
            num_alignments = 0
            nms = 0
            pmdsover2 = 0
            pmdsover4 = 0

        # now for the current read
        seq = read.query_sequence
        readlength = len(seq)
        cigar = read.cigarstring
        md = read.get_tag('MD')
        nms += read.get_tag('NM')
        pmd = float(read.get_tag('DS'))
        if(pmd>2):
            pmdsover2 += 1
        if(pmd>4):
            pmdsover4 += 1
        phred = read.query_qualities
        flagsum = read.flag
        num_alignments += 1 

        # get the mismatch table for this read, don't bother rerunning the internal pmd calculation though
        subs = mismatch_table(seq,cigar,md,flagsum,phred) 
        for sub in subs:
            key = "".join(str(sub))
            if key in currentsubdict:
                currentsubdict[key] +=1
            else:
                currentsubdict[key] = 1

        for key in subs:
            if key[2] == 0:
                print("Something is wrong. Printing the problem read. Give the following info to Bianca please:")
                print(f"Substitution: {key}, Read: {read.query_name}, md: {md}, cigar: {cigar}")

        if stop:
            break

        # quick catch for the starting read
        if oldreadname == "":
            oldreadname = readname

    print("Gathered substitution and kmer data for " + str(len(node_data)) + " taxonomic nodes. Now processing to compute damage and k-mer complexity...\n ")

    bamfile.close() 
    lcafile.close()

    return node_data

def format_subs(subs, nreads, stranded):
    formatted_subs = []
    other_subs = {}
    
    for key, value in subs.items():
        # Extract position and check if it is within the range -15 to 15
        parts = key.strip("[]").replace("'", "").split(", ")
        pos = int(parts[2])
        if -15 <= pos <= 15:
            substitution = parts[0] + parts[1]
            formatted_key = "".join(parts)
            formatted_value = round(value / nreads, 3)
            
            if substitution == 'CT': # keep any c>t
                formatted_subs.append((pos, f"{formatted_key}:{formatted_value}"))
            elif substitution == 'GA' and stranded == "ds": # keep g>a only if you are double stranded
                formatted_subs.append((pos, f"{formatted_key}:{formatted_value}"))
            else:
                if pos not in other_subs:
                    other_subs[pos] = 0.0
                other_subs[pos] += formatted_value
    
    # Add the summarized 'O' substitutions to the list
    for pos, value in other_subs.items():
        formatted_subs.append((pos, f"O{pos}:{round(value, 3)}"))
    
    # Sort the formatted_subs based on the specified order
    formatted_subs.sort(key=lambda x: (x[0] > 0, (x[0])))
    
    # Return the formatted substitutions as a string
    return " ".join(sub[1] for sub in formatted_subs)

def parse_and_write_node_data(nodedata, stats_path, subs_path, k, stranded):
    # parses a dictionary where keys are node tax ids, and entries are total_reads, meanlength, total_alignments, subs and kmers

    statsfile = open(stats_path, 'w', newline='')
    subsfile = open(subs_path, 'w', newline='')
    header = ['TaxNodeID', 'TaxName', 'TotalReads', 'PMDsover2', 'PMDSover4', 'Damaged+1', 'Damaged-1', 'TotalAlignments', 'MeanLength', 'ANI', 'PerReadKmerGI', 'ReadSetKmerGI', 'AvgGC', 'taxpath'] # ...divergence... damage... 
    writer = csv.writer(statsfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
    subswriter = csv.writer(subsfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
    writer.writerow(header)

    rows = []
    subsrows = {}

    for node in nodedata:
        tn = nodedata[node]
        
        # get formatted subs
        fsubs = format_subs(tn['subs'], tn['total_reads'], stranded)

        # get kmer metrics
        kg = calculate_kmer_complexity(tn['kmer_count'], k) # read set kmer gini index: between 0 and 1. 0 is good, 1 is bad
        taxname = tn['tax_path'].split(",")[0].split(":")[1]

        row = [int(node), taxname, tn['total_reads'], round(tn['pmdsover2'], 2), round(tn['pmdsover4'], 2), tn['dp1'], tn['dm1'], tn['total_alignments'], 
               round(tn['meanlength'], 4), round(tn['ani'], 5), round(tn['avgperreadgini'], 4), round(kg, 4), round(tn['avggc'], 3), tn['tax_path']] 
        rows.append(row)

        subsrows[int(node)] = [int(node), taxname, fsubs]
    
    # Sort rows by total reads
    rows.sort(key=lambda x: x[2], reverse=True)

    # Write sorted rows to stats file
    for row in rows:
        writer.writerow(row)
    
    # Write sorted subsrows based on the order of sorted rows 
    for row in rows:
        subswriter.writerow(subsrows[row[0]])

    statsfile.close()
    subsfile.close()

    print("Wrote final stats and subs files. Done!")

def main(in_lca, in_bam, out_lca, out_bam, out_stats, out_subs, stranded, mincount, k, upto, lcaheaderlines, minsim, other_lca_keywords):

    # STEP 1: Reduce the size of the lca and bam files by removing things in both which:
    #  - don't meet your tax threshold 
    #  - don't meet your min node read count 
    #  - don't appear in the lca file (just got filtered from crap evenness of coverage or couldn't get assigned due to some taxonomy issue)
    # While we're at it, add a PMD tag to each read in the shortened bam file. 

    # Moving things around a little bit. Step 1.1 is now to rewrite just the lca file in two passes.
    write_shortened_lca(in_lca,out_lca,upto,mincount,lcaheaderlines,other_lca_keywords)
    write_shortened_bam(in_bam,out_lca,out_bam,stranded,lcaheaderlines,minsim,other_lca_keywords) 

    # STEP 2: calculate substitution table, kmers, etc per node. write temp file. compute damage and k-mer read set complexity. write subs and stats files. 
    nodedata = gather_subs_and_kmers(out_bam, out_lca, k = k, upto = upto,stranded = stranded, lcaheaderlines = lcaheaderlines)
    parse_and_write_node_data(nodedata,out_stats,out_subs,k,stranded)  


if __name__ == "__main__":
    
    # Initialize
    parser = argparse.ArgumentParser(
        description="BamDam processes LCA and bam files for ancient environmental DNA. Please ensure your input bam is sorted in the same order as your input lca file, and that every read in your LCA file is also \
          in your bam file (the opposite need not be true). This will be true by default if you use fresh ngsLCA output. Also, please provide a ngsLCA lca file, not a metaDMG lca file (they have different formats).",
        epilog="\
          BamDam was written by Bianca De Sanctis in 2024. For assistance please contact bddesanctis@gmail.com.")
    
    # Mandatory arguments
    parser.add_argument("--in_lca", type=str, required=True, help="Path to the original (sorted) LCA file")
    parser.add_argument("--in_bam", type=str, required=True, help="Path to the original (sorted) BAM file")
    parser.add_argument("--out_lca", type=str, required=True, help="Path to the short output LCA file")
    parser.add_argument("--out_bam", type=str, required=True, help="Path to the short output BAM file")
    parser.add_argument("--out_stats", type=str, required=True, help="Path to the output stats file")
    parser.add_argument("--out_subs", type=str, required=True, help="Path to the output subs file")
    parser.add_argument("--stranded", type=str, required=True, help="Either ss for single stranded or ds for double stranded")

    # Optional arguments with defaults
    parser.add_argument("--mincount", type=int, default=5, help="Minimum read count to keep a node (default: 5)")
    parser.add_argument("--k", type=int, default=5, help="Value of k for kmer complexity calculations (default: 5)")
    parser.add_argument("--upto", type=str, default="family", help="Keep nodes up to and including this tax threshold, use root to disable (default: family)")
    parser.add_argument("--lcaheaderlines", type=int, default=0, help="Number of header lines in input LCA file (default: 0)")
    parser.add_argument("--minsim", type=float, default=0.95, help="Minimum similarity to reference to keep a read (default: 0.95)")
    parser.add_argument("--keep_keywords", type=str, nargs='+', default=[], help="Other keyword(s) in LCA file for filtering to keep, e.g. Eukaryota (default: none)")
    parser.add_argument("--exclude_keywords", type = str, nargs='+', default=["Hominidae"], help="Other keyword(s) in LCA file for filtering to delete (default Hominidae)- NOT IMPLEMENTED YET")
    parser.add_argument("--exclude_keyword_file", type = str, default="", help="Text file containing list of keywords or tax paths to remove in LCA file, one on each line- NOT IMPLEMENTED YET")

    if '--help' in sys.argv or '-h' in sys.argv:
        parser.print_help()
        sys.exit()
    
    args = parser.parse_args()

    # Validation checks
    if args.stranded not in ["ss", "ds"]:
        parser.error(f"Invalid value for stranded: {args.stranded}. Must be 'ss' or 'ds'.")
    if not isinstance(args.mincount, int):
        parser.error(f"Invalid integer value for mincount: {args.mincount}")
    if not isinstance(args.k, int) or args.k > 10:
        parser.error(f"Invalid integer value for k: {args.k} (max 10)")
    if not re.match("^[a-z]+$", args.upto):
        parser.error(f"Invalid value for upto: {args.upto}. Must be a string of only lowercase letters.")
    if not isinstance(args.lcaheaderlines, int):
        parser.error(f"Invalid integer value for lcaheaderlines: {args.lcaheaderlines}")
    if not isinstance(args.minsim, float):
        parser.error(f"Invalid float value for minsim: {args.minsim}")
    if not os.path.exists(args.in_lca):
        parser.error(f"Input LCA path does not exist: {args.in_lca}")
    if not os.path.exists(args.in_bam):
        parser.error(f"Input BAM path does not exist: {args.in_bam}")
    
    # Print message before calling main
    print("Hello! You are running BamDam with the following arguments:")
    print(f"in_lca: {args.in_lca}")
    print(f"in_bam: {args.in_bam}")
    print(f"out_lca: {args.out_lca}")
    print(f"out_bam: {args.out_bam}")
    print(f"out_stats: {args.out_stats}")
    print(f"out_subs: {args.out_subs}")
    print(f"stranded: {args.stranded}")
    print(f"mincount: {args.mincount}")
    print(f"k: {args.k}")
    print(f"upto: {args.upto}")
    print(f"lcaheaderlines: {args.lcaheaderlines}")
    print(f"minsim: {args.minsim}")
    print(f"other_lca_keywords: {args.other_lca_keywords}")
 
    main(
        args.in_lca, args.in_bam, args.out_lca, args.out_bam, 
        args.out_stats, args.out_subs, args.stranded, 
        args.mincount, args.k, args.upto, args.lcaheaderlines, args.minsim, args.other_lca_keywords
    )




###### Comment section for Bianca. 

# Sometimes ngslca spits out duplicate lines for no reason. I don't know why. Get rid of them before running this!! Like this
# awk '!seen[$0]++' input_lca > output_lca

# EXAMPLE
''' 
python BamDam4.py \
    --in_lca "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7008961409-LV7005366316-LV3005888478_S32_onlyeukaryotfamilies.score96.lca" \
    --in_bam "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7008961409-LV7005366316-LV3005888478_S32.sorted.bam" \
    --out_lca "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7008961409-LV7005366316-LV3005888478_S32_shortened.lca" \
    --out_bam "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7008961409-LV7005366316-LV3005888478_S32_shortened.bam" \
    --out_stats "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7008961409-LV7005366316-LV3005888478_S32_allstats.txt" \
    --out_subs "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7008961409-LV7005366316-LV3005888478_S32_allsubs.txt" \
    --stranded ds \
    --mincount 5 \
    --k 5 \
    --upto "family" \
    --lcaheaderlines 0 \
    --minsim 0.96 \
    --other_lca_keywords Eukaryota
'''

# 
in_lca = "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.lca"
out_lca = "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.shortened.lca"
in_bam = "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.sort.bam"
out_bam = "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.shortened.bam"
out_stats = "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.stats.txt"
out_subs = "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.subs.txt"
stranded = "ds"
mincount = "5"
upto = "family"
minsim = 0.95
other_lca_keywords = ["Eukaryota"]
lcaheaderlines = 2

# and then 
# ./PlotDamage.R -f "/Users/bianca/Dropbox/Documents/academic/postdoc_durbin/betadmg/data/LV7001882898-LV7005224373-CGG3-012108.subs.txt" -t "48230" -o "dryas.png" -s "ds"

### BIANCA STILL TO DO ###
# - write a plotting function in R that runs on (a) the subs file to plot damage and (b) the lca + bam to plot pmd distribution
# - write another plotting function for read length distribution
# - write a python function, separately from this, to quickly extract reads from a bam or fq file assigned to a specific genus (see below) 
# - consider writing a function to add some metadmg parameters to the output stats file 
# - consider writing an option to NOT include nodes under ones with mincount if those underneath nodes don't meet mincount themselves
    # (you could just pass an option keepunder = False or something like this)
# - implement progress bars w tqdm
# - write a function to combine multiple stats and subs files from multiple different files, like from ellesmere 
    # then you can combine the "read-pulling-out function" with this and make a full-on damage distribution for different taxa in ellesmere! 
# - consider outputting a metadmg-like Dfit once you understand what the hell it is doing 
# - add a logger
# - implement remove keywords from lca like Homonidae , maybe take in a whole txt file list of keywords from controls lcas?
# - consider "damtools" 
# - consider rewriting bamfilter 
# - rory: are more strands ending in purines (A vs G) than anything else? could check this! https://www.pnas.org/doi/abs/10.1073/pnas.0704665104 
# - output ani without the damage substitutions . could make it comparable b/w ss and ds by cutting the read in half 
# - clean up the keyword 

# here is a bash line to get all reads associated to a certain tax. but it is slow and opens the whole file and i could do it faster.
'''
wanna see all reads associated with a given tax node? do:
grep "Alnus" LV7008961409-LV7005366316-LV3005888478_S32_shortened.lca | \
awk '{print $1}' | \
awk '{gsub(/^M_/, ""); split($0, fields, ":"); print fields[1]":"fields[2]":"fields[3]":"fields[4]":"fields[5]":"fields[6]":"fields[7]}' | \
grep -F -f - <(samtools view LV7008961409-LV7005366316-LV3005888478_S32_shortened.bam) | awk {'print $9'} | uniq -c | sort
i should write a lil function to do this faster later. grep is slow and will freak out if it has to open a 1tb bam file. 
''' 




