
# bamdam compute

import sys 
import re
import csv
import pysam
import math
import argparse
import os
import hyperloglog
import subprocess
import random
from functools import lru_cache

try: # optional library only needed for plotting 
    import matplotlib.pyplot as plt
    matplotlib_imported = True
except:
    matplotlib_imported = False

# local
from bamdam import utils
from bamdam import alignment_utils

def run_compute(args):
    lca_file_type = utils.find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("Error: It looks like you're trying to run bamdam compute with a metaDMG-style lca file. Please use an ngsLCA-style lca file. Bamdam shrink can automatically convert it for you.")
        sys.exit()
    nodedata, pmds_in_bam = gather_subs_and_kmers(args.in_bam, args.in_lca, kn=args.k, upto=args.upto, mode=args.mode, stranded=args.stranded, show_progress=args.show_progress, cpg_split=args.udg)
    parse_and_write_node_data(nodedata, args.out_tsv, args.out_subs, args.stranded, pmds_in_bam, args.mode, cpg_split=args.udg) 
    if args.plotdupdust:
        plotdupdust(args.out_tsv, args.plotdupdust) 

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
                rep_kmers.append(utils.get_rep_kmer(kmer))
                total_kmers +=1 
    else:
        print(f"Warning: One of your reads is shorter than k.")
    return rep_kmers, total_kmers

def gather_subs_and_kmers(bamfile_path, lcafile_path, kn, upto, mode, stranded, show_progress = False, cpg_split = False):
    print("\nGathering substitution and kmer metrics per node...")
    # this function is organized in a potentially unintuitive way. it uses a bunch of nested loops to pop between the bam and lca files line by line. 
    # it matches up bam read names and lca read names and aggregates some things per alignment, some per read, and some per node, the last of which are added into a large structure node_data.
    # altogether this uses very little ram

    lcaheaderlines = 0
    with open(lcafile_path, 'r') as lcafile:
        for lcaline in lcafile:
            if "no rank" in lcaline and "#" not in lcaline:
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
    best_alignments = []
    best_alignment_score = -1000
    refgc = 0

    # do we want a progress bar?
    progress_bar = None
    if show_progress:
        try: # optional library for progress bars 
            from tqdm import tqdm
            tqdm_imported = True
        except ImportError:
            tqdm_imported = False
        if not tqdm_imported: # you can't have one after all
            print("The library tqdm is not available, so progress bars will not be shown. This will not impact performance.")
        else:
            totallcalines = utils.line_count(lcafile_path)
            # should be super fast compared to anything else; worth it to initiate a progress bar if you want it. 
            progress_bar = tqdm(total=totallcalines-lcaheaderlines, unit='lines', disable=not tqdm_imported) if tqdm_imported else None
            if progress_bar:
                update_interval = 100 # this is arbitrary of course. could also be e.g. (totallcalines // 1000) 

    for _ in range(lcaheaderlines +1):
        currentlcaline = next(lcafile) 

    for alignment in bamfile:
        readname = alignment.query_name
        # find out if it's a new read. if so, you just finished the last read, so do a bunch of stuff for it. 
        # the first read will skip this if statement because of the second condition
        if readname != oldreadname and oldreadname != "":

            if mode == 1:
                # compute on the best alignment here, now that we're done with that read name
                if len(best_alignments) > 1:
                    # pick one
                    best_alignment = random.choice(best_alignments)
                elif len(best_alignments) == 1:
                    best_alignment = best_alignments[0]
                else:
                    print("Error: unexpected behavior. Can't find a best alignment for this read.")
                    sys.exit(1)
                seq = best_alignment.query_sequence
                readlength = len(seq)
                cigar = best_alignment.cigarstring
                md = best_alignment.get_tag('MD')
                nms = best_alignment.get_tag('NM')
                if are_pmds_in_the_bam: 
                    try:
                        pmd = float(best_alignment.get_tag('DS'))
                    except KeyError:
                        pmd = 0
                    if(pmd>2):
                        pmdsover2 += 1 
                    if(pmd>4):
                        pmdsover4 += 1
                flagsum = best_alignment.flag
                # go and get the mismatch table for this best alignment
                subs, matches, reference_seq = alignment_utils.mismatch_table(seq,cigar,md,flagsum) 
                refgc = (reference_seq.count('C') + reference_seq.count('G')) / len(reference_seq)
                allsubs = subs + matches
                for sub in allsubs:
                    key = tuple(sub)
                    if key in currentsubdict:
                        currentsubdict[key] +=1
                    else:
                        currentsubdict[key] = 1
                # other modes will already have that info stored, as it was being tracked over all alignments as we ran through them
                # okay. now add CpG-annotated entries if flag is set
                if cpg_split:
                    # this is kind of a slow way to do this particular thing but whatever
                    for sub in allsubs:
                        from_base, to_base, pos = sub
                        if from_base == 'C':
                            next_pos_exists = False
                            is_cpg = False
                            for other_sub in allsubs:
                                if other_sub[2] == pos + 1:
                                    next_pos_exists = True
                                    if other_sub[0] == 'G':
                                        is_cpg = True
                                    break
                            if next_pos_exists:
                                if is_cpg:
                                    cpg_key = ('CpG', to_base, pos)
                                else:
                                    cpg_key = ('nonCpG', to_base, pos)
                                if cpg_key in currentsubdict:
                                    currentsubdict[cpg_key] += 1
                                else:
                                    currentsubdict[cpg_key] = 1 

            # do k-mer things for this read
            dust = utils.calculate_dust(seq)
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
                            'ani': 0, 'avgdust' : 0, 'avgreadgc': 0, 'avgrefgc': 0,  'tax_path' : "", 'subs': {}, 'damagerelevantreadsp1': 0, 'damagerelevantreadsm1': 0,
                            'dp1' : 0, 'dm1' : 0,  'hll': hyperloglog.HyperLogLog(0.01), 'totalkmers' : 0, 'unaggregatedreads' : 0}
                    else:
                        node_data[node] = {'total_reads': 0,'meanlength': 0, 'total_alignments': 0, 
                            'ani': 0, 'avgdust' : 0, 'avgreadgc': 0, 'avgrefgc': 0, 'tax_path' : "", 'subs': {}, 'damagerelevantreadsp1': 0, 'damagerelevantreadsm1': 0,
                            'dp1' : 0, 'dm1' : 0, 'hll': hyperloglog.HyperLogLog(0.01), 'totalkmers' : 0, 'unaggregatedreads' : 0}
                        
                # now populate/update it        
                node_data[node]['meanlength'] = ((node_data[node]['meanlength'] * node_data[node]['total_reads']) + readlength) / (node_data[node]['total_reads'] + 1)

                if dust != -1:
                    node_data[node]['avgdust'] = ( (node_data[node]['avgdust'] * node_data[node]['total_reads']) + dust) / (node_data[node]['total_reads'] + 1)
                else:
                    readswithNs +=1 

                # in theory, ani = (1 - num_mismatches) / readlength 
                # but you will have multiple values because each read has multiple alignments, so they need to be dealt with accordingly
                if mode == 1:
                    # nms will be made from 1 (the best) alignment. easy to transform to ani
                    ani_for_this_read = (readlength - nms) / readlength
                    node_data[node]['ani'] = (ani_for_this_read + node_data[node]['ani'] * node_data[node]['total_reads']) / (node_data[node]['total_reads'] + 1)
                if mode == 2:
                    # nms will be sum of num_alignments alignments. average over alignments first, then get ani
                    ani_for_this_read = (readlength - (nms/num_alignments))/readlength 
                    node_data[node]['ani'] = (ani_for_this_read + node_data[node]['ani'] * node_data[node]['total_reads']) / (node_data[node]['total_reads'] + 1)
                if mode == 3:
                    # nms will be sum of num_alignments alignments
                    # logic: say i have num_alignments = 3 and i kept all the nms separate, so there are three of them, and therefore three values for ani
                    # then i want to get the summed ani, so sum_ani =( (readlength -nms[1] /readlength + (readlength -nms[2]) /readlength  + (readlength -nms[3])/readlength )
                    # or (num_alignments - nms/(readlength*num_alignments) ), or
                    ani_for_this_read =  (readlength*num_alignments - nms) / readlength
                    node_data[node]['ani'] = (ani_for_this_read + node_data[node]['ani'] * node_data[node]['total_alignments']) / (node_data[node]['total_alignments'] + num_alignments)
                
                gc_content_for_this_read = (seq.count('C') + seq.count('G')) / readlength
                node_data[node]['avgreadgc'] = ((node_data[node]['avgreadgc'] * node_data[node]['total_reads']) + gc_content_for_this_read) / (node_data[node]['total_reads'] + 1)
                if mode == 1:
                    # there will be 1 contributor to refgc
                    node_data[node]['avgrefgc'] = ((node_data[node]['avgrefgc'] * node_data[node]['total_reads']) + refgc) / (node_data[node]['total_reads'] + 1)
                if mode == 2:
                    # there will be num_alignments contributors to refgc, but we want to average it 
                    node_data[node]['avgrefgc'] = ((node_data[node]['avgrefgc'] * node_data[node]['total_reads']) + (refgc/num_alignments) ) / (node_data[node]['total_reads'] + 1)
                if mode == 3:
                    # there will be num_alignments contributors to refgc
                    node_data[node]['avgrefgc'] = ((node_data[node]['avgrefgc'] * node_data[node]['total_alignments']) + refgc) / (node_data[node]['total_alignments'] + num_alignments)

                # only now we can increment total alignments
                node_data[node]['total_alignments'] += num_alignments

                # you can do pmds with the modes too
                # if mode == 1, then you can set num_alignments = 1
                if are_pmds_in_the_bam:
                    if mode == 2:
                        node_data[node]['pmdsover2'] += pmdsover2 / num_alignments
                        node_data[node]['pmdsover4'] += pmdsover4 / num_alignments
                    if mode == 1 or mode == 3:
                        node_data[node]['pmdsover2'] += pmdsover2 
                        node_data[node]['pmdsover4'] += pmdsover4 
                        # we fix the denominator later in the mode = 3 case (in the next function)

                # update hyperloglogs
                for kmer in rep_kmers:
                    node_data[node]['hll'].add(kmer)
                node_data[node]['totalkmers'] += total_kmers

                other_sub_count = 0 # useless for now
                if currentsubdict:
                    for sub, count in currentsubdict.items():
                        if not ((sub[0] == 'C' and sub[1] == 'T') or (sub[0] == 'G' and sub[1] == 'A')):
                            other_sub_count += count # don't include c>t or g>a in any case, regardless of library 
                        if sub in node_data[node]['subs']: 
                            if mode == 1 or mode == 3:
                                node_data[node]['subs'][sub] += count
                            if mode == 2:
                                node_data[node]['subs'][sub] += count / num_alignments
                        else:
                            if mode == 1 or mode == 3:
                                node_data[node]['subs'][sub] = count
                            if mode == 2:
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
            # outside of the loop, we should update the unaggregated read count only for the actual node we hit, not all the ones on top
            node_data[nodestodumpinto[0]]['unaggregatedreads'] += 1

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
            best_alignments = []
            best_alignment_score = -1000
            nms = 0
            refgc = 0
            if are_pmds_in_the_bam:
                pmdsover2 = 0
                pmdsover4 = 0

        # now for the current alignment.
        # the following might change for different alignments of the same read: 
        num_alignments += 1 

        if mode == 1:
            try:
                alignment_score = alignment.get_tag("AS")
            except KeyError:
                print("Error: There are not AS (alignment score) tags in your bam file. Bowtie2 produces them by default, but other aligners might not. These tags are needed to determine the best alignment in compute mode 1. You can either run a different compute mode (we suggest mode 2), or annotate your bam and try mode 1 again.")
                sys.exit(1)

            # much less computing needed, hold it all off until the end when we know the best alignment
            if best_alignments == []:
                # nothing tracked yet, initialize
                best_alignments = [alignment]
                best_alignment_score = alignment_score
            else:
                if alignment_score > best_alignment_score:
                    best_alignments = [alignment]
                    best_alignment_score = alignment_score
                elif alignment_score == best_alignment_score:
                    best_alignments.append(alignment)
                else:
                    continue  # onto the next alignment! until you are done with the read
                # then pick a random best alignment at the end of the read and compute on it

        if mode == 2 or mode == 3:
            # we are computing on every alignment
            seq = alignment.query_sequence
            readlength = len(seq)
            cigar = alignment.cigarstring
            md = alignment.get_tag('MD')
            nms += alignment.get_tag('NM')
            if are_pmds_in_the_bam: 
                try:
                    pmd = float(alignment.get_tag('DS'))
                except KeyError:
                    pmd = 0
                if(pmd>2):
                    pmdsover2 += 1 
                if(pmd>4):
                    pmdsover4 += 1
            flagsum = alignment.flag
            # go and get the mismatch table for this read if the name/md/cigar/flagsum is different to before (this is expensive, so there is a catch to avoid it when possible)
            if (readname != oldreadname) or (cigar != oldcigar) or (md != oldmd) or (flagsum != oldflagsum):
                subs, matches, refseq = alignment_utils.mismatch_table(seq,cigar,md,flagsum) 
                oldcigar = cigar; oldmd = md; oldflagsum = flagsum
            # now refgc
            refgc += (refseq.count('C') + refseq.count('G')) / len(refseq)  # sum of refgc percent over all alignments 
            allsubs = subs + matches
            for sub in allsubs:
                key = tuple(sub)
                if key in currentsubdict:
                    currentsubdict[key] +=1
                else:
                    currentsubdict[key] = 1
            # add CpG-annotated entries if flag is set
            if cpg_split:
                for sub in allsubs:
                    from_base, to_base, pos = sub
                    if from_base == 'C':
                        next_pos_exists = False
                        is_cpg = False
                        for other_sub in allsubs:
                            if other_sub[2] == pos + 1:
                                next_pos_exists = True
                                if other_sub[0] == 'G':
                                    is_cpg = True
                                break
                        if next_pos_exists:
                            if is_cpg:
                                cpg_key = ('CpG', to_base, pos)
                            else:
                                cpg_key = ('nonCpG', to_base, pos)
                            if cpg_key in currentsubdict:
                                currentsubdict[cpg_key] += 1
                            else:
                                currentsubdict[cpg_key] = 1

        # quick catch for the starting read; check if the first read (and then presumably the whole bam) has a pmd score 
        if oldreadname == "":
            oldreadname = readname
            try:
                pmd = float(alignment.get_tag('DS'))
            except KeyError:
                are_pmds_in_the_bam = False

    if progress_bar:
        progress_bar.update(totallcalines - lcaheaderlines - progress_bar.n)
        progress_bar.close()

    bamfile.close() 
    lcafile.close()

    if lcalinesskipped > 0:
        print("\nWarning: " + str(lcalinesskipped) + " reads in the input LCA file did not appear in the input bam file and so were not used. This may have happened if the minimum similarity used in bamdam shrink was more stringent than that used in ngsLCA. This will not affect output statistics, except that these reads will not be included.")
        
    if readswithNs >0 :
        print("\nSkipped k-mer counting and DUST score computation for " + str(readswithNs) + " reads with non-ACGT characters. \n" +
            "If all reads assigned to a taxonomic node have non-ACGT characters, both the mean DUST and duplicity will be 0 for that node.")

    print("\nGathered substitution and kmer data for " + str(len(node_data)) + " taxonomic nodes. Now sorting and writing output files... ")

    return node_data, are_pmds_in_the_bam

def format_subs(subs, nreads, cpg_split):
    formatted_subs = []
    
    for key, value in subs.items():
        # extract position and check if it is within the range -15 to 15 (more is unnecessary for damage)
        # easy to remove the condition if you want to write all of the subs though!
        parts = [str(key[0]), str(key[1]), str(key[2])]
        pos = int(parts[2])
        if (-15 <= pos <= 15) and (parts[0] in {'A', 'C', 'T', 'G'}) and (parts[1] in {'A', 'C', 'T', 'G'}):
            formatted_key = "".join(parts)
            formatted_value = round(value / nreads, 3)  
            formatted_subs.append((pos, f"{formatted_key}:{formatted_value}"))
        if cpg_split:
            # don't forget about the cpg ones, if you need those. 
            is_cpg_entry = parts[0] in {'CpG', 'nonCpG'}
            if (-15 <= pos <= 15) and is_cpg_entry and (parts[1] in {'A', 'C', 'T', 'G'}):
                formatted_key = "".join(parts)
                formatted_value = round(value / nreads, 3)  
                formatted_subs.append((pos, f"{formatted_key}:{formatted_value}"))

    formatted_subs.sort(key=lambda x: (x[0] > 0, (x[0])))

    return " ".join(sub[1] for sub in formatted_subs)




def calculate_node_damage(subs, stranded, cpg_split):

    ctp1 = 0  # C>T at 5' position 1
    ctm1 = 0  # C>T at 3' position -1 for ss
    gam1 = 0  # G>A at 3' position -1 for ds
    c_p1 = 0  # total C at 5' position 1
    c_m1 = 0  # total C at 3' position -1 for ss
    g_m1 = 0  # total G at 3' position -1 for ds

    if cpg_split:
        cpg_ctp1 = 0
        cpg_c_p1 = 0
        noncpg_ctp1 = 0
        noncpg_c_p1 = 0

    for key, count in subs.items():
        from_base, to_base, pos = str(key[0]), str(key[1]), int(key[2])

        # from_base, to_base, pos = parts[0], parts[1], int(parts[2])

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

        # cpg-specific damage calculation (5' only), ds vs ss doesn't matter
        # note 3' is not splittable because we don't know if a G comes after a C on the 3'
        # we could print -2 instead of the 3' -1, which is what pmdtools does, but i don't think it's so informative
        # so i will just add the split 5' c-to-t for udg mode
        if cpg_split:
            if from_base == 'CpG' and pos == 1:
                if to_base == 'T':
                    cpg_ctp1 += count
                cpg_c_p1 += count
            
            if from_base == 'nonCpG' and pos == 1:
                if to_base == 'T':
                    noncpg_ctp1 += count
                noncpg_c_p1 += count


    dp1 = ctp1 / c_p1 if c_p1 > 0 else 0
    dm1 = (ctm1 / c_m1 if c_m1 > 0 else 0) if stranded == "ss" else (gam1 / g_m1 if g_m1 > 0 else 0)

    if cpg_split:
        cpg_dp1 = cpg_ctp1 / cpg_c_p1 if cpg_c_p1 > 0 else 0
        noncpg_dp1 = noncpg_ctp1 / noncpg_c_p1 if noncpg_c_p1 > 0 else 0
        return dp1, dm1, cpg_dp1, noncpg_dp1
    else:
        return dp1, dm1
def parse_and_write_node_data(nodedata, tsv_path, subs_path, stranded, pmds_in_bam, mode, cpg_split = False):
    # parses a dictionary where keys are node tax ids, and entries are total_reads, meanlength, total_alignments, etc 
    
    statsfile = open(tsv_path, 'w', newline='')
    subsfile = open(subs_path, 'w', newline='')
    
    if pmds_in_bam:
        if cpg_split:
            header = ['TaxNodeID', 'TaxName', 'TotalReads','Duplicity', 'MeanDust','Damage+1', 'Damage-1', 'Damage+1_CpG', 'Damage+1_nonCpG',
                'MeanLength', 'ANI','AvgReadGC','AvgRefGC', 'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments', 'UnaggregatedReads', 
                'PMDsover2', 'PMDSover4','taxpath']
        else:
            header = ['TaxNodeID', 'TaxName', 'TotalReads','Duplicity', 'MeanDust','Damage+1', 'Damage-1','MeanLength', 'ANI','AvgReadGC','AvgRefGC',
                'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments', 'UnaggregatedReads', 'PMDsover2', 'PMDSover4','taxpath']
    else:
        if cpg_split:
            header = ['TaxNodeID', 'TaxName', 'TotalReads','Duplicity', 'MeanDust','Damage+1', 'Damage-1', 'Damage+1_CpG', 'Damage+1_nonCpG',
                'MeanLength', 'ANI','AvgReadGC','AvgRefGC', 'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments', 'UnaggregatedReads', 'taxpath']
        else:
            header = ['TaxNodeID', 'TaxName', 'TotalReads','Duplicity', 'MeanDust','Damage+1', 'Damage-1','MeanLength', 'ANI','AvgReadGC','AvgRefGC',
                'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments', 'UnaggregatedReads', 'taxpath']
    
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
        if mode == 1 or mode == 2: 
            # if mode is 1, then subs is based off only the top alignment(s)
            # if mode is 2, then subs is already averaged over all alignments for each read
            # in both cases, you now want to weight it by total_reads
            fsubs = format_subs(tn['subs'], tn['total_reads'],cpg_split)
        if mode == 3:
            # if mode is 3, then subs is currently summed over all alignments 
            fsubs = format_subs(tn['subs'], tn['total_alignments'],cpg_split)

        # number of unique k-mers approximated by the hyperloglog algorithm 
        numuniquekmers = len(tn['hll'])
        if numuniquekmers > 0:
            duplicity = tn['totalkmers'] / numuniquekmers
            uniquekmersperread = numuniquekmers / tn['total_reads']     
            # it used to be the below (ratio duplicated kmers):  
            #  (tn['totalkmers'] - numuniquekmers) / tn['totalkmers']
            # but a reviewer requested we instead do UniqKmersPerRead, so that's what it is now
            # in any case you can get all of this from the existing output in whatever ratio you'd like.
        else:
            duplicity = 0
            uniquekmersperread = 0
        if uniquekmersperread < 0 :
            # you can get tiny negative numbers from essentially zero duplicity and error 
            # (there is up to 1% error in the uniq kmer counting)
            uniquekmersperread = 0

        taxname = tn['tax_path'].split(";")[0].split(":")[1]

        # only now calculate damage per node from the subs dict (see output subs file)
        # do not use formatted subs ; these should be raw numbers: how many READS for this taxa have these matches/mismatches?
        # (possibly this is avg'd over all the alignments per read, so maybe not an integer) 
        # dp1, dm1  = calculate_node_damage(tn['subs'], stranded, cpg_split)
        if cpg_split:
            dp1, dm1, cpg_dp1, noncpg_dp1 = calculate_node_damage(tn['subs'], stranded, cpg_split)
        else:
            dp1, dm1 = calculate_node_damage(tn['subs'], stranded, cpg_split)

        # write , number of columns depends on if you wanted pmds or not
        if pmds_in_bam:
            if mode == 1 or mode == 2:
                pd2 = tn['pmdsover2'] / tn['total_reads']
                pd4 = tn['pmdsover4'] / tn['total_reads']
            if mode == 3:
                pd2 = tn['pmdsover2'] / tn['total_alignments']
                pd4 = tn['pmdsover4'] / tn['total_alignments']
            if cpg_split:
                row = [int(node), taxname, tn['total_reads'], round(duplicity, 3), round(tn['avgdust'], 2), 
                    round(dp1, 4), round(dm1, 4), round(cpg_dp1, 4), round(noncpg_dp1, 4),
                    round(tn['meanlength'], 2), round(tn['ani'], 4), round(tn['avgreadgc'], 3), round(tn['avgrefgc'], 3), 
                    numuniquekmers, round(uniquekmersperread, 3), tn['total_alignments'],
                    round(pd2, 3), round(pd4, 3), tn['unaggregatedreads'], tn['tax_path']]
            else:
                row = [int(node), taxname, tn['total_reads'], round(duplicity, 3), round(tn['avgdust'], 2), 
                    round(dp1, 4), round(dm1, 4), 
                    round(tn['meanlength'], 2),  round(tn['ani'], 4), round(tn['avgreadgc'], 3),round(tn['avgrefgc'], 3), numuniquekmers, round(uniquekmersperread,3),tn['total_alignments'],
                    round(pd2, 3), round(pd4, 3), tn['unaggregatedreads'], tn['tax_path']]
        else:
            if cpg_split:
                row = [int(node), taxname, tn['total_reads'], round(duplicity, 2), round(tn['avgdust'], 2), 
                    round(dp1, 4), round(dm1, 4), round(cpg_dp1, 4), round(noncpg_dp1, 4),
                    round(tn['meanlength'], 2), round(tn['ani'], 4), round(tn['avgreadgc'], 3), round(tn['avgrefgc'], 3),
                    numuniquekmers, round(uniquekmersperread, 3), tn['total_alignments'], tn['unaggregatedreads'], tn['tax_path']]
            else:
                row = [int(node), taxname, tn['total_reads'], round(duplicity, 2), round(tn['avgdust'], 2), 
                    round(dp1, 4), round(dm1, 4), 
                    round(tn['meanlength'], 2), round(tn['ani'], 4), round(tn['avgreadgc'], 3), round(tn['avgrefgc'], 3),
                    numuniquekmers, round(uniquekmersperread, 3), tn['total_alignments'], tn['unaggregatedreads'], tn['tax_path']]
        
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
    


def plotdupdust(tsv_path, plot_path):
    # makes an optional duplicity-dust plot at the genus level. by reviewer request after we put something similar in the manuscript.
    
    if not matplotlib_imported:
        print("Error: couldn't import matplotlib.")
        sys.exit()

    genera = []
    duplicities = []
    dusts = []
    damages = []

    with open(tsv_path, 'r') as f:
        header = next(f).strip().split('\t')
        
        # find column indices
        try:
            duplicity_idx = header.index('Duplicity')
            dust_idx = header.index('MeanDust')
            damage_idx = header.index('Damage+1')
            reads_idx = header.index('TotalReads')
            taxpath_idx = header.index('taxpath')
        except ValueError as e:
            print(f"Error: Could not find required column in header: {e}")
            return
        
        for line in f:
            cols = line.split('\t')
            if len(cols) <= max(duplicity_idx, dust_idx, damage_idx, reads_idx, taxpath_idx):
                continue
            
            taxpath = cols[taxpath_idx]
            reads = int(cols[reads_idx])
            duplicity = float(cols[duplicity_idx])
            dust = float(cols[dust_idx])
            damage = float(cols[damage_idx])
            
            if 'genus' in taxpath and 'species' not in taxpath and 'subgenus' not in taxpath and reads > 50:
                genera.append(taxpath)
                duplicities.append(duplicity)
                dusts.append(dust)
                damages.append(damage)
    
    if not duplicities or not dusts:
        print("No data matching criteria (genus level with >50 reads).")
        return

    # calculate limits with 5% padding
    dup_min, dup_max = min(duplicities), max(duplicities)
    dust_min, dust_max = min(dusts), max(dusts)
    dup_range = dup_max - dup_min
    dust_range = dust_max - dust_min

    plt.figure(figsize=(8,6))
    scatter = plt.scatter(dusts, duplicities, c=damages, cmap='viridis', edgecolor='k')

    plt.colorbar(scatter, label='Damage')
    plt.xlabel('Dust score')
    plt.ylabel('Duplicity')
    plt.title(f'Duplicity vs dust, {os.path.basename(tsv_path)}')

    plt.xlim(dust_min - 0.05 * dust_range, dust_max + 0.05 * dust_range)
    plt.ylim(dup_min - 0.05 * dup_range, dup_max + 0.05 * dup_range)

    plt.tight_layout()
    
    file_extension = os.path.splitext(plot_path)[1].lower()
    if file_extension == '.pdf':
        plt.savefig(plot_path, format='pdf')
    else:
        if file_extension != '.png':
            print("Warning: Invalid plot file suffix. Your plot file is being saved in png format with the filename you requested.")
        plt.savefig(plot_path)
    
    plt.close()

