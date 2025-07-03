
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
from functools import lru_cache

try: # optional library for progress bars 
    from tqdm import tqdm
    tqdm_imported = True
except ImportError:
    tqdm_imported = False

try: # optional library only needed for plotting 
    import matplotlib.pyplot as plt
    matplotlib_imported = True
except:
    matplotlib_imported = False

# local
import utils
import alignment_utils

def run_compute(args):
    lca_file_type = utils.find_lca_type(args.in_lca)
    if(lca_file_type == "metadmg"):
        print("Error: It looks like you're trying to run bamdam compute with a metaDMG-style lca file. Please use an ngsLCA-style lca file. Bamdam shrink can automatically convert it for you.")
        sys.exit()
    nodedata, pmds_in_bam = gather_subs_and_kmers(args.in_bam, args.in_lca, kn=args.k, upto=args.upto, stranded=args.stranded)
    parse_and_write_node_data(nodedata, args.out_tsv, args.out_subs, args.stranded, pmds_in_bam) 
    if args.plotdupdust:
        plotdupdust(args.out_tsv, args.plotdupdust) 


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
                rep_kmers.append(utils.get_rep_kmer(kmer))
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
        totallcalines = utils.line_count(lcafile_path)
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
            subs, matches, refseq = alignment_utils.mismatch_table(seq,cigar,md,flagsum) 
            oldcigar = cigar; oldmd = md; oldflagsum = flagsum

        allsubs = subs + matches
        for sub in allsubs:
            key = tuple(sub)
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
        parts = [str(key[0]), str(key[1]), str(key[2])]
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
        next(f)  # skip header
        for line in f:
            cols = line.split('\t')
            if len(cols) < 6:
                continue
            taxpath = cols[14]
            reads = int(cols[2])
            duplicity = float(cols[3])
            dust = float(cols[4])
            damage = float(cols[5])

            if 'genus' in taxpath and 'species' not in taxpath and 'subgenus' not in taxpath and reads > 50:
                genera.append(taxpath)
                duplicities.append(duplicity)
                dusts.append(dust)
                damages.append(damage)
    
    if not duplicities or not dusts:
        print("No data matching criteria.")
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
    plt.savefig(plot_path)
    plt.close()



