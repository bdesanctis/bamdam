# bamdam plotdamage

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

try: # optional library only needed for plotting 
    import matplotlib.pyplot as plt
    matplotlib_imported = True
except:
    matplotlib_imported = False

def run_plotdamage(args):
    if args.udg:
        make_damage_plot_udg(args.in_subs_list, args.in_subs, args.tax, args.outplot, args.ymax)
    else:
        make_damage_plot(args.in_subs_list, args.in_subs, args.tax, args.outplot, args.ymax)

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
        try:
            from_base, to_base, pos = mutation[:1], mutation[1:2], int(mutation[2:])
        except:
            print("Error in splitting mutations from subs file. If you used bamdam compute --udg, you need to set the --udg flag in plotdamage as well.")
            sys.exit(1)
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

# udg section
# udg parsing logic is sufficiently different that it was easier for me to write separate functions, instead of a ton of if statements in the above function.
# right now it only does c>t and not g>a. but that should be enough to analyze damage patterns. 

def calculate_damage_for_plot_udg(items):

    # so here we will have cpg > t, noncpg > t, total cpg, total noncpg, and other again. 

    ctp_cpg_5prime = {i: 0 for i in range(1, 16)}
    ctp_noncpg_5prime = {i: 0 for i in range(1, 16)}
    total_cpg_5prime = {i: 0 for i in range(1, 16)}
    total_noncpg_5prime = {i: 0 for i in range(1, 16)}

    ctp_cpg_3prime = {i: 0 for i in range(1, 16)}
    ctp_noncpg_3prime = {i: 0 for i in range(1, 16)}
    total_cpg_3prime = {i: 0 for i in range(1, 16)}
    total_noncpg_3prime = {i: 0 for i in range(1, 16)}

    other_5prime = {i: 0 for i in range(1, 16)}  # other mismatches at 5' positions (not including C>T)
    other_3prime = {i: 0 for i in range(1, 16)}  # other mismatches at 3' positions (not including C>T)
    total_nondeam_5prime = {i: 0 for i in range(1, 16)} #  denominator part : all non C>T matches or mismatches at each position
    total_nondeam_3prime = {i: 0 for i in range(1, 16)} 

    # skip CX mutations and only look at CpGX and nonCpGX 

    for item in items:
        mutation, proportion = item.split(":")
        count = float(proportion)

        # get the from_base right
        if mutation.startswith("nonCpG"):
            from_base = 'nonCpG'
            to_base = mutation[6]; pos = int(mutation[7:])
        if mutation.startswith("CpG"):
            from_base = 'CpG'
            to_base = mutation[3]; pos = int(mutation[4:])
        elif mutation.startswith("C"):
            # this is exactly the case you dont want : it starts with C but not CpG.
            continue
        if mutation.startswith("A") or mutation.startswith("G") or mutation.startswith("T"):
            from_base = mutation[0]; to_base = mutation[1]; pos = int(mutation[2:])

        # 5'
        if 1 <= pos <= 15:
            if from_base == 'nonCpG':
                if to_base == 'T':
                    ctp_cpg_5prime[pos] += count  # nonCpG>T at 5' positions 1 to 15
                total_cpg_5prime[pos] += count  
            if from_base == 'CpG':
                if to_base == 'T':
                    ctp_noncpg_5prime[pos] += count  # CpG>T at 5' positions 1 to 15
                total_noncpg_5prime[pos] += count  
            if from_base == to_base:
                total_nondeam_5prime[pos] += count
            if from_base != to_base and from_base != "nonCpG" and from_base != "CpG":
                other_5prime[pos] += count
                total_nondeam_5prime[pos] += count

        # 3' 
        elif -15 <= pos <= -1:
            abs_pos = abs(pos)
            if from_base == 'nonCpG':
                if to_base == 'T':
                    ctp_cpg_3prime[abs_pos] += count  # nonCpG>T at 5' positions 1 to 15
                total_cpg_3prime[abs_pos] += count  
            if from_base == 'CpG':
                if to_base == 'T':
                    ctp_noncpg_3prime[abs_pos] += count  # CpG>T at 5' positions 1 to 15
                total_noncpg_3prime[abs_pos] += count  
            if from_base == to_base:
                total_nondeam_3prime[abs_pos] += count
            if from_base != to_base and from_base != "nonCpG" and from_base != "CpG":
                other_3prime[abs_pos] += count
                total_nondeam_3prime[abs_pos] += count

    proportion_cpgt_5prime = {i: ctp_cpg_5prime[i] / total_cpg_5prime[i] if total_cpg_5prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_noncpgt_5prime = {i: ctp_noncpg_5prime[i] / total_noncpg_5prime[i] if total_noncpg_5prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_cpgt_3prime = {i: ctp_cpg_3prime[i] / total_cpg_3prime[i] if total_cpg_3prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_noncpgt_3prime = {i: ctp_noncpg_3prime[i] / total_noncpg_3prime[i] if total_noncpg_3prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_other_5prime = {i: other_5prime[i] / total_nondeam_5prime[i] if total_nondeam_5prime[i] > 0 else 0 for i in range(1, 16)}
    proportion_other_3prime = {i: other_3prime[i] / total_nondeam_3prime[i] if total_nondeam_3prime[i] > 0 else 0 for i in range(1, 16)}

    return proportion_cpgt_5prime, proportion_noncpgt_5prime, proportion_cpgt_3prime, proportion_noncpgt_3prime, proportion_other_5prime, proportion_other_3prime
    

def make_damage_plot_udg(in_subs_list, in_subs, tax, plotfile, ymax=0):

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

    values_all_files = {key: [] for key in ['Other', 'CpGT', 'nonCpGT']}  # store data for all files

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

        cpgt_5, noncpgt_5, cpgt_3, noncpgt_3, other_5prime, other_3prime = calculate_damage_for_plot_udg(data_items)

        values = {'Other': [0] * 30, 'CpGT': [0] * 30, 'nonCpGT': [0] * 30}
        for i in range(1, 16):
            # 5' 
            values['CpGT'][positions.index(i)] = cpgt_5[i]
            values['nonCpGT'][positions.index(i)] = noncpgt_5[i]
            values['Other'][positions.index(i)] = other_5prime[i]

            # 3' 
            values['CpGT'][positions.index(-i)] = cpgt_3[i]
            values['nonCpGT'][positions.index(-i)] = noncpgt_3[i]
            values['Other'][positions.index(-i)] = other_3prime[i]

        for key in values:
            values_all_files[key].append(values[key])

    if not any(values_all_files['Other']):
        print("Warning: No valid rows found for the given tax ID.")
        return

    if ymax == 0 or ymax == "0":
        max_y = min(1.0, max(max(item for sublist in values_all_files[key] for item in sublist) for key in ['Other', 'CpGT', 'nonCpGT']) * 1.2)
    else:
        max_y = float(ymax)

    # do the plotting
    plt.figure(figsize=(10, 5))  
    color_palette = {'Other': '#009E73', 'CpGT': '#F8766D', 'nonCpGT': '#56B4E9'}  # colorblind-friendly palette

    # 5' 
    ax1 = plt.subplot(1, 2, 1)
    handles1 = []  # store the line handles for ax1
    labels1 = ['Other', 'C to T (CpG context)', 'C to T (non CpG context)']  # labels for legend, one entry for each damage type
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






