# bamdam combine

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


def run_combine(args):
    # parse real quick
    input_files = []
    if args.in_tsv:
        input_files = args.in_tsv
    elif args.in_tsv_list:
        with open(args.in_tsv_list, 'r') as file:
            input_files = [line.strip() for line in file if line.strip()]   
    parsed_data = {}
    headers = None
    
    for file_path in input_files:
        sample_name = file_path.split('/')[-1].replace('.tsv', '')  
        parsed_data[sample_name] = []
        with open(file_path, 'r') as file:
            for i, line in enumerate(file):
                fields = line.strip().split('\t')
                if i == 0: # skip the header line
                    if headers is None:
                        headers = fields
                    continue
                parsed_data[sample_name].append(fields)
    
    tsvs_to_matrix(parsed_data, headers, args.out_tsv, args.include, args.minreads)


def tsvs_to_matrix(parsed_data, headers, output_file, include=['none'], minreads=50):
    # for combine. pretty straightforward, read them all in and construct the matrix by pulling values out as needed
    # "meandamage" output (both for +1 and -1) is actually weighted by number of reads
    # otherwise everything else is a straightforward unweighted average 

    # by default this will include the tax name, total reads, damage, duplicity and dust. 

    col_indices = {}
    required_cols = ['TaxNodeID','TaxName','TotalReads','Duplicity','MeanDust','Damage+1','taxpath']
    available_cols = ['Damage-1', 'MeanLength','ANI', 'AvgReadGC', 'AvgRefGC', 'UniqueKmers', 'UniqKmersPerRead', 'TotalAlignments', 'UnaggregatedReads']
    potential_udg_cols = ['Damage+1_CpG', 'Damage+1_nonCpG']

    # make it so potential_udg_cols are options for input, and included in "all" if they're present in the file,
    # but you don't print a warning if they're not there unless the user explicitly asks for them. 

    for col in required_cols:
        try:
            col_indices[col] = headers.index(col)
        except ValueError:
            print(f"Error: Required column '{col}' not found in input files. Are all your input files bamdam compute output with their original headers?")
            sys.exit(1)
    
    for col in available_cols + potential_udg_cols:
        if col in headers:
            col_indices[col] = headers.index(col)

    if include == ['none']:
        include_cols = ['Duplicity', 'MeanDust', 'Damage+1']  # default set from required_cols
    elif 'all' in include:
        include_cols = ['Duplicity', 'MeanDust', 'Damage+1'] + [col for col in available_cols if col in col_indices]
        for udg_col in potential_udg_cols:
            if udg_col in headers:
                include_cols.append(udg_col)
        # search for potential_udg_cols in the actual file(s), if they're there, add them to include_cols. 
    else:
        include_cols = ['Duplicity', 'MeanDust', 'Damage+1']  # start with defaults
        for col in include:
            if col in col_indices and col not in include_cols:
                include_cols.append(col)
            elif col not in col_indices:
                print(f"Error: Requested column '{col}' not found in input files.")
                sys.exit(1)

    # make sure everything the user wants is available in the input files
    for col in include_cols:
        try:
            col_indices[col] = headers.index(col)
        except ValueError:
            print(f"Warning: Column '{col}' not found in input files. Are all your input files bamdam compute output with their original headers?")


    tax_data = {}
    for sample_name, records in parsed_data.items():
        if sample_name.endswith('.tsv'):
            sample_name = sample_name.replace('.tsv', '')

        print(f"Processing sample: {sample_name} with {len(records)} records.")
        for record in records:
            taxpath = record[col_indices['taxpath']]
            tax_node_id = record[col_indices['TaxNodeID']]
            tax_name = record[col_indices['TaxName']]
            
            if tax_node_id not in tax_data:
                tax_data[tax_node_id] = {
                    "TaxNodeID": tax_node_id,
                    "TaxName": tax_name,
                    "TaxPath": taxpath,
                    "samples": {},
                    "TotalReads": 0,
                    "WeightSum": 0,
                }
                # Initialize weighted sums for included columns
                for col in include_cols:
                    tax_data[tax_node_id][f"Weighted{col}"] = 0

            try:
                total_reads = int(record[col_indices['TotalReads']])
                sample_data = {"TotalReads": total_reads}
                
                for col in include_cols:
                    try:
                        sample_data[col] = float(record[col_indices[col]])
                    except (ValueError, IndexError):
                        sample_data[col] = None
                        
            except (ValueError, IndexError):
                print("Error: At least one of the input files does not appear to be an output file from bamdam compute.")
                sys.exit(1)
            
            tax_data[tax_node_id]['samples'][sample_name] = sample_data
            tax_data[tax_node_id]['TotalReads'] += total_reads
            tax_data[tax_node_id]["WeightSum"] += total_reads
            
            for col in include_cols:
                if sample_data[col] is not None:
                    tax_data[tax_node_id][f"Weighted{col}"] += sample_data[col] * total_reads

    for tax_id, data in tax_data.items():
        for col in include_cols:
            if data['WeightSum'] > 0:
                data[f'Avg{col}'] = data[f'Weighted{col}'] / data['WeightSum']
            else:
                data[f'Avg{col}'] = 'NA'

    tax_data = {tax_id: data for tax_id, data in tax_data.items() if data['TotalReads'] >= minreads}
    sorted_tax_data = sorted(tax_data.items(), key=lambda x: x[1]['TotalReads'], reverse=True)

    with open(output_file, "w") as outfile:
        header = ['TaxNodeID', 'TaxName', 'TotalReads']
        
        for col in include_cols:
            header.append(f'Avg{col}')
        
        for sample_name in sorted(parsed_data.keys()):
            if sample_name.endswith('.tsv'):
                sample_name = sample_name.replace(".tsv", "")
            header.append(f"{sample_name}_TotalReads")
            for col in include_cols:
                header.append(f"{sample_name}_{col}")
        
        header.append("taxpath")
        outfile.write("\t".join(header) + "\n")

        for tax_id, data in sorted_tax_data:
            row = [data['TaxNodeID'], data['TaxName'], str(data['TotalReads'])]
            
            for col in include_cols:
                avg_val = data[f'Avg{col}']
                if isinstance(avg_val, float):
                    row.append(str(round(avg_val, 4)))
                else:
                    row.append('NA')

            for sample_name in sorted(parsed_data.keys()):
                if sample_name.endswith('.tsv'):
                    sample_name = sample_name.replace(".tsv", "")
                sample_data = data['samples'].get(sample_name, {})
                
                row.append(str(sample_data.get("TotalReads", 0)))
                
                for col in include_cols:
                    val = sample_data.get(col)
                    if val is not None:
                        row.append(str(val))
                    else:
                        row.append('NA')
            
            row.append(data["TaxPath"])
            outfile.write('\t'.join(row) + '\n')

