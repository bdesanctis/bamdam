
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


