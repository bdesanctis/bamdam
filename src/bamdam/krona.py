
# bamdam krona

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


def krona(args):
    make_krona_xml(args.in_tsv, args.in_tsv_list, args.out_xml, args.minreads, args.maxdamage)


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

