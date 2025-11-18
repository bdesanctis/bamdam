
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

# makes krona plot for import into ktImportXML

def run_krona(args):
    make_krona_xml(args.in_tsv, args.in_tsv_list, args.out_xml, args.minreads, args.maxdamage, args.aggregate_to)

def make_krona_xml(in_tsv, in_tsv_files, out_xml, minreads, maxdamage, aggregateto):
    # create an xml text file to be loaded into kronatools for visualization. adds damage as a colour option.

    level_order = ["domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"] 
    if aggregateto not in level_order:
        print(f"Warning: Potentially unsafe value for aggregateto, {aggregateto}. Will try to continue anyway.")

    # did we get one or more files? check they exist then parse the input style 
    input_files = []
    if in_tsv:
        input_files = in_tsv
    # catch if we got a list of tsvs accidentally with the wrong flag
    if len(input_files) == 1:
        with open(input_files[0], 'r') as f:
            first_line = f.readline().strip()
            if '.tsv' in first_line:
                print(f"Error: It looks like {input_files[0]} is a list of TSV files, not a TSV file itself. If this is the case, use --in_tsv_list instead of --in_tsv.")
                exit(-1)
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
                print(header)
                print(f"Error: It looks like your input tsv file {file} was not generated from bamdam compute or bamdam combine, or perhaps you edited the header.")
                continue
            aggregate_marker = f":{aggregateto}"
            if not any(aggregate_marker in line for line in lines):
                print(f"Warning: File {file} does not contain any taxa at the '{aggregateto}' level. This file will be excluded. If you keep seeing this message, consider changing the --aggregate_to flag.")
                continue
            fields = lines[1].strip().split('\t')
            try:
                reads = int(fields[2])
                if reads < minreads:
                    print(f"Warning: File {file} only contains taxa with at most ({reads}) reads, which is below minreads ({minreads}). This file will be excluded.")
                    continue
            except (ValueError, IndexError):
                print(f"Error: File {file} has an unexpected format and will be excluded.")
                continue
        if file in valid_files:
            print(f"File {file} listed more than once in the input, but will only be included once.")
        if file not in valid_files: # avoid accidental duplicates
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

        taxpath_col = next((i for i, col in enumerate(header) if col.lower() == "taxpath"), -1)        
        if taxpath_col == -1:
            print(f"Error: There doesn't seem to be a tax path column in your header.")
            exit(-1)

        max_reads = int(lines[1].strip().split('\t')[2])  # tsvs are ordered; top line is max reads 
        sample_max_reads[sample_name] = max_reads

        howmanylinesfailed = 0
        for line in lines[1:]:  # skip the header
            fields = line.strip().split('\t')
            reads = float(fields[2]);
            if reads < minreads:
                break # you're done with this file
            dup = float(fields[3]); dust = float(fields[4]); damage = float(fields[5]); lengths = float(fields[7])
            taxpath = fields[-1].strip('"').strip("'") # 16 used to be -1; could replace this with taxpath_col anyway
            taxpathsplit = taxpath.split(";")
            fullnode = taxpathsplit[0].split(":")
            taxid = fullnode[0]
            taxname = fullnode[1]
            level = fullnode[2]
            just_added = False
            if f":{aggregateto};" not in taxpath:
                howmanylinesfailed +=1
                just_added = True
                if howmanylinesfailed == 100 and just_added:
                    print(f"Warning: More than 100 lines do not contain the string ':{aggregateto};'. These lines are being omitted. This may not be a problem if you intend to exclude these lines, but on the other hand you may be about to get an error because of it. Here is an example of one such tax path:")
                    print(taxpath)
                continue # skip this line, i don't need to deal with this rare weird case.
            if reads >= minreads:
                # first of all, is the node in the tree or not? 
                # if not, put it in. if yes, add this sample to the existing list for this node. 
                if taxid not in tree: # add a new node to the tree structure if it's not already there
                    tree[taxid] = {
                        'taxname': taxname,
                        'samples': {sample_name: {'reads': reads, 'duplicity': dup, 'dust': dust, 'damage': damage, 'length': lengths}},
                        'taxpath': taxpath,
                        'children': set(),
                        'node_filled' : True # this tracks which nodes we need to aggregate into later
                    }
                else: # otherwise just add this sample info in
                    if sample_name not in tree[taxid]['samples']:
                        tree[taxid]['samples'][sample_name] = {'reads': reads, 'duplicity': dup, 'dust': dust, 'damage': damage, 'length': lengths}
                        tree[taxid]['node_filled'] = True
                    else:  # this should not ever happen anymore
                        print(f"Warning: It looks like the file {sample_name} has two lines for the same tax id {taxid}.")
                        continue
                # second of all, are its ancestors (up to aggregateto) in the tree or not?
                # if not we can initialize them (then fill the unfilled ones by averaging over children values at the end)
                # similarly to the node, create ancestor nodes if needed, or add this sample to the ancestor node if it's not already there
                if level != aggregateto:
                    # if it was aggregateto, you just created that node, so you're good
                    # otherwise let's run through its ancestors to aggregateto level 
                    # all the nodes upto "upto" (from bamdam compute) will get filled in naturally
                    # and then everything between upto and aggregateto will get filled in afterwards but initialized here, 
                    #   with attached information about which samples they need to be present for in 'node_in_input'
                        current_taxid = taxid
                        for i in range(1,len(taxpathsplit)): # you did the 0 one above, because it carried information in this line
                            parent = taxpathsplit[i].split(":")
                            parent_taxid = parent[0]; parent_taxname = parent[1]; parent_level = parent[2]
                            if parent_taxid not in tree:
                                # initialize it
                                tree[parent_taxid] = {'taxname': parent_taxname, 'samples': {}, 'taxpath': ", ".join(taxpathsplit[1:]), 'children': set(), 'node_filled': False}
                                # will change node_filled if it appears later in any actual input file
                            # link the child node if it's not already linked
                            if current_taxid not in tree[parent_taxid]['children']: # do i need to check this? it's a set
                                tree[parent_taxid]['children'].add(current_taxid)
                            # then add empty information for this sample if it's not already there
                            # nb: it could already be there! rare but possible (and fine)
                            current_taxid = parent_taxid # update for the next ancestor up
                            if parent_level == aggregateto:
                                break

        # also get max damage of all the taxa for this sample  
        sample_max_damage[sample_name] = max(
            node['samples'][sample_name]['damage']
            for node in tree.values()
            if sample_name in node['samples'] and node['samples'][sample_name]['damage'] is not None
        )

    # aggregation time
    nodes_to_fill = [taxid for taxid in tree if not tree[taxid]['node_filled']]
    iterations = 0 # just in case as an infinite recursion catch
    while len(nodes_to_fill) > 0 and iterations < 100: # this will loop by however many tax levels are in between upto and aggregateto but it might have a lot of different clades? 
        for node in nodes_to_fill:
            # we can only fill it in this round if all of its children have been filled
            # have all the children been filled? this is true if none of the children have node_filled = False
            ready_to_fill = not any([tree[child]['node_filled'] is False for child in tree[node]['children']])
            # i guess you could add a check for no children, that would be a serious error lol
            if ready_to_fill:
                who_are_we_filling = tree[node]['taxname']
                # ok, fill it! i guess one sample at a time
                # but you don't have any sample info yet, so you have to go collect it from the children
                # like, which sample(s) does this node belong to? 
                # maybe i just pop through the children, aggregate to the samples per entry, then divide after
                total_reads, dupsum, dustsum, damsum, lengthsum = {}, {}, {}, {}, {}  
                for child in tree[node]['children']:
                    for sample_name in tree[child]['samples']:
                        sample_data = tree[child]['samples'][sample_name]
                        reads = sample_data['reads']
                        total_reads[sample_name] = total_reads.get(sample_name, 0) + reads
                        dupsum[sample_name] = dupsum.get(sample_name, 0) + sample_data['duplicity'] * reads
                        dustsum[sample_name] = dustsum.get(sample_name, 0) + sample_data['dust'] * reads
                        damsum[sample_name] = damsum.get(sample_name, 0) + sample_data['damage'] * reads
                        lengthsum[sample_name] = lengthsum.get(sample_name, 0) + sample_data['length'] * reads
                for sample_name in total_reads:
                    tree[node]['samples'][sample_name] = {'reads': total_reads[sample_name], 
                                                            'duplicity': dupsum[sample_name]  / total_reads[sample_name], 
                                                            'dust': dustsum[sample_name] / total_reads[sample_name], 
                                                            'damage': damsum[sample_name] / total_reads[sample_name], 
                                                            'length': lengthsum[sample_name] / total_reads[sample_name]}
                tree[node]['node_filled'] = True
        # reset
        iterations += 1 
        nodes_to_fill = [taxid for taxid in tree if not tree[taxid]['node_filled']]

    # a very handy debug point, i will leave it here in case anyone is ever changing this code and needs it (needs pprint package):
#    from pprint import pprint
#    treefile = "tree.tmp"
#    def write_tree_to_file(tree, filename):
#        with open(filename, "w") as f:
#            pprint(tree, stream=f)  # print the tree in a readable format
#        print(f"Tree structure written to file: {filename}")
#    write_tree_to_file(tree, treefile)

    # now build the xml manually by iterating through the tree 

    # now build a summary of everything if there is more than one thing; 
    # just iterate through existing nodes to do this then treat the summary like a bonus sample 
    all_children = {child for node in tree.values() for child in node['children']}
    tree_items = list(tree.items())
    if len(input_files) > 1:

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

    # calculate all the root reads for everything else
    for sample_name in sample_names:
        sample_root_reads = sum(
                node['samples'][sample_name]['reads']
                for taxid, node in tree_items
                if sample_name in node['samples'] and taxid not in all_children
            )
        sample_reads_at_root[sample_name] = sample_root_reads

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
    # override:
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
                node_str += f"<val>{round(value,3)}</val>"
            node_str += f"</{attr}>\n"

        node_str += f"{indent}\t<taxid><val>{node_id}</val></taxid>\n" # add the tax id too 

        for child_id in node_data["children"]:
            node_str += add_node_to_xml(child_id, tree, indent_level + 1)
        
        node_str += f"{indent}</node>\n"
        return node_str

    with open(out_xml, "w") as file:
        file.write('<krona>\n')
        file.write("\t<attributes magnitude=\"reads\">\n")  # Added defaultColorAttribute
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




