#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Variant distribution table.

Description: Script to generate tables containing the distribution of different
    information from the variants. Outputs 6 files:
        - known/novel variant distribution
        - known/novel variant distribution per gene
        - deltaRSCU table
        - conservation table

Created on: 2021-11-11
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import re
import os
import time

start_time = time.time()


def check_type(transcript_line, gene_name, existence, pos):
    """
    Function that identifies the consequence and assign the variant for the
    correct dictionary and gene.

    :param str transcript_line: csq line from the vcf.
    :param str gene_name: gene symbol.
    :param str existence: dbSNP id for known or empty for novel variants.
    :param str pos: position of the variant in the format chr:position
    :return None
    """

    # Set a flag for breaking the loop after finding a consequence.
    found = False

    # Loop through the consequence dictionary keys.
    for dict_key in consequence_dict:
        if not found:

            # Split the consequences that have deletion/insertion or
            # acceptor/donor.
            for key in dict_key.split('/'):

                # If the variant consequence is found, save in the dictionary.
                if key in transcript_line:
                    gene_consequence[pos] = (gene_name, dict_key, existence)

                    found = True

    # If the loop ends and no consequence was found, we set to other, since it
    # won't interest us.
    if not found:
        gene_consequence[pos] = (gene_name, 'other', existence)


def conservation_value(transcript_line, pos, existence):
    """
    Function to retrieve the conservation value together with the consequence.

    :param str transcript_line: csq line from the vcf.
    :param str pos: position of the variant in the format chr:position.
    :param str existence: dbSNP id for known or empty for novel variants.
    :return: None
    """

    # Set a flag for breaking the loop after finding the conservation values.
    found = False

    # Isolate the PhyloP and Gerp values from the fixed csq.
    conserv = csq.split(',')[0].split('|')[6:]

    # Loop through the consequence dictionary keys.
    for dict_key in consequence_dict:
        if not found:
            # Split the consequences that have deletion/insertion or
            # acceptor/donor.
            for key in dict_key.split('/'):
                if key in transcript_line:
                    conserv_dict[pos] = (dict_key, conserv[0],
                                         conserv[1], existence)
                    found = True
                    break

    # If the loop ends and no consequence was found, we set to other, since it
    # won't interest us.
    if not found:
        conserv_dict[pos] = ('other', conserv[0],
                             conserv[1], existence)


# Retrieve the targeted genes in the sequencing.
print("Loading the targeted genes...")
genes_file = '/home/ar7343bo-s/whole_gene_list.txt'
targeted_genes = set()
with open(genes_file, 'r') as gene_list:
    for line in gene_list:
        targeted_genes.add(line.strip())
print("Done!")

# Create a dictionary with the interesting consequences, the consequences are
# ordered in importance order, so in the loop, the first consequence found will
# be counted.
consequence_dict = {'stop_gained': 0,
                    'frameshift': 0,
                    'splice_acceptor/donor': 0,
                    'missense': 0,
                    'inframe_deletion/insertion': 0,
                    'splice_region': 0,
                    'synonymous_variant': 0,
                    '5_prime_UTR': 0,
                    '3_prime_UTR': 0,
                    'intron': 0,
                    'intergenic': 0,
                    'other': 0
                    }

# Dictionaries to store the information retrieved.
gene_consequence = {}
rscu_dict = {}
conserv_dict = {}

# Prepare the files that will be used, removing anything that is not a vcf.
files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)


print("Reading the files...")
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    with open(files_directory + file, 'r') as vcf_file:
        for line in vcf_file:

            # Ignore header lines.
            if not line.startswith("#"):

                # Check if the variant passed the filters.
                filt = line.split()[6]
                if filt != 'PASS':
                    continue

                # Join the chromosome number and the position for
                # uniqueness.
                chrom = line.split('\t')[0]
                position = line.split('\t')[1]
                position = f'{chrom}:{position}'

                # Isolate the csq from the line and retrieve the dbSNP id.
                csq = re.search(r'CSQ=(\S+)', line).group(1)
                exist = csq.split('|')[0]

                # Get the delta RSCU value if it is a synonymous variant.
                if 'synonymous_variant' in csq:
                    # Captures digits separated by a dot, the digits are
                    # after 'synonymous_variant|codon/codon|encode|' and
                    # before '|'. The '-' sign is optional.
                    rscu = re.search(
                        r'synonymous_variant\S*?\|'
                        r'\w+/\w+\|.*?\|(-?\d+\.\d+)\|',
                        csq).group(1)
                    rscu_dict[position] = (rscu, exist)

                # Call the conservation parser function.
                conservation_value(csq, position, exist)

                # Only get the consequence for variants that are in the
                # targeted genes.
                for gene in targeted_genes:
                    if gene in csq:
                        check_type(csq, gene, exist, position)

    # Raise the file count.
    file_count += 1
print("\nDone!")

# Create copies of the consequence dict to count the number of known and novel.
variant_types_known = consequence_dict.copy()
variant_types_novel = consequence_dict.copy()
known_gene_dict = {}
novel_gene_dict = {}

print("\nPreparing the output files...")
# Loop through the consequence dictionary.
for position in gene_consequence:
    # The gene name, consequence and dbSNP id are in a tuple.
    gene, consequence, exist = gene_consequence[position]

    # Count the known variants in their respective dictionaries.
    if exist:
        variant_types_known[consequence] += 1
        if gene not in known_gene_dict:
            known_gene_dict[gene] = consequence_dict.copy()
        known_gene_dict[gene][consequence] += 1

    # Count the novel variants in their respective dictionaries.
    else:
        variant_types_novel[consequence] += 1
        if gene not in novel_gene_dict:
            novel_gene_dict[gene] = consequence_dict.copy()
        novel_gene_dict[gene][consequence] += 1
print("Done!")

print("\nWriting the files...")
# Create the files
with open('known_variant_distribution.txt', 'w') as known, \
        open('novel_variant_distribution.txt', 'w') as novel:
    for i in variant_types_known:
        print(f'{i}\t{variant_types_known[i]}', file=known)
        print(f'{i}\t{variant_types_novel[i]}', file=novel)

with open('genes_known_variant_distribution.txt', 'w') as gene_known, \
        open('genes_novel_variant_distribution.txt', 'w') as gene_novel:
    for i in sorted(known_gene_dict):
        for j in known_gene_dict[i]:
            print(f'{i}\t{j}\t{known_gene_dict[i][j]}', file=gene_known)
    for i in sorted(novel_gene_dict):
        for j in novel_gene_dict[i]:
            print(f'{i}\t{j}\t{novel_gene_dict[i][j]}', file=gene_novel)

with open('rscu_table.txt', 'w') as rscu_table:
    for i in rscu_dict:
        if rscu_dict[i][1]:
            print(f'{rscu_dict[i][0]}\tknown', file=rscu_table)
        else:
            print(f'{rscu_dict[i][0]}\tnovel', file=rscu_table)

with open('conservation_table.txt', 'w') as conserv_table:
    print("Exist\tConsequence\tPhyloP\tGERP", file=conserv_table)
    for i in conserv_dict:
        if conserv_dict[i][3]:
            print(f'known\t{conserv_dict[i][0]}\t{conserv_dict[i][1]}\t'
                  f'{conserv_dict[i][2]}', file=conserv_table)
        else:
            print(f'novel\t{conserv_dict[i][0]}\t{conserv_dict[i][1]}\t'
                  f'{conserv_dict[i][2]}', file=conserv_table)

print("Done!")
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
