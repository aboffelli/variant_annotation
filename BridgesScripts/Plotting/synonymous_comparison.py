#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Synonymous comparison and AF calculation

Description: Script to find the synonymous variants from the SWEA synonymous
    table in the BRIDGES samples. The script will also calculate the allele
    frequency for these variants.
    A file containing the unique gene symbols found in the BRIDGES files will
    also be created to find out which genes were used in the variant call.

Created on: 2022-02-16
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import re
import glob
import time

start_time = time.time()

# Test
# list_of_files = glob.glob("/Users/student/Box/Notes/TestData/Bridges/**/*.vcf",
#                           recursive=True)
# for file in list_of_files.copy():
#     if '/encode_vep' not in file:
#         list_of_files.remove(file)

# Original files
list_of_files = glob.glob("Annotation/Control/**/*.vcf",
                          recursive=True)
for file in list_of_files.copy():
    if '/vep' not in file:
        list_of_files.remove(file)

# Retrieve the variant positions from the synonymous table
swea_synonymous = {}
unique_genes = set()
with open('/home/ar7343bo-s/SWEA/synonymous_table.txt', 'r') as syn_table:
    header = syn_table.readline()
    for line in syn_table:
        split_line = line.strip().split('\t')
        position = split_line[0]
        swea_synonymous[position] = split_line

bridges_af = {}
bridges_af_perc = {}

file_count = 1
for file in list_of_files:
    print(file_count)

    sample_name = file.split('/')[-1].lstrip('encode_vep_')
    print(sample_name)
    with open(file, 'r') as vcf:
        for vcfline in vcf:
            if not vcfline.startswith('#'):
                split_line = vcfline.split('\t')
                position = ':'.join(split_line[0:2]).lstrip('chr')

                csq = re.search(r'CSQ=(\S+)', vcfline).group(1)
                csq = csq.split(',')

                for transcript in csq:
                    gene = transcript.split('|')[1]
                    if gene:
                        unique_genes.add(gene)

                if position in swea_synonymous:
                    allele = vcfline.split('\t')[-1][0:3]
                    if position not in bridges_af:
                        bridges_af[position] = {}
                    if allele == '0/1':  # only one allele
                        allele = 0.5
                    else:  # 1/1 two alleles
                        allele = 1
                    bridges_af[position][sample_name] = allele

    file_count += 1
    
for position in bridges_af:
    bridges_af_perc[position] = 0
    # Sum all the frequencies from the samples.
    for sample in bridges_af[position]:
        bridges_af_perc[position] += bridges_af[position][sample]

    # Divide the sum by the number of samples and store it in the dictionary.
    bridges_af_perc[position] = \
        f'{bridges_af_perc[position] / len(bridges_af[position]):.3f} ' \
        f'({len(bridges_af[position])} samples)'


new_header = header.split('\t')
new_header.insert(6, 'AF_BRIDGES')
new_header = '#' + '\t'.join(new_header).strip()

with open('new_synonymous_table.txt', 'w') as outfile:
    print(new_header, file=outfile)
    for position in swea_synonymous:
        if position in bridges_af_perc:
            new_line = swea_synonymous[position]
            new_line.insert(6, bridges_af_perc[position])

        else:
            new_line = swea_synonymous[position]
            new_line.insert(6, 'NA')
        print('\t'.join(new_line), file=outfile)

with open('bridges_gene_list.txt', 'w') as gene_list:
    for gene_name in sorted(list(unique_genes)):
        print(gene_name, file=gene_list)
