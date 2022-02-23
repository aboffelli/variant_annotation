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
#     if '/vep' not in file:
#         list_of_files.remove(file)

# Original files
list_of_files = glob.glob("Annotation/**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '/vep' not in file:
        list_of_files.remove(file)

# Retrieve the variant positions from the synonymous table
swea_synonymous = {}
unique_genes = set()

print('Loading synonymous_table...')
with open('/home/ar7343bo-s/SWEA/synonymous_table.txt', 'r') as syn_table:
    header = syn_table.readline()
    for line in syn_table:
        split_line = line.strip().split('\t')
        position = ''.join(split_line[0:2])
        swea_synonymous[position] = split_line
print('Done.\nReading files...')

bridges_af = {'Control': {}, 'Samples': {}}
bridges_af_perc = {'Control': {}, 'Samples': {}}

file_count = 1
for file in list_of_files:
    print(file_count)
    file_type = file.split('/')[1]
    sample_name = file.split('/')[-1].lstrip('encode_vep_')
    with open(file, 'r') as vcf:
        for vcfline in vcf:
            if not vcfline.startswith('#'):
                split_line = vcfline.split('\t')
                position = ':'.join(split_line[0:2]).lstrip('chr')
                alt_base = vcfline.split('\t')[4]
                position += alt_base

                csq = re.search(r'CSQ=(\S+)', vcfline).group(1)
                csq = csq.split(',')

                for transcript in csq:
                    gene = transcript.split('|')[1]
                    if gene:
                        unique_genes.add(gene)

                if position in swea_synonymous:
                    allele = vcfline.split('\t')[-1][0:3]
                    if position not in bridges_af[file_type]:
                        bridges_af[file_type][position] = {}
                    if allele == '1/1':  # homozygous - two alleles
                        allele = 1
                    else:  # heterozygous - one allele (0/1, 1/0 or 1/2)
                        allele = 0.5
                    bridges_af[file_type][position][sample_name] = allele

    file_count += 1
    
for dict_type in bridges_af:
    for position in bridges_af[dict_type]:
        bridges_af_perc[dict_type][position] = 0
    # Sum all the frequencies from the samples.
        for sample in bridges_af[dict_type][position]:
            bridges_af_perc[dict_type][position] += \
                bridges_af[dict_type][position][sample]

        # Divide the sum by the number of samples and store it in the dictionary.
        bridges_af_perc[dict_type][position] = \
            f"""{bridges_af_perc[dict_type][position] /
                 len(bridges_af[dict_type][position]):.3f} ({
            len(bridges_af[dict_type][position])} samples)"""


new_header = header.split('\t')
new_header.insert(7, 'AF_BRIDGES_Control')
new_header.insert(8, 'AF_BRIDGES_Samples')
new_header = '\t'.join(new_header).strip()

with open('bridges_synonymous_table.txt', 'w') as outfile:
    print(new_header, file=outfile)
    for dict_type in bridges_af_perc:
        if dict_type == 'Control':
            insert_position = 7
        else:
            insert_position = 8
        for position in swea_synonymous:
            if position in bridges_af_perc[dict_type]:
                new_line = swea_synonymous[position]
                new_line.insert(insert_position,
                                bridges_af_perc[dict_type][position])

            else:
                new_line = swea_synonymous[position]
                new_line.insert(insert_position, 'NA')

            if dict_type == 'Samples':
                print('\t'.join(new_line), file=outfile)

with open('bridges_gene_list.txt',
          'w') as gene_list:
    for gene_name in sorted(list(unique_genes)):
        print(gene_name, file=gene_list)

# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
