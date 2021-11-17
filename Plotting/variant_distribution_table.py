#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2021-11-11
Author: Arthur Boffelli Castro
"""
import re
import os


def check_type(transcript_line, gene_name, existence, pos):
    """
    Function that identifies the consequence and assign the variant for the
        correct dictionary and gene.

    :param transcript_line: csq line from the vcf.
    :param gene_name: gene symbol
    :param existence: rs code
    :param pos: position of the variant
    """

    found = False

    for dict_key in consequence_dict:
        for key in dict_key.split('/'):
            if key in transcript_line:
                gene_consequence[pos] = (gene_name, dict_key, existence)
                found = True
        if found:
            break
    else:
        gene_consequence[pos] = (gene_name, 'other', existence)


def conservation_value(transcript_line, pos, existence):
    found = False
    conserv = csq.split(',')[0].split('|')[6:]
    for dict_key in consequence_dict:
        for key in dict_key.split('/'):
            if key in transcript_line:
                conserv_dict[pos] = (dict_key, conserv[0],
                                     conserv[1], existence)
                found = True
                break
        if found:
            break
    else:
        conserv_dict[pos] = ('other', conserv[0],
                             conserv[1], existence)

# Test
# genes_file = '/Users/student/Documents/whole_gene_list.txt'
genes_file = '/home/ar7343bo-s/whole_gene_list.txt'
targeted_genes = set()
with open(genes_file, 'r') as gene_list:
    for line in gene_list:
        targeted_genes.add(line.strip())

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

gene_consequence = {}

rscu_dict = {}
conserv_dict = {}

# Test
# files_directory = '/Users/student/Box/Notes/TestData/CustomAnnotation/'
files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)
# list_of_files = ['test.vcf']
for file in list_of_files:
    if '.vcf' in file:
        with open(files_directory + file, 'r') as vcf_file:
            for line in vcf_file:
                if not line.startswith("#"):
                    filt = line.split()[6]
                    if filt != 'PASS':
                        continue
                    position = line.split()[1]
                    csq = re.search(r'CSQ=(\S+)', line).group(1)
                    exist = csq.split('|')[0]
                    if 'synonymous_variant' in csq:
                        # Captures digits separated by a dot, the digits are
                        # after 'synonymous_variant|codon/codon|encode|' and
                        # before '|'. The '-' sign is optional.
                        rscu = re.search(
                            r'synonymous_variant\S*?\|'
                            r'\w+/\w+\|.*?\|(-?\d+\.\d+)\|',
                            csq).group(1)
                        rscu_dict[position] = (rscu, exist)

                    conservation_value(csq, position, exist)
                    for gene in targeted_genes:
                        if gene in csq:
                            check_type(csq, gene, exist, position)

variant_types_known = consequence_dict.copy()
variant_types_novel = consequence_dict.copy()
known_gene_dict = {}
novel_gene_dict = {}

for position in gene_consequence:
    gene = gene_consequence[position][0]
    consequence = gene_consequence[position][1]
    exist = gene_consequence[position][2]
    if exist:
        variant_types_known[consequence] += 1
        if gene not in known_gene_dict:
            known_gene_dict[gene] = consequence_dict.copy()
        known_gene_dict[gene][consequence] += 1
    else:
        variant_types_novel[consequence] += 1
        if gene not in novel_gene_dict:
            novel_gene_dict[gene] = consequence_dict.copy()
        novel_gene_dict[gene][consequence] += 1

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

# Test
# for i in variant_types_known:
#     print(f'{i}\t{variant_types_known[i]}')
#     print(f'{i}\t{variant_types_novel[i]}')

# for i in sorted(known_gene_dict):
#     for j in known_gene_dict[i]:
#         print(f'{i}\t{j}\t{known_gene_dict[i][j]}')
# for i in sorted(novel_gene_dict):
#     for j in novel_gene_dict[i]:
#         print(f'{i}\t{j}\t{novel_gene_dict[i][j]}')

# for i in rscu_dict:
#     if rscu_dict[i][2]:
#         print(f'{rscu_dict[i][0]}\tknown\t{rscu_dict[i][1]}')
#     else:
#         print(f'{rscu_dict[i][0]}\tnovel\t{rscu_dict[i][1]}')

# print("Gene\tExist\tConsequence\tPhyloP\tGERP")
# for i in conserv_dict:
#     if conserv_dict[i][4]:
#         print(f'{conserv_dict[i][0]}\tknown\t{conserv_dict[i][1]}\t'
#               f'{conserv_dict[i][2]}\t{conserv_dict[i][3]}')
#     else:
#         print(f'{conserv_dict[i][0]}\tnovel\t{conserv_dict[i][1]}\t'
#               f'{conserv_dict[i][2]}\t{conserv_dict[i][3]}')
