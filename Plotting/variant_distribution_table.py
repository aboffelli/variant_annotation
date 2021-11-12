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
    global known_gene_consequence, novel_gene_consequence
    """
    Function to add the position and consequence for each gene.

    :param consequence: string containing the consequence of the variant
    :return:
    """
    found = False
    for dict_key in consequence_dict:
        for key in dict_key.split('/'):
            if key in transcript_line:
                if existence:
                    known_gene_consequence[gene_name][pos] = dict_key

                else:
                    novel_gene_consequence[gene_name][pos] = dict_key

                found = True
        if found:
            break
    else:
        print(transcript_line)
        if existence:
            known_gene_consequence[gene_name][pos] = 'other'

        else:
            novel_gene_consequence[gene_name][pos] = 'other'


def create_gene_dict(gene_name):
    """
    Function to populate the gene_dict.
    :param gene_name:
    :return: 
    """
    global known_gene_consequence, novel_gene_consequence
    if gene:
        if gene_name not in known_gene_consequence:
            known_gene_consequence[gene_name] = {}
            novel_gene_consequence[gene_name] = {}


genes_file = '/home/ar7343bo-s/whole_gene_list.txt'
# genes_file = r'C:\Users\Arthu\Documents\CustomAnnotation\whole_gene_list.txt'
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

known_gene_consequence = {}
novel_gene_consequence = {}

files_directory = "CustomAnnotation/"
# files_directory = r'C:\Users\Arthu\Documents\CustomAnnotation\\'
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)
# list_of_files = ['test.vcf']
for file in list_of_files.copy():
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
                    for gene in targeted_genes:
                        if gene in csq:
                            create_gene_dict(gene)
                            check_type(csq, gene, exist, position)

variant_types_known = consequence_dict.copy()
variant_types_novel = consequence_dict.copy()
known_gene_dict = {}
novel_gene_dict = {}

for gene in known_gene_consequence:
    for position in known_gene_consequence[gene]:
        consequence = known_gene_consequence[gene][position]
        variant_types_known[consequence] += 1
        if gene not in known_gene_dict:
            known_gene_dict[gene] = consequence_dict.copy()
        known_gene_dict[gene][consequence] += 1

for gene in novel_gene_consequence:
    if novel_gene_consequence[gene]:
        for position in novel_gene_consequence[gene]:
            if novel_gene_consequence[gene]:
                consequence = novel_gene_consequence[gene][position]
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
