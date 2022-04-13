#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description: Script that creates tables with ClinVar information to plot.

Created on: 2021-12-16
Author: Arthur Boffelli Castro
"""
# TODO: Add comments

import os
import re
import time

start_time = time.time()


def is_pathogenic(vcfline, pos):
    if "ClinVar" in vcf_line:
        clinvar_info = re.search(r'ClinVar=(\S*)', vcfline).group(1)
        ref = clinvar_info.split('|')[-2]
        alt = clinvar_info.split('|')[-1]
        gene = clinvar_info.split('|')[4]
        pos = f'{pos}_{ref}_{alt}_{gene}'
        if 'Pathogenic' in clinvar_info or "Likely_pathogenic" in clinvar_info:
            if 'Glycogen_storage_disease_due_to_glucose-6-' \
               'phosphatase_deficiency_type_IA' not in clinvar_info:
                patho_count[pos] = 'Pathogenic/Likely pathogenic related to ' \
                                   'breast cancer'
                pathogenic_samples(sample_name, clinvar_info, pos)
                most_common_variant(pos)
            else:
                patho_count[pos] = 'Pathogenic/Likely pathogenic ' \
                                   'not related to breast cancer'

        else:
            patho_count[pos] = 'Non pathogenic'
        type_of_clinical_significance(clinvar_info, pos)

    else:
        patho_count[pos] = 'No ClinVar info'
        if 'synonymous_variant' in vcf_line:
            if pos not in synonymous_variants:
                synonymous_variants[pos] = 'No ClinVar info'


def type_of_clinical_significance(clinvar, pos):
    clinvar_info = clinvar.split('|')
    type_dict[pos] = clinvar_info[6]
    if 'synonymous_variant' in vcf_line:
        if pos not in synonymous_variants:
            synonymous_variants[pos] = clinvar_info[6]


def pathogenic_samples(sample, clinvar, pos):
    if sample not in samples_pathogenic:
        samples_pathogenic[sample] = {}
    if pos not in samples_pathogenic[sample]:
        samples_pathogenic[sample][pos] = clinvar.split('|')[6]


def most_common_variant(pos):
    if pos not in most_common:
        most_common[pos] = 1
    else:
        most_common[pos] += 1


patho_count = {}
type_dict = {}
samples_pathogenic = {}
most_common = {}
synonymous_variants = {}

files_directory = "ClinVar/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

file_count = 1
for file in list_of_files:
    print(file_count, end='\r')
    sample_name = re.search(r'filtered_(\S+)_VEP', file).group(1)
    with open(files_directory + file, 'r') as vcf_file:
        for vcf_line in vcf_file:
            if not vcf_line.startswith('#'):
                filt = vcf_line.split('\t')[6]
                if filt != 'PASS':
                    continue
                chrom = vcf_line.split('\t')[0]
                position = vcf_line.split('\t')[1]
                position = f'{chrom}:{position}'
                is_pathogenic(vcf_line, position)

        file_count += 1

# Print the tables
out_dir = 'ClinVarTables/'
# Pathogenic count
with open(out_dir + 'pathogenic_count.txt', 'w') as outfile:
    for key in patho_count:
        print(f'{key}\t{patho_count[key]}', file=outfile)

# Clinical type
with open(out_dir + 'clinical_type.txt', 'w') as outfile:
    for key in type_dict:
        print(f'{key}\t{type_dict[key]}', file=outfile)

# Samples pathogenic
with open(out_dir + 'samples_pathogenic.txt', 'w') as outfile:
    for sample in samples_pathogenic:
        for key in samples_pathogenic[sample]:
            info = key.split(':')
            chrom = info[0]
            pos, ref, alt, gene = info[1].split('_')
            print(f'{sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t'
                  f'{samples_pathogenic[sample][key]}',
                  file=outfile)

# Most common pathogenic
with open(out_dir + 'most_common_pathogenic_var.txt', 'w') as outfile:
    for key in most_common:
        print(f'{key}\t{most_common[key]}', file=outfile)


# Synonymous variants
with open(out_dir+'synonymous_variants.txt', 'w') as outfile:
    for key in synonymous_variants:
        print(f'{key}\t{synonymous_variants[key]}', file=outfile)

print('Run time: {:.2f} seconds'.format(time.time() - start_time))
