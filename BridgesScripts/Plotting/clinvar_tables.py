#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description: Script that creates tables with ClinVar information to plot.

Created on: 2021-12-16
Author: Arthur Boffelli Castro
"""
# TODO: Add comments

import glob
import re
import os
import time

start_time = time.time()


def is_pathogenic(vcfline, pos):
    """

    :param vcfline:
    :param pos:
    :return:
    """
    if "ClinVar" in vcf_line:
        clinvar_info = re.search(r'ClinVar=(\S*)', vcfline).group(1)
        ref = clinvar_info.split('|')[-2]
        alt = clinvar_info.split('|')[-1]
        gene = clinvar_info.split('|')[4]
        pos = f'{pos}_{ref}_{alt}_{gene}'
        if 'Pathogenic' in clinvar_info or "Likely_pathogenic" in clinvar_info:
            if 'Glycogen_storage_disease_due_to_glucose-6-' \
               'phosphatase_deficiency_type_IA' not in clinvar_info:
                patho_count[fam_hist][pos] = 'Pathogenic/Likely pathogenic related to ' \
                                   'breast cancer'
                pathogenic_samples(sample_name, clinvar_info, pos)
                most_common_variant(pos)
            else:
                patho_count[fam_hist][pos] = 'Pathogenic/Likely pathogenic ' \
                                   'not related to breast cancer'

        else:
            patho_count[fam_hist][pos] = 'Non pathogenic'
        type_of_clinical_significance(clinvar_info, pos)

    else:
        patho_count[fam_hist][pos] = 'No ClinVar info'
        if 'synonymous_variant' in vcf_line:
            if pos not in synonymous_variants[fam_hist]:
                synonymous_variants[fam_hist][pos] = 'No ClinVar info'


def type_of_clinical_significance(clinvar, pos):
    """

    :param clinvar:
    :param pos:
    :return:
    """
    clinvar_info = clinvar.split('|')
    type_dict[fam_hist][pos] = clinvar_info[6]
    if 'synonymous_variant' in vcf_line:
        if pos not in synonymous_variants[fam_hist]:
            synonymous_variants[fam_hist][pos] = clinvar_info[6]


def pathogenic_samples(sample, clinvar, pos):
    """

    :param sample:
    :param clinvar:
    :param pos:
    :return:
    """
    if sample not in samples_pathogenic[fam_hist]:
        samples_pathogenic[fam_hist][sample] = {}
    if pos not in samples_pathogenic[fam_hist][sample]:
        samples_pathogenic[fam_hist][sample][pos] = clinvar.split('|')[6]


def most_common_variant(pos):
    """

    :param pos:
    :return:
    """
    if pos not in most_common[fam_hist]:
        most_common[fam_hist][pos] = 1
    else:
        most_common[fam_hist][pos] += 1


list_of_files = glob.glob("./**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

file_type = os.path.abspath(os.getcwd()).split('/')[-1].lower()

with open(f"family_bc_history_{file_type}.txt", 'r') as family_samples:
    family = set()
    for line in family_samples:
        sample = re.search(r"vep_(\S+)", line).group(1)
        family.add(sample.strip())

patho_count = {"Family_hist": {},
               "No_family_hist": {}}

type_dict = {"Family_hist": {},
             "No_family_hist": {}}

samples_pathogenic = {"Family_hist": {},
                      "No_family_hist": {}}

most_common = {"Family_hist": {},
               "No_family_hist": {}}

synonymous_variants = {"Family_hist": {},
                       "No_family_hist": {}}

file_count = 1
for file in list_of_files:
    print(file_count)
    sample_name = re.search(r'vep_(\S+)\.raw', file).group(1)
    if sample_name in family:
        fam_hist = 'Family_hist'
    else:
        fam_hist = "No_family_hist"

    with open(file, 'r') as vcf_file:
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

# Create an output directory if it doesn't exists.
out_dir = 'ClinVarTables/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Print the tables to the files.

# Pathogenic count
with open(out_dir + 'pathogenic_count.txt', 'w') as outfile:
    for history in patho_count:
        for key in patho_count[history]:
            print(f'{history}\t{key}\t{patho_count[history][key]}', file=outfile)

# Clinical type
with open(out_dir + 'clinical_type.txt', 'w') as outfile:
    for history in type_dict:
        for key in type_dict[history]:
            print(f'{history}\t{key}\t{type_dict[history][key]}', file=outfile)

# Samples pathogenic
with open(out_dir + 'samples_pathogenic.txt', 'w') as outfile:
    for history in samples_pathogenic:
        for sample in samples_pathogenic[history]:
            for key in samples_pathogenic[history][sample]:
                info = key.split(':')
                chrom = info[0]
                pos, ref, alt, gene = info[1].split('_')
                print(f'{history}\t{sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t'
                      f'{gene}\t{samples_pathogenic[history][sample][key]}',
                      file=outfile)

# Most common pathogenic
with open(out_dir + 'most_common_pathogenic_var.txt', 'w') as outfile:
    for history in most_common:
        for key in most_common[history]:
            print(f'{history}\t{key}\t{most_common[history][key]}', file=outfile)


# Synonymous variants
with open(out_dir+'synonymous_variants.txt', 'w') as outfile:
    for history in synonymous_variants:
        for key in synonymous_variants[history]:
            print(f'{history}\t{key}\t{synonymous_variants[history][key]}', file=outfile)

print('Run time: {:.2f} seconds'.format(time.time() - start_time))
