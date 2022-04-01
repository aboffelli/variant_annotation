#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Clinvar Tables.

Description: Script that parse the ClinVar annotation in the vcfs, and creates
    5 tables containing the parsed information.
    - Variants with pathogenic report.
    - Clinical relevance of all unique variants.
    - Samples that contain pathogenic variant.
    - Number of occurrences of each pathogenic variant.
    - Only synonymous variants and their clinical relevance.

Created on: 2021-12-16
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import glob
import re
import os
import time

start_time = time.time()


def is_pathogenic(vcfline, pos):
    """
    Function to isolate the ClinVar section of the vcf line, check if it
    contains a pathogenic report and add to the right dictionaries.

    :param vcfline: string containing the whole vcf line.
    :param pos: string with position of the variant in format chr:positionAlt.
    :return: None
    """

    # Check is there is a ClinVar section on the line
    if "ClinVar" in vcf_line:

        # Isolate the ClinVar section.
        clinvar_info = re.search(r'ClinVar=(\S*)', vcfline).group(1)

        # Retrieve the reference and alternative bases, and the gene affected.
        ref = clinvar_info.split('|')[-2]
        alt = clinvar_info.split('|')[-1]
        gene = clinvar_info.split('|')[4]
        # Update the position by adding the bases and gene.
        pos = f'{pos}_{ref}_{alt}_{gene}'

        # Check if the variant is reported as pathogenic.
        if 'Pathogenic' in clinvar_info or "Likely_pathogenic" in clinvar_info:

            if ('Glycogen_storage_disease_due_to_glucose-6-'
               'phosphatase_deficiency_type_IA') not in clinvar_info:
                # Store the pathogenic count for pathogenic variants that are
                # cancer related and call the functions to store in the
                # respective dictionaries.
                patho_count[fam_hist][pos] = ('Pathogenic/Likely pathogenic ' 
                                              'related to breast cancer')
                pathogenic_samples(sample_name, clinvar_info, pos)
                most_common_variant(pos)

            # The Glycogen report is not related to cancer, so we can handle it
            # as not related to breast cancer.
            else:
                patho_count[fam_hist][pos] = ('Pathogenic/Likely pathogenic '
                                              'not related to breast cancer')

        # If the variant does not contain a pathogenic variant, store it as
        # "Non pathogenic".
        else:
            patho_count[fam_hist][pos] = 'Non pathogenic'

        # Get the clinical significance for all variants that have the ClinVar
        # info.
        type_of_clinical_significance(clinvar_info, pos)

    # If the variant does not contain a ClinVar section, add to the dictionaries
    # as "No ClinVar info" in both pathogenic count and synonymous dictionary,
    # if the variant is synonymous.
    else:
        patho_count[fam_hist][pos] = 'No ClinVar info'
        if 'synonymous_variant' in vcf_line:
            if pos not in synonymous_variants[fam_hist]:
                synonymous_variants[fam_hist][pos] = 'No ClinVar info'


def type_of_clinical_significance(clinvar, pos):
    """
    Function to isolate the clinical significance from all variants and adds
    the variant and significance to the type dictionary. If the variant is
    synonymous also adds to the synonymous variants dictionary.

    :param clinvar: string containing ClinVar section of the vcf line.
    :param pos: string with position of the variant in format
    chr:position_ref_alt_gene.
    :return: None
    """
    clin_sig = clinvar.split('|')[6]
    type_dict[fam_hist][pos] = clin_sig
    if 'synonymous_variant' in vcf_line:
        if pos not in synonymous_variants[fam_hist]:
            synonymous_variants[fam_hist][pos] = clin_sig


def pathogenic_samples(sample, clinvar, pos):
    """
    Function that isolates the clinical significance from the variants
    containing pathogenic significance and adds the variant to their respective
    sample in the pathogenic dictionary.

    :param sample: string containing the sample name
    :param clinvar: string containing ClinVar section of the vcf line.
    :param pos: string with position of the variant in format
    chr:position_ref_alt_gene.
    :return: None
    """
    if sample not in samples_pathogenic[fam_hist]:
        samples_pathogenic[fam_hist][sample] = {}
    if pos not in samples_pathogenic[fam_hist][sample]:
        samples_pathogenic[fam_hist][sample][pos] = clinvar.split('|')[6]


def most_common_variant(pos):
    """
    Function that counts the most common pathogenic variants. If the variant
    is not present in the dictionary, the key is created. Otherwise, the
    value is summed.

    :param pos: string with position of the variant in format
    chr:position_ref_alt_gene.
    :return: None
    """

    if pos not in most_common[fam_hist]:
        most_common[fam_hist][pos] = 1
    else:
        most_common[fam_hist][pos] += 1


# Prepare the VCF files that will be used.
list_of_files = glob.glob("./**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

# Use the current path to decide if the samples are cases or controls.
file_type = os.path.abspath(os.getcwd()).split('/')[-1].lower()

# Load the family history file and save the sample names in  a set.
with open(f"family_history_{file_type}.txt", 'r') as family_samples:
    family = {}
    for line in family_samples:
        # Retrieve the sample name from the file name.
        sample, category = line.strip().split('\t')
        # Remove the new line in the end.
        family[sample] = category

# Initiate all the dictionaries with two divisions, one for samples with family
# history and one for samples without family history.
patho_count = {}
type_dict = {}
samples_pathogenic = {}
most_common = {}
synonymous_variants = {}


file_count = 1
for file in list_of_files:
    # Print the file number to the screen.
    print(file_count)

    # Get the sample name from the file name.
    sample_name = re.search(r'vep_(\S+)\.raw', file).group(1)

    # Check if the sample has family history or not and create the categories
    # in all dictionaries.
    if sample_name in family:
        fam_hist = family[sample_name]
    else:
        fam_hist = "No_family_hist"

    for d in [patho_count, type_dict, samples_pathogenic,
              most_common, synonymous_variants]:
        if fam_hist not in d:
            d[fam_hist] = {}

    # Open the file
    with open(file, 'r') as vcf_file:
        for vcf_line in vcf_file:

            # Ignore the header lines.
            if not vcf_line.startswith('#'):
                # Ignore the variants that did not pass the filters.
                filt = vcf_line.split('\t')[6]
                if filt == 'PASS':
                    # Join the chromosome number and the position for
                    # uniqueness.
                    position = ":".join(vcf_line.split('\t')[0:2])
                    # Call the main function.
                    is_pathogenic(vcf_line, position)

        # Raise the file count.
        file_count += 1

# Create an output directory if it does not exist.
out_dir = 'ClinVarTables/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Print the results of each dictionary to the files.

# Pathogenic variants count
with open(out_dir + 'pathogenic_count.txt', 'w') as outfile:
    # Family history and No family history.
    for history in patho_count:
        # Variant.
        for key in patho_count[history]:
            print(f'{history}\t{key}\t{patho_count[history][key]}',
                  file=outfile)

# Clinical type for all variants
with open(out_dir + 'clinical_type.txt', 'w') as outfile:
    # Family history and No family history.
    for history in type_dict:
        # Variant.
        for key in type_dict[history]:
            print(f'{history}\t{key}\t{type_dict[history][key]}', file=outfile)

# Samples that contain pathogenic variants
with open(out_dir + 'samples_pathogenic.txt', 'w') as outfile:
    # Family history and No family history.
    for history in samples_pathogenic:
        # Sample.
        for sample in samples_pathogenic[history]:
            # Variant.
            for key in samples_pathogenic[history][sample]:
                # Separate the variant information.
                info = key.split(':')
                chrom = info[0]
                pos, ref, alt, gene = info[1].split('_')

                # Print each in a different column.
                print(f'{history}\t{sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t'
                      f'{gene}\t{samples_pathogenic[history][sample][key]}',
                      file=outfile)

# Most common pathogenic variants
with open(out_dir + 'most_common_pathogenic_var.txt', 'w') as outfile:
    # Family history and No family history.
    for history in most_common:
        # Variant.
        for key in most_common[history]:
            print(f'{history}\t{key}\t{most_common[history][key]}',
                  file=outfile)


# Synonymous variants
with open(out_dir+'synonymous_variants.txt', 'w') as outfile:
    # Family history and No family history.
    for history in synonymous_variants:
        # Variant.
        for key in synonymous_variants[history]:
            print(f'{history}\t{key}\t{synonymous_variants[history][key]}',
                  file=outfile)

print('Run time: {:.2f} seconds'.format(time.time() - start_time))
