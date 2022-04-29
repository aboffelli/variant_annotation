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

import os
import re
import time

start_time = time.time()


def is_pathogenic(vcfline, pos):
    """
    Function to isolate the ClinVar section of the vcf line, check if it
    contains a pathogenic report and add to the right dictionaries.

    :param str vcfline: string containing the whole vcf line.
    :param str pos: string with position of the variant in format chr:positionAlt.
    :return: None
    """

    # Check if there is a ClinVar section in the vcf line.
    if "ClinVar" in vcf_line:

        # Isolate the ClinVar section.
        clinvar_info = re.search(r'ClinVar=(\S*)', vcfline).group(1)

        # Retrieve the reference and alternative bases and the gene name.
        ref = clinvar_info.split('|')[-2]
        alt = clinvar_info.split('|')[-1]
        gene = clinvar_info.split('|')[4]

        # Join the position, bases and gene together for uniqueness.
        pos = f'{pos}_{ref}_{alt}_{gene}'

        # Check if there is a pathogenic or likely phatogenic classification.
        if 'Pathogenic' in clinvar_info or "Likely_pathogenic" in clinvar_info:

            # Separate the Glycogen storage disease from the cancer.
            if ('Glycogen_storage_disease_due_to_glucose-6-'
               'phosphatase_deficiency_type_IA') not in clinvar_info:

                # Add the pathogenic result to the counting dictionary.
                patho_count[pos] = ('Pathogenic/Likely pathogenic related to '
                                    'breast cancer')

                # Call the functions
                pathogenic_samples(sample_name, clinvar_info, pos)
                most_common_variant(pos)

            # Add the result as "not related to breast cancer" for the Glycogen
            else:
                patho_count[pos] = ('Pathogenic/Likely pathogenic '
                                    'not related to breast cancer')

        # If there is no pathogenic classification, add the result as "non
        # pathogenic".
        else:
            patho_count[pos] = 'Non pathogenic'

        # Call the function for type of clinical significance for all the
        # variants, independently if it is pathogenic or not.
        type_of_clinical_significance(clinvar_info, pos)

    # If there is no Clinvar info, set the result as "No Clinvar info"
    else:
        patho_count[pos] = 'No ClinVar info'

        # Add the result also in the synonymous dictionary.
        if 'synonymous_variant' in vcf_line:
            if pos not in synonymous_variants:
                synonymous_variants[pos] = 'No ClinVar info'


def type_of_clinical_significance(clinvar, pos):
    """
    Function to isolate the clinical significance from all variants and adds
    the variant and significance to the type dictionary. If the variant is
    synonymous also adds to the synonymous variants dictionary.

    :param str clinvar: string containing ClinVar section of the vcf line.
    :param str pos: string with position of the variant in format
    chr:position_ref_alt_gene.
    :return: None
    """

    clinvar_info = clinvar.split('|')
    type_dict[pos] = clinvar_info[6]

    if 'synonymous_variant' in vcf_line:
        if pos not in synonymous_variants:
            synonymous_variants[pos] = clinvar_info[6]


def pathogenic_samples(sample, clinvar, pos):
    """
    Function that isolates the clinical significance from the variants
    containing pathogenic significance and adds the variant to their respective
    sample in the pathogenic dictionary.

    :param str sample: string containing the sample name
    :param str clinvar: string containing ClinVar section of the vcf line.
    :param str pos: string with position of the variant in format
    chr:position_ref_alt_gene.
    :return: None
    """

    if sample not in samples_pathogenic:
        samples_pathogenic[sample] = {}

    if pos not in samples_pathogenic[sample]:
        samples_pathogenic[sample][pos] = clinvar.split('|')[6]


def most_common_variant(pos):
    """
    Function that counts the most common pathogenic variants. If the variant
    is not present in the dictionary, the key is created. Otherwise, the
    value is summed.

    :param str pos: string with position of the variant in format
    chr:position_ref_alt_gene.
    :return: None
    """

    if pos not in most_common:
        most_common[pos] = 1

    else:
        most_common[pos] += 1


# Initiate all the dictionaries.
patho_count = {}
type_dict = {}
samples_pathogenic = {}
most_common = {}
synonymous_variants = {}

# Prepare the files that will be used, removing anything that is not a vcf.
files_directory = "ClinVar/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

print("Reading files...")
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    # Get the sample name from the file name.
    sample_name = re.search(r'filtered_(\S+)_VEP', file).group(1)

    with open(files_directory + file, 'r') as vcf_file:
        for vcf_line in vcf_file:

            # Ignore header files.
            if not vcf_line.startswith('#'):

                # Check if the variant passed the filters.
                filt = vcf_line.split('\t')[6]
                if filt != 'PASS':
                    continue

                # Join the chromosome number and the position for
                # uniqueness.
                chrom = vcf_line.split('\t')[0]
                position = vcf_line.split('\t')[1]
                position = f'{chrom}:{position}'

                # Call the main function.
                is_pathogenic(vcf_line, position)

        # Raise the file count.
        file_count += 1
print("\nDone!")

print("\nWriting the output files...")
# Create an output directory if it does not exist.
out_dir = 'ClinVarTables/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Print the results of each dictionary to the files.

# Pathogenic variants count.
with open(out_dir + 'pathogenic_count.txt', 'w') as outfile:
    for key in patho_count:
        print(f'{key}\t{patho_count[key]}', file=outfile)

# Clinical type for all variants.
with open(out_dir + 'clinical_type.txt', 'w') as outfile:
    for key in type_dict:
        print(f'{key}\t{type_dict[key]}', file=outfile)

# Samples that contain pathogenic variants.
with open(out_dir + 'samples_pathogenic.txt', 'w') as outfile:
    for sample in samples_pathogenic:
        for key in samples_pathogenic[sample]:
            info = key.split(':')
            chrom = info[0]
            pos, ref, alt, gene = info[1].split('_')
            print(f'{sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t'
                  f'{samples_pathogenic[sample][key]}',
                  file=outfile)

# Most common pathogenic variants.
with open(out_dir + 'most_common_pathogenic_var.txt', 'w') as outfile:
    for key in most_common:
        print(f'{key}\t{most_common[key]}', file=outfile)


# Synonymous variants.
with open(out_dir+'synonymous_variants.txt', 'w') as outfile:
    for key in synonymous_variants:
        print(f'{key}\t{synonymous_variants[key]}', file=outfile)

print("Done!")

print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
