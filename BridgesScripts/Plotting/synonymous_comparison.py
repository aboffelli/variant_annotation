#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Synonymous comparison and AF calculation

Description: Script to find the synonymous variants from the SWEA synonymous
    table in the BRIDGES samples. The script will also calculate the allele
    frequency for these variants.
    A file containing the unique gene symbols found in the BRIDGES vcfs will
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
list_of_files = glob.glob("raidset/FilteredClinVar/**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '/filtered' not in file:
        list_of_files.remove(file)

print("Loading family history cases...")
# Load the sample names that have family history.
with open("raidset/FilteredClinVar/Cases/family_history_cases.txt") as family:
    # Divide into all type of family history and only family history with less
    # than 50 years old.
    fam_hist_samples = set()
    fam_hist_50 = set()
    for line in family:
        # The first column of the file contains the sample name.
        fam_hist_samples.add(line.split('\t')[0])

        # Search for the less than 50 years breast or ovarian cancer code.
        if re.search(r"50[BO]C", line):
            fam_hist_50.add(line.split('\t')[0])
print("Done!")

# Retrieve the variant positions from the synonymous table and save it in a
# dictionary.
swea_synonymous = {}

# Set for storing the gene symbols annotated in BRIDGES.
# unique_genes = set()

print('\nLoading synonymous_table...')
with open('/home/ar7343bo-s/SWEA/new_synonymous_table.txt', 'r') as syn_table:
    # Remove the header from the loop.
    header = syn_table.readline()

    for line in syn_table:
        split_line = line.strip().split('\t')

        # Join the alternative base to the position to handle different
        # alternative bases in the same position.
        position = ''.join(split_line[0:2])

        # Save the whole line in the dictionary
        swea_synonymous[position] = split_line

print('Done!\n\nReading files...')

# Start two dictionaries separating Controls and Cases, one for storing the
# genotype value for each variant, and the other to calculate the percentage.
bridges_af = {'Controls': {}, 'Cases': {}}
bridges_af_perc = {'Controls': {}, 'Cases': {}}
bridges_fam_perc = {}
bridges_fam_perc_50 = {}

# File count that will be printed in the screen.
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    # Retrieve the type (Controls/Cases) from the file path.
    file_type = file.split('/')[2]

    # Isolate the sample name.
    sample = re.search(r'vep_(\S+).raw', file).group(1)
    with open(file, 'r') as vcf:
        for vcfline in vcf:

            # Ignore the comment lines.
            if not vcfline.startswith('#'):
                split_line = vcfline.split('\t')

                # If the variant did not pass the filters ignore it.
                filt = split_line[6]
                if filt != 'PASS':
                    continue

                # Get the position as chr:position and remove the chr since the
                # SWEA nomenclature does not use it.
                position = ':'.join(split_line[0:2]).lstrip('chr')

                # Get the alternative base and add it to the end of the
                # position.
                alt_base = vcfline.split('\t')[4]
                position += alt_base

                # # Isolate the csq section.
                # csq = re.search(r'CSQ=(\S+)', vcfline).group(1)
                #
                # # Split the transcripts.
                # csq = csq.split(',')
                #
                # # Loop through the transcripts ignoring the first one, that is
                # # the fixed csq.
                # for transcript in csq[1:]:
                #     # Get the gene symbol.
                #     gene = transcript.split('|')[1]
                #
                #     # If not empty, add to the set.
                #     if gene:
                #         unique_genes.add(gene)

                if position in swea_synonymous:
                    # Get the genotype information in the last column of the
                    # line.
                    allele = vcfline.split('\t')[-1][0:3]

                    # Add the position to the bridges dictionary, containing
                    # a dictionary if it is the first time we see this variant.
                    if position not in bridges_af[file_type]:
                        bridges_af[file_type][position] = {}

                    # Transform the genotype into numbers.
                    if allele == '1/1':  # homozygous - two alleles
                        allele = 1
                    else:  # heterozygous - one allele (0/1, 1/0 or 1/2)
                        allele = 0.5

                    # Store in the respective file type, position, and sample.
                    bridges_af[file_type][position][sample] = allele

    # Raise the file count.
    file_count += 1
print("\nDone!")
print("\nPreparing the output file...")
# Loop through the bridges dictionary to calculate the percentage of
for dict_type in bridges_af:

    # The number of samples is different depending on the type.
    if dict_type == "Cases":
        total_num = 60239

    else:  # Controls
        total_num = 53306

    for position in bridges_af[dict_type]:

        # Start the value for the current variant
        bridges_af_perc[dict_type][position] = 0
        bridges_fam_perc[position] = 0
        bridges_fam_perc_50[position] = 0
        fam_perc_count = 0
        fam_perc_50_count = 0

        # Sum all the frequencies from the samples.
        for sample in bridges_af[dict_type][position]:
            bridges_af_perc[dict_type][
                position] += bridges_af[dict_type][position][sample]

            if sample in fam_hist_samples:
                bridges_fam_perc[position
                                 ] += bridges_af[dict_type][position][sample]
                fam_perc_count += 1

            if sample in fam_hist_50:
                bridges_fam_perc_50[position
                ] += bridges_af[dict_type][position][sample]
                fam_perc_50_count += 1

        # Divide the sum by the total number of samples and store it in
        # the dictionary. Add the number of samples that had the variant.
        bridges_af_perc[dict_type][position] = (
            f"""{bridges_af_perc[dict_type][position] /
                 total_num:.6f} ({
            len(bridges_af[dict_type][position])} samples)""")

        bridges_fam_perc[position] = (
            f"""{bridges_fam_perc[position] /
                 len(fam_hist_samples):.6f} ({fam_perc_count} samples)""")

        bridges_fam_perc_50[position] = (
            f"""{bridges_fam_perc_50[position] /
                 len(fam_hist_samples):.6f} ({fam_perc_50_count} samples)""")

# Add the new columns to the header of the table and join it together again.
new_header = header.split('\t')
new_header.insert(7, 'AF_BRIDGES_Controls')
new_header.insert(8, 'AF_BRIDGES_Cases')
new_header.insert(9, "AF_BRIDGES_Cases_Family_Hist")
new_header.insert(10, "AF_BRIDGES_Cases_Family_Hist_50")
new_header = '\t'.join(new_header).strip()
print("Done!")
print("\nWriting the output file...")
# Write everything to a new file.
with open('bridges_synonymous_table.txt', 'w') as outfile:
    print(new_header, file=outfile)

    # Set the position to  insert the information based on the type.
    for dict_type in bridges_af_perc:
        if dict_type == 'Controls':
            insert_position = 7
        else:
            insert_position = 8

        # Loop through the swea dictionary.
        for position in swea_synonymous:
            # If the position is present in the bridges dictionary, insert the
            # the value in the right column.
            if position in bridges_af_perc[dict_type]:
                new_line = swea_synonymous[position]
                new_line.insert(insert_position,
                                bridges_af_perc[dict_type][position])

            # If the position is not in the bridges dictionary, add NA in the '
            # respective column.
            else:
                new_line = swea_synonymous[position]
                new_line.insert(insert_position, 'NA')

            # If we are on Cases already, we can join and print to the file,
            # otherwise we still need to include the samples numbers.
            if dict_type == 'Cases':
                if position in bridges_fam_perc:
                    new_line.insert(9, bridges_fam_perc[position])
                else:
                    new_line.insert(9, "NA")

                if position in bridges_fam_perc_50:
                    new_line.insert(10, bridges_fam_perc_50[position])
                else:
                    new_line.insert(10, bridges_fam_perc_50[position])

                print('\t'.join(new_line), file=outfile)

# # Print all the unique genes in a file.
# with open('bridges_gene_list.txt',
#           'w') as gene_list:
#     for gene_name in sorted(list(unique_genes)):
#         print(gene_name, file=gene_list)

# Print the run time.
print("Done!")
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
