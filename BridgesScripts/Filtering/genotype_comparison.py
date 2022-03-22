#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Genotype comparison

Description: File that uses the genotype csv file provided with the BRIDGES
    files to sanity check the vcf calls.

Created on: 2022-03-17
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import time
import re
import glob

start_time = time.time()

# Prepare the files that will be used.
test_files = glob.glob(r"/Users/student/Box/Notes/TestData/Bridges/"
                       r"FilteredClinVar/Cases/**/*.vcf", recursive=True)

# Start an empty dictionary to store the sample names from the test set that we
# are using.
variants_dict = {}
for x in test_files:
    # Add a new key for every sample, that contains a dictionary.
    variants_dict[re.search(r'vep_(\S+).raw', x).group(1)] = {}


# Load the genotype file
print("Loading the genotype file...")
with open('/Users/student/Box/Notes/TestData/Bridges/'
          'concept_699_kvist_bridges_genotypes.csv') as genotype_table:
    # Store the header for future reference of the variants.
    header = genotype_table.readline().strip().split(',')

    # For every line, check the sample name and only store the ones that are
    # in the test set. Use the header index to store the value for each variant.
    for line in genotype_table:
        line = line.strip().split(',')
        if line[0] in variants_dict:
            # Since we are skipping the index 0 (sample name) we start enumerate
            # in 1 instead of 0, so we match the right position in the header.
            for index, variant in enumerate(line[1:], 1):
                variants_dict[line[0]][header[index]] = variant
print("Done!")

# From here on there are different versions that were used to check the
# variants. This script is not totally necessary for the analysis.
check = {}
unique_vars = set()

for file in test_files:
    with open(file, 'r') as vcf:

        # Isolate the sample name from the file path.
        sample = re.search(r'vep_(\S+).raw', file).group(1)

        # Create the key for the sample and start some counting variables.
        check[sample] = []
        same = 0
        overlap = 0
        total = 0
        for line in vcf:
            # Ignore the header lines
            if not line.startswith("#"):
                line = line.strip()

                # Remove the CSQ section.
                no_csq = re.sub(r";CSQ=\S+", "", line)
                split_line = line.split('\t')

                total += 1

                # Get the variant in the same format of the genotype file
                # (chr_position_ref_alt)
                var = '_'.join(map(split_line.__getitem__, [0, 1, 3, 4]))

                # Isolate the allele fraction number.
                af = float(re.search(r";AF=([^A-Z;]+)", line).group(1))

                # Transform to the respective genotype number.
                if af >= 0.8:
                    vcf_value = 2
                elif 0.2 < af < 0.8:
                    vcf_value = 1
                else:
                    vcf_value = 0

                # Only add to the dictionary if the vcf value is 1 or 2, and
                # the value from the genotype file is 0.
                if var in variants_dict[sample]:
                    if vcf_value > 0 and variants_dict[sample][var] == '0':
                        check[sample].append(var)

for sample in check:
    # Make a copy of the dictionary and remove the current sample.
    new_dict = check.copy()
    new_dict.pop(sample)

    # Transform all variants from all samples into a set.
    var_list = set([x for sublist in new_dict.values() for x in sublist])

    # Start some counting variables to calculate the percentage.
    count = 0
    total = 0

    # If the variant is present in the set, meaning it was present in at least
    # one more sample, raise the count.
    for var in check[sample]:
        if var in var_list:
            count += 1
        # Total raises independently if it is in the set or not.
        total += 1

    # Print the percentage of variants that are present in at least two samples.
    print(f"{round(count/total*100, 4)}% are present in multiple samples")


        #         if var in variants_dict[sample]:
        #             overlap += 1
        #
        #             print(var, af, vcf_value,
        #                   variants_dict[sample][var], no_csq,
        #                   sep=',', file=outtable)
        #             if vcf_value == variants_dict[sample][var]:
        #                 same += 1
        #         else:
        #             print(var, af, vcf_value, 'NA', no_csq,
        #                   sep=',', file=outtable)
        #
        # print(f"""\n##{sample}: {same}/{overlap}/{total
        # } {round(same*100/total, 2)}%/{round(overlap*100/total, 2)
        # }% variants had the same number""", file=outtable)


# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
