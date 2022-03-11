#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Filter Fixer

Description: Script to fix the filtering options for BRIDGES vcf files, based on
    the values specified in the paper
    https://www.nejm.org/doi/10.1056/NEJMoa1913948.

Created on: 2022-03-11
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import time
import re
import glob

start_time = time.time()

# Prepare the vcf files that will be used.
# list_of_files = glob.glob("raidset/ClinVar/**/*.vcf", recursive=True)
# for file in list_of_files.copy():
#     # Remove any files that are not the last version of the vcfs.
#     if '/clinvar' not in file:
#         list_of_files.remove(file)

# Test data
list_of_files = glob.glob(
    "/Users/student/Box/Notes/TestData/Bridges/Clinvar/**/*.vcf",
    recursive=True)

# Base values to the threshold
base_values = {
    "QUAL": 30,
    "AF": 0.2,
    "MQ": 60,
    "NM": 2.0
}

filter_set = {"PASS", "q20", "Q10", "NM5.25", "f0.05"}

file_count = 1
for file in list_of_files:
    # Print the file number to the screen.
    print(file_count)
    # Change the path and file names for the output file
    new_file = re.sub(r'ClinVar(\S*/)(clinvar)',
                      r'FixedClinVar\1fixed_\2', file).lstrip('raidset/')
    with open(file, 'r') as vcf:  # , open(new_file, 'w') as outvcf:
        for vcf_line in vcf:
            # Header lines
            if vcf_line.startswith("#"):
                # TODO: Update the header.
                pass

            # Variant lines
            else:
                # We can keep the variants that did not pass the filters, since
                # the original values are higher than the ones being used.
                filt = vcf_line.split('\t')[6]

                # Check if any of the filters are the ones that we are changing.
                if all(i not in filter_set for i in filt.split(';')):
                    continue


# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
