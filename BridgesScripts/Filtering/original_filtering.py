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
list_of_files = glob.glob("raidset/ClinVar/**/*.vcf", recursive=True)
for file in list_of_files.copy():
    # Remove any files that are not the last version of the vcfs.
    if '/clinvar' not in file:
        list_of_files.remove(file)

# Test data
# list_of_files = glob.glob(
#     "/Users/student/Box/Notes/TestData/Bridges/ClinVar/**/*.vcf",
#     recursive=True)

# Base values to the threshold
base_values = {
    "AF": 0.2,
    "QUAL": 30,
    "MQ": 60,
    "AFxDP": 7.5
}

# New names for each filter
filter_names = {
    "QUAL": "q30",
    "AF": "f0.2",
    "MQ": "Q60",
    "NM": "NM2.0",
    "AFxDP": "fd7.5"
}

# Old filters names that we want to change
filter_set = {"PASS", "q20", "Q10", "NM5.25", "f0.05"}

# Connection from the old filter names to the filters.
old_filters = {
    "q20": "QUAL",
    "Q10": "MQ",
    "f0.05": "AF",
}

file_count = 1
for file in list_of_files:
    # Print the file number to the screen.
    print(file_count)
    # Change the path and file names for the output file
    new_file = re.sub(r'ClinVar(\S*/)(clinvar)',
                      r'FilteredClinVar\1filtered_\2', file)

    with open(file, 'r') as vcf, open(new_file, 'w') as outvcf:
        for vcf_line in vcf:
            # Header lines
            if vcf_line.startswith("#"):
                new_line = ''
                if '##FILTER' in vcf_line:
                    filt_name = re.search(r"ID=(\w+\.?\w+)", vcf_line).group(1)
                    if filt_name in old_filters:
                        filt_type = old_filters[filt_name]
                        new_line = re.sub(filt_name, filter_names[filt_type],
                                          vcf_line)
                        new_line = re.sub(r'\d\.?\d+(?=">)',
                                          str(base_values[filt_type]), new_line)

                    elif filt_name == "NM5.25":
                        new_line = re.sub(r"NM5\.25", "NM2.0", vcf_line)
                        new_line = re.sub(r'>=\s\d\.\d+(?=,)', '> 2.0',
                                          new_line)

                if new_line:
                    vcf_line = new_line
                outvcf.write(vcf_line)

            # Variant lines
            else:
                split_line = vcf_line.split('\t')
                # We can keep the variants that did not pass the filters, since
                # the original values are higher than the ones being used.
                filt = split_line[6]
                filt = filt.split(';')
                # If none of the filters is either PASS or one of the filter
                # that we want to change, leave the line as it is skip the line.
                if all(i not in filter_set for i in filt):
                    outvcf.write(vcf_line)
                    continue
                values = re.search(
                    r";DP=([^A-Z;]+)\S*;AF=([^A-Z;]+)\S*;QUAL=([^A-Z;]+)\S*"
                    r";MQ=([^A-Z;]+)\S*;NM=([^A-Z;]+)", vcf_line)

                dp = float(values.group(1))
                af = float(values.group(2))
                qual = float(values.group(3))
                mq = float(values.group(4))
                nm = float(values.group(5))

                # Values that must be less than
                variant_values = [af, qual, mq, af*dp]

                for index, key in enumerate(base_values):
                    if variant_values[index] < base_values[key]:
                        filt.append(filter_names[key])
                if nm > 2.0:
                    filt.append(filter_names['NM'])
                filt = set(filt) - filter_set
                filt = ';'.join(list(filt))
                if not filt:
                    filt = 'PASS'
                split_line[6] = filt
                new_line = '\t'.join(split_line)
                outvcf.write(new_line)

    file_count += 1


# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
