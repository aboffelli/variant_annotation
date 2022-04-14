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


# The base values that we are using for the filter are:
#   QUAL < 30
#   AF < 2.0
#   MQ < 60
#   NM > 2.0
#   AFxDP < 7.5


# Base values to the threshold, we exclude the NM from this dictionary, since it
# is the only one that is set to bigger than.
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

# Connection from the old filter names to the value names.
old_filters = {
    "q20": "QUAL",
    "Q10": "MQ",
    "f0.05": "AF",
}

file_count = 1
for file in list_of_files:
    # Print the file number to the screen.
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)
    # Change the path and file names for the output file.
    new_file = re.sub(r'ClinVar(\S*/)(clinvar)',
                      r'FilteredClinVar\1filtered_\2', file)

    with open(file, 'r') as vcf, open(new_file, 'w') as outvcf:
        for vcf_line in vcf:

            # For the header lines, we need to update the filter names.
            if vcf_line.startswith("#"):

                # Set an empty new line.
                new_line = ''
                # Get the filter header lines
                if '##FILTER' in vcf_line:
                    # Isolate the filter name by finding the ID=.
                    filt_name = re.search(r"ID=(\w+\.?\w+)", vcf_line).group(1)

                    # Check if the filter is one of the ones that we want to
                    # change
                    if filt_name in old_filters:
                        # Retrieve the value name
                        filt_type = old_filters[filt_name]

                        # Change the filter name for the new one.
                        new_line = re.sub(filt_name, filter_names[filt_type],
                                          vcf_line)

                        # Change the base value for the filter.
                        new_line = re.sub(r'\d\.?\d+(?=">)',
                                          str(base_values[filt_type]), new_line)

                    # Since NM is not in the dictionary, we need to account for
                    # it.
                    elif filt_name == "NM5.25":
                        # Change the line accordingly.
                        new_line = re.sub(r"NM5\.25", "NM2.0", vcf_line)
                        new_line = re.sub(r'>=\s\d\.\d+(?=,)', '> 2.0',
                                          new_line)

                # If new line is not empty, overwrite the vcf line by the new
                # line.
                if new_line:
                    vcf_line = new_line

                # Add the line to the file, changed or not.
                outvcf.write(vcf_line)

            # Variant lines
            else:
                split_line = vcf_line.split('\t')

                # Since we are not changing all the filters, first we have to
                # check if the variant was filtered by one of the filters that
                # we are changing.
                filt = split_line[6]
                filt = filt.split(';')

                # If none of the filters is either PASS or one of the filter
                # that we want to change, leave the line as it is skip the line.
                if all(i not in filter_set for i in filt):
                    outvcf.write(vcf_line)
                    continue

                # Retrieve all the values that we want from the vcf line.
                values = re.search(
                    r";DP=([^A-Z;]+)\S*;AF=([^A-Z;]+)\S*;QUAL=([^A-Z;]+)\S*"
                    r";MQ=([^A-Z;]+)\S*;NM=([^A-Z;]+)", vcf_line)

                # Transform all of them into floats.
                dp = float(values.group(1))
                af = float(values.group(2))
                qual = float(values.group(3))
                mq = float(values.group(4))
                nm = float(values.group(5))

                # Store the values that are in the dictionary in the right
                # order for the loop.
                variant_values = [af, qual, mq, af*dp]

                # Loop through the keys of the dictionary and check if the
                # values are less than the base value.
                for index, key in enumerate(base_values):
                    if variant_values[index] < base_values[key]:
                        # Add the filter name if the variant did not pass.
                        filt.append(filter_names[key])

                # Check if NM is higher than 2 and add the filter name if yes.
                if nm > 2.0:
                    filt.append(filter_names['NM'])

                # Remove the old filter names from the list of filters and join.
                filt = set(filt) - filter_set
                filt = ';'.join(list(filt))

                # If filt is empty at this point, it means that the variant
                # passed all the filters. So we can assign PASS to it.
                if not filt:
                    filt = 'PASS'

                # Update the filter in the line and join the line.
                split_line[6] = filt
                new_line = '\t'.join(split_line)

                # Write the new line in the new file.
                outvcf.write(new_line)

    # Raise the file count.
    file_count += 1


# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
