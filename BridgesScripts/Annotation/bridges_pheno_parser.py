#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Phenotype file parser.

Description: Script to parse the phenotype tables containing information from
    family history of Breast Cancer and Ovarian Cancer. The script will output
    files with the path to the samples containing history of breast cancer
    and/or ovarian cancer and a summary of the total number of patients for
    each category in stdout.

Created on: 2022-02-25 
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import time

start_time = time.time()

# Path to the phenotype files provided by BRIDGES.
pheno_file = '/home/ar7343bo-s/BRIDGES/concept_699_kvist_pheno_v13.txt'
control_pheno_file = pheno_file.rstrip('.txt') + '_controls.txt'

# Create a dictionary to store the info for each sample dividing Cases and
# Controls.
samples_dict = {'Cases': {},
                'Controls': {}
                }


# Open the files.
with open(pheno_file, 'r') as pheno, open(control_pheno_file, 'r') as control:
    # Store file objects in a tuple, so we can loop through both of them.
    files = (pheno, control)

    # Read the header lines from both files
    header1 = pheno.readline()
    header2 = control.readline()

    # Loop through both files
    for file in files:
        for line in file:
            # Split the line into a list
            line = line.split('\t')

            # The control file has one column less than the case, so we fix that
            # by inserting an object at position 1.
            if file == control:
                line.insert(1, '')

            # Create a categories list.
            category_list = ["FirstBC", "SecondBC", "FirstOC", "SecondOC",
                             "50BC", "50OC"]

            # Store all the information needed.
            sample = line[0]
            sample_list = list(map(line.__getitem__, [22, 24, 26, 28, 30, 31]))

            # Transform the items in integers, so zeros will be False.
            sample_list = [int(x) for x in sample_list]

            # Set the sample type to be used in the dictionaries.
            if file == pheno:
                sample_type = 'Cases'
            else:
                sample_type = 'Controls'

            # Using enumerate we get the index of the list, if it is True (not
            # 0) and not 888 (which means missing data) we add the sample name
            # to the dictionary and the respective category.
            for index, category in enumerate(category_list):
                if sample_list[index] and sample_list[index] != 888:
                    if sample not in samples_dict[sample_type]:
                        samples_dict[sample_type][sample] = [category]
                    else:
                        samples_dict[sample_type][sample].append(category)


# Save the sample name and all the categories for that sample in a file.
for sample_type in samples_dict:
    # Create two files, one for cases and one for controls.
    with open(f'''family_history_{
    sample_type.lower()}.txt''', 'w') as family_out:
        for sample in samples_dict[sample_type]:
            print(f"{sample}\t{'/'.join(samples_dict[sample_type][sample])}",
                  file=family_out)

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
