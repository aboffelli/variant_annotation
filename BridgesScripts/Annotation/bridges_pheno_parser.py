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

# Create a dictionary to store the info for breast cancer dividing Cases and
# Controls
bc_dict = {
    'Cases': {
        'FirstBC': set(),
        'SecondBC': set(),
        'BC_50': set()
    },
    'Controls': {
        'FirstBC': set(),
        'SecondBC': set(),
        'BC_50': set()
    }
}

# Create a dictionary to store the info for ovarian cancer dividing Cases and
# Controls
oc_dict = {
    'Cases': {
        'FirstOC': set(),
        'SecondOC': set(),
        'OC_50': set()
    },
    "Controls": {
        'FirstOC': set(),
        'SecondOC': set(),
        'OC_50': set()
    }

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

            # Store all the information needed.
            sample = line[0]
            first_bc = line[22]
            second_bc = line[24]
            first_oc = line[26]
            second_oc = line[28]
            bc_50 = line[30]
            oc_50 = line[31]

            # Store the information in two different lists (breast cancer and
            # ovarian cancer) and transform the contents in integers.
            bc_list = [int(first_bc), int(second_bc), int(bc_50)]
            oc_list = [int(first_oc), int(second_oc), int(oc_50)]

            # Set the sample type to be used in the dictionaries.
            if file == pheno:
                sample_type = 'Cases'
            else:
                sample_type = 'Controls'

            # Using enumerate we get the index of the list, if it is True (not
            # 0) and not 888 (which means missing data) we add the sample name
            # to the set in the respective dictionary and category.
            for index, category in enumerate(bc_dict[sample_type]):
                if bc_list[index] and bc_list[index] != 888:
                    bc_dict[sample_type][category].add(sample)
            for index, category in enumerate(oc_dict[sample_type]):
                if oc_list[index] and oc_list[index] != 888:
                    oc_dict[sample_type][category].add(sample)


# To retrieve the total numbers to print to stdout, we loop through the sample
# types in one of the dictionaries.
for sample_type in bc_dict:
    # Start empty sets to store the combinations of sets.
    bc_intersect = set()
    oc_intersect = set()

    # Print with sample type is being used.
    print(sample_type + ':')

    # Loop through all the categories in the dictionary.
    for key in bc_dict[sample_type]:
        # Replace BC to OC for the Ovarian Cancer.
        oc_key = key.replace("B", "O")

        # Print the length of the set for each category.
        print(f'{key}: {len(bc_dict[sample_type][key])} patients')
        print(f'''{oc_key}: {len(oc_dict[sample_type][oc_key])} patients''')

        # Join the sets to keep only the unique patients.
        bc_intersect = bc_intersect.union(bc_dict[sample_type][key])
        oc_intersect = oc_intersect.union(oc_dict[sample_type][oc_key])

    # Print the total of unique patients for BC and OC.
    print(f"""\nTotal unique patients with family history breast cancer: {
    len(bc_intersect)}""")
    print(f"""Total unique patients with family history ovarian cancer: {
    len(oc_intersect)}\n""")

    # Use the joint sets to print the file path to the files.
    with open(f'''family_bc_history_{
            sample_type.lower()}.txt''', 'w') as bc_out, open(
            f'''family_oc_history_{sample_type.lower()}.txt''', 'w') as oc_out:
        for sample in bc_intersect:

            # Different path for intermediate panel.
            if sample.startswith("Inter"):
                if sample_type == 'Controls':
                    sample = f'''intermediate_panel/ctrl{
                    "_".join(sample.split("_")[0:2])
                    }/clinvar_custom_encode_vep_{sample}'''
                else:
                    sample = f'''intermediate_panel/{
                    "_".join(sample.split("_")[0:2])
                    }/clinvar_custom_encode_vep_{sample}'''

            # Path for bridges panel
            else:
                directory = f"""{sample.lower().split('_')[0]}{
                sample.lower().split('_')[2].lstrip('lib')}"""
                sample = f'''bridges_panel1_run1/{directory
                }/clinvar_custom_encode_vep_{sample}'''
            print(sample, file=bc_out)

        # Do the same for the OC set.
        for sample in oc_intersect:
            if sample.startswith("Inter"):
                if sample_type == 'Controls':
                    sample = f'''intermediate_panel/ctrl{
                    "_".join(sample.split("_")[0:2])
                    }/clinvar_custom_encode_vep_{sample}'''
                else:
                    sample = f'''intermediate_panel/{
                    "_".join(sample.split("_")[0:2])
                    }/clinvar_custom_encode_vep_{sample}'''
            else:
                directory = f"""{sample.lower().split('_')[0]}{
                sample.lower().split('_')[2].lstrip('lib')}"""
                sample = f'''bridges_panel1_run1/{directory
                }/clinvar_custom_encode_vep_{sample}'''
            print(sample, file=oc_out)


# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
