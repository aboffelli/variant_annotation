#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-02-25 
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""

import time

start_time = time.time()

# Path to the phenotype files provided by BRIDGES.
pheno_file = 'concept_699_kvist_pheno_v13.txt'
control_pheno_file = pheno_file.rstrip('.txt') + '_controls.txt'

# Create a dictionary to store the info for the cases
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

# Create a dictionary to store the info for the controls
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

with open(pheno_file, 'r') as pheno, open(control_pheno_file, 'r') as control:
    files = (pheno, control)
    i = 0
    header1 = pheno.readline()

    header2 = control.readline().split()
    header2.insert(1, '')
    for file in files:
        for line in file:
            line = line.split('\t')
            if file == control:
                line.insert(1, '')
            sample = line[0]
            first_bc = line[22]
            second_bc = line[24]
            first_oc = line[26]
            second_oc = line[28]
            bc_50 = line[30]
            oc_50 = line[31]
            bc_list = [int(first_bc), int(second_bc), int(bc_50)]
            oc_list = [int(first_oc), int(second_oc), int(oc_50)]

            if i == 0:
                sample_type = 'Cases'
            else:
                sample_type = 'Controls'

            for index, category in enumerate(bc_dict[sample_type]):
                if bc_list[index] and bc_list[index] != 888:
                    bc_dict[sample_type][category].add(sample)
            for index, category in enumerate(oc_dict[sample_type]):
                if oc_list[index] and oc_list[index] != 888:
                    oc_dict[sample_type][category].add(sample)

        else:
            i += 1


for sample_type in bc_dict:
    bc_intersect = set()
    oc_intersect = set()
    print(sample_type + ':')
    for key in bc_dict[sample_type]:
        oc_key = key.replace("B", "O")
        print(f'{key}: {len(bc_dict[sample_type][key])} patients')
        print(f'''{oc_key}: {len(oc_dict[sample_type][oc_key])} patients''')
        bc_intersect = bc_intersect.union(bc_dict[sample_type][key])
        oc_intersect = oc_intersect.union(oc_dict[sample_type][oc_key])
    print(f"""\nTotal unique patients with family history breast cancer: {
    len(bc_intersect)}""")
    print(f"""Total unique patients with family history ovarian cancer: {
    len(oc_intersect)}\n""")

    with open(f'''family_bc_history_{
            sample_type.lower()}.txt''', 'w') as bc_out, open(
            f'''family_oc_history_{sample_type.lower()}.txt''', 'w') as oc_out:
        for sample in bc_intersect:
            if sample.startswith("Inter"):
                sample = f'''intermediate_panel/{
                "_".join(sample.split("_")[0:2])
                }/clinvar_custom_encode_vep_{sample}'''
            else:
                directory = f"""{sample.lower().split('_')[0]}{
                sample.lower().split('_')[2].lstrip('lib')}"""
                sample = f'''bridges_panel1_run1/{directory
                }/clinvar_custom_encode_vep_{sample}'''
            print(sample, file=bc_out)
        for sample in oc_intersect:
            if sample.startswith("Inter"):
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
