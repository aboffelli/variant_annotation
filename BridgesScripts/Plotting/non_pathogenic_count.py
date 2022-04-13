#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-01-06 
Author: Arthur Boffelli Castro
"""

import os
import re
import time

start_time = time.time()

# Search for possible causing variants in samples that do not have a
# pathogenic variants.
# Those are stop_gained (non-sense), splice_acceptor, splice_donor, or
# frameshift


def identify_variant(file_name):
    variant_dictionary = {}
    with open(files_directory + file_name, 'r') as vcf_file:
        for vcf_line in vcf_file:
            if not vcf_line.startswith("#"):
                filt = vcf_line.split('\t')[6]
                if filt != 'PASS':
                    continue
                for vtype in variant_types:
                    if vtype in vcf_line:
                        chrom = vcf_line.split('\t')[0]
                        position = vcf_line.split('\t')[1]
                        position = f'{chrom}:{position}'
                        variant_dictionary[position] = vtype
                        break
    return variant_dictionary


samples_dictionary = {}
variant_types = {'stop_gained', 'splice_acceptor',
                 'splice_donor', 'frameshift'}

# Check what the new directory is called.
files_directory = r"SamplesWithoutPathogenic/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r')
    sample_name = re.search(r'filtered_(\S+)_VEP', file).group(1)
    samples_dictionary[sample_name] = identify_variant(file)
    file_count += 1

with open('non_patho_other_variants.txt', 'w') as outfile:
    for sample in samples_dictionary:
        for variant in samples_dictionary[sample]:
            var_type = samples_dictionary[sample][variant]
            print(f'{sample}\t{variant}\t{var_type}', file=outfile)

print('Run time: {:.2f} seconds'.format(time.time() - start_time))
