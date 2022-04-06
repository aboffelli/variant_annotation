#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description: Script that creates two input files for PLINK.
    One tab delimited MAP file containing the information for each variant, i.e.
    chromosome number, snp ID, and position.
    One tab delimited PED file containing the information for each sample, i.e.
    sample ID, sex, phenotype, and genotype.

Created on: 2022-04-06
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import time
import re
import glob

start_time = time.time()

# Dictionary to store the variants information for the map file.
map_dict = {}

# Dictionary with the sample information for the ped file.
sample_dict = {}

# Prepare the VCF files that will be used.
list_of_files = glob.glob("./**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

var_num = 1
for file in list_of_files:
    # Get the file type and sample name for the ped file.
    file_type = file.split('/')[1]
    sample_name = re.search(r"vep_(\S+)\.raw", file).group(1)

    with open(file, 'r') as vcf_file:
        for vcfline in vcf_file:
            # Exclude header lines.
            if not vcfline.startswith("#"):
                split_line = vcfline.split('\t')

                # Exclude variants that did not pass the filters.
                filt = split_line[6]
                if filt != "PASS":
                    continue

                # Retrieve chr number and position for the map file.
                chr = split_line[0][3:]
                position = split_line[1]

                # Retrieve the dbSNP id and assign a name for the variant if
                # there isn't one. In this case we're using varN.
                rs = re.search(r"CSQ=(?:(\w+)|\|)", vcfline).group(1)
                if not rs:
                    rs = 'var' + str(var_num)
                    var_num += 1

                # Save the map info in the map dictionary.
                if rs not in map_dict:
                    map_dict[rs] = (chr, position)

# For printing the map file, first prepare the lines into a list.
var_list = []
for rs in map_dict:
    chr, pos = map_dict[rs]
    var_list.append(f"{chr}\t{rs}\t{pos}")
# Then sort the list based on the chromosome number, and print to the file.
for item in sorted(var_list, key=lambda x: int(x.split('\t')[0])):
    print(item)


# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
