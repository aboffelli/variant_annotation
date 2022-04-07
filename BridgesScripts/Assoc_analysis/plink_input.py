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

# Dictionary to save the ref base and alt base for all variants.
variant_dict = {}


# Prepare the VCF files that will be used.
list_of_files = glob.glob("./**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

var_num = 1
for file in list_of_files:
    # Get the phenotype, sample name and sex for the ped file.
    # The phenotype is Case or Control, for now.
    if file.split('/')[1] == "Controls":
        phenotype = "1"
    else:
        phenotype = "2"
    sample_name = re.search(r"vep_(\S+)\.raw", file).group(1)
    sex = '2'
    sample_dict[sample_name] = [sex, phenotype, {}]

    with open(file, 'r') as vcf_file:
        for vcfline in vcf_file:
            # Exclude header lines.
            if not vcfline.startswith("#"):
                split_line = vcfline.split('\t')

                # Exclude variants that did not pass the filters.
                filt = split_line[6]
                if filt != "PASS":
                    continue

                # Get the reference and alternative base.
                ref_base = split_line[3].upper().strip()
                alt_base = split_line[4].upper().strip()

                rs = re.search(r"CSQ=(?:(\w+)|\|)", vcfline).group(1)

                if rs == "rs555607708":
                    rs = "positive_control"
                    ref_base = ref_base[1]

                if len(ref_base + alt_base) != 2:
                    continue

                # Retrieve chr number and position for the map file.
                chrom = split_line[0][3:]
                if chrom == 'X':
                    chrom = '23'
                elif chrom == 'Y':
                    chrom = '24'
                position = split_line[1]

                # Retrieve the dbSNP id and assign a name for the variant if
                # there isn't one. In this case we're using varN.

                if not rs:
                    rs = 'var' + str(var_num)
                    var_num += 1

                # Save the map info in the map dictionary.
                if rs not in map_dict:
                    map_dict[rs] = (chrom, position)

                if rs not in variant_dict:
                    variant_dict[rs] = (ref_base, alt_base)

                # Allele fraction to differentiate between hetero and
                # homozygous.
                af = float(re.search(r"AF=(\d\.\d+)", vcfline).group(1))

                if af >= 0.8:
                    sample_dict[sample_name][2][rs] = [alt_base, alt_base]
                elif 0.2 <= af < 0.8:
                    sample_dict[sample_name][2][rs] = [ref_base, alt_base]


# For printing the map file, first prepare the lines into a list.
map_list = []
for rs in map_dict:
    chrom, pos = map_dict[rs]
    map_list.append(f"{chrom}\t{rs}\t{pos}")

# Then sort the list based on the chromosome number and position, and print to
# the file.
map_list = sorted(map_list, key=lambda x: (int(x.split('\t')[0]),
                                           int(x.split('\t')[2])))

with open("test.map", "w") as outfile:
    for item in map_list:
        print(item, file=outfile)


ped_dict = {}
for sample in sample_dict:
    ped_dict[sample] = sample_dict[sample][0:2]
    for var in map_list:
        var = var.split('\t')[1]
        if var in sample_dict[sample][2]:
            ped_dict[sample].extend(sample_dict[sample][2][var])
        else:
            ped_dict[sample].extend(variant_dict[var][0] * 2)

with open("test.ped", 'w') as outfile:
    for sample in ped_dict:
        line = '\t'.join(ped_dict[sample])
        print(f"{sample}\t{line}", file=outfile)


# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
