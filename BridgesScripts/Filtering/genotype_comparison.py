#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-03-17
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""

# TODO: Add comments.

import time
import re
import glob

start_time = time.time()

test_files = glob.glob(r"/Users/student/Box/Notes/TestData/Bridges/"
                       r"FilteredClinVar/Cases/**/*.vcf", recursive=True)
variants_dict = {}
for x in test_files:
    variants_dict[re.search(r'vep_(\S+).raw', x).group(1)] = {}


print("Loading the genotype file...")
with open('/Users/student/Box/Notes/TestData/Bridges/'
          'concept_699_kvist_bridges_genotypes.csv') as genotype_table:
    header = genotype_table.readline().strip().split(',')

    for line in genotype_table:
        line = line.strip().split(',')
        if line[0] in variants_dict:
            for index, variant in enumerate(line[1:], 1):
                variants_dict[line[0]][header[index]] = variant
print("Done!")

check = {}
unique_vars = set()
for file in test_files:
    with open(file, 'r') as vcf:

        sample = re.search(r'vep_(\S+).raw', file).group(1)
        check[sample] = []
        same = 0
        overlap = 0
        total = 0
        for line in vcf:
            if not line.startswith("#"):
                line = line.strip()
                no_csq = re.sub(r";CSQ=\S+", "", line)
                split_line = line.split('\t')

                total += 1
                var = '_'.join(map(split_line.__getitem__, [0, 1, 3, 4]))
                af = float(re.search(r";AF=([^A-Z;]+)", line).group(1))
                if af >= 0.8:
                    vcf_value = 2
                elif 0.2 < af < 0.8:
                    vcf_value = 1
                else:
                    vcf_value = 0

                if var in variants_dict[sample]:
                    if vcf_value > 0 and variants_dict[sample][var] == '0':
                        check[sample].append(var)

for sample in check:
    new_dict = check.copy()
    new_dict.pop(sample)
    var_list = set([x for sublist in new_dict.values() for x in sublist])

    count = 0
    total = 0
    for var in check[sample]:
        if var in var_list:
            count += 1
        total += 1
    print(f"{round(count/total*100, 4)}% are present in all samples")


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
