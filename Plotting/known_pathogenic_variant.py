#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2021-12-07 
Author: Arthur Boffelli Castro
"""

import re

anders_file = '/home/ar7343bo-s/SWEA_export_GATKv3.8_HC_till_HP.txt'
files_directory = 'CustomAnnotationComplete/'
with open(anders_file, 'r') as anders, open('known_pathogenic_var.txt',
                                            'w') as outfile:
    header = anders.readline()
    for line in anders:
        line = line.split('\t')

        file_name = 'custom_gene_name_encode_edited_fs_filtered_' +\
                    line[3].strip('.vcf') + '_VEP.vcf'

        for index in [10, 19, 28]:
            if line[index]:
                if int(line[index]) >= 4:
                    chrom, pos = line[index - 4].split('-')[:2]

                    with open(files_directory + file_name, 'r') as vcf:
                        for var in vcf:
                            target_var = re.search(chrom + r'\t' + pos, var)
                            if target_var:
                                print(f'{line[3]}\t{var.strip()}',
                                      file=outfile)
