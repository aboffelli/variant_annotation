#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2021-12-10
Author: Arthur Boffelli Castro
"""
import os
import time

start_time = time.time()
files_directory = r"/Users/student/Box/Notes/TestData/CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)
original_variant = {}
with open(r'/Users/student/Box/Arthur/variant_summary_GRCh37.txt',
          'r') as variant:
    for line in variant:
        chrom = line.split('\t')[18]
        if chrom not in original_variant:
            original_variant[chrom] = {line}
        else:
            original_variant[chrom].add(line)

file_count = 1

new_info = '##INFO=<ID=ClinVar,Number=.,Type=String,Description=' \
           '"ClinVar annotation from python script clinvar.py. Format: ' \
           'AlleleID|Type|Name|GeneID|GeneSymbol|HGNC_ID|' \
           'ClinicalSignificance|ClinSigSimple|LastEvaluated|' \
           'RS#(dbSNP)|nsv/esv(dbVar)|RCVaccession|PhenotypeIDS|' \
           'PhenotypeList|Origin|OriginSimple|Assembly|ChromosomeAccession|' \
           'Chromosome|Start|Stop|ReferenceAllele|AlternateAllele|' \
           'Cytogenetic|ReviewStatus|NumberSubmitters|Guidelines|' \
           'TestedInGTR|OtherIDs|SubmitterCategories|VariationID|' \
           'PositionVCF|ReferenceAlleleVCF|AlternateAlleleVCF">\n'
for file in list_of_files:
    file_start = time.time()
    variant = original_variant.copy()
    print(file_count)
    with open(files_directory + file, 'r') as vcf, \
            open(files_directory + 'clinvar_'+file, 'w') as outvcf:
        for line in vcf:
            if not line.startswith("#"):
                line = line.split('\t')
                chrom = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]

                search = f'{pos}\t{ref}\t{alt}'
                for new_line in original_variant[chrom]:
                    if search in new_line:
                        new_line = new_line.replace('|', ',')
                        new_line = new_line.replace(' ', '_')
                        new_line = '|'.join(new_line.split('\t'))
                        line[7] = line[7] + ';ClinVar=' + new_line.strip()
                        # print('\t'.join(line))
                        break
                outvcf.write('\t'.join(line))
            else:
                if line.startswith('#CHROM'):
                    outvcf.write(new_info)
                outvcf.write(line)
        file_count += 1
        print('{} seconds file'.format(time.time()-file_start))
print('{} seconds total'.format(time.time()-start_time))
# 769.9284782409668 seconds
# 723.2530450820923 seconds
# 536.6957507133484 seconds
# TODO: Script too slow, maybe annotate just the pathogenic variants
