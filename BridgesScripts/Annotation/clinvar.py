#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: clinvar.py

Description: Script that annotates the clinical information of the variant
    based on the file variant_summary.txt obtained from ClinVar at
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.

Created on: 2021-12-10
Author: Arthur Boffelli Castro
"""
import glob
import re
import time

start_time = time.time()

# Prepare the VCF files that will be run
list_of_files = glob.glob("raidset/CustomAnnotation/Samples/**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '/custom' not in file:
        list_of_files.remove(file)

# Create a dictionary for the ClinVar file info.
clinvar_dict = {}
with open(r'/home/ar7343bo-s/variant_summary_GRCh37.txt',
          'r') as variant:
    print("Loading ClinVar info...")
    for line in variant:
        # # Get only the pathogenic/likely pathogenic variants.
        # if "pathogenic" in vcf_line.lower():
        # Use the chromosome number as keys and a nested dictionary.
        chrom = 'chr' + line.split('\t')[18]
        position = line.split('\t')[31]
        if chrom not in clinvar_dict:
            clinvar_dict[chrom] = {}
        # Add the position as keys in the nested dictionary, and line as a set.
        if position not in clinvar_dict[chrom]:
            clinvar_dict[chrom][position] = {line}
        else:
            clinvar_dict[chrom][position].add(line)
    print('Done.')

# New line to add in the vcf header.
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

file_count = 1
for file in list_of_files:
    # Print the file number to the screen.
    print(file_count)
    # Change the path and file names for the output file
    new_file = re.sub(r'CustomAnnotation(\S*/)(custom)',
                      r'ClinVar\1clinvar_\2', file).lstrip('raidset/')
    with open(file, 'r') as vcf, open(new_file, 'w') as outvcf:
        for vcf_line in vcf:
            # Variant lines
            if not vcf_line.startswith("#"):
                vcf_line = vcf_line.split('\t')
                # Get chromosome number, position and ref and alt bases.
                chrom = vcf_line[0]
                pos = vcf_line[1]
                ref = vcf_line[3]
                alt = vcf_line[4]

                # The search in the ClinVar file is based on the position, ref
                # and alt bases.
                search = f'{pos}\t{ref}\t{alt}'
                # Search only in the respective chromosome set.
                if pos in clinvar_dict[chrom]:
                    for clinvar_line in clinvar_dict[chrom][pos]:
                        if search in clinvar_line:
                            # Change the '|' in the clinvarline to ',', and the
                            # spaces to '_'. Join the ClinVar line with '|'.
                            clinvar_line = clinvar_line.replace('|', ',')
                            clinvar_line = clinvar_line.replace(' ', '_')
                            clinvar_line = '|'.join(clinvar_line.split('\t'))
                            # Add the clinvar line to the info section of the vcf
                            # line
                            vcf_line[7] = vcf_line[7] + ';ClinVar=' + \
                                clinvar_line.strip()
                            # print('\t'.join(line))
                            break
                outvcf.write('\t'.join(vcf_line))
            # Header lines
            else:
                if vcf_line.startswith('#CHROM'):
                    # Add the new line in the header, just before the last
                    # header line.
                    outvcf.write(new_info)
                outvcf.write(vcf_line)
        file_count += 1

print('Run time: {:.2f} seconds'.format(time.time() - start_time))
