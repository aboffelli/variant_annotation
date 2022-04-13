#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-02-10
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import os
import re
import time

start_time = time.time()


def variant_parser(vcfline, bases):
    """
    Function that identifies a synonymous variant and retrieves the information
        needed for the table.

    :param bases: List of tuples containing position and nucleotides
    :param vcfline: line from the vcf file
    :return: tuple containing the DbSNP ID, gene affected, the allele
        frequency, conservation value PhyloP, conservation value GERP,
        delta rscu, ese, ess, and encode information.
    """
    # Check if variant is SNV.
    for base in bases:
        if len(base[1]) == 2:
            new_pos = base[0]
            # Isolate only the annotation info, using regex.
            csq = re.search(r'CSQ=(\S+?),(\S+?)(;ClinVar|\s)', vcfline)
            # The fixed csq contains information that are the same for all
            # transcripts
            fixed_csq = csq.group(1).split('|')
            known = fixed_csq[0].split('&')[0]
            af = fixed_csq[3]
            phylop = fixed_csq[6]
            gerp = fixed_csq[7]

            # It is possible that the same variant affects different genes, so we
            # loop through all transcripts and retrieve the gene name and
            # consequence and join it after.
            all_transcripts = csq.group(0).strip(';ClinVar\t').split(',')[1:]

            # We use sets to keep just the unique ones
            genes = set()
            consequence = set()

            for transcript in all_transcripts:
                transcript = transcript.split('|')
                genes.add(transcript[1])
                consequence.add(transcript[6])

            # Remove the genes that are not in the screening genes list by
            # intersecting the sets.
            gene = genes & seq_genes
            gene = ';'.join(gene)
            consequence = ';'.join(consequence)

            # Retrieve the pathogenicity relevance.
            pathogenicity = ''

            if "ClinVar" in vcf_line:
                clinvar_info = re.search(r'ClinVar=(\S*)', vcfline).group(1)
                pathogenicity = clinvar_info.split('|')[6]

            # Store all the results.
            result = (known, gene, af, phylop, gerp, consequence, pathogenicity)

            # Store the result for the variant in the dictionary.
            if new_pos not in variant_table:
                variant_table[new_pos] = result

            # Retrieve the allele frequency to calculate later, assign the value
            # to the dictionary.
            if new_pos not in swea_af:
                swea_af[new_pos] = {}
            if allele == '1/1':  # homozygous, two alleles
                allele_perc = 1
            else:  # heterozygous, one allele (0/1, 1/0 or 1/2)
                allele_perc = 0.5
            swea_af[new_pos][sample_name] = allele_perc


# Load the screened genes in a set.
seq_genes = set()
with open('../whole_gene_list.txt', 'r') as gene_list:
    for line in gene_list:
        seq_genes.add(line.strip())

# Put all files in a list, removing anything that is not a vcf file.
files_directory = "FixedESEESS/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

# Initiate the three dictionaries.
variant_table = {}
swea_af = {}
swea_af_perc = {}

# File count that will be printed in the screen.
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r')
    # Get the sample name from the file name using regex.
    sample_name = re.search(r'filtered_(\S+)_VEP', file).group(1)

    with open(files_directory + file, 'r') as vcf_file:
        for vcf_line in vcf_file:
            # Ignore the header lines.
            if not vcf_line.startswith('#'):
                # Ignore variants that did not pass the filters.
                filt = vcf_line.split('\t')[6]
                if filt != 'PASS':
                    continue

                # Retrieve the chromosome number and the position of the
                # variant and save it together as chr:position.
                chrom = vcf_line.split('\t')[0]
                position = vcf_line.split('\t')[1]
                position = f'{chrom}:{position}'

                # Account for variants with two alternative bases.
                allele = vcf_line.split('\t')[-1][:3]
                ref_base = vcf_line.split('\t')[3]
                alt_bases = vcf_line.split('\t')[4].split(',')

                # Save the position(s) in a list adding the alternative base in
                # the end.
                if '2' in allele:
                    nucleotides = [(position + alt_bases[0],
                                    ref_base + alt_bases[0]),
                                   (position + alt_bases[1],
                                    ref_base + alt_bases[1])
                                   ]
                else:
                    nucleotides = [(position + alt_bases[0],
                                    ref_base + alt_bases[0])]

                # Call the synonymous parser function
                variant_parser(vcf_line, nucleotides)

        # Raise the file count.
        file_count += 1

# Calculate the allele frequency of each variant.
for position in swea_af:
    swea_af_perc[position] = 0
    # Sum all the frequencies from the samples.
    for sample in swea_af[position]:
        swea_af_perc[position] += swea_af[position][sample]

    # Divide the sum by the number of samples and store it in the dictionary.
    swea_af_perc[position] = \
        f'{swea_af_perc[position] / 3403:.6f} ' \
        f'({len(swea_af[position])} samples)'

# Save all the information in a text file tab delimited.
with open('complete_table.txt', 'w') as outfile:
    # Add the header
    print("#Variant\tAlternative_base\tDbSNP_ID\tGene\tAF_SweGen\tAF_SWEA\t"
          "PhyloP\tGERP\tConsequence\tPathogenicity",
          file=outfile)

    for variant in sorted(variant_table):
        # Insert the SWEA allele frequency in the result list.
        result_line = list(variant_table[variant])
        result_line.insert(3, swea_af_perc[variant])

        # Change any empty space to NA
        for j in range(len(result_line)):
            if not result_line[j]:
                result_line[j] = 'NA'

        # Insert the line into the file.
        print("{}\t{}\t{}".format(variant[:-1], variant[-1],
                                  '\t'.join(result_line)),
              file=outfile)


# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
