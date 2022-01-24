#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Synonymous variants' visualization table.

Description: Script to generate a table containing relevant information
to explore the synonymous variants in samples that do not have any pathogenic
variant reported.

Created on: 2022-01-19
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""

import os
import re
import time

start_time = time.time()


def synonymous_parser(vcfline):
    """
    Function that identifies a synonymous variant and retrieves the information
        needed for the table.

    :param vcfline: line from the vcf file
    :return: tuple containing the DbSNP ID, gene affected, the allele
        frequency, conservation value PhyloP, conservation value GERP,
        delta rscu, ese, ess, and encode information.
    """
    # Check if variant is synonymous.
    if 'synonymous_variant' in vcfline:
        # Isolate only the annotation info, using regex.
        csq = re.search(r'CSQ=(\S+?),(\S+?)(;ClinVar|\s)', vcfline)
        # The fixed csq contains information that are the same for all
        # transcripts
        fixed_csq = csq.group(1).split('|')
        known = fixed_csq[0].split('&')[0]
        af = fixed_csq[3]
        phylop = fixed_csq[6]
        gerp = fixed_csq[7]

        # The delta rscu and the encode are the same for all transcripts, so
        # they can be retrieved from the first transcript.
        first_transcript = csq.group(2).split('|')
        rscu = first_transcript[9]
        encode_info = first_transcript[8].split('&')
        encode = []
        # Since the encode information can contain more than one protein, in a
        # loop, get all proteins in a list and join them together.
        for prot in encode_info:
            encode.append(prot.split(':')[0])
        encode = ';'.join(encode)

        # It is possible that the same variant affects different genes, so we
        # loop through all transcripts and retrieve the gene name and join it
        # after.
        all_transcripts = csq.group(0).strip(';ClinVar\t').split(',')[1:]
        gene = [first_transcript[1]]

        # The ese and ess can be present or not depending on the transcript,
        # however, they will be the same for all transcripts. So. after we
        # retrieve a value we can stop looking for it.
        ese = ''
        ess = ''
        ese_ess = ''

        for transcript in all_transcripts:
            transcript = transcript.split('|')
            if transcript[1] not in gene:
                gene.append(transcript[1])
            if not ese_ess:
                # There is a function to evaluate the ese and ess.
                ese_ess = ese_ess_parser(transcript)

        # If any ese or ess was found separate them into different variables.
        if ese_ess:
            ese, ess = ese_ess
        # Join the genes.
        gene = ';'.join(gene)

        # Retrieve the allele frequency to calculate later, assign the value to
        # the dictionary.
        allele = vcfline.split('\t')[-1][0:3]
        if position not in swea_af:
            swea_af[position] = {}
        if allele == '0/1':  # only one allele
            allele = 0.5
        else:  # 1/1 two alleles
            allele = 1
        swea_af[position][sample_name] = allele
        return known, gene, af, phylop, gerp, rscu, ese, ess, encode


def ese_ess_parser(transcript):
    """
    Function to evaluate if there was a change in the ese and ess hexamers,
        that can be defined as:
            alter - if the hexamer or number of hexamers changed.
            create - if no hexamer was found in the reference but is found in
                the altered.
            remove - if one or more hexamer is found in the reference but no
                hexamer is found in the altered.

    :param transcript: list containing the information of the transcript line.
    :return: tuple containing the results for ese and ess (alter, remove, or
        create).
    """
    ese_ess = {'ese': '', 'ess': ''}
    # ESE is found in the index 10 and 11, ESS is found in the index 12 and 13.
    for i in range(10, 14, 2):
        if i == 10:
            es_type = 'ese'
        else:
            es_type = 'ess'

        # If reference and alteration differ.
        if transcript[i] != transcript[i + 1]:
            reference = transcript[i]
            alteration = transcript[i + 1]

            # Check if both have values
            if reference and alteration:
                # Create a list for all the hexamers
                ref_list = reference.split(';')
                alt_list = alteration.split(';')

                # If the sizes differ, there was an alteration
                if len(ref_list) != len(alt_list):
                    ese_ess[es_type] = 'alter'

                # If the sizes are the same, compare all hexamers to see if
                # there is a difference.
                else:
                    for hexamer in ref_list:
                        if hexamer not in alt_list:
                            ese_ess[es_type] = 'alter'

            # If the reference is empty and the alteration is not.
            elif not reference and alteration:
                ese_ess[es_type] = 'create'

            # If the alteration is empty and the reference is not.
            else:  # if reference and not alteration.
                ese_ess[es_type] = 'remove'

    # Return the values only if one of them is not empty. This makes possible
    # to check the other transcripts if both are empty.
    ese, ess = ese_ess.values()
    if ese or ess:
        return ese, ess


# Put all files in a list, removing anything that is not a vcf file.
files_directory = "SamplesWithoutPathogenic/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

# Initiate the three dictionaries.
synonymous_table = {}
swea_af = {}
swea_af_perc = {}

# File count that will be printed in the screen.
file_count = 1
for file in list_of_files:
    print(file_count)
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

                # Call the synonymous parser function
                result = synonymous_parser(vcf_line)

                # Store the result for the variant in the dictionary.
                if result:
                    if position not in synonymous_table:
                        synonymous_table[position] = result

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
        f'{swea_af_perc[position] / len(swea_af[position]):.3f} ' \
        f'({len(swea_af[position])} samples)'

# Save all the information in a text file tab delimited.
with open('synonymous_table.txt', 'w') as outfile:
    # Add the header
    print("Variant\tDbSNP ID\tGene\tAF SweGen\tAF SWEA\tPhyloP\t"
          "GERP\tRSCU\tESE\tESS\tRBP",
          file=outfile)

    for variant in synonymous_table:
        # Insert the SWEA allele frequency in the result list.
        result = list(synonymous_table[variant])
        result.insert(3, swea_af_perc[variant])

        # Change any empty space to NA
        for j in range(len(result)):
            if not result[j]:
                result[j] = 'NA'

        # Insert the line into the file.
        print("{}\t{}".format(variant, '\t'.join(result)),
              file=outfile)

# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))
