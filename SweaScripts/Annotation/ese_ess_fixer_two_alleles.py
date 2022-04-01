#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-02-21
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""
# TODO: Add comments

import os
import re
import time

start_time = time.time()


def two_alleles(vcf_line, alt_bases_list):

    csq = re.search(r'CSQ=(\S+?),(\S+?)(;ClinVar|\s)', vcf_line)
    transcripts = csq.group(2).split(',')
    for index, transcript in enumerate(transcripts):
        transcript = transcript.split('|')
        exon = transcript[4]
        strand = transcript[3]
        if exon:
            ref_motifs_ese = set()
            ref_motifs_ess = set()
            alt_motifs_ese = set()
            alt_motifs_ess = set()
            flanking_seq = re.search(r"FSEQ=(\S[A-Z]+\[.+\/.+\][A-Z]+)",
                                     vcfline).group(1)
            # Retrieve the flanking sequences
            seq = re.search(r'(\w{5})\[(\S+)\/\S+\](\w{5})',
                            flanking_seq)
            ref_seq = ''.join(seq.group(1, 2, 3))
            for base in alt_bases_list:
                alt_seq = f'{seq.group(1)}{base}{seq.group(3)}'

                # Check if the variant is located in the reverse strand
                # and get the reverse complement.
                if strand == '-1':
                    ref_seq = reverse_complement(ref_seq)
                    alt_seq = reverse_complement(alt_seq)

                # Check all possible hexamers in the reference sequence
                # and if they are present in the the ESE and ESS sets.
                for i in range(len(ref_seq) - 5):
                    hexamer = ref_seq[i:i + 6]
                    if hexamer in ese_set:
                        ref_motifs_ese.add(hexamer)
                    if hexamer in ess_set:
                        ref_motifs_ess.add(hexamer)

                # Check all possible hexamers in the altered sequence
                # and if they are present in the the ESE and ESS sets.
                for i in range(len(alt_seq) - 5):
                    hexamer = alt_seq[i:i + 6]
                    if hexamer in ese_set:
                        alt_motifs_ese.add(hexamer)
                    if hexamer in ess_set:
                        alt_motifs_ess.add(hexamer)

            transcript[10] = ';'.join(ref_motifs_ese)
            transcript[11] = ';'.join(alt_motifs_ese)
            transcript[12] = ';'.join(ref_motifs_ess)
            transcript[13] = ';'.join(alt_motifs_ess)
        transcripts[index] = '|'.join(transcript)
    transcripts = ','.join(transcripts)
    new_vcfline = vcf_line.split(csq.group(2))
    new_vcfline.insert(1, transcripts)
    new_vcfline = ''.join(new_vcfline)
    return new_vcfline


def reverse_complement(sequence):
    """
    Function that creates the reverse complement of a DNA sequence.

    :param sequence: sequence of DNA
    :return: reverse complement
    """
    # Complement dictionary
    complement_dictionary = {'A': 'T', 'C': 'G',
                             'T': 'A', 'G': 'C'}
    # Reverse the sequence and get the complementary base from the dictionary.
    rev_comp = ''.join(
        complement_dictionary.get(base, base) for base in reversed(sequence))
    return rev_comp


# Load the ESS and ESE hexamers from the tables
ese_file = '/home/ar7343bo-s/Resources/RESCUE-ESE_hexamers_200703.txt'
ess_file = '/home/ar7343bo-s/Resources/ESS_hexamers_200824.txt'
ese_set = set()
ess_set = set()

with open(ese_file) as ese, open(ess_file) as ess:
    for line in ese:
        ese_set.add(line.strip())

    for line in ess:
        ess_set.add(line.strip())


# Put all files in a list, removing anything that is not a vcf file.
files_directory = r"ClinvarMMSplice/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)


file_count = 1
for file in list_of_files:
    print(file_count)

    with open(files_directory + file, 'r') as vcf, open('FixedESEESS/'+file,
                                                        'w') as outfile:
        for vcfline in vcf:
            if not vcfline.startswith('#'):
                split_line = vcfline.split('\t')
                allele = split_line[-1][0:3]

                if '2' in allele:
                    alt_bases = split_line[4].split(',')
                    is_snv = False
                    base_sum_1 = len(split_line[3] + alt_bases[0])
                    base_sum_2 = len(split_line[3] + alt_bases[1])
                    if base_sum_1 == 2 and base_sum_2 == 2:
                        is_snv = True
                    if is_snv:
                        new_line = two_alleles(vcfline, alt_bases)
                        print(new_line, file=outfile, end='')
                else:
                    print(vcfline, file=outfile, end='')
            else:
                print(vcfline, file=outfile, end='')

    file_count += 1
