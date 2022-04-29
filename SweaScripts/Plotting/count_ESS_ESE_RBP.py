#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: ESE/ESS/RBP tables

Description: Script that evaluates and counts the ESE/ESS alterations and
    overlapping on RNA-binding proteins binding sites.

Created on: 2021-11-17
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import os
import re
import time

start_time = time.time()


def es_counter(csq_line, exist):
    """
    Function to count the alterations in ESE/ESS hexamers.
    :param str csq_line: csq section of the vcf line.
    :param str exist: dbSNP id of known variants or empty for novel variants.
    :return: None
    """

    # Define if the variant is known or novel.
    if exist:
        exist = 'known'
    else:
        exist = 'novel'

    # Isolate the transcripts into a list.
    transcript_line = csq_line.split(',')[1:]

    # Set a flag ese/ess.
    found = False

    # Loop through the transcripts.
    for transcript in transcript_line:
        transcript = transcript.split('|')

        # Retrieve the ese/ess information.
        ese_ref, ese_alt, ess_ref, ess_alt = transcript[10:14]

        # Check if there is an alteration in the ESE.
        ese_conseq = comparison(ese_ref, ese_alt)

        # If yes, add to the count and set found to True
        if ese_conseq:
            ese_count[ese_conseq][exist] += 1
            found = True

        # Check if there is an alteration in the ESS.
        ess_conseq = comparison(ess_ref, ess_alt)

        # If yes, add to the count and set found to True
        if ess_conseq:
            ess_count[ess_conseq][exist] += 1
            found = True

        # If any alteration was found, break the loop, since all transcripts
        # will have the same alteration.
        if found:
            break


def comparison(reference, alternative):
    """
    Function that evaluates the alteration of the ESE/ESS hexamer into
    alteration, creation or removal.
    :param str reference: one or more hexamers with the reference base
    :param str alternative: one or more hexamers with the alternative base
    :return: string with the type of alteration.
    """

    result = ''

    # If the reference and alternative differ.
    if reference != alternative:

        # Alteration: if both have ESEs but they are different
        if reference and alternative:

            # Account for more than one hexamer
            ref_list = reference.split(';')
            alt_list = alternative.split(';')

            # If there is more or less hexamers, set as alter.
            if len(ref_list) != len(alt_list):
                result = 'alter'

            # If there is no difference in quantity, check the hexamers one by
            # one to see if they differ.
            else:
                for hexamer in ref_list:
                    if hexamer not in alt_list:
                        result = 'alter'

        # If only the alternative has hexamers, it is being created.
        elif not reference and alternative:
            result = 'create'

        # If only the reference has hexamers, it it being removed.
        else:
            result = 'remove'

    return result


def rbp_info(csq_line, pos, exist):
    """
    Function to
    :param str csq_line: csq section of the vcf line.
    :param str pos: position of the variant in the format chr:position.
    :param str exist: dbSNP id of known variants or empty for novel variants.
    :return: None.
    """

    # Define if the variant is known or novel.
    if exist:
        exist = 'known'
    else:
        exist = 'novel'

    # Only run if the position is not in the dictionary already.
    if pos not in rbp_number:

        # Isolate the transcripts into a list.
        transcript_line = csq_line.split(',')[1:]
        unique_rbps = set()

        # Loop through the transcripts.
        for transcript in transcript_line:
            transcript = transcript.split('|')

            # Isolate the rbp section of the transcript.
            rbp = transcript[8]

            # If not empty.
            if rbp:

                # Account for multiple proteins in the same transcript, they are
                # separate by &.
                rbp = rbp.split('&')
                for prot in rbp:

                    # Store only the protein name in the set.
                    prot = prot.split(':')[0]
                    unique_rbps.add(prot)

        # If the set is not empty.
        if unique_rbps:

            # Raise the count for existing rbp.
            rbp_count[exist]['yes'] += 1

            # Store the number of proteins for the variant and if it's known or
            # novel.
            rbp_number[pos] = (len(unique_rbps), exist)

            # Store the proteins for the variant and if it's known or novel.
            rbp_prot[pos] = (unique_rbps, exist)

        # If the set is empty, raise the count for not rbp.
        else:
            rbp_count[exist]['no'] += 1


# Test
# files_directory = "/Users/student/Box/Notes/TestData/CustomAnnotation/"
# files_directory = r"C:\Users\Arthu\Documents\TestData\CustomAnnotation\\"

# Prepare the files that will be used, removing anything that is not a vcf file.
files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

# Start all the dictionaries needed.
pos_set = set()

rbp_count = {'known': {'yes': 0, 'no': 0},
             'novel': {'yes': 0, 'no': 0}
             }
rbp_prot = {}
rbp_number = {}

ese_count = {'alter': {'known': 0, 'novel': 0},
             'create': {'known': 0, 'novel': 0},
             'remove': {'known': 0, 'novel': 0}
             }

ess_count = {'alter': {'known': 0, 'novel': 0},
             'create': {'known': 0, 'novel': 0},
             'remove': {'known': 0, 'novel': 0}
             }

print("Reading files...")
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    with open(files_directory + file, 'r') as vcf_file:

        for line in vcf_file:
            # Ignore the header lines.
            if not line.startswith("#"):

                # Check if the variant passed the filters
                filt = line.split('\t')[6]
                if filt != 'PASS':
                    continue

                # Join the chr and position.
                chrom = line.split('\t')[0]
                position = line.split('\t')[1]
                position = f'{chrom}:{position}'

                # Isolate the csq section.
                csq = re.search(r'CSQ=(\S+)', line).group(1)

                # Get the dbSNP id.
                existence = csq.split('|')[0]

                # Run the functions only once for the variants, since it will be
                # the same information.
                if position not in pos_set:
                    pos_set.add(position)
                    es_counter(csq, existence)
                    rbp_info(csq, position, existence)

    # Raise the file count.
    file_count += 1
print('\nDone!')

print("\nWriting the output files...")
# Write everything in the output files
with open('ese_count.txt', 'w') as ese_out, \
        open('ess_count.txt', 'w') as ess_out:
    for key in ese_count:
        for var_type in ese_count[key]:
            print(f"{key}\t{var_type}\t{ese_count[key][var_type]}",
                  file=ese_out)

    for key in ess_count:
        for var_type in ess_count[key]:
            print(f"{key}\t{var_type}\t{ess_count[key][var_type]}",
                  file=ess_out)

with open('rbp_count.txt', 'w') as rbpc:
    for key in rbp_count:
        for key2 in rbp_count[key]:
            print(f"{key}\t{key2}\t{rbp_count[key][key2]}", file=rbpc)


with open('rbp_protein_frequency.txt', 'w') as rbpp:
    for key in rbp_prot:
        for prot in rbp_prot[key][0]:
            print(f'{prot}\t{rbp_prot[key][1]}', file=rbpp)


with open('rbp_variant_frequency.txt', 'w') as rbpv:
    for key in rbp_number:
        print(f'{key}\t{rbp_number[key][0]}\t{rbp_number[key][1]}', file=rbpv)

print("Done!")
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
