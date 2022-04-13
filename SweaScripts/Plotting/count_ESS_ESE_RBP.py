#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2021-11-17
Author: Arthur Boffelli Castro
"""
# TODO: add comments.

import os
import re


def es_counter(csq_line, exist):
    if exist:
        exist = 'known'
    else:
        exist = 'novel'

    transcript_line = csq_line.split(',')[1:]
    found = False
    for transcript in transcript_line:
        transcript = transcript.split('|')
        ese_ref, ese_alt, ess_ref, ess_alt = transcript[10:14]

        # ESE check
        ese_conseq = comparison(ese_ref, ese_alt)
        if ese_conseq:
            ese_count[ese_conseq][exist] += 1
            found = True

        # ESS check
        ess_conseq = comparison(ess_ref, ess_alt)
        if ess_conseq:
            ess_count[ess_conseq][exist] += 1
            found = True

        if found:
            break


def comparison(reference, alteration):
    result = ''
    if reference != alteration:
        # Alteration: if both have ESEs but they are different
        if reference and alteration:
            ref_list = reference.split(';')
            alt_list = alteration.split(';')
            if len(ref_list) != len(alt_list):
                result = 'alter'
            else:
                for hexamer in ref_list:
                    if hexamer not in alt_list:
                        result = 'alter'
        elif not reference and alteration:
            result = 'create'

        else:  # if reference and not alteration.
            result = 'remove'

    return result


def rbp_info(csq_line, pos, exist):
    if exist:
        exist = 'known'
    else:
        exist = 'novel'

    if pos not in rbp_number:
        # Check if there is an rbp
        transcript_line = csq_line.split(',')[1:]
        unique_rbps = set()
        for transcript in transcript_line:
            transcript = transcript.split('|')
            rbp = transcript[8]
            if rbp:

                rbp = rbp.split('&')
                for prot in rbp:
                    prot = prot.split(':')[0]
                    unique_rbps.add(prot)

        if unique_rbps:
            rbp_count[exist]['yes'] += 1
            rbp_number[pos] = (len(unique_rbps), exist)
            rbp_prot[pos] = (unique_rbps, exist)
        else:
            rbp_count[exist]['no'] += 1


# Test
# files_directory = "/Users/student/Box/Notes/TestData/CustomAnnotation/"
# files_directory = r"C:\Users\Arthu\Documents\TestData\CustomAnnotation\\"
files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

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

count = 1
for file in list_of_files:
    print(count, end='\r')
    with open(files_directory + file, 'r') as vcf_file:
        for line in vcf_file:
            if not line.startswith("#"):
                filt = line.split('\t')[6]
                if filt != 'PASS':
                    continue
                chrom = line.split('\t')[0]
                position = line.split('\t')[1]
                position = f'{chrom}:{position}'
                csq = re.search(r'CSQ=(\S+)', line).group(1)
                existence = csq.split('|')[0]
                if position not in pos_set:
                    pos_set.add(position)
                    es_counter(csq, existence)
                    rbp_info(csq, position, existence)
    count += 1
print('Files are over')


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
