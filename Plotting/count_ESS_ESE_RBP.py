#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2021-11-17
Author: Arthur Boffelli Castro
"""

import os
import re


def es_counter(csq_line, pos, exist):
    if exist:
        exist = 'known'
    else:
        exist = 'novel'
    if pos not in set.union(set(ese_dict), set(ess_dict)):
        transcript_line = csq_line.split(',')[1:]
        found = False
        for transcript in transcript_line:
            transcript = transcript.split('|')
            ese_ref, ese_alt, ess_ref, ess_alt = transcript[10:14]

            # ESE check
            ese_conseq = comparison(ese_ref, ese_alt)
            if ese_conseq:
                ese_dict[pos] = (ese_conseq, exist)
                found = True

            # ESS check
            ess_conseq = comparison(ess_ref, ess_alt)
            if ess_conseq:
                ess_dict[pos] = (ess_conseq, exist)
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

    if pos not in rbp_dict:
        # Check if there is an rbp
        transcript_line = csq_line.split(',')[1:]
        found = False
        for transcript in transcript_line:
            transcript = transcript.split('|')
            rbp = transcript[8]
            if rbp:
                rbp = rbp.split('&')
                rbps = []
                for prot in rbp:
                    prot = prot.split(':')[0]
                    rbps.append(prot)
                rbp_dict[pos] = (rbps, exist)

# Test
# files_directory = r"C:\Users\Arthu\Documents\CustomAnnotation\\"
files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

ese_dict = {}
ess_dict = {}
rbp_dict = {}

for file in list_of_files:
    with open(files_directory + file, 'r') as vcf_file:
        for line in vcf_file:
            if not line.startswith("#"):
                position = line.split()[1]
                csq = re.search(r'CSQ=(\S+)', line).group(1)
                existence = csq.split('|')[0]
                es_counter(csq, position, existence)
                #rbp_info(csq, position, existence)

ese_count = {'alter': {'known': 0, 'novel': 0},
             'create': {'known': 0, 'novel': 0},
             'remove': {'known': 0, 'novel': 0}
             }
ess_count = ese_count.copy()

with open('ese_count.txt', 'w') as ese_out, \
        open('ess_count.txt', 'w') as ess_out:
    for position in ese_dict:
        info = ese_dict[position]
        ese_count[info[0]][info[1]] += 1
    for key in ese_count:
        for var_type in ese_count[key]:
            print(f"{key}\t{var_type}\t{ese_count[key][var_type]}",
                  file=ese_out)

    for position in ess_dict:
        info = ess_dict[position]
        ess_count[info[0]][info[1]] += 1
    for key in ess_count:
        for var_type in ess_count[key]:
            print(f"{key}\t{var_type}\t{ess_count[key][var_type]}",
                  file=ess_out)

# print(rbp_dict)
# TODO:
#  RBP:
#   Check if there are a RBP (just yes or no)
#   Get the protein name, to sort by frequency (unique prot in the variants).
#   Count the number of RBPs in each variant (also unique prot).
