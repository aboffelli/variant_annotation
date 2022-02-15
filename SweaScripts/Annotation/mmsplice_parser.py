#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-02-10
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""

import os
import re
import time

start_time = time.time()

# Put all files in a list, removing anything that is not a vcf file.
files_directory = "/Users/student/Box/Notes/TestData/ClinVar/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file or 'mmsplice' not in file:
        list_of_files.remove(file)

for vcf_file in list_of_files:
    clinvar_file = 'clinvar' + vcf_file.replace('mmsplice',                                                     '')
    with open(files_directory + vcf_file, 'r') as mmsplice, open(
            files_directory + clinvar_file, 'r') as clinvar, open(
        files_directory + 'mmsplice_' + clinvar_file, 'w') as outfile:
        mmsplice_dict = {}
        header_lines = []
        for mmsplice_line in mmsplice:
            if not mmsplice_line.startswith("#"):
                mm_position = '\t'.join(mmsplice_line.split('\t')[0:2])
                mmsplice_dict[mm_position] = {}
                mmsplice_info = re.search(r"MMSplice=(\S+)",
                                          mmsplice_line).group(1).split(',')
                for mmsplice_transc in mmsplice_info:
                    transc_name = mmsplice_transc.split('|')[0]
                    delta_logit = mmsplice_transc.split('|')[1]
                    if delta_logit:
                        delta_logit = str(round(float(
                            mmsplice_transc.split('|')[1]), 4))

                    mmsplice_dict[mm_position][transc_name] = delta_logit

        for vcfline in clinvar:
            vcfline = vcfline.strip()
            if vcfline.startswith('#'):
                if vcfline.startswith("##INFO=<ID=CSQ"):
                    vcfline = vcfline.rstrip('">')
                    vcfline += '|mmsplice_delta_logit_psi">'
                    print(vcfline, file=outfile)

                elif vcfline.startswith("#CHROM"):
                    print("##mmsplice_delta_logit_psi=delta logit psi score of "
                          "variant", file=outfile)
                    print(vcfline, file=outfile)

                else:
                    print(vcfline, file=outfile)

            else:
                position = '\t'.join(vcfline.split('\t')[0:2])
                csq = re.search(r"CSQ=(\S+;Cl|\S+)", vcfline
                                ). group(1).strip(';Cl')

                transcripts = csq.split(',')[1:]

                for index, transcript in enumerate(transcripts):
                    transcript = transcript.split('|')
                    if transcript[2]:
                        mmsplice_value = mmsplice_dict[position][transcript[2]]
                    else:
                        mmsplice_value = ''
                    transcript.append(mmsplice_value)
                    transcript = '|'.join(transcript)
                    transcripts[index] = transcript

                transcripts = ','.join(transcripts)
                to_replace = re.search(r"CSQ=\S+?,(\S+;Cl|\S+)", vcfline
                                       ). group(1).strip(';Cl')
                vcfline = vcfline.split(to_replace)
                vcfline.insert(1, transcripts)
                vcfline = ''.join(vcfline)
                print(vcfline, file=outfile)

# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))