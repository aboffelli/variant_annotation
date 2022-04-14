#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MMSplice annotation parser

Description: Since the MMSplice annotation is a slow process, it was running
    in the background while the workflow continued. This script takes the
    MMSplice annotation and adds it to the more updated files.

Created on: 2022-02-10
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import os
import re
import time

start_time = time.time()

# Put all files with MMSplice annotation in a list, removing anything that is
# not a vcf file.
files_directory = "MMSplice/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

# Set the directory for the files that we will add the MMSplice value.
clinvar_directory = 'ClinVar/'

file_count = 1
# Loop through the MMSplice files
for vcf_file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)
    # Find the same file in the ClinVar directory.
    clinvar_file = 'clinvar' + vcf_file.replace('mmsplice', '')

    with open(files_directory + vcf_file, 'r') as mmsplice, open(
            clinvar_directory + clinvar_file, 'r') as clinvar, open(
            'ClinvarMMSplice/mmsplice_' + clinvar_file, 'w') as outfile:

        # Firstly we load the MMSplice information in a dictionary.
        mmsplice_dict = {}
        for mmsplice_line in mmsplice:

            # We are only interested in the variant lines.
            if not mmsplice_line.startswith("#"):
                # Set the position to work as a dictionary key.
                mm_position = '\t'.join(mmsplice_line.split('\t')[0:2])

                # For each position, create a dictionary
                mmsplice_dict[mm_position] = {}

                # Retrieve the MMSplice annotation in the line and split in the
                # transcripts.
                mmsplice_info = re.search(r"MMSplice=(\S+)",
                                          mmsplice_line).group(1).split(',')

                # For each transcript, save the transcript name and the delta
                # logit psi value.
                for mmsplice_transc in mmsplice_info:
                    transc_name = mmsplice_transc.split('|')[0]
                    delta_logit = mmsplice_transc.split('|')[1]

                    # If there is a delta logit psi value, round it to 4
                    # numbers.
                    if delta_logit:
                        delta_logit = str(round(float(
                            mmsplice_transc.split('|')[1]), 4))

                    # The information is then added in the nested dictionary.
                    # For each position, we have the transcript names as keys
                    # and delta logit psi as values.
                    mmsplice_dict[mm_position][transc_name] = delta_logit

        # Now we loop throug the clinvar file.
        for vcfline in clinvar:
            vcfline = vcfline.strip()

            # Header lines.
            if vcfline.startswith('#'):

                # Update the CSQ info.
                if vcfline.startswith("##INFO=<ID=CSQ"):
                    vcfline = vcfline.rstrip('">')
                    vcfline += '|mmsplice_delta_logit_psi">'
                    print(vcfline, file=outfile)

                # Add the line explaining the delta logit psi as the last line
                # before the header for the variants.
                elif vcfline.startswith("#CHROM"):
                    print("##mmsplice_delta_logit_psi=delta logit psi score of "
                          "variant", file=outfile)
                    print(vcfline, file=outfile)

                # Print all the other header lines.
                else:
                    print(vcfline, file=outfile)

            # Variant lines.
            else:
                # Store the position of the variant
                position = '\t'.join(vcfline.split('\t')[0:2])

                # Retrieve the CSQ section of the line.
                csq = re.search(r"CSQ=(\S+;Cl|\S+)", vcfline
                                ).group(1).strip(';Cl')

                # Isolate the transcripts, removing the fixed CSQ.
                transcripts = csq.split(',')[1:]

                # Loop through the transcripts.
                for index, transcript in enumerate(transcripts):
                    transcript = transcript.split('|')
                    # Get the respective value in the mmsplice dictionary using
                    # the transcript name
                    if transcript[2]:
                        mmsplice_value = mmsplice_dict[position][transcript[2]]

                    # If there is no transcript name, set the mmsplice to empty.
                    else:
                        mmsplice_value = ''

                    # Add the value to the end of the list, join the transcript
                    # and put the new transcript in the right position in the
                    # transcripts list.
                    transcript.append(mmsplice_value)
                    transcript = '|'.join(transcript)
                    transcripts[index] = transcript

                # After all transcripts were changed, join the transcripts list.
                transcripts = ','.join(transcripts)

                # Here we isolate the part that we want to change, which is only
                # the transcripts part.
                to_replace = re.search(r"CSQ=\S+?,(\S+;Cl|\S+)", vcfline
                                       ).group(1).strip(';Cl')

                # Then we split the line in two, using the part that we are
                # replacing. This removes the transcripts section.
                vcfline = vcfline.split(to_replace)

                # We insert the new transcripts in the right position.
                vcfline.insert(1, transcripts)

                # Join the line and print in the file.
                vcfline = ''.join(vcfline)
                print(vcfline, file=outfile)

    file_count += 1

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
