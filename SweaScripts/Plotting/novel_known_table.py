#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-04-20
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""

import time
import glob
import re

start_time = time.time()

# Prepare the files that are being read.
list_of_files = glob.glob("ClinvarMMSplice/*.vcf", recursive=True)
for file in list_of_files.copy():
    # Remove any files that are not vcfs.
    if '.vcf' not in file:
        list_of_files.remove(file)

var_dict = {"Known": set(),
            "Novel": set()}

file_count = 1
print("Reading files...")
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    with open(file, 'r') as vcf:
        for vcfline in vcf:
            # Ignore the header lines.
            if not vcfline.startswith("#"):
                split_line = vcfline.split('\t')

                # Exclude variants that did not pass the filters.
                filt = split_line[6]
                if filt != 'PASS':
                    continue

                # Get the chromosome and position of the variant
                pos = ':'.join(split_line[0:2])

                # Isolate only the annotation info, using regex.
                csq = re.search(r'CSQ=(\S+?),\S+?(;ClinVar|\s)', vcfline)

                # The fixed csq contains information that are the same for all
                # transcripts
                fixed_csq = csq.group(1).split('|')

                # Retrieve the dbSNP ID
                known = fixed_csq[0].split('&')[0]

                # Set the type of variant (Known/Novel).
                if not known:
                    var = 'Novel'
                else:
                    var = "Known"

                # Add the position to the dictionary in the respective file
                # type and variant.
                var_dict[var].add(pos)

    file_count += 1
print("\nDone!")

# Print the file.
print("Writing the table...")
with open('novel_known_count.txt', 'w') as outtable:
    for var_type in var_dict:
        print(f'''{var_type}\t{len(var_dict[var_type])}''', file=outtable)
print("Done!")

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
