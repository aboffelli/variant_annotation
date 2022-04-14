#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description:

Created on: 2022-03-14
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/
"""

import time
import glob
import re
import sys

start_time = time.time()

# Prepare the vcf files that will be used.
list_of_files = glob.glob(f"{sys.argv[1]}/**/*.vcf", recursive=True)
for file in list_of_files.copy():
    # Remove any files that are not vcfs.
    if '.vcf' not in file:
        list_of_files.remove(file)

# Test data
# list_of_files = glob.glob(
#     "/Users/student/Box/Notes/TestData/Bridges/Clinvar/**/*.vcf",
#     recursive=True)

# Create a dictionary to store all the information that we need, we use sets to
# keep unique variants.
var_dict = {"Controls": {'Known': set(),
                         "Novel": set()},
            "Cases": {"Known": set(),
                      "Novel": set()}
            }

af_dict = {"Controls": {'Known': [],
                        "Novel": [],
                        "Pathogenic": []},
           "Cases": {"Known": [],
                     "Novel": [],
                     "Pathogenic": []}
           }
file_count = 1
print("Reading files")
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    # Define if it is a Control or Case file based on the path.
    file_type = file.split('/')[1]
    with open(file, 'r') as vcf:
        for vcfline in vcf:
            if not vcfline.startswith('#'):
                split_line = vcfline.split('\t')
                filt = split_line[6]
                if filt != 'PASS':
                    continue
                # Get the chromosome and position of the variant
                pos = ':'.join(split_line[0:2])

                # Isolate only the annotation info, using regex.
                csq = re.search(r'CSQ=(\S+?),(\S+?)(;ClinVar|\s)', vcfline)

                # The fixed csq contains information that are the same for all
                # transcripts
                fixed_csq = csq.group(1).split('|')
                # Retrieve the dbSNP section
                known = fixed_csq[0].split('&')[0]

                clinvar = csq.group(3)

                if not known:
                    var = 'Novel'
                else:
                    var = "Known"

                var_dict[file_type][var].add(pos)

                if "Pathogenic" in vcfline or "Likely_pathogenic" in vcfline:
                    var = "Pathogenic"

                allele_fraction = re.search(r";AF=([^A-Z;]+)", vcfline).group(1)
                af_dict[file_type][var].append(allele_fraction)

    file_count += 1

print("\nDone!")

with open('novel_known_count.txt',
          'w') as outtable:
    for file_type in var_dict:
        for var_type in var_dict[file_type]:
            print(f'''{file_type}\t{var_type}\t{
            len(var_dict[file_type][var_type])}''', file=outtable)

with open('table.txt',
          'w') as outtable:
    for file_type in af_dict:
        for var_type in af_dict[file_type]:
            for af in af_dict[file_type][var_type]:
                print(f'''{file_type}\t{var_type}\t{af}''', file=outtable)

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
