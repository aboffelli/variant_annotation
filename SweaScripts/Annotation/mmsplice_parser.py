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
            files_directory + clinvar_file, 'r') as clinvar:
        pass

# Print the run time.
print('Run time: {:.2f} seconds'.format(time.time() - start_time))