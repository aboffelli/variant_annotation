#!~/miniconda3/bin/python3
"""
Title: SOR3 to SOR4

Description: Script created to change the filter for the SOR value without
    rerunning the whole pipeline. The initial value was > 3, and the value was
    changed to > 4.

    Script recognize SOR3 in the filter column, then check the SOR value. If it
    is less than 4, the program checks if the filter line contains only SOR3
    and change it to PASS. If there are more filters, the program only removes
    the SOR3 from the list. If the number is higher than 4, only replace SOR3
    to SOR4.

Created on: 2021-10-15
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""
import os
import re
import time

start_time = time.time()

# Retrieve a list of files in the Annotation directory.
list_of_files = os.listdir("Annotation/")
# Remove any file that is not a vcf file.
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    # For each file of the list, create an edited file.
    with open("Annotation/"+file, 'r') as original, open(
            'Annotation/Edited/edited_' + file, 'w') as edited:
        for line in original:
            if not line.startswith('#'):  # Not header lines
                # store the filter column
                split_line = line.split('\t')
                filter_col = split_line[6]
                # Check if the filter contains SOR3
                if 'SOR3' in filter_col:
                    # Get the INFO column to retrieve the SOR value.
                    info = split_line[7]
                    # Match one or more digits followed by a dot and one or
                    # more digits, after SOR3=. Only retrieve the digits.
                    match = re.match("SOR=(\d+\.\d+)", info).group(1)
                    # Transform it to float and check if it's less than 4.
                    sor_value = float(match)
                    if sor_value < 4:
                        # Change to PASS if SOR3 is the only filter.
                        if filter_col == 'SOR3':
                            split_line[6] = 'PASS'
                        # Only remove SOR3 if more filters are found.
                        else:
                            split_line[6] = filter_col.replace(';SOR3', '')

                    else:  # value > 4
                        # Only replace SOR3 to SOR4
                        split_line[6] = filter_col.replace('SOR3', 'SOR4')

                # Transform everything into string and join all parts.
                line = '\t'.join([str(x) for x in split_line])
                # Write each line in the edited file. 
                edited.write(line)
            else:  # Header line
                # Finds the line containing the SOR3 Description and change it
                # to SOR4.
                if "SOR3" in line:
                    line = '##FILTER=<ID=SOR4,Description="SOR > 4.0">\n'
                # write all the header lines in the edited file.
                edited.write(line)

    file_count += 1

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
