#!~/miniconda3/bin/python3

import os
import re

list_of_files = os.listdir("Annotation/")
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

for file in list_of_files:
    with open("Annotation/"+file, 'r') as original, open('Annotation/Edited/edited_' + file, 'w') as edited:
        for line in original:
            if not line.startswith('#'):
                split_line = line.split('\t')
                filter_col = split_line[6]
                if 'SOR3' in filter_col:
                    info = split_line[7]
                    match = re.match("\S*SOR=(\d+.\d+)", info)
                    print(match.group(1))
                    sor_value = float(match.group(1))
                    if sor_value < 4:
                        if filter_col == 'SOR3':
                            split_line[6] = 'PASS'
                        else:
                            split_line[6] = filter_col.replace(';SOR3', '')
                    else:
                        if filter_col == 'SOR3':
                            split_line[6] = 'SOR4'
                        else:
                            split_line[6] = filter_col.replace('SOR3', 'SOR4')
                line = '\t'.join([str(x) for x in split_line])
                edited.write(line)
            else:  # Header line
                if "SOR3" in line:
                    line = '##FILTER=<ID=SOR4,Description="SOR > 4.0">\n'
                edited.write(line)
