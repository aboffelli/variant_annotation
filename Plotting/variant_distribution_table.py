import re
import os

variant_types_known = {}
variant_types_novel = {}

files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' in file:
        with open(files_directory+file,'r') as vcf_file:
            for line in vcf_file:
                if not line.startswith("#"):
                    filt = line.split()[6]
                    if filt != 'PASS':
                        continue
                    csq = re.search(r'CSQ=(\S+)', line).group(1)
                    exist = csq.split('|')[0]
                    transcripts = csq.split(',')[1:]
                    for transcript in transcripts:
                        conseq = transcript.split('|')[4]
                        if exist:
                            if conseq not in variant_types_known:
                                variant_types_known[conseq] = 1
                            else:
                                variant_types_known[conseq] += 1
                        else:
                            if conseq not in variant_types_novel:
                                variant_types_novel[conseq] = 1
                            else:
                                variant_types_novel[conseq] += 1


with open('known_variant_distribution.txt', 'w') as known, \
        open('novel_variant_distribution.txt', 'w') as novel:
    for i in variant_types_known:
        print(i + '\t', variant_types_known[i], file=known)
    for j in variant_types_novel:
        print(j + '\t', variant_types_novel[j], file=novel)
