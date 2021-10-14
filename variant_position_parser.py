#!~/miniconda3/bin/python3

"""
Script to parse, compare and identify in which gene each variant is located.
"""

import os
import sys

flanking_regions = {'BRCA1': 150000, 'BRCA2': 150000, 'CDH1': 50000,
                    'CHEK2': 50000, 'PTEN': 50000, 'STK11': 50000,
                    'TP53': 50000, 'ATM': 20000, 'BARD1': 20000,
                    'BRIP1': 20000, 'CDKN2A': 20000, 'MRE11': 20000,
                    'NBN': 20000, 'PALB2': 20000, 'RAD50': 20000,
                    'RAD51C': 20000, 'RAD51D': 20000
                    }

# Retrieve the positions of the genes that where sequenced
if len(sys.argv) >= 3:
    gff = sys.argv[1]
    whole_genes = sys.argv[2]

    gl = list()
    gene_dict = dict()
    with open(gff, 'r') as resource, open(whole_genes, 'r') as gene_list:
        for line_split in gene_list:
            gl.append(line_split.strip())
        for gff_line in resource:
            if not gff_line.startswith("#"):
                gff_line = gff_line.strip().split('\t')
                if gff_line[2] == "gene":
                    gene_name = gff_line[8].split(';')[3].strip('gene_name=')
                    for gene in gl:
                        if gene == gene_name:
                            chromosome = gff_line[0].strip("chr")
                            if gene in flanking_regions:
                                fr = flanking_regions[gene]
                            else:
                                fr = 0
                            if chromosome not in gene_dict:
                                gene_dict[chromosome] = {gene_name: (int(gff_line[3])-fr, int(gff_line[4])+fr)}
                                break

                            else:
                                gene_dict[chromosome][gene_name] = (int(gff_line[3])-fr, int(gff_line[4])+fr)
                                break

########################################################################################################################
    # Create the new directories
    list_of_files = os.listdir("Annotation/")
    for file in list_of_files.copy():
        if '.vcf' not in file:
            list_of_files.remove(file)
    new_directories = ['Annotation/KnownOnFail/',
                       'Annotation/KnownOffFail/',
                       'Annotation/KnownOnPass/',
                       'Annotation/KnownOffPass/',
                       'Annotation/NovelOnFail/',
                       'Annotation/NovelOffFail/',
                       'Annotation/NovelOnPass/',
                       'Annotation/NovelOffPass/']

    for directory in new_directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

########################################################################################################################
    # Start a counter for filters
    filters = {"QD2S": 0, "MQ40": 0, "FS60": 0,
               "SOR3": 0, "MQRS-12.5": 0, "RPRS-8": 0,
               "QD2I": 0, "FS200": 0, "SOR10": 0, "RPRS-20": 0
               }

########################################################################################################################
    # Open count files and create the headers
    with open('type_comparison.txt', 'w') as out, open("known_on_target_fail.txt", 'w') as k_on_fail, \
            open("known_off_target_fail.txt", "w") as k_off_fail, open("novel_on_target_fail.txt", 'w') as n_on_fail, \
            open("novel_off_target_fail.txt", "w") as n_off_fail:

        # Header for type_comparison.txt
        out.write("Sample\tKnownOnTargetPass\tKnownOnTargetFail\tKnownOffTargetPass\tKnownOffTargetFail\t"
                  "NovelOnTargetPass\tNovelOnTargetFail\tNovelOffTargetPass\tNovelOffTargetFail")

        # Header for fail output files
        k_on_fail.write("Sample")
        k_off_fail.write("Sample")
        n_on_fail.write("Sample")
        n_off_fail.write("Sample")

        for key in filters:
            k_on_fail.write("\t" + key)
            k_off_fail.write("\t" + key)
            n_on_fail.write("\t" + key)
            n_off_fail.write("\t" + key)

        # Loop through the vcf files
        for file in list_of_files:
            # For each file reset the variables
            sample = file.lstrip('fs_filtered_').rstrip(".vcf")
            k_on_fail_filters = filters.copy()
            k_off_fail_filters = filters.copy()
            n_on_fail_filters = filters.copy()
            n_off_fail_filters = filters.copy()
            counter = {"k_on_pass": 0, "k_on_fail": 0, "k_off_pass": 0, "k_off_fail": 0,
                       "n_on_pass": 0, "n_on_fail": 0, "n_off_pass": 0, "n_off_fail": 0}

            with open('Annotation/'+file, 'r') as vcf, \
                    open('Annotation/KnownOnFail/k_on_fail_'+file, 'w') as k_on_fail_out, \
                    open('Annotation/KnownOffFail/k_off_fail_'+file, 'w') as k_off_fail_out, \
                    open('Annotation/KnownOnPass/k_on_pass_'+file, 'w') as k_on_pass_out, \
                    open('Annotation/KnownOffPass/k_off_pass_'+file, 'w') as k_off_pass_out, \
                    open('Annotation/NovelOnFail/n_on_fail_'+file, 'w') as n_on_fail_out, \
                    open('Annotation/NovelOffFail/n_off_fail_'+file, 'w') as n_off_fail_out, \
                    open('Annotation/NovelOnPass/n_on_pass_'+file, 'w') as n_on_pass_out, \
                    open('Annotation/NovelOffPass/n_off_pass_'+file, 'w') as n_off_pass_out:

                for line in vcf:
                    # Get the information we need of each line.
                    if not line.startswith("#"):
                        line_split = line.split()
                        chrom = line_split[0]
                        pos = line_split[1]
                        filter_col = line_split[6]
                        exist = line_split[7].split('|')[1]

                        # On or Off target
                        if chrom in gene_dict:
                            for gene in gene_dict[chrom]:
                                start = gene_dict[chrom][gene][0]
                                end = gene_dict[chrom][gene][1]
                                if start < int(pos) < end:  # On target
                                    if exist:  # Known
                                        if filter_col == 'PASS':  # Pass
                                            counter["k_on_pass"] += 1
                                            k_on_pass_out.write(line)
                                            break
                                        else:  # Fail
                                            counter["k_on_fail"] += 1
                                            k_on_fail_out.write(line)
                                            for f in filter_col.split(";"):
                                                k_on_fail_filters[f] += 1
                                            break
                                    else:  # Novel
                                        if filter_col == 'PASS':  # Pass
                                            counter["n_on_pass"] += 1
                                            n_on_pass_out.write(line)
                                            break
                                        else:  # Fail
                                            counter["n_on_fail"] += 1
                                            n_on_fail_out.write(line)
                                            for f in filter_col.split(";"):
                                                n_on_fail_filters[f] += 1
                                            break
                            else:  # Off target
                                if exist:  # Known
                                    if filter_col == 'PASS':  # Pass
                                        counter["k_off_pass"] += 1
                                        k_off_pass_out.write(line)
                                    else:  # Fail
                                        counter["k_off_fail"] += 1
                                        k_off_fail_out.write(line)
                                        for f in filter_col.split(";"):
                                            k_off_fail_filters[f] += 1
                                else:  # Novel
                                    if filter_col == 'PASS':  # Pass
                                        counter["n_off_pass"] += 1
                                        n_off_pass_out.write(line)
                                    else:  # Fail
                                        counter["n_off_fail"] += 1
                                        n_off_fail_out.write(line)
                                        for f in filter_col.split(";"):
                                            n_off_fail_filters[f] += 1

                        else:  # Chromosome that is not in the dict.
                            if exist:  # Known
                                if filter_col == 'PASS':  # Pass
                                    counter["k_off_pass"] += 1
                                    k_off_pass_out.write(line)
                                else:  # Fail
                                    counter["k_off_fail"] += 1
                                    k_off_fail_out.write(line)
                                    for f in filter_col.split(";"):
                                        k_off_fail_filters[f] += 1
                            else:  # Novel
                                if filter_col == 'PASS':  # Pass
                                    counter["n_off_pass"] += 1
                                    n_off_pass_out.write(line)

                                else:  # Fail
                                    counter["n_off_fail"] += 1
                                    n_off_fail_out.write(line)
                                    for f in filter_col.split(";"):
                                        n_off_fail_filters[f] += 1

            # Print the values on file type_comparison.txt
            out.write("\n" + sample)
            for key in counter:
                out.write("\t" + str(counter[key]))

            # Print the values on fail output files
            k_on_fail.write("\n" + sample)
            k_off_fail.write("\n" + sample)
            n_on_fail.write("\n" + sample)
            n_off_fail.write("\n" + sample)
            for key in filters:
                k_on_fail.write("\t" + str(k_on_fail_filters[key]))
                k_off_fail.write("\t" + str(k_off_fail_filters[key]))
                n_on_fail.write("\t" + str(n_on_fail_filters[key]))
                n_off_fail.write("\t" + str(n_off_fail_filters[key]))

else:
    print("Usage:\npython variant_position_parser.py file.gff gene_list.txt")
