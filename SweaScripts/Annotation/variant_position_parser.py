#!~/miniconda3/bin/python3

"""
Title: Variant Position Parser

Description: Script to parse, compare and identify if the variants are located
    in the targeted genes.

Created on: 2021-10-04
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import os
import sys
import re
import time

start_time = time.time()

# Create a dictionary with the flanking region size used in the sequencing for
# the
flanking_regions = {'BRCA1': 150000, 'BRCA2': 150000, 'CDH1': 50000,
                    'CHEK2': 50000, 'PTEN': 50000, 'STK11': 50000,
                    'TP53': 50000, 'ATM': 20000, 'BARD1': 20000,
                    'BRIP1': 20000, 'CDKN2A': 20000, 'MRE11': 20000,
                    'NBN': 20000, 'PALB2': 20000, 'RAD50': 20000,
                    'RAD51C': 20000, 'RAD51D': 20000
                    }

# Retrieve the positions of the genes that where sequenced. The first argument
# in the command line should be the gff file with the positions of the genes,
# and the second argument should be a file with all the genes targeted in the
# sequencing listed.
if len(sys.argv) >= 3:
    gff = sys.argv[1]
    whole_genes = sys.argv[2]

    # Initiate a list to store the list of genes
    gl = list()
    gene_dict = dict()
    with open(gff, 'r') as resource, open(whole_genes, 'r') as gene_list:

        # Load the list genes from the file.
        for line_split in gene_list:
            gl.append(line_split.strip())

        # Loop through the gff file
        for gff_line in resource:

            # Ignore comment lines.
            if not gff_line.startswith("#"):
                gff_line = gff_line.strip().split('\t')

                # Get only the gene lines.
                if gff_line[2] == "gene":

                    # The gene name is located in the 9th column, after
                    # gene_name.
                    gene_name = re.search(r"gene_name=(\w+)",
                                          gff_line[8]).group(1)

                    # Loop through the list of targeted genes.
                    for gene in gl:

                        # If the same gene.
                        if gene == gene_name:

                            # Get the chromosome number.
                            chromosome = gff_line[0].strip("chr")

                            # Check if it is one of the genes that were
                            # sequenced with flanking regions.
                            if gene in flanking_regions:
                                # Store the size of the flanking region.
                                fr = flanking_regions[gene]
                            else:
                                # Set it to zero if there is no flanking region.
                                fr = 0

                            # Create a new key in the dictionary if the
                            # chromosome is not there yet.
                            if chromosome not in gene_dict:
                                # For each chromosome, create a dictionary with
                                # the start and end of each gene, adding the
                                # flanking regions.
                                gene_dict[chromosome] = {gene_name: (
                                    int(gff_line[3]) - fr,
                                    int(gff_line[4]) + fr)}

                                break

                            # If the chromosome already exists in the
                            # dictionary.
                            else:
                                # Add the new gene with the respective start and
                                # end positions, adding the flanking regions.
                                gene_dict[chromosome][gene_name] = (
                                    int(gff_line[3]) - fr,
                                    int(gff_line[4]) + fr)
                                break

    # Prepare the files that will be used, removing anything that is not a vcf.
    list_of_files = os.listdir("Annotation/Edited")
    for file in list_of_files.copy():
        if '.vcf' not in file:
            list_of_files.remove(file)

    # Create the new directories for the output files.
    new_directories = ['Annotation/Edited/KnownOnFail/',
                       'Annotation/Edited/KnownOffFail/',
                       'Annotation/Edited/KnownOnPass/',
                       'Annotation/Edited/KnownOffPass/',
                       'Annotation/Edited/NovelOnFail/',
                       'Annotation/Edited/NovelOffFail/',
                       'Annotation/Edited/NovelOnPass/',
                       'Annotation/Edited/NovelOffPass/']

    for directory in new_directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

    # Start a counter for each filter
    filters = {"QD2S": 0, "MQ40": 0, "FS60": 0,
               "SOR4": 0, "MQRS-12.5": 0, "RPRS-8": 0,
               "QD2I": 0, "FS200": 0, "SOR10": 0, "RPRS-20": 0
               }

    # Open the count files and create the headers
    with open('type_comparison.txt', 'w') as out, open(
            "known_on_target_fail.txt", 'w') as k_on_fail, open(
            "known_off_target_fail.txt", "w") as k_off_fail, open(
            "novel_on_target_fail.txt", 'w') as n_on_fail, open(
            "novel_off_target_fail.txt", "w") as n_off_fail:

        # Header for type_comparison.txt
        out.write("Sample\tKnownOnTargetPass\tKnownOnTargetFail\t"
                  "KnownOffTargetPass\tKnownOffTargetFail\tNovelOnTargetPass\t"
                  "NovelOnTargetFail\tNovelOffTargetPass\tNovelOffTargetFail")

        # Header for fail output files
        k_on_fail.write("Sample")
        k_off_fail.write("Sample")
        n_on_fail.write("Sample")
        n_off_fail.write("Sample")

        # Add the filter names in the header.
        for key in filters:
            k_on_fail.write("\t" + key)
            k_off_fail.write("\t" + key)
            n_on_fail.write("\t" + key)
            n_off_fail.write("\t" + key)

        # Loop through the vcf files
        for file in list_of_files:

            # For each file reset all the counters.
            sample = file.lstrip('fs_filtered_').rstrip(".vcf")
            k_on_fail_filters = filters.copy()
            k_off_fail_filters = filters.copy()
            n_on_fail_filters = filters.copy()
            n_off_fail_filters = filters.copy()
            counter = {"k_on_pass": 0, "k_on_fail": 0, "k_off_pass": 0,
                       "k_off_fail": 0, "n_on_pass": 0, "n_on_fail": 0,
                       "n_off_pass": 0, "n_off_fail": 0}

            # Open all the output files.
            with open('Annotation/Edited/' + file, 'r') as vcf, open(
                'Annotation/Edited/KnownOnFail/k_on_fail_' + file,
                    'w') as k_on_fail_out, open(
                'Annotation/Edited/KnownOffFail/k_off_fail_' + file,
                    'w') as k_off_fail_out, open(
                'Annotation/Edited/KnownOnPass/k_on_pass_' + file,
                    'w') as k_on_pass_out, open(
                'Annotation/Edited/KnownOffPass/k_off_pass_' + file,
                    'w') as k_off_pass_out, open(
                'Annotation/Edited/NovelOnFail/n_on_fail_' + file,
                    'w') as n_on_fail_out, open(
                'Annotation/Edited/NovelOffFail/n_off_fail_' + file,
                    'w') as n_off_fail_out, open(
                'Annotation/Edited/NovelOnPass/n_on_pass_' + file,
                    'w') as n_on_pass_out, open(
                'Annotation/Edited/NovelOffPass/n_off_pass_' + file,
                    'w') as n_off_pass_out:

                for line in vcf:

                    # Get the information we need from each variant.
                    if not line.startswith("#"):
                        line_split = line.split()
                        chrom = line_split[0]
                        pos = line_split[1]
                        filter_col = line_split[6]
                        exist = line_split[7].split('|')[1]

                        # Check if the variant is on or off the targeted areas.
                        if chrom in gene_dict:

                            # Loop through all the genes in our dictionary, and
                            # store the start an end positions in variables.
                            for gene in gene_dict[chrom]:
                                start = gene_dict[chrom][gene][0]
                                end = gene_dict[chrom][gene][1]

                                # If the variant position is between the start
                                # and end of a gene, it is on target.
                                if start <= int(pos) <= end:

                                    # Check if it is a known variant.
                                    if exist:

                                        # Check if the variant passed the
                                        # filters and add to the respective
                                        # counter.
                                        if filter_col == 'PASS':
                                            counter["k_on_pass"] += 1
                                            k_on_pass_out.write(line)
                                            break
                                        else:  # Fail
                                            counter["k_on_fail"] += 1
                                            k_on_fail_out.write(line)

                                            # The variant may fail more than one
                                            # filter, so make sure we count all
                                            # of them.
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

                            # If the variant is not between the start and end,
                            # it is off target.
                            else:

                                # Check if it is a known variant.
                                if exist:

                                    # Check if the variant passed the
                                    # filters and add to the respective
                                    # counter.
                                    if filter_col == 'PASS':  # Pass
                                        counter["k_off_pass"] += 1
                                        k_off_pass_out.write(line)

                                    else:  # Fail
                                        counter["k_off_fail"] += 1
                                        k_off_fail_out.write(line)

                                        # The variant may fail more than one
                                        # filter, so make sure we count all
                                        # of them.
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

                        # If the chromosome is not in the dictionary, it is also
                        # off target.
                        else:

                            # Check if it is a known variant.
                            if exist:

                                # Check if the variant passed the
                                # filters and add to the respective
                                # counter.
                                if filter_col == 'PASS':  # Pass
                                    counter["k_off_pass"] += 1
                                    k_off_pass_out.write(line)

                                else:  # Fail
                                    counter["k_off_fail"] += 1
                                    k_off_fail_out.write(line)

                                    # The variant may fail more than one
                                    # filter, so make sure we count all
                                    # of them.
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

            # Print the values on file type_comparison.txt before going to the
            # next vcf file.
            out.write("\n" + sample)
            for key in counter:
                out.write("\t" + str(counter[key]))

            # Print the values on fail output files. Start the line with the
            # sample name, and add the count.
            k_on_fail.write("\n" + sample)
            k_off_fail.write("\n" + sample)
            n_on_fail.write("\n" + sample)
            n_off_fail.write("\n" + sample)
            for key in filters:
                k_on_fail.write("\t" + str(k_on_fail_filters[key]))
                k_off_fail.write("\t" + str(k_off_fail_filters[key]))
                n_on_fail.write("\t" + str(n_on_fail_filters[key]))
                n_off_fail.write("\t" + str(n_off_fail_filters[key]))

    print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))

else:
    print("Usage:\npython variant_position_parser.py file.gff gene_list.txt")
