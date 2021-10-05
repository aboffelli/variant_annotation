#!~/miniconda3/bin/python3

"""
Script to parse and compare and identify in which gene each variant is located.
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

filters = {"QD2S": 0, "MQ40": 0, "FS60": 0,
           "SOR3": 0, "MQRS-12.5": 0, "RPRS-8": 0,
           "QD2I": 0, "FS200": 0, "SOR10": 0, "RPRS-20": 0
           }


# Retrieve the positions of the genes that where sequenced
if len(sys.argv) >= 3:
    gff = sys.argv[1]
    whole_genes = sys.argv[2]

    gl = list()
    gene_dict = dict()
    with open(gff, 'r') as resource, open(whole_genes, 'r') as gene_list:
        for line in gene_list:
            gl.append(line.strip())
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

    list_of_files = os.listdir("Annotation/")

    with open('number_var.txt', 'w') as out, open("on_target_fail.txt", 'w') as on, open("off_target_fail.txt",
                                                                                         "w") as off:
        # Header for number_var.txt
        print("Sample\tOnTargetPass\tOnTargetFail\tOffTargetPass\tOffTargetFail", file=out)
        # Header for on_target_fail.txt and off_target_fail.txt
        print("Sample", file=on, end="\t")
        print("Sample", file=off, end="\t")
        for key in filters:
            print(key, file=on, end="\t")
            print(key, file=off, end="\t")
        print("", file=on)
        print("", file=off)

        for file in list_of_files:
            sample = file.lstrip('fs_filtered_').rstrip(".vcf")
            on_fail_filters = filters.copy()
            off_fail_filters = filters.copy()
            on_pass = 0
            on_fail = 0
            off_pass = 0
            off_fail = 0
            with open('Annotation/'+file, 'r') as vcf:
                for line in vcf:
                    if not line.startswith("#"):
                        line = line.split()
                        chrom = line[0]
                        pos = line[1]
                        filter_col = line[6]
                        if chrom in gene_dict:
                            for gene in gene_dict[chrom]:
                                start = gene_dict[chrom][gene][0]
                                end = gene_dict[chrom][gene][1]
                                if start < int(pos) < end:
                                    if filter_col == 'PASS':
                                        on_pass += 1
                                        break
                                    else:  # Not PASS
                                        on_fail += 1
                                        for f in filter_col.split(";"):
                                            on_fail_filters[f] += 1
                                        break
                            else:  # Not inside any gene
                                if filter_col == 'PASS':
                                    off_pass += 1
                                else:  # FAIL
                                    off_fail += 1
                                    for f in filter_col.split(";"):
                                        off_fail_filters[f] += 1
                        else:  # Chromosome that is not in the dict.
                            if filter_col == 'PASS':
                                off_pass += 1
                            else:  # FAIL
                                off_fail += 1
                                for f in filter_col.split(";"):
                                    off_fail_filters[f] += 1

            # Print the values on file number_var.txt
            print("{}\t{}\t{}\t{}\t{}".format(sample, on_pass, on_fail, off_pass, off_fail), file=out)

            # Print the values on on_target_fail.txt and off_target_fail.txt
            print(sample, file=on, end="\t")
            print(sample, file=off, end="\t")
            for key in filters:
                print(on_fail_filters[key], file=on, end="\t")
                print(off_fail_filters[key], file=off, end="\t")
            print("", file=on)
            print("", file=off)

else:
    print("Usage:\npython variant_position_parser.py file.gff gene_list.txt")
