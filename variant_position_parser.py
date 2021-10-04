#!~/miniconda3/bin/python3

"""
Script to parse and compare and identify in which gene each variant is located.
"""

import os

# Retrieve the positions of whole genes the genes that where sequenced
gff = '/home/ar7343bo-s/Resources/gencode.v38lift37.annotation.gff3'
whole_genes = '/home/ar7343bo-s/whole_gene_list.txt'

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
                        if chromosome not in gene_dict:
                            gene_dict[chromosome] = {gene_name: (int(gff_line[3]), int(gff_line[4]))}
                            break
                        else:
                            gene_dict[chromosome][gene_name] = (int(gff_line[3]), int(gff_line[4]))
                            break

list_of_files = os.listdir("/home/ar7343bo-s/VcfTest/Annotation")
# list_of_files = os.listdir("/home/ar7343bo-s/GATKv3.8_HC/Annotation")

with open('number_var.txt', 'w') as out:
    print("Sample\tOnTargetPass\tOnTargetFail\tOffTargetPass\tOffTargetFail", file=out)
    for file in list_of_files:
        sample = file.lstrip('fs_filtered_').rstrip(".vcf")
        on_pass = 0
        on_fail = 0
        off_pass = 0
        off_fail = 0
        with open(file, 'r') as vcf:
            for line in vcf:
                if not line.startswith("#"):
                    line = line.split()
                    chrom = line[0]
                    pos = line[2]
                    filter_col = line[6]
                    if chrom in gene_dict:
                        for gene in gene_dict[chrom]:
                            start = gene_dict[chrom][gene][0]
                            end = gene_dict[chrom][gene][1]
                            if start < int(pos) < end:
                                if filter_col == 'PASS':
                                    on_pass += 1
                                else:  # Not PASS
                                    on_fail += 1
                            else:  # Not inside the gene
                                if filter_col == 'PASS':
                                    off_pass += 1
                                else:  # FAIL
                                    off_fail += 1
                    else:  # Chromosome that is not in the dict.
                        if filter_col == 'PASS':
                            off_pass += 1
                        else:  # FAIL
                            off_fail += 1

        print("{}\t{}\t{}\t{}\t{}".format(sample, on_pass, on_fail, off_pass, off_fail), file=out)
