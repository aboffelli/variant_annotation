import re
import os


def check_type(transcript_line, gene_name, existence):
    global variant_types_known, variant_types_novel
    global known_gene_dict, novel_gene_dict
    """
    Function to add the the count for the right consequence.

    :param consequence: string containing the consequence of the variant
    :return:
    """
    for dict_key in variant_types_known:
        for key in dict_key.split('/'):
            if key in transcript_line:
                if existence:
                    variant_types_known[dict_key] += 1
                    known_gene_dict[gene_name][dict_key] += 1

                else:
                    variant_types_novel[dict_key] += 1
                    novel_gene_dict[gene_name][dict_key] += 1
                return True


def create_gene_dict(gene_name):
    """
    Function to populate the gene_dict.
    :param gene_name:
    :return: 
    """
    global known_gene_dict, novel_gene_dict
    if gene:
        if gene_name not in known_gene_dict:
            known_gene_dict[gene_name] = variant_types_known.copy()
            novel_gene_dict = known_gene_dict.copy()


targeted_genes = set()
with open('/home/ar7343bo-s/whole_gene_list.txt', 'r') as gene_list:
    for line in gene_list:
        targeted_genes.add(line.strip())

variant_types_known = {'stop_gained': 0,
                       'frameshift': 0,
                       'splice_acceptor/donor': 0,
                       'missense': 0,
                       'inframe_deletion/insertion': 0,
                       'splice_region': 0,
                       'synonymous_variant': 0,
                       '5_prime_UTR': 0,
                       '3_prime_UTR': 0,
                       'intron': 0,
                       'intergenic': 0,
                       'other': 0
                       }
variant_types_novel = variant_types_known.copy()

known_gene_dict = {}
novel_gene_dict = {}

files_directory = "CustomAnnotation/"
list_of_files = os.listdir(files_directory)
for file in list_of_files.copy():
    if '.vcf' in file:
        with open(files_directory + file, 'r') as vcf_file:
            for line in vcf_file:
                if not line.startswith("#"):
                    filt = line.split()[6]
                    if filt != 'PASS':
                        continue
                    csq = re.search(r'CSQ=(\S+)', line).group(1)
                    exist = csq.split('|')[0]
                    for gene in targeted_genes:
                        if gene in csq:
                            create_gene_dict(gene)
                            check_type(csq, gene, exist)


# Create the files
with open('known_variant_distribution.txt', 'w') as known, \
        open('novel_variant_distribution.txt', 'w') as novel:
    for i in variant_types_known:
        print(f'{i}\t{variant_types_known[i]}', file=known)
        print(f'{i}\t{variant_types_novel[i]}', file=novel)

with open('genes_known_variant_distribution.txt', 'w') as gene_known, \
        open('genes_novel_variant_distribution.txt', 'w') as gene_novel:
    for i in known_gene_dict:
        for j in known_gene_dict[i]:
            print(f'{i}\t{j}\t{known_gene_dict[i][j]}', file=gene_known)
            print(f'{i}\t{j}\t{novel_gene_dict[i][j]}', file=gene_novel)
