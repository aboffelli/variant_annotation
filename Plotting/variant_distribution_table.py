import re
import os


def check_type(dict_key, consequence):
    global variant_types_known, variant_types_novel
    """
    Function to add the the count for the right consequence.

    :param consequence: string containing the consequence of the variant
    :return:
    """
    for key in dict_key.split('/'):
        if key in consequence:
            if exist:
                variant_types_known[dict_key] += 1
            else:
                variant_types_novel[dict_key] += 1
            return True


def create_gene_dict(gene_name):
    """
    Function to populate the gene_dict.
    :param gene_name:
    :return: 
    """
    global known_gene_dict, novel_gene_dict
    if gene_name not in known_gene_dict:
        known_gene_dict[gene_name] = variant_types_known.copy()
        novel_gene_dict[gene_name] = known_gene_dict.copy()


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
                    transcripts = csq.split(',')[1:]
                    for transcript in transcripts:
                        gene = transcript.split('|')[1]
                        create_gene_dict(gene)
                        conseq = transcript.split('|')[6]
                        for keys in variant_types_known:
                            if check_type(keys, conseq):
                                break
                        else:
                            if exist:
                                variant_types_known['other'] += 1
                            else:
                                variant_types_novel['other'] += 1

# TODO: Create table of distribution of genes.

# Create the files
with open('known_variant_distribution.txt', 'w') as known, \
        open('novel_variant_distribution.txt', 'w') as novel:
    for i in variant_types_known:
        print(i + '\t' + variant_types_known[i], file=known)
        print(i + '\t' + variant_types_novel[i], file=novel)

with open('genes_known_variant_distribution.txt', 'w') as gene_known, \
        open('genes_novel_variant_distribution.txt', 'w') as gene_novel:
    for i in known_gene_dict:
        for j in known_gene_dict[i]:
            print(f'{i}\t{j}\t{known_gene_dict[i][j]}', file=gene_known)
            print(f'{i}\t{j}\t{novel_gene_dict[i][j]}', file=gene_novel)
