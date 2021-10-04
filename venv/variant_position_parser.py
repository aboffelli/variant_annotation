#!~/miniconda3/bin/python3

"""
Script to parse and compare and identify in which gene each variant is located.
"""
home_dir = '/Users/student/'


# Retrieve the positions of whole genes the genes that where sequenced
gff = home_dir + 'Documents/gencode.v38lift37.annotation.gff3'
whole_genes = home_dir + 'Documents/whole_gene_list.txt'

gl = list()
gene_dict = dict()
with open(gff, 'r') as resource, open(whole_genes, 'r') as gene_list:
    for line in gene_list:
        gl.append(line.strip())
    for inf in resource:
        if not inf.startswith("#"):
            inf = inf.strip().split('\t')
            if inf[2] == "gene":
                gene_name = inf[8].split(';')[3].strip('gene_name=')
                for gene in gl:
                    if gene == gene_name:
                        if inf[0] not in gene_dict:
                            gene_dict[inf[0]] = [(gene_name, int(inf[3]), int(inf[4]))]
                            break
                        else:
                            gene_dict[inf[0]].append((gene_name, int(inf[3]), int(inf[4])))
                            break

print(gene_dict)
for key in gene_dict:
    print(key, gene_dict[key])
