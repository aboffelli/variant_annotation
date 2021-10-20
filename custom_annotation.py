import re


def reverse_complement(sequence):
    """
    Function that creates the reverse complement of a DNA sequence.

    :param sequence: sequence of DNA
    :return: reverse complement
    """
    rev_comp = ''.join(complement_dictionary.get(base, base) for base in reversed(sequence))
    return rev_comp


# Complemnt dictionary
complement_dictionary = {'A': 'T', 'C': 'G',
                         'T': 'A', 'G': 'C'}
# Load RSCU dictionary
rscu_dictionary = {"TTT": 0.92, "TTC": 1.08, "TTA": 0.46, "TTG": 0.77, "CTT": 0.80, "CTC": 1.16, "CTA": 0.43,
                   "CTG": 2.38, "ATT": 1.09, "ATC": 1.39, "ATA": 0.52, "ATG": 1.00, "GTT": 0.73, "GTC": 0.95,
                   "GTA": 0.47, "GTG": 1.85, "TCT": 1.12, "TCC": 1.30, "TCA": 0.92, "TCG": 0.33, "CCT": 1.15,
                   "CCC": 1.30, "CCA": 1.09, "CCG": 0.46, "ACT": 1.00, "ACC": 1.41, "ACA": 1.15, "ACG": 0.45,
                   "GCT": 1.05, "GCC": 1.61, "GCA": 0.91, "GCG": 0.44, "TAT": 0.89, "TAC": 1.11, "TAA": 0.84,
                   "TAG": 0.66, "CAT": 0.84, "CAC": 1.16, "CAA": 0.53, "CAG": 1.47, "AAT": 0.95, "AAC": 1.05,
                   "AAA": 0.88, "AAG": 1.12, "GAT": 0.93, "GAC": 1.07, "GAA": 0.85, "GAG": 1.15, "TGT": 0.92,
                   "TGC": 1.08, "TGA": 1.50, "TGG": 1.00, "CGT": 0.48, "CGC": 1.11, "CGA": 0.65, "CGG": 1.23,
                   "AGT": 0.90, "AGC": 1.44, "AGA": 1.26, "AGG": 1.27, "GGT": 0.64, "GGC": 1.35, "GGA": 1.00,
                   "GGG": 1.01}

# Load the ESS and ESE hexamers from the tables
ese_table = '/Users/student/Documents/ExampleData/RESCUE-ESE_hexamers_200703.txt'
ess_table = '/Users/student/Documents/ExampleData/ESS_hexamers_200824.txt'

ese_set = set()
ess_set = set()

with open(ese_table) as ese, open(ess_table) as ess:
    for line in ese:
        ese_set.add(line.strip())

    for line in ess:
        ess_set.add(line.strip())

# TODO: add this line to the headers, after '##INFO=<ID=ExcessHet'
#   ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
# TODO: Change flanking sequence (FS) to something else, including in the header.

# TODO: Change CSQ line in the header '##INFO=<ID=CSQ'


# Only look for ESEs and ESSs for exonic variants
#   Reverse comp if strand -1

# Round PhyloP and GERP to four numbers

# Calc deltaRSCU for synonymous variants

# Write the new line to the new vcf


with open('/Users/student/Documents/test.vcf') as vcf:
    for line in vcf:
        if line.startswith("#"):
            pass
        else:
            split_line = line.split()
            line_info = split_line[7]
            flanking_seq = re.search("FS=(\S[A-Z]+\[.+\/.+\][A-Z]+)", line_info).group(1)
            transcripts = re.search('CSQ=(\S.+)', line_info).group(1).split(',')
            for transcript in transcripts.copy():
                csq = transcript.split('|')
                transcripts.append('a')
                for x in enumerate(csq): print(x)