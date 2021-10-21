import re
import platform


def reverse_complement(sequence):
    """
    Function that creates the reverse complement of a DNA sequence.

    :param sequence: sequence of DNA
    :return: reverse complement
    """
    rev_comp = ''.join(
        complement_dictionary.get(base, base) for base in reversed(sequence))
    return rev_comp


if platform.system() == 'Windows':
    test_file = 'C:/Users/Arthu/Box/test.vcf'
    ese_file = 'C:/Users/Arthu/Box/RESCUE-ESE_hexamers_200703.txt'
    ess_file = 'C:/Users/Arthu/Box/ESS_hexamers_200824.txt'

else:
    test_file = '/Users/student/Documents/test.vcf'
    ese_file = '/Users/student/Documents/ExampleData/RESCUE-ESE_hexamers_200703.txt'
    ess_file = '/Users/student/Documents/ExampleData/ESS_hexamers_200824.txt'


# Complement dictionary
complement_dictionary = {'A': 'T', 'C': 'G',
                         'T': 'A', 'G': 'C',
                         '-': '-'}
# Load RSCU dictionary
rscu_dictionary = {"TTT": 0.92, "TTC": 1.08, "TTA": 0.46, "TTG": 0.77,
                   "CTT": 0.80, "CTC": 1.16, "CTA": 0.43,
                   "CTG": 2.38, "ATT": 1.09, "ATC": 1.39, "ATA": 0.52,
                   "ATG": 1.00, "GTT": 0.73, "GTC": 0.95,
                   "GTA": 0.47, "GTG": 1.85, "TCT": 1.12, "TCC": 1.30,
                   "TCA": 0.92, "TCG": 0.33, "CCT": 1.15,
                   "CCC": 1.30, "CCA": 1.09, "CCG": 0.46, "ACT": 1.00,
                   "ACC": 1.41, "ACA": 1.15, "ACG": 0.45,
                   "GCT": 1.05, "GCC": 1.61, "GCA": 0.91, "GCG": 0.44,
                   "TAT": 0.89, "TAC": 1.11, "TAA": 0.84,
                   "TAG": 0.66, "CAT": 0.84, "CAC": 1.16, "CAA": 0.53,
                   "CAG": 1.47, "AAT": 0.95, "AAC": 1.05,
                   "AAA": 0.88, "AAG": 1.12, "GAT": 0.93, "GAC": 1.07,
                   "GAA": 0.85, "GAG": 1.15, "TGT": 0.92,
                   "TGC": 1.08, "TGA": 1.50, "TGG": 1.00, "CGT": 0.48,
                   "CGC": 1.11, "CGA": 0.65, "CGG": 1.23,
                   "AGT": 0.90, "AGC": 1.44, "AGA": 1.26, "AGG": 1.27,
                   "GGT": 0.64, "GGC": 1.35, "GGA": 1.00,
                   "GGG": 1.01}

# Load the ESS and ESE hexamers from the tables
ese_set = set()
ess_set = set()

with open(ese_file) as ese, open(ess_file) as ess:
    for line in ese:
        ese_set.add(line.strip())

    for line in ess:
        ess_set.add(line.strip())


# TODO: Change flanking sequence (FS) to something else, including in the header.
with open(test_file) as vcf:
    for line in vcf:
        # Do header stuff
        if line.startswith("#"):
            # TODO: add this line to the headers, after '##INFO=<ID=ExcessHet'
            #   ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using
            #   Fisher's exact test to detect strand bias">

            # TODO: Change CSQ line in the header '##INFO=<ID=CSQ'
            #   New CSQ Feature|Existing_variant|STRAND|EXON|INTRON|Consequence|Codons|AF|
            #   EUR_AF|SweGen_AF|gnomAD_AF|gnomAD_NFE_AF|PhyloP|GERP|delta_RSCU|ESEs_REF|
            #   ESEs_ALT|ESSs_REF|ESSs_ALT"
            pass

        # Non-header stuff
        else:
            split_line = line.split()
            line_info = split_line[7]
            flanking_seq = re.search(r"FS=(\S[A-Z]+\[.+\/.+\][A-Z]+)", line_info
                                     ).group(1)
            transcripts = re.search(r'CSQ=(\S.+)', line_info).group(1).split(
                ',')

            for index, transcript in enumerate(transcripts.copy()):
                csq = transcript.split('|')
                strand = csq[2]
                exon = csq[3]
                consequence = csq[5]
                codon = csq[6].upper()
                delta_rscu = ''
                ref_motifs_ese = []
                ref_motifs_ess = []
                alt_motifs_ese = []
                alt_motifs_ess = []

                # Round PhyloP and GERP to four numbers
                if csq[12]:
                    csq[12] = str(round(float(csq[12]), 4))
                if csq[13]:
                    csq[13] = str(round(float(csq[13]), 4))

                # Calculate the deltaRSCU for synonymous variants:
                if 'synonymous_variant' in consequence:
                    ref_codon = codon[:3]
                    alt_codon = codon[-3:]
                    delta_rscu = round(rscu_dictionary[alt_codon]
                                       - rscu_dictionary[ref_codon], 2)

                # Check for ESE and ESS in exonic variants
                if exon:
                    # get surrounding sequence
                    # TODO: Ask how to treat indels
                    seq = re.search('(\w{5})\[(\S+)\/(\S+)\](\w{5})',
                                    flanking_seq)
                    ref_seq = ''.join(seq.group(1, 2, 4))
                    alt_seq = ''.join(seq.group(1, 3, 4))

                    # Get reverse complement if reverse strand
                    if strand == '-1':
                        ref_seq = reverse_complement(ref_seq)
                        alt_seq = reverse_complement(alt_seq)

                    # get hexamers matches for ref ESE and ESS
                    for i in range(len(ref_seq) - 5):
                        hexamer = ref_seq[i:i+6]
                        if hexamer in ese_set:
                            ref_motifs_ese.append(hexamer)
                        if hexamer in ess_set:
                            ref_motifs_ess.append(hexamer)

                    # get hexamers matches for alt ESE and ESS
                    for i in range(len(alt_seq) - 5):
                        hexamer = alt_seq[i:i+6]
                        if hexamer in ese_set:
                            alt_motifs_ese.append(hexamer)
                        if hexamer in ess_set:
                            alt_motifs_ess.append(hexamer)

                # Append the new info to csq
                csq.append(str(delta_rscu))
                csq.append(','.join(ref_motifs_ese).strip(','))
                csq.append(','.join(alt_motifs_ese).strip(','))
                csq.append(','.join(ref_motifs_ess).strip(','))
                csq.append(','.join(alt_motifs_ess).strip(','))
                csq = '|'.join(csq)
                transcripts[index] = csq

            # Join the line together again
            transcripts = ','.join(transcripts)
            split_line[7] = re.sub(r'CSQ=(\S.+)', 'CSQ=' + transcripts,
                                   line_info)
            print(line.strip())
            line = '\t'.join(split_line)
            print(line + '\n')
            # TODO: Write the line to the output
