#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:

Description: Script to parse the significant synonymous variants from the plink
    output file, and generate a table containing all the variant information.

Created on: 2022-04-25
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import time
import re
import glob
import sys

start_time = time.time()


def synonymous_parser(vcfline, posit):
    """
    Function that identifies a synonymous variant and retrieves the information
        needed for the table.

    :param vcfline: line from the vcf file
    :param posit: list containing the chromosome, position of the
                     variant and alternative base.
    :return: tuple containing the DbSNP ID, gene affected, the allele
        frequency, conservation value PhyloP, conservation value GERP,
        delta rscu, ese, ess, and encode information.
    """
    if posit not in synonymous_table:
        # Isolate only the annotation info, using regex.
        csq = re.search(r'CSQ=(\S+?),(\S+?)(;ClinVar|\s)', vcfline)
        # The fixed csq contains information that are the same for all
        # transcripts
        fixed_csq = csq.group(1).split('|')
        known = fixed_csq[0].split('&')[0]
        af = fixed_csq[3]
        phylop = fixed_csq[6]
        gerp = fixed_csq[7]

        # If the position has different alternative bases we loop through both.

        # Retrieve the transcript that contain the synonymous variant.
        synonymous_transcript = re.search(
            r'synonymous_variant\S*?\|-?\d.\d+\|', csq.group(2))

        # Since it's possible that the same position contains a different base that
        # is not a synonymous variant, we first check if we actually have a
        # synonymous transcript.
        if synonymous_transcript:
            synonymous_transcript = synonymous_transcript.group(0)

            # Remove variants with NMD transcript and missense variants.
            if 'NMD_transcript_variant' in synonymous_transcript:
                return False
            if 'missense_variant|' in csq.group(0):
                return False

            # Get the codon, rscu and encode information from the right
            # transcript (in cases of two alternative bases, there is one
            # transcript for each base.
            synonymous_transcript = synonymous_transcript.split('|')
            codon = synonymous_transcript[1]
            rscu = synonymous_transcript[3]
            encode_info = synonymous_transcript[2].split('&')
            encode = []

            # Since the encode information can contain more than one protein,
            # use a loop to get all proteins in a list and join them together.
            for prot in encode_info:
                encode.append(prot.split(':')[0])
            encode = ';'.join(encode)

            # It is possible that the same variant affects different genes, so
            # we loop through all transcripts and store all the gene names
            # and consequence in a set to keep them unique to join it after.
            all_transcripts = csq.group(0).strip(';ClinVar\t').split(',')[1:]

            genes = set()
            consequence = set()

            # The ese and ess can be present or not depending on the transcript,
            # however, they will be the same for all transcripts. So. after we
            # retrieve a value we can stop looking for it.
            ese = ''
            ess = ''
            ese_ess = ''

            for transcript in all_transcripts:
                transcript = transcript.split('|')
                genes.add(transcript[1])
                consequence.add(transcript[6])
                if 'synonymous' in transcript[6]:
                    if not ese_ess:
                        # There is a function to evaluate the ese and ess.
                        ese_ess = ese_ess_parser(transcript)

            # If any ese or ess was found separate them into different
            # variables.
            if ese_ess:
                ese, ess = ese_ess

            genes = genes & seq_genes
            genes = ';'.join(genes)
            consequence = ';'.join(sorted(list(consequence)))

            # Store all the info in a tuple.
            func_result = (known, genes, codon, af, phylop, gerp, rscu,
                           ese, ess, encode, consequence)

            # Store the result in the dictionary with the position of the
            # variant.
            synonymous_table[posit] = func_result

            # Add the consequence to a quality check file.
            if posit not in qc:
                qc[file_type][posit] = consequence

    # Calculate the allele frequency and store the result in the af
    # dictionary.
    if posit not in bridges_af[file_type]:
        bridges_af[file_type][posit] = {}
    if allele == '1/1':  # homozygous - two alleles
        allele_value = 1
        bridges_af[file_type][posit][sample] = allele_value

    else:  # heterozygous - one allele (0/1, 1/0 or 1/2)
        allele_value = 0.5
        bridges_af[file_type][posit][sample] = allele_value


def ese_ess_parser(transcript):
    """
    Function to evaluate if there was a change in the ese and ess hexamers,
        that can be defined as:
            alter - if the hexamer or number of hexamers changed.
            create - if no hexamer was found in the reference but is found in
                the altered.
            remove - if one or more hexamer is found in the reference but no
                hexamer is found in the altered.

    :param transcript: list containing the information of the transcript line.
    :return: tuple containing the results for ese and ess (alter, remove, or
        create).
    """
    ese_ess = {'ese': '', 'ess': ''}
    # ESE is found in the index 10 and 11, ESS is found in the index 12 and 13.
    for i in range(10, 14, 2):
        if i == 10:
            es_type = 'ese'
        else:
            es_type = 'ess'

        # If reference and alteration differ.
        if transcript[i] != transcript[i + 1]:
            reference = transcript[i]
            alteration = transcript[i + 1]

            # Check if both have values
            if reference and alteration:
                # Create a list for all the hexamers
                ref_list = reference.split(';')
                alt_list = alteration.split(';')

                # If the sizes differ, there was an alteration
                if len(ref_list) != len(alt_list):
                    ese_ess[es_type] = 'alter'

                # If the sizes are the same, compare all hexamers to see if
                # there is a difference.
                else:
                    for hexamer in ref_list:
                        if hexamer not in alt_list:
                            ese_ess[es_type] = 'alter'

            # If the reference is empty and the alteration is not.
            elif not reference and alteration:
                ese_ess[es_type] = 'create'

            # If the alteration is empty and the reference is not.
            else:  # if reference and not alteration.
                ese_ess[es_type] = 'remove'

    # Return the values only if one of them is not empty. This makes possible
    # to check the other transcripts if both are empty.
    ese, ess = ese_ess.values()
    if ese or ess:
        return ese, ess


# Retrieve the synonymous variants from the plink output. Since we have
# different possible outputs from plink, we can use sys to be able to use the
# same script.
try:
    plink_file_name = sys.argv[1]

    # Retrieve all the variants position in the respective map file
    # Use the plink file to get the map file name.
    map_file_name = re.search(r"plink_(\S+)\.assoc",
                              plink_file_name).group(1) + ".map"

    print("Reading the map file...")
    with open(map_file_name, 'r') as map_file:
        variant_pos = {}
        # Search the variant in the map file and save the positions in the
        # variant_pos dictionary in the format chr:position.
        for m_line in map_file:
            map_line = m_line.strip().split('\t')
            variant_pos[map_line[1]] = f"{map_line[0]}:{map_line[2]}"
        print("Done!")

    print("\nRetrieving the significant variants...")
    # Go through the plink file and get only the significant variants.
    with open(plink_file_name, 'r') as plink_file:
        significant_variants = {}

        # Ignore the header.
        header = plink_file.readline()

        # Since the plink file is separated with unpatterned single spaces and
        # not tabs(YES! WHO DOES THAT?), we have to fix it first.
        for p_line in plink_file:
            # Change all one or more spaces to tabs.
            plink_line = re.sub(r'\s+', '\t', p_line.strip())
            plink_line = plink_line.split('\t')

            # Since the plink output has different types of correction, choose
            # the column containing the wanted correction.
            p_value = plink_line[8]

            # If we have a significant value, store it in the dictionary.
            if (p_value == 'INF') or (float(p_value) < 0.05):
                snv_id = plink_line[1]
                position = variant_pos[snv_id]
                significant_variants[position] = snv_id
        print("Done!")

# Exception for missing argument
except IndexError:
    print("Missing file name.")
    exit()

# Exception for file exists.
except FileNotFoundError:
    print("Association file not found.")
    exit()

# Exception for map file.
except AttributeError:
    print("Map file not found.")
    exit()

# Initiate the three dictionaries.
synonymous_table = {}

bridges_af = {"Cases": {},
              "Controls": {}}

bridges_af_perc = {"Cases": {},
                   "Controls": {}}

# QC checking
qc = {"Cases": {},
      "Controls": {}}

seq_genes = set()
# Load the list of genes used in the sequencing.
print('\nLoading targeted genes...')
with open("/home/ar7343bo-s/Resources/BRIDGES_fluidigm_juno_panel_primers.txt",
          'r') as primers:
    # Remove the header.
    header = primers.readline()

    for line in primers:
        # Get gene from the line.
        exon = line.split('\t')[4]
        gene = re.search(r"^(\w+?)_", exon).group(1)
        seq_genes.add(gene)
print("Done!")

# Prepare the VCF files that will be used.
print("\nReading files...")
list_of_files = glob.glob("./**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    # Get the type (Case or Control) from the file path.
    file_type = re.search(r"/(\w+)/(?:bridges|inter)", file).group(1)

    # Isolate the sample name.
    sample = re.search(r'vep_(\S+).raw', file).group(1)

    with open(file, 'r') as vcf:
        for vcfline in vcf:

            # Ignore the comment lines.
            if not vcfline.startswith('#'):
                split_line = vcfline.split('\t')

                # If the variant did not pass the filters ignore it.
                filt = split_line[6]
                if filt != 'PASS':
                    continue

                # Get the position as chr:position and remove the chr since the
                # map file does not contain it.
                position = ':'.join(split_line[0:2]).lstrip('chr')

                if position in significant_variants:
                    # Get the alternative base and add it to the end of the
                    # position.
                    alt_base = vcfline.split('\t')[4]
                    position += alt_base
                    allele = vcfline.split('\t')[-1][0:3]

                    synonymous_parser(vcfline, position)

    # Raise the file count
    file_count += 1
print("\nDone!")

print("\nCalculating the allele frequency...")
# Loop through the bridges dictionary to calculate the percentage of
for dict_type in bridges_af:

    # The number of samples is different depending on the type.
    if dict_type == "Cases":
        total_num = 60239

    else:  # Controls
        total_num = 53306

    for position in bridges_af[dict_type]:

        # Start the value for the current variant
        bridges_af_perc[dict_type][position] = 0

        # Sum all the frequencies from the samples.
        for sample in bridges_af[dict_type][position]:
            bridges_af_perc[dict_type][
                position] += bridges_af[dict_type][position][sample]

        # Divide the sum by the total number of samples and store it in
        # the dictionary. Add the number of samples that had the variant.
        bridges_af_perc[dict_type][position] = (
            f"""{bridges_af_perc[dict_type][position] /
                 total_num:.6f} ({
            len(bridges_af[dict_type][position])} samples)""")
print("Done!")

# Save all the information in a text file tab delimited.
print("\nWriting the output file...")
output_name = re.search(r"bridges_(\w+)", plink_file_name).group(1)
with open(f'{output_name}_table.txt', 'w') as outfile:
    # Add the header
    print("#Variant\tAlternative_base\tDbSNP_ID\tGene\tCodon_(ref/alt)\t"
          "AF_SweGen\tAF_BRIDGES_Control\tAF_BRIDGES_Cases\tPhyloP\tGERP\t"
          "deltaRSCU\tESE\tESS\tRBP\tConsequence",
          file=outfile)

    for variant in sorted(synonymous_table,
                          key=lambda x: (int(x.split(':')[0]),
                                         int(x.split(':')[1][:-1]))):
        # Insert the Bridges allele frequencies in the result list. Two minus
        # the header position since we don't have the position in this list.
        result = list(synonymous_table[variant])
        if variant in bridges_af_perc["Controls"]:
            result.insert(4, bridges_af_perc["Controls"][variant])
        else:
            result.insert(4, 'NA')

        if variant in bridges_af_perc["Cases"]:
            result.insert(5, bridges_af_perc["Cases"][variant])
        else:
            result.insert(5, 'NA')

        # Change any empty space to NA
        for j in range(len(result)):
            if not result[j]:
                result[j] = 'NA'

        # Insert the line into the file. Separate the base letter from the
        # position.
        print("{}\t{}\t{}".format(variant[:-1], variant[-1], '\t'.join(result)),
              file=outfile)


# Create the qc file.
with open('qc_file.txt', 'w') as qc_file:
    for file_type in qc:
        for pos in qc[file_type]:
            print(file_type, pos, qc[file_type][pos], file=qc_file, sep='\t')
print("Done!")

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
