#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Plink input files.

Description: Script that creates two input files for PLINK.
    One tab delimited MAP file containing the information for each variant, i.e.
    chromosome number, snp ID, and position.
    One tab delimited PED file containing the information for each sample, i.e.
    sample ID, sex, phenotype, and genotype.

Created on: 2022-04-06
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation
"""

import time
import re
import glob
import sys

start_time = time.time()


def qc_check(vcf_line):
    """
    Function to select the variants that will continue for the analysis. The
    quality is based on the filter, if it is not and indel, and if it is inside
    the primers target.
    :param vcf_line: string containing the variant vcf line.
    :return: True if the variant passed, False if the variant did not pass.
    """
    global fail_count, indel_count, off_target_count

    # Get all information that we need from the line.
    vcf_split = vcf_line.split('\t')
    chrom = vcf_split[0][3:]
    position = vcf_split[1]
    filt = vcf_split[6]
    bases = ''.join(vcf_split[3:5])

    # Create the passing variable.
    qc_pass = True

    # The positive control is a deletion, so we remove one of the bases to
    # pass the qc.
    if rs == "positive_control":
        bases = bases[1:]

    # Test for the filter and snv.
    if filt != "PASS":
        fail_count += 1
        qc_pass = False

    if len(bases) != 2:
        indel_count += 1
        qc_pass = False

    # Test for on target variant.
    target = False
    # Since the same chromosome has several primers, we have to check all of
    # them.
    for primer in primer_positions[chrom]:
        if primer[0] <= position <= primer[1]:
            target = True

    if not target:
        off_target_count += 1
        qc_pass = False

    return qc_pass


# Dictionary containing the primer positions.
primer_positions = {}

# Dictionary to store the variants information for the map file.
map_dict = {}

# Dictionary with the sample information for the ped file.
sample_dict = {}

# Dictionary to save the ref base and alt base for all variants.
variant_dict = {}

# Dictionary to save the number of samples that containing each variant.
variant_samples = {}

# Create a dictionary for novel variants
novel_variants_id = {}

# Variables to show the number of variants removed in the end.
fail_count = 0
indel_count = 0
off_target_count = 0
removed_variants = 0
n_samples = 0
no_fam_samples = 0
no_synonymous_count = 0
pathogenic_count = 0

usage = ("""Usage: plink_input.py [N] [--family] [--synonymous] [-h]

    N: positional argument integer with the cut-off number for the number of \
samples containing a variant
    --family: Boolean to remove samples that do not have family history
    --synonymous: Boolean to remove non-synonymous variants
    -h --help: Prints this help message and exits.""")

# Number to exclude variants based on the number of samples.

try:
    cut_off = int(sys.argv[1])

except IndexError:
    cut_off = 0

except ValueError:
    print(usage)
    exit()

# Check the arguments passed in the command line.
if "--help" in sys.argv or "-h" in sys.argv:
    print(usage)
    exit()

if '--family' in sys.argv:
    fam_hist = True
else:
    fam_hist = False

if '--synonymous' in sys.argv:
    synonymous = True
else:
    synonymous = False

if "--pathogenic" in sys.argv:
    pathogenic = True
    pathogenic_set = set()
else:
    pathogenic = False

# Load the primer positions from the file.
print('Loading the primer positions...')
with open("/home/ar7343bo-s/Resources/BRIDGES_fluidigm_juno_panel_primers.txt",
          'r') as primers:
    # Remove the header.
    header = primers.readline()

    for line in primers:
        # Get the chromosome, start, and end of the primers.
        p_chrom, p_start, p_end = line.strip().split('\t')[-3:]
        if p_chrom not in primer_positions:
            primer_positions[p_chrom] = []
        # Save the starts and ends as tuples, so we can access them later.
        primer_positions[p_chrom].append((p_start, p_end))
print("Done!")

# If the family history flag is activated, load the sample names from the file,
# and store it in a set.
if fam_hist:
    fam_hist_samples = set()
    print("\nLoading cases samples with family history...")
    with open("/home/ar7343bo-s/BRIDGES/raidset/FilteredClinVar/Cases/"
              "family_history_cases.txt", 'r') as family:
        for line in family:
            fam_hist_samples.add(line.split('\t')[0])
    print("Done!")

# Prepare the VCF files that will be used.
print("\nReading files...")
list_of_files = glob.glob("./**/*.vcf", recursive=True)
for file in list_of_files.copy():
    if '.vcf' not in file:
        list_of_files.remove(file)

# The novel variant will be named as varN. Every time a novel variant is found
# we add 1 to this variable.
var_num = 1
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)
    # Get the phenotype, sample name and sex for the ped file.
    # The phenotype is Case or Control, for now.
    if file.split('/')[1] == "Controls":
        phenotype = "1"
    else:
        phenotype = "2"
    sample_name = re.search(r"vep_(\S+)\.raw", file).group(1)
    sex = '2'

    # If the family history flag is activated, ignore the cases that are not in
    # the set containing the samples with family history.
    if fam_hist:
        if phenotype == '2' and sample_name not in fam_hist_samples:
            no_fam_samples += 1
            file_count += 1
            continue

    sample_dict[sample_name] = [sex, phenotype, {}]
    with open(file, 'r') as vcf_file:
        for vcfline in vcf_file:
            # Exclude header lines.
            if not vcfline.startswith("#"):
                split_line = vcfline.split('\t')

                # If the synonymous flag is activated, remove variants that are
                # not synonymous variants.
                if synonymous:
                    if 'synonymous_variant' not in vcfline:
                        no_synonymous_count += 1
                        continue

                # Get the reference and alternative base.
                ref_base = split_line[3].upper().strip()
                alt_base = split_line[4].upper().strip()

                # Get the dbSNP id.
                rs = re.search(r"CSQ=(?:(\w+)|\|)", vcfline).group(1)

                # Use one of the known variants as a positive control. We are
                # using the 1100delC deletion as the positive control, since we
                # know that this variant is related to breast cancer.
                if rs == "rs555607708":
                    rs = "positive_control"
                    ref_base = ref_base[1]

                # Check if the variant passes in all criteria.
                if qc_check(vcfline):

                    # Check if the variant is pathogenic after the qc, since we
                    # want to remove samples that have real pathogenic variants.
                    if pathogenic:
                        if "Pathogenic" in vcfline \
                                or "Likely_pathogenic" in vcfline:
                            pathogenic_set.add(sample_name)
                            continue

                    # Retrieve chr number and position for the map file.
                    chrom = split_line[0][3:]
                    if chrom == 'X':
                        chrom = '23'
                    elif chrom == 'Y':
                        chrom = '24'
                    position = split_line[1]

                    # Assign a name for the variant if it is a novel variant.
                    # In this case we're using varN.
                    if not rs:
                        chrom_pos = chrom + ':' + position
                        if chrom_pos not in novel_variants_id:
                            novel_variants_id[chrom_pos] = 'var' + str(var_num)
                            var_num += 1
                        rs = novel_variants_id[chrom_pos]

                    # Save the map info in the map dictionary.
                    if rs not in map_dict:
                        map_dict[rs] = (chrom, position)

                    if rs not in variant_dict:
                        variant_dict[rs] = (ref_base, alt_base)

                    # Allele fraction to differentiate between hetero and
                    # homozygous.
                    af = float(re.search(r"AF=(\d\.\d+)", vcfline).group(1))

                    # If allele fraction is 0.8 or higher, consider it
                    # homozygous alt, save the alternative base twice.
                    if af >= 0.8:
                        sample_dict[sample_name][2][rs] = [alt_base, alt_base]

                    # If allele fraction is between 0.2 and 0.8, consider it
                    # heterozygous, save the reference and alternative bases.
                    elif 0.2 <= af < 0.8:
                        sample_dict[sample_name][2][rs] = [ref_base, alt_base]

                    # Add to the dictionary for counting the number of samples.
                    if rs not in variant_samples:
                        variant_samples[rs] = set()
                    variant_samples[rs].add(sample_name)

                else:   # If qc_check returned False.
                    removed_variants += 1
    file_count += 1
print("\nDone!")

with open("number_samples.txt", "w") as out:
    for variant in variant_samples:
        print(f"{variant}\t{len(variant_samples[variant])}", file=out)

print("\nPreparing the map file...")
# For printing the map file, first prepare the lines into a list.
map_list = []

# Set a cut off value for the minimum number of patients with the variant.
print(f"Excluding variants that are present in less then {cut_off} samples...")

# Append the map lines already in the right format if the number of samples
# is enough.
for rs in map_dict:
    if len(variant_samples[rs]) > cut_off:
        chrom, pos = map_dict[rs]
        map_list.append(f"{chrom}\t{rs}\t{pos}")

    else:   # Variants that failed the minimum number of samples.
        removed_variants += 1
        n_samples += 1

# Then sort the list based on the chromosome number and position, and print to
# the file.
map_list = sorted(map_list, key=lambda x: (int(x.split('\t')[0]),
                                           int(x.split('\t')[2])))
print("Done!")

# Set the output file name.
output = f"bridges_filt{cut_off}"
if fam_hist:
    output = f"{output}_family"

if synonymous:
    output = f"{output}_synonymous"

if pathogenic:
    output = f"{output}_pathogenic"

# Write the .map file.
print("\nWriting the .map file...")
with open(f"{output}.map", "w") as outfile:
    for item in map_list:
        print(item, file=outfile)
print("Done!")

# Now we can use the map_list to get the same order of variants in both files.
print("\nWriting the .ped file...")
# Write all the information in the ped file.
with open(f"{output}.ped", 'w') as outfile:
    for sample in sample_dict:
        if sample in pathogenic_set:
            pathogenic_count += 1
            continue
        # Get the sex and phenotype for the sample.
        outfile.write(f"{sample}\t" + '\t'.join(sample_dict[sample][0:2]))

        # Using the order of the map_list add the bases for each variant in the
        # sample line.
        for var in map_list:
            # Use the variant id.
            var = var.split('\t')[1]

            # Get the bases for this particular variant if they are in the
            # dictionary.
            if var in sample_dict[sample][2]:
                outfile.write("\t" + '\t'.join(sample_dict[sample][2][var]))

            # If the sample does not have that variant, we assume that it is
            # homozygous reference, so we append the reference base twice.
            else:
                outfile.write("\t" + '\t'.join(variant_dict[var][0] * 2))
        outfile.write('\n')
print("Done!")

# Print the number of variants that were removed during the process on the
# screen.
print(f"\nA total of {removed_variants} variants were removed.")
print(f"{fail_count} variants did not pass the filters.")
print(f"{indel_count} variants were indels.")
print(f"{off_target_count} variants were off the primer targets.")
print(f"{n_samples} variants were found only in less than {cut_off} samples.")
if fam_hist:
    print(f"{no_fam_samples} samples did not have family history.")
if synonymous:
    print(f"{no_synonymous_count} variants were not synonymous.")
if pathogenic:
    print(f"{pathogenic_count} samples had a pathogenic variant.")

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
