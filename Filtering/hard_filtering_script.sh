#!/bin/bash

# Script to perform the Hard Filtering using GATK. The files must be divided since the values are different for indels and SNVs. The filter is applied to both and the files are merged together again.

# Create the new directories to store the results.
mkdir HardFiltering
cd HardFiltering
mkdir SNPs
mkdir Indels

# Separate the indels from the snps using the option SelectVariants from GATK.
for file in ../Data/*.vcf; do
  gatk SelectVariants -V $file -select-type SNP -O SNPs/snps_${file##*/} 2>>SNPs/snps_error;
  gatk SelectVariants -V $file -select-type INDEL -select-type MIXED -O Indels/indels_${file##*/} 2>>Indels/indels_error;
done

# Filter SNVs using the values sugested by the GATK best practices. With the exception of SOR, which we use > 4.0 instead of 3.0. Save the results in a new directory.
mkdir FilteredSNPs
for file in SNPs/snps_*.vcf; do
  gatk VariantFiltration -V $file \
  -filter "QD < 2.0 " --filter-name QD2S \
  -filter "MQ < 40.0" --filter-name MQ40 \
  -filter "FS > 60.0" --filter-name FS60 \
  -filter "SOR > 4.0" --filter-name SOR4 \
  -filter "MQRankSum < -12.5" --filter-name MQRS-12.5 \
  -filter "ReadPosRankSum < -8.0" --filter-name RPRS-8 \
  -O FilteredSNPs/filtered_${file##*/} 2>>FilteredSNPs/snps_filter_error;
 done

# Filter Indels using the values sugested by the GATK best practices. Save the results in a new directory.
mkdir FilteredIndels
for file in Indels/indels_*.vcf; do
  gatk VariantFiltration -V $file \
  -filter "QD < 2.0 " --filter-name QD2I  \
  -filter "FS > 200.0" --filter-name FS200 \
  -filter "SOR > 10.0" --filter-name SOR10 \
  -filter "ReadPosRankSum < -20.0" --filter-name RPRS-20 \
  -O FilteredIndels/filtered_${file##*/} 2>>FilteredIndels/indels_filter_error;
done

# Merge the Indels and SNVs together and store it in a new directory.
mkdir Merged
for file in FilteredIndels/*.vcf; do
  fname=${file#*indels_}; gatk MergeVcfs -I $file -I FilteredSNPs/filtered_snps_$fname -O Merged/filtered_$fname 2>>Merged/merge_error;
done
