#!/bin/bash

# Remember to activate conda splice first
mkdir MMSplice
for file in CustomAnnotation/*.vcf
do
   vep -i $file -o MMSplice/mmsplice_$(basename "$file") \
   --vcf \
   --format vcf \
   --cache \
   --force \
   --assembly GRCh37 \
   --plugin MMSplice \
   --vcf_info_field MMSplice \
   --fasta ~/.vep/homo_sapiens/104_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
   --fields "Feature","mmsplice_delta_logit_psi" \
   --keep_csq \
   --no_stats;
 done
