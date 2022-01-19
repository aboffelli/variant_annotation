#!/bin/bash

# Script to perform the annotation in the vcf files using VEP. Installation of all necessary software is detailed in the README file.

# REMEMBER to activate conda env e_vep.

# Annotate the files using VEP and store in a new directory.
# Custom annotation to retrieve allele frequency from SweGen database.
# Custom annotation for conservation values from PhyloP and GERP.
mkdir Annotation
rsync -av -f"+ */" -f"- *" Data/ Annotation/

find Data/ -type f | while read file
do
     path_name=${file%/*}
     vep -i $file -o Annotation/${path_name#Data/}/vep_$(basename "$file") \
     --cache \
     --vcf \
     --offline \
     --force \
     --assembly GRCh37 \
     --distance 0 \
     --custom ~/Resources/SweGen/anon-SweGen_STR_NSPHS_1000samples_freq_hg19.vcf.gz,SweGen,vcf,exact,0,AF \
     --af \
     --af_1kg \
     --af_gnomad \
     --check_existing \
     --custom ~/Resources/All_hg19_RS.bw,GERP,bigwig,exact \
     --custom ~/Resources/hg19.100way.phyloP100way.bw,PhyloP,bigwig,exact \
     --fields "Gene","SYMBOL","Feature","Existing_variation","STRAND","EXON","INTRON","Consequence","Codons","AF","EUR_AF","SweGen_AF","gnomAD_AF","gnomAD_NFE_AF","PhyloP","GERP", \
     --no_stats
 done
