#!/bin/bash

# Remember to activate conda e_vep first

# Get the flanking sequences of the SelectVariants
mkdir FsVcf
for file in HardFiltering/Merged/*.vcf;
do
    fill-fs -r /home/ar7343bo-s/.vep/homo_sapiens/104_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz -l 15 $file > FsVcf/fs_$(basename "$file");
done


# VEP command
mkdir Annotation
for file in FsVcf/*.vcf
do
     vep -i $file -o Annotation/$(basename "$file" .vcf)_VEP.vcf \
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
     --fields "Feature","Existing_variation","STRAND","EXON","INTRON","Consequence","Codons","AF","EUR_AF","SweGen_AF","gnomAD_AF","gnomAD_NFE_AF","PhyloP","GERP", \
     --no_stats
 done
