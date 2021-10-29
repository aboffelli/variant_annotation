#!/bin/bash

mkdir BindingSites

for f in *.vcf;
do
    vep -i $f -o BindingSites/binding_sites_$(basename "$f") \
    --vcf \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ~/Resources/EncodeFiles/encode_rna_binding_1.bed.gz,Encode,bed,overlap \
    --fields "Feature","Encode" \
    --keep_csq \
    --vcf_info_field Encode \
    --no_stats;
done
