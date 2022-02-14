#!/bin/bash

# Script to add the Encode annotation of overlapping RNA binding protein binding sites on the vcfs.

# Create a new directory to store the annotated files and copy the directory
# structure from Annotation/.
mkdir EncodeAnnotation
rsync -av -f"+ */" -f"- *" Annotation/ EncodeAnnotation/

find Annotation/ -type f | while read file;
do
    path_name=${file%/*}
    vep -i $file -o EncodeAnnotation/${path_name#Annotation/}/encode_$(basename "$file") \
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
