mkdir GeneName

for f in *.vcf; do
    vep -i $f -o GeneName/gene_name_$(basename "$f") \
    --cache \
    --vcf \
    --offline \
    --assembly GRCh37 \
    --no_stats \
    --fields "Gene","SYMBOL" \
    --force \
    --keep_csq \
    --vcf_info_field Gene;
done
