#!/bin/bash

# Create a tab delimited table with the number of variants removed with each filter for each file.

echo -e "Sample\tQD2S\tMQ40\tFS60\tSOR3\tMQRS-12.5\tRPRS-8\tQD2I\tFS200\tSOR10\tRPRS-20" > filters_tab.txt

for file in HardFiltering/Merged/*.vcf; do
    sample=$(basename "$file" .vcf);
    line=$(
        grep -v "^#" $file | grep -c "QD2S";
        grep -v "^#" $file | grep -c "MQ40";
        grep -v "^#" $file | grep -c "FS60";
        grep -v "^#" $file | grep -c "SOR3";
        grep -v "^#" $file | grep -c "MQRS-12.5";
        grep -v "^#" $file | grep -c "RPRS-8";
        grep -v "^#" $file | grep -c "QD2I";
        grep -v "^#" $file | grep -c "FS200";
        grep -v "^#" $file | grep -c "SOR10";
        grep -v "^#" $file | grep -c "RPRS-20");
    echo ${sample} $line | tr ' ' '\t' >> filters_tab.txt;
done

# Create a table with the total number of variants that Passed or not with the hard filtering and machine learning approach.

echo -e "Filtering_Type\tPass\tFail" > hf-ml_comparison.txt

for file in HardFiltering/Merged/*.vcf; do
    line=$(
        grep -v "#" $file | grep -c "PASS";
        grep -v "#" $file | grep -v "PASS" | wc -l);
    echo "HF" $line | tr " " "\t" >> hf-ml_comparison.txt;
done

for file in VQSR/Filtered/*.vcf; do
    line=$(
        grep -v "#" $file | grep -c "PASS";
        grep -v "#" $file | grep -v "PASS" | wc -l);
    echo "ML" $line | tr " " "\t" >> hf-ml_comparison.txt;
done
