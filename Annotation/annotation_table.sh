#!/bin/bash

# Get how many novel and known variants were filtered.
for dir in *Fail/; do
    filename=${dir##/};
    echo -e "Sample\tQD2S\tMQ40\tFS60\tSOR4\tMQRS-12.5\tRPRS-8\tQD2I\tFS200\tSOR10\tRPRS-20" > ${filename%/}_filters_tab.txt

    for file in $dir/*.vcf; do
        sample=$(basename "$file" .vcf);
        line=$(
            grep -v "^#" $file | grep -c "QD2S";
            grep -v "^#" $file | grep -c "MQ40";
            grep -v "^#" $file | grep -c "FS60";
            grep -v "^#" $file | grep -c "SOR4";
            grep -v "^#" $file | grep -c "MQRS-12.5";
            grep -v "^#" $file | grep -c "RPRS-8";
            grep -v "^#" $file | grep -c "QD2I";
            grep -v "^#" $file | grep -c "FS200";
            grep -v "^#" $file | grep -c "SOR10";
            grep -v "^#" $file | grep -c "RPRS-20");
        echo ${sample} $line | tr ' ' '\t' >> ${filename%/}_filters_tab.txt;
    done
done


# Get the number of variants per consequence.

# Save all the unique consequences from all files
cons_list=$(grep -v "#" Annotation/*.vcf | cut -d "|" -f6 | sort -u | tr "\n" " ")
# turn the string into an array
read -r -a list <<< $cons_list
echo "Sample" $cons_list | tr ' ' '\t' > consequence_list.txt

for f in Annotation/*.vcf; do
    sample=$(basename "$f" .vcf)
    line=$(
        for cons in $cons_list; do
            grep -v "#" $f | grep $cons | wc -l;
        done);
    echo $sample $line | tr ' ' '\t' >> consequence_list.txt
done
