#!/bin/bash

# Get how many novel and known variants were filtered.
echo -e "Sample\tKnown_Pass\tKnown_Fail\tNovel_Pass\tNovel_Fail" > novel_vs_know_filter.txt

for f in Annotation/*.vcf; do
    sample=$(basename "$f" .vcf);
    line=$(
        grep -v "#" $f | grep "PASS" | cut -d "|" -f2 | grep "rs" | wc -l;
        grep -v "#" $f | grep -v "PASS" | cut -d "|" -f2 | grep "rs" | wc -l;
        grep -v "#" $f | grep "PASS" | cut -d "|" -f2 | grep -v "rs" | wc -l;
        grep -v "#" $f | grep -v "PASS" | cut -d "|" -f2 | grep -v "rs" | wc -l);
    echo $sample $line | tr " " "\t" >> novel_vs_know_filter.txt;
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
