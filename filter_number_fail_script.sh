#!/bin/bash

for dir in KnownO*/; do
    filename=${dir##/};
    grep "SOR3" $dir/*.vcf | grep -oP "SOR=.+?;" | tr -d "SOR=;" | while read number; do
        echo ${filename%/} $number | tr " " "\t" >> ${filename%/}_sor3_values_table.txt;
    done
done

for dir in NovelO*/; do
    filename=${dir##/};
    grep "SOR3" $dir/*.vcf | grep -oP "SOR=.+?;" | tr -d "SOR=;" | while read number; do
        echo ${filename%/} $number | tr " " "\t" >> ${filename%/}_sor3_values_table.txt;
    done
done
