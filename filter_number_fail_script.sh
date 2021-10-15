#!/bin/bash

for dir in KnownO*/; do
    filename=${dir##/};
    grep "QD2S" $dir/*.vcf | grep -oP "QD=.+?;" | tr -d "QD=;" | while read number; do
        echo ${filename%/} $number | tr " " "\t" >> ${filename%/}_qd2s_values_table.txt;
    done
done

for dir in NovelO*/; do
    filename=${dir##/};
    grep "QD2S" $dir/*.vcf | grep -oP "QD=.+?;" | tr -d "QD=;" | while read number; do
        echo ${filename%/} $number | tr " " "\t" >> ${filename%/}_qd2s_values_table.txt;
    done
done
