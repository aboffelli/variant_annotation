#!/bin/bash

# Script to retrieve the number of variants that were caught by a specific filter.

# Get the value of the filter specified, regex matches one or more of any character except new lines until a ';', as few characters as possible.
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
