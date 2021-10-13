#!/bin/bash

for dir in *Fail; do
    filename=${dir##/};
    line=$(
        grep "SOR3" $dir/*.vcf | grep -oP "SOR=.+?;" | tr -d "SOR=;"| tr "\n" "\t"| less -N)
        echo ${filename} $line | tr " " "\t" >> sor3_values_table.txt;
done
