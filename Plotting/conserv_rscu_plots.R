## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Mon Nov 15 15:11:11 2021
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description:
##  Script to plot the distribution of the conservation and delta-rscu values. 
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## ----------------------------------------------------------------------------- 

library(ggplot2)
library(gridExtra)
options(scipen = 100)

setwd("~/Box/Notes/Tables")

rscu_table <- read.table('rscu_table.txt', sep='\t')
rscu_table$V2 <- factor(rscu_table$V2)
conserv_table <- read.table('conservation_table.txt', sep='\t', header=T)
conserv_table[c('Gene', 'Exist', 'Consequence')] <- lapply(
    conserv_tablec('Gene', 'Exist', 'Consequence'), factor)


ggplot(data=rscu_table, aes(x=V1, fill=V2)) +
    geom_histogram(alpha=0.3, col='black', bins = 30) +
    geom_freqpoly(aes(col=V2)) +
    theme_classic()

ggplot(data=rscu_table, aes(y=V1, x=V2, fill=V2)) +
    geom_violin() +
    labs(x='', y='delta-rscu') +
    theme_classic()

