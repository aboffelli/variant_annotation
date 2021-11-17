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
conserv_table <- conserv_table[conserv_table$Consequence!='inframe_deletion/insertion',]
conserv_table[c('Exist', 'Consequence')] <- lapply(
    conserv_table[c('Exist', 'Consequence')], factor)



rscu_hist <- ggplot(data=rscu_table, aes(x=V1, fill=V2)) +
    geom_histogram(alpha=0.3, col='black', bins = 30) +
    geom_freqpoly(aes(col=V2)) +
    labs(x="Delta-RSCU", y='Count') +
    theme_classic()
ggsave("Plots/rscu_histogram.pdf", plot=rscu_hist, width=25, height = 20, 
       units = 'cm')

rscu_violin <- ggplot(data=rscu_table, aes(y=V1, x=V2, fill=V2)) +
    geom_violin() +
    labs(x='', y='Delta-RSCU') +
    theme_classic()
ggsave("Plots/rscu_violin.pdf", plot=rscu_violin, width=25, height = 20, 
       units = 'cm')

conserv_PhyloP <- ggplot(data=conserv_table, 
                         aes(y=PhyloP, x=Consequence, fill=Exist)) +
    geom_violin() +
    labs(title = 'PhyloP Distribution') +
    theme_classic()
ggsave("Plots/conservation_phylop.pdf", plot = conserv_PhyloP, width=35, 
       height = 20, units = 'cm')

conserv_GERP <- ggplot(data=conserv_table, 
                       aes(y=GERP, x=Consequence, fill=Exist)) +
    geom_violin() +
    labs(title = 'GERP Distribution') +
    theme_classic()
ggsave("Plots/conservation_gerp.pdf", plot=conserv_GERP, width = 35, 
       height = 20, units = 'cm')

    