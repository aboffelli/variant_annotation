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

# TODO: Add comments
library(ggplot2)
library(gridExtra)
options(scipen = 100)

setwd("~/Box/Notes/Tables/SWEA/RscuConserv")

rscu_table <- read.table('rscu_table_SWEA.txt', sep='\t')
rscu_table$V2 <- factor(rscu_table$V2)
conserv_table <- read.table('conservation_table_SWEA.txt', sep='\t', header=T)
conserv_table <- conserv_table[conserv_table$Consequence!='inframe_deletion/insertion',]
conserv_table[c('Exist', 'Consequence')] <- lapply(
    conserv_table[c('Exist', 'Consequence')], factor)

counts <- setNames(data.frame(table(conserv_table$Consequence[conserv_table$Exist=='known'])), c("Cons", "FreqK"))
counts$FreqN <- data.frame(table(conserv_table$Consequence[conserv_table$Exist=='novel']))[,2]
counts$Freq <- paste0(counts$FreqK, '/', counts$FreqN)

rscu_hist <- ggplot(data=rscu_table, aes(x=V1, fill=V2)) +
    geom_histogram(alpha=0.3, col='black', bins = 30) +
    geom_freqpoly(aes(col=V2)) +
    labs(x="Delta-RSCU", y='Count') +
    scale_fill_discrete(name='Variant type') +
    scale_color_discrete(name="Variant type") +
    theme_classic() +
    theme(text = element_text(size=20))
ggsave("Plots/rscu_histogram_SWEA.pdf", plot=rscu_hist, width=25, height = 20, 
       units = 'cm')

rscu_violin <- ggplot(data=rscu_table, aes(y=V1, x=V2, fill=V2)) +
    geom_violin() +
    labs(x='', y='Delta-RSCU') +
    theme_classic() +
    scale_fill_discrete(name='Variant type') +
    scale_x_discrete(labels=element_blank()) +
    theme(text = element_text(size=20))
ggsave("Plots/rscu_violin_SWEA.pdf", plot=rscu_violin, width=25, height = 20, 
       units = 'cm')

conserv_PhyloP <- ggplot(data=conserv_table, 
                         aes(y=PhyloP, x=Consequence, fill=Exist)) +
    geom_violin() +
    labs(title = 'PhyloP Distribution', y='PhyloP Score') +
    scale_x_discrete(labels=paste0(counts$Cons, '\n', counts$Freq),
                     guide = guide_axis(n.dodge=2)) +
    scale_fill_discrete(name='Variant type') +
    theme_classic() +
    theme(text = element_text(size=20))
ggsave("Plots/conservation_phylop_SWEA.pdf", plot = conserv_PhyloP, width=35, 
       height = 20, units = 'cm')

conserv_GERP <- ggplot(data=conserv_table, 
                       aes(y=GERP, x=Consequence, fill=Exist)) +
    geom_violin() +
    labs(title = 'GERP Distribution', y='GERP Score') +
    scale_x_discrete(labels=paste0(counts$Cons, '\n', counts$Freq),
                     guide = guide_axis(n.dodge=2)) +
    scale_fill_discrete(name='Variant type') +
    theme_classic() +
    theme(text = element_text(size=20))
ggsave("Plots/conservation_gerp_SWEA.pdf", plot=conserv_GERP, width = 35, 
       height = 20, units = 'cm')
