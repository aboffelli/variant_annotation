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

library(tidyverse)
library(ggplot2)
library(gridExtra)
options(scipen = 100)

setwd("C:/Users/Arthu/Box/Notes/Tables/SWEA/RscuConserv")


## Delta RSCU plots ----
## Two plots showing the distribution of delta RSCU between the known and novel
## synonymous variants. 

# Load the tables
rscu_table <- read.table('rscu_table_SWEA.txt', sep='\t')
rscu_table$V2 <- factor(rscu_table$V2)

# Create a histogram.
rscu_hist <- ggplot(data=rscu_table, aes(x=V1, fill=V2)) +
    geom_histogram(alpha=0.3, col='black', bins = 30) +
    geom_freqpoly(aes(col=V2)) +
    labs(x="Delta-RSCU", y='Count') +
    scale_fill_discrete(name='Variant type') +
    scale_color_discrete(name="Variant type") +
    theme_classic() +
    theme(text = element_text(size=20))

# Save in pdf and png formats.
ggsave("Plots/rscu_histogram_SWEA.pdf", plot=rscu_hist, width=25, height = 20, 
       units = 'cm')
ggsave("Plots/rscu_histogram_SWEA.png", plot=rscu_hist, width=25, height = 20, 
       units = 'cm')

# Create a violin plot.
rscu_violin <- ggplot(data=rscu_table, aes(y=V1, x=V2, fill=V2)) +
    geom_violin() +
    labs(x='', y='??RSCU', title="Distribution of ??RSCU") +
    theme_classic() +
    scale_fill_brewer(palette="Greys") +
    scale_x_discrete(labels=c("Known variants", "Novel variants")) +
    theme(text = element_text(size=20)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    guides(fill="none")
rscu_violin

# Save in pdf and png formats.
ggsave("Plots/rscu_violin_SWEA.pdf", plot=rscu_violin, width=25, height = 20, 
       units = 'cm')
ggsave("Plots/rscu_violin_SWEA.png", plot=rscu_violin, width=25, height = 20, 
       units = 'cm')

## Conservation plots ----
## Distribution of the conservation values between the variants divided as 
## consequence type. One plot for PhyloP and one for GERP. 


# Load the tables.
conserv_table <- read.table('conservation_table_SWEA.txt', sep='\t', header=T)

# Remove indels.
conserv_table <- conserv_table %>% 
    filter(Consequence!='inframe_deletion/insertion')

# Count the number of variants for the novel and known.
counts <- setNames(data.frame(
    table(conserv_table$Consequence[conserv_table$Exist=='known'])), 
    c("Cons", "FreqK"))
counts$FreqN <- data.frame(
    table(conserv_table$Consequence[conserv_table$Exist=='novel']))[,2]

# Join the known and novel values.
counts$Freq <- paste0(counts$FreqK, '/', counts$FreqN)


# Create the violin plot for PhyloP.
conserv_PhyloP <- ggplot(data=conserv_table, 
                         aes(y=PhyloP, x=Consequence, fill=Exist)) +
    geom_violin() +
    labs(title = 'phyloP distribution for SNVs', y='PhyloP Score') +
    scale_x_discrete(labels=paste0(counts$Cons, '\nn= ', counts$Freq),
                     guide = guide_axis(n.dodge=2)) +
    scale_fill_brewer(name='Variant type', palette="Greys") +
    theme_classic() +
    theme(text = element_text(size=18)) +
    geom_hline(yintercept = 0, linetype="dashed")

conserv_PhyloP
ggsave("Plots/conservation_phylop_SWEA.pdf", plot = conserv_PhyloP, width=35, 
       height = 20, units = 'cm')
ggsave("Plots/conservation_phylop_SWEA.png", plot = conserv_PhyloP, width=35, 
       height = 20, units = 'cm')

# Create the violin plot for GERP.
conserv_GERP <- ggplot(data=conserv_table, 
                       aes(y=GERP, x=Consequence, fill=Exist)) +
    geom_violin() +
    labs(title = 'GERP Distribution', y='GERP Score') +
    scale_x_discrete(labels=paste0(counts$Cons, '\n', counts$Freq),
                     guide = guide_axis(n.dodge=2)) +
    scale_fill_discrete(name='Variant type') +
    theme_classic() +
    theme(text = element_text(size=20)) +
    geom_hline(yintercept = 0, linetype="dashed")
ggsave("Plots/conservation_gerp_SWEA.pdf", plot=conserv_GERP, width = 35, 
       height = 20, units = 'cm')
ggsave("Plots/conservation_gerp_SWEA.png", plot=conserv_GERP, width = 35, 
       height = 20, units = 'cm')