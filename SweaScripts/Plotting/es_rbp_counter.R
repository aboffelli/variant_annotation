## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: 2021-11-19
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description: Script to plot bar plots for the ESE and ESS alteration counts 
##  and bar plots for the counts of existence of overlap in RBP binding sites,
##   the count of proteins involved in the overlapping, and the frequency of 
##   number of RBPs by variant.
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## -----------------------------------------------------------------------------

library(ggrepel)
library(tidyverse)

es_plot <- function(table, name) {
    ##--------------------------------------------------------------------------
    ## Function to create a bar plot using ggplot2.
    ##--------------------------------------------------------------------------
    x <- ggplot(data=table, aes(x=X2, y=X3, fill=X1)) +
        geom_bar(stat='identity', col='black') +
        theme_classic() +
        labs(x='Variant Type', y='Percentage', title=name) +
        scale_fill_manual(name='Consequence',
                          values=c('black', "gray30", "gray50")) +
        geom_label(aes(y = label_y, label = paste0(Perc, '%')), vjust=1.1, 
                  colour = "white")
        
        
    return(x)
}

setwd("~/Box/Notes/Tables/SWEA/EseEssRbp")

## ESE/ESS ------------------------------------------------------------
## ESE/ESS alterations count bar plots.

ese_count <- read_tsv('ese_count_SWEA.txt', col_names = F) %>%
    arrange(X2, rev(X1)) %>% 
    # Add the percentage column and the position for the labels.
    group_by(X2) %>% 
    # Add the percentage and position of the flags.
    mutate(Perc = round(X3/sum(X3) *100, 2),
           label_y = cumsum(X3)) %>% 
    add_column(Type="ESE")
    

ess_count <- read_tsv('ess_count_SWEA.txt', col_names = F) %>% 
    arrange(X2, rev(X1)) %>% 
    # Add the percentage column and the position for the labels.
    group_by(X2) %>% 
    # Add the percentage and position of the flags.
    mutate(Perc = round(X3/sum(X3) *100, 2),
           label_y = cumsum(X3)) %>% 
    add_column(Type="ESS")

count <- ese_count %>% full_join(ess_count)

# Create the plot and save them as pdfs.
count_plot <- ggplot(data=count, aes(x=X2, y=X3, fill=X1)) +
    geom_bar(stat='identity', col='black') +
    theme_classic() +
    labs(x='Variant Type', y='Count', title="ESE/ESS count") +
    scale_fill_manual(name='Consequence',
                      values=c('black', "gray30", "gray50")) +
    geom_text(aes(y = label_y, label = paste0(Perc, '%')), vjust=1.1, 
               colour = "white", check_overlap = T) +
    facet_wrap(~ Type)
count_plot
ggsave("Plots/ese_ess_count.pdf", count_plot, width=30, height=20, units='cm')
ggsave("Plots/ese_ess_count.png", count_plot, width=30, height=20, units='cm')

## RBP existence -------------------------------------------------------
## Bar plot with the percentage of existence of overlapping RBP binding sites 
## in the variants.

rbp_count <- read_tsv('rbp_count_SWEA.txt', col_names = F) %>%
    arrange(X1, rev(X2)) %>% 
    group_by(X1) %>% 
    mutate(Perc = round(X3/sum(X3) *100, 2),
           label_y = cumsum(Perc))

# Create the plot and save it as a pdf.
rbp_count_plot <- ggplot(data=rbp_count, aes(x=X1, y=Perc, fill=X2))+
    geom_bar(stat='identity') +
    theme_classic() +
    labs(x='Variant type', y='Percentage', 
         title='RNA Binding Protein binding site overlap') +
    scale_fill_manual(name='RBP site presence', values=c('gray60', 'gray20')) +
    geom_label(aes(y=label_y, label=paste0(Perc, '%', "\n", "n = ", X3)), 
               vjust=0.5, colour="white")
rbp_count_plot

ggsave('Plots/rbp_count_SWEA.pdf', rbp_count_plot)
ggsave('Plots/rbp_count_SWEA.png', rbp_count_plot)

## RBP proteins ----------------------------------------------------------
## Plot for the frequency of proteins with overlapping RBP binding sites.

rbp_protein <- read.table('rbp_protein_frequency_SWEA.txt', sep = '\t')

# Isolate only the protein name from the first column.
rbp_protein$V1 <- sapply(strsplit(rbp_protein$V1, split='_', fixed=T), '[[', 1)

# Create the plot and save it as a pdf.
rbp_protein_plot <- ggplot(data=as.data.frame(table(rbp_protein)), 
                           aes(x=reorder(V1, -Freq), y=Freq, fill=V2)) +
    geom_bar(stat='identity') +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    theme_classic() +
    theme(axis.text.y=element_text(size=5)) +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20')) +
    labs(x='Protein Name', y='Frequency', title= "Frequency of RBP") + 
    coord_flip()

rbp_protein_plot
ggsave('Plots/rbp_protein_SWEA.pdf', rbp_protein_plot)

## Number of RBPs per variant ---------------------------------------------
## Bar plot for the frequency of number of overlapping RBP binding sites per 
## variant.

rbp_variant <- read.table('rbp_variant_frequency_SWEA.txt', sep = '\t')

# Count the frequency of the number of RBPs.
rbp_variant <- as.data.frame(table(rbp_variant[,-1]))

# Create the plot and save it as a pdf.
rbp_variant_plot <- ggplot(data=rbp_variant, aes(x=V2, y= Freq, fill=V3)) +
    geom_bar(stat='identity') + 
    theme_classic() +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20')) +
    labs(x='Number of RBPs', y='Frequency', 
    title='Frequency of number of RBP on each variant')

ggsave('Plots/rbp_variant_SWEA.pdf', rbp_variant_plot)
