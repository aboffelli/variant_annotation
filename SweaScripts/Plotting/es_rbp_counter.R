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


library(ggplot2)
library(dplyr)

es_plot <- function(table, name) {
    ##--------------------------------------------------------------------------
    ## Function to create a bar plot using ggplot2.
    ##--------------------------------------------------------------------------
    x <- ggplot(data=table, aes(x=V2, y=V3, fill=V1)) +
        geom_bar(stat='identity') +
        theme_classic() +
        labs(x='Variant Type', y='Count', title=name) +
        scale_fill_manual(name='Consequence',
                          values=c('gray60', 'gray20', 'gray2')) +
        geom_text(aes(y = label_y, label = Percentage), vjust=1.2, colour = "white", size=4)
        
        
    return(x)
}

setwd("~/Box/Notes/Tables/SWEA/EseEssRbp")

##------------------------------------------------------------------------------
## ESE/ESS alterations count bar plots.

ese_count <- read.table('ese_count_SWEA.txt', sep = '\t')
ess_count <- read.table('ess_count_SWEA.txt', sep = '\t')

# Add the percentage column and the position for the labels in both tables.
ese_count <-  cbind(ese_count, 
                    Percentage=paste0(round(ese_count$V3/sum(ese_count$V3)*100, 
                                            2),'%')) %>%
    arrange(V2, rev(V1)) %>%
    group_by(V2) %>%
    mutate(label_y=cumsum(V3))

ess_count <-  cbind(ess_count, 
                    Percentage=paste0(round(ess_count$V3/sum(ess_count$V3)*100,
                                            2),'%')) %>%
    arrange(V2, rev(V1)) %>%
    group_by(V2) %>%
    mutate(label_y=cumsum(V3))

# Create the plots and save them as pdfs.
ese_count_plot <- es_plot(ese_count, 'ESE Count')
ggsave('Plots/ese_count_SWEA.pdf', ese_count_plot)
ess_count_plot <- es_plot(ess_count, 'ESS Count')
ggsave('Plots/ess_count_SWEA.pdf', ess_count_plot)


##------------------------------------------------------------------------------
## Bar plot with the percentage of existence of overlapping RBP binding sites 
## in the variants.

rbp_count <- read.table('rbp_count_SWEA.txt', sep = '\t')

# Create the plot and save it as a pdf.
rbp_count_plot <- ggplot(data=rbp_count, aes(x=V2, y=V3, fill=V1))+
    geom_bar(stat='identity') +
    theme_classic() +
    labs(x='RNA Binding Protein Site', y='Count', 
         title='RNA Binding Protein Site Presence') +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20'))

ggsave('Plots/rbp_count_SWEA.pdf', rbp_count_plot)

##------------------------------------------------------------------------------
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

ggsave('Plots/rbp_protein_SWEA.pdf', rbp_protein_plot)

##------------------------------------------------------------------------------
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
