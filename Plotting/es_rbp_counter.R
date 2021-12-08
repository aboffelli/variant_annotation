## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: 2021-11-19
##
## GitHub: https://github.com/aboffelli/
##
## Description:
##    
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## -----------------------------------------------------------------------------
library(ggplot2)

es_plot <- function(table, name) {
    x <- ggplot(data=table, aes(x=V1, y=V3, fill=V2)) +
        geom_bar(stat='identity') +
        theme_classic() +
        labs(x='Consequence', y='Count', title=name) +
        scale_fill_manual(name='Variant type', values=c('gray60', 'gray20'))
    return(x)
}


setwd("C:/Users/Arthu/Box/Notes/Tables/EseEssRbp")


ese_count <- read.table('ese_count.txt', sep = '\t')
ess_count <- read.table('ess_count.txt', sep = '\t')

rbp_count <- read.table('rbp_count.txt', sep = '\t')
rbp_protein <- read.table('rbp_protein_frequency.txt', sep = '\t')
rbp_variant <- read.table('rbp_variant_frequency.txt', sep = '\t')

ese_count_plot <- es_plot(ese_count, 'ESE Count')
ggsave('ese_count.pdf', ese_count_plot)
ess_count_plot <- es_plot(ess_count, 'ESS Count')
ggsave('ess_count.pdf', ess_count_plot)

rbp_count_plot <- ggplot(data=rbp_count, aes(x=V2, y=V3, fill=V1))+
    geom_bar(stat='identity') +
    theme_classic() +
    labs(x='RNA Binding Protein Site', y='Count', 
         title='RNA Binding Protein Site Presence') +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20'))
ggsave('rbp_count.pdf', rbp_count_plot)


rbp_protein_plot <- ggplot(data=as.data.frame(table(rbp_protein)), 
                           aes(x=reorder(V1, -Freq), y=Freq, fill=V2)) +
    geom_bar(stat='identity') +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90, size=5)) +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20')) +
    labs(x='Protein Name', y='Frequency', title= "Frequency of RBP")
ggsave('rbp_protein.pdf', rbp_protein_plot)

rbp_variant <- as.data.frame(table(rbp_variant[,-1]))
rbp_variant_plot <- ggplot(data=rbp_variant, aes(x=V2, y= Freq, fill=V3)) +
    geom_bar(stat='identity') + 
    theme_classic() +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20')) +
    labs(x='Number of RBPs', y='Frequency', 
    title='Frequency of number of RBP on each variant')
ggsave('rbp_variant.pdf', rbp_variant_plot)
