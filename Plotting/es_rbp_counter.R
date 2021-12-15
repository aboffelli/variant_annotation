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
library(dplyr)

es_plot <- function(table, name) {
    x <- ggplot(data=table, aes(x=V2, y=V3, fill=V1)) +
        geom_bar(stat='identity') +
        theme_classic() +
        labs(x='Variant Type', y='Count', title=name) +
        scale_fill_manual(name='Consequence',
                          values=c('gray60', 'gray20', 'gray2')) +
        geom_text(aes(y = label_y, label = Percentage), vjust=1.2, colour = "white", size=4)
        
        
    return(x)
}

setwd("C:/Users/Arthu/Box/Notes/Tables/EseEssRbp")


ese_count <- read.table('ese_count.txt', sep = '\t')
ess_count <- read.table('ess_count.txt', sep = '\t')

ese_count <-  cbind(ese_count, 
                    Percentage=paste0(round(ese_count$V3/sum(ese_count$V3)*100, 2),'%')) %>%
    arrange(V2, rev(V1)) %>%
    group_by(V2) %>%
    mutate(label_y=cumsum(V3))

ess_count <-  cbind(ess_count, 
                    Percentage=paste0(round(ess_count$V3/sum(ess_count$V3)*100, 2),'%')) %>%
    arrange(V2, rev(V1)) %>%
    group_by(V2) %>%
    mutate(label_y=cumsum(V3))

rbp_count <- read.table('rbp_count.txt', sep = '\t')
rbp_protein <- read.table('rbp_protein_frequency.txt', sep = '\t')
rbp_variant <- read.table('rbp_variant_frequency.txt', sep = '\t')

rbp_protein$V1 <- sapply(strsplit(rbp_protein$V1, split='_', fixed=T), '[[', 1)

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
    theme(axis.text.y=element_text(size=5)) +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20')) +
    labs(x='Protein Name', y='Frequency', title= "Frequency of RBP") + 
    coord_flip()
ggsave('rbp_protein.pdf', rbp_protein_plot)

rbp_variant <- as.data.frame(table(rbp_variant[,-1]))
rbp_variant_plot <- ggplot(data=rbp_variant, aes(x=V2, y= Freq, fill=V3)) +
    geom_bar(stat='identity') + 
    theme_classic() +
    scale_fill_manual(name='Variant type', values=c('gray60', 'gray20')) +
    labs(x='Number of RBPs', y='Frequency', 
    title='Frequency of number of RBP on each variant')
ggsave('rbp_variant.pdf', rbp_variant_plot)
