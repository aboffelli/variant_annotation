## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: 2022-01-21
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description: Script to plot an interactive complex heatmap containing the 
##  synonymous variants from samples that do not have a pathogenic variant 
##  reported in ClinVar.
##    
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## -----------------------------------------------------------------------------

library(ggplot2)
library(plotly)
library(tidyverse)
library(iheatmapr)
options(scipen = 100)

setwd("~/Box/Notes/Tables/SweaSynVar")

synonymous <- read.table('synonymous_table.txt', header=T, sep='\t')
synonymous <- synonymous[order(synonymous$Gene),]

synonymous_table <- synonymous
synonymous_table[!is.na(synonymous_table$RBP),]$RBP<- 'Presence'



main.heatmap <-synonymous_table[,5:6] %>% 
    separate(AF_SWEA, c('AF_SWEA', NA), "([(])")
main.heatmap$AF_SWEA <- as.numeric(main.heatmap$AF_SWEA)

swea_presence <- separate(data=synonymous_table, 
                                 col=AF_SWEA, 
                                 into= c(NA, 'AF_SWEA'), 
                                 sep="([(])")$AF_SWEA

rownames(main.heatmap) <- paste(synonymous_table$Variant, 
                                paste('Codon (ref/alt):', 
                                      synonymous_table$Codon_.ref.alt.), 
                                paste('dbSNP ID:', synonymous_table$DbSNP_ID),
                                paste('SWEA presence:', 
                                      substr(swea_presence, 1, 
                                             nchar(swea_presence)-1)),
                                paste('RBP proteins:', synonymous$RBP),
                                sep='<br>')


main.heatmap <- as.matrix(main.heatmap)
ese_ess <- synonymous_table[,10:12]
phylop <- as.matrix(synonymous_table[,7])
colnames(phylop) <- 'PhyloP'
gerp <- as.matrix(synonymous_table[,8])
colnames(gerp) <- 'GERP'
rscu <- as.matrix(synonymous_table[,9])
colnames(rscu) <- 'deltaRSCU'

genes <- synonymous_table[,3]
genes[is.na(genes)] <- 'Not a targeted gene'

the_plot <- main_heatmap(main.heatmap, name='Allele frequency', 
             colors= c('lightblue', 'darkblue'), tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>%
    add_col_labels() %>% 
    add_row_annotation(data.frame('ESE' = ese_ess$ESE,
                       'ESS' = ese_ess$ESS,
                       'RBP' = ese_ess$RBP
                       ),
                       colors = list('ESE' = c('yellow',
                                               'orange',
                                               'black'),
                                     'ESS' = c('yellow',
                                               'orange',
                                               'black'),
                                     'RBP' = c('forestgreen', 'white'))) %>%
    add_row_groups(genes, title= 'Gene', side='left', show_colorbar=F) %>% 
    add_main_heatmap(phylop, name='PhyloP', colors='RdBu', 
                     tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>%
    add_col_labels() %>% 
    add_main_heatmap(gerp, name='GERP', colors='RdBu', 
                     tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>%
    add_col_labels() %>% 
    add_main_heatmap(rscu, name='deltaRSCU', colors='RdBu', 
                     tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>% 
    add_col_labels()

the_plot

the_plot %>% 
    save_iheatmap("synonymous_heatmap_int.html")
the_plot %>% 
    save_iheatmap('~/Box/Arthur/SweaSynVar/synonymous_heatmap_int.html')
