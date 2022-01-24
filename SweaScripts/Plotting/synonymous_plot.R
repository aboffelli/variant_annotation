## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: 2022-01-21
##
## GitHub: https://github.com/aboffelli/
##
## Description: Script to plot the synonymous variants for exploration.
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
options(scipen = 100)

setwd("~/Box/Arthur/SweaSynVar")

synonymous_table <- read.table('synonymous_table.txt', header=T, sep='\t')

af.swegen <- 0.1
af.swea <- 0.9
rscu <- 0
phylop <- 2
gerp <- 0


synonymous_table[!is.na(synonymous_table$AF.SweGen) & synonymous_table$AF.SweGen > af.swegen,]$AF.SweGen <- NA
synonymous_table[!is.na(synonymous_table$AF.SweGen) & synonymous_table$AF.SweGen < af.swegen,]$AF.SweGen <- 'AF SweGen'

synonymous_table[!is.na(synonymous_table$AF.SWEA) & synonymous_table$AF.SWEA > af.swea,]$AF.SWEA <- NA
synonymous_table[!is.na(synonymous_table$AF.SWEA) & synonymous_table$AF.SWEA < af.swea,]$AF.SWEA <- 'AF SWEA'

synonymous_table[!is.na(synonymous_table$RSCU) & synonymous_table$RSCU > rscu,]$RSCU <- NA
synonymous_table[!is.na(synonymous_table$RSCU) & synonymous_table$RSCU < rscu,]$RSCU <- 'RSCU'

synonymous_table[!is.na(synonymous_table$PhyloP) & synonymous_table$PhyloP < phylop,]$PhyloP <- NA
synonymous_table[!is.na(synonymous_table$PhyloP) & synonymous_table$PhyloP > phylop,]$PhyloP <- 'PhyloP'

synonymous_table[!is.na(synonymous_table$GERP) & synonymous_table$GERP <= gerp,]$GERP <- NA
synonymous_table[!is.na(synonymous_table$GERP) & synonymous_table$GERP > gerp,]$GERP <- 'GERP'

synonymous_table[!is.na(synonymous_table$ESE),]$ESE <- 'ESE'
synonymous_table[!is.na(synonymous_table$ESS),]$ESS <- 'ESS'

synonymous_table[!is.na(synonymous_table$RBP),]$RBP<- 'RBP'

synonymous_table_filtered <- synonymous_table[,-2]

synonymous_to_plot <- data.frame()

for (i in 1:nrow(synonymous_table_filtered)) { 
    for (j in 3:ncol(synonymous_table_filtered)) {
        synonymous_to_plot <- rbind(synonymous_to_plot, c(synonymous_table_filtered[i,1], synonymous_table_filtered[i,j], synonymous_table_filtered[i,2], unite(synonymous_table_filtered[i,], "text", colnames(synonymous_table_filtered[i,3:ncol(synonymous_table_filtered)]), sep=',', na.rm=T)$text))
}}
colnames(synonymous_to_plot) <- c('Variant', 'Type', "Gene", "Text")
synonymous_to_plot <- synonymous_to_plot[!is.na(synonymous_to_plot$Type),]
synonymous_to_plot <- drop_na(synonymous_to_plot)

plot <- ggplot(data=synonymous_to_plot, aes(y=Variant, x=Type, fill='', text= Text, gene=Gene)) +
    geom_tile() +
    theme(legend.position='none')

ggplotly(plot, tooltip = c('y','gene', 'text')) %>% layout(hovermode="y unified")

