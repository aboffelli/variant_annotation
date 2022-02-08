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
## Notes: This script uses the packages tidyverse and iheatmapr, make sure 
## that they are installed before running the script.
## 
## -----------------------------------------------------------------------------


library(tidyverse)
library(iheatmapr)
options(scipen = 100)

setwd("~/Box/Notes/Tables/SweaSynVar")


synonymous <- read.table('synonymous_table.txt', header=T, sep='\t')

# Reorder the table based on the Gene column, other columns can be added inside
# the order function.
synonymous <- synonymous[order(synonymous$Gene),]

synonymous_table <- synonymous
# Transform all the values that are not NA from the RBP column.
synonymous_table[!is.na(synonymous_table$RBP),]$RBP<- 'Presence'


# Separate the two columns with AF, and isolate just the value from SWEA.
main.heatmap <-synonymous_table[,5:6] %>% 
    separate(AF_SWEA, c('AF_SWEA', NA), "([(])")
# Transform into numeric characters.
main.heatmap$AF_SWEA <- as.numeric(main.heatmap$AF_SWEA)

# Create a vector with the number of samples in SWEA.
swea_presence <- separate(data=synonymous_table, 
                                 col=AF_SWEA, 
                                 into= c(NA, 'AF_SWEA'), 
                                 sep="([(])")$AF_SWEA

# Add all the information the we want to appear in the hoverinfo in the row 
# names separated by line breaks.
rownames(main.heatmap) <- paste(synonymous_table$Variant, 
                                paste('Codon (ref/alt):', 
                                      synonymous_table$Codon_.ref.alt.), 
                                paste('dbSNP ID:', synonymous_table$DbSNP_ID),
                                paste('SWEA presence:', 
                                      substr(swea_presence, 1, 
                                             nchar(swea_presence)-1)),
                                paste('RBP proteins:', synonymous$RBP),
                                # Plotly uses HTML tags.
                                sep='<br>')


# Transform the data frames into matrix to use in iheatmapr. Only the heatmaps 
# need to be matrixes, the annotations can be data frames.
main.heatmap <- as.matrix(main.heatmap)
ese_ess <- synonymous_table[,10:12]
phylop <- as.matrix(synonymous_table[,7])
# Rename the columns for the column label in the plot.
colnames(phylop) <- 'PhyloP'
gerp <- as.matrix(synonymous_table[,8])
colnames(gerp) <- 'GERP'
rscu <- as.matrix(synonymous_table[,9])
colnames(rscu) <- 'deltaRSCU'

# Assign a value for NAs in the gene column.
genes <- synonymous_table[,3]
genes[is.na(genes)] <- 'Not a targeted gene'


# Create the plot
the_plot <- main_heatmap(main.heatmap, name='Allele frequency', 
             colors= c('lightblue', 'darkblue'), 
             # Change the hover info 'Row' to 'Variant'
             tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>%
    # Add the column label for the AF heatmap
    add_col_labels() %>%
    
    # Add the three annotation columns for ESE, ESS and RBP
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
    
    # Add the Gene annotation in the left side
    add_row_groups(genes, title= 'Gene', side='left', show_colorbar=F) %>% 
    
    # Add the PhyloP heatmap and the column labels
    add_main_heatmap(phylop, name='PhyloP', colors='RdBu', 
                     tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>%
    add_col_labels() %>% 
    
    # Add the GERP heatmap and the column label
    add_main_heatmap(gerp, name='GERP', colors='RdBu', 
                     tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>%
    add_col_labels() %>% 
    
    # Add the RSCU heatmap and the column label
    add_main_heatmap(rscu, name='deltaRSCU', colors='RdBu', 
                     tooltip=setup_tooltip_options(prepend_row='Variant: ')) %>% 
    add_col_labels()

# Show the plot in the Viewer tab
the_plot

# Save the HTML plot.
the_plot %>% 
    save_iheatmap("synonymous_heatmap_int.html")

# To save as a static plot the package webshot is necessary, after installing 
# it, just change the extension of the file name to pdf/png/jpeg.