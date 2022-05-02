## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Mon Mar 14 11:03:05 2022
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description: Plots for quality control of the BRIDGES files. One plot 
##  containing the percentage of known and novel variants, and one plot with
##  the distribution of the allele fraction for all variants. 
##
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## ----------------------------------------------------------------------------- 

library(ggrepel)
library(tidyverse)
library(gridExtra)

setwd("~/Box/Notes/Tables/BRIDGES/KnownNovel")

pie_chart <- function(file_table, plot_name) {
    ## -------------------------------------------------------------------------
    ## Function to create a pie chart using ggplot2.
    ## -------------------------------------------------------------------------
    
    # Rename the first and second column to avoid errors
    colnames(file_table)[c(1,2)] <- c('V1', "V2")
    x <- ggplot(data=file_table, aes(x='', y=Perc, fill=paste(V1,V2))) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        
        # Label flags with the percentage.
        geom_label_repel(data = file_table,
                         aes(y = pos, label = paste0(Perc, "%")),
                         size = 4.5, nudge_x = 0.6, show.legend = FALSE) +
        labs(title=plot_name)
    return(x)
}


## -----------------------------------------------------------------------------
## Plot for the percentage of known and novel variant before and after the 
## second filtration.

known_novel <- read.table('Tables/Tables/novel_known_count_BRIDGES.txt', sep='\t') %>% 
    arrange(V1)
filt_known_novel <- read.table('Tables/filtered_novel_known_count_BRIDGES.txt', 
                               sep='\t') %>% 
    arrange(V1)

# Add a percentage column.
known_novel <- transform(known_novel, 
                         Perc = ave(V3, V1, FUN = function(x) round(x/sum(x), 
                                                                    2)*100))
filt_known_novel <- transform(filt_known_novel, Perc = ave(V3, V1, FUN = function(x) round(x/sum(x), 2)*100))

# Add the positions for the flags in the pie chart.
known_novel <- known_novel %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

filt_known_novel <- filt_known_novel %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

# Create the two pie charts, join the plots and save it as pdf.
known_novel_plot <- pie_chart(known_novel, "Known vs novel percentage")
filt_known_novel_plot <- pie_chart(filt_known_novel, "Known vs novel percentage after filtration")
grid.arrange(known_novel_plot, filt_known_novel_plot, ncol=2)
ggsave('Plots/Plots/known_vs_novel_piechart_BRIDGES.pdf ', 
       arrangeGrob(known_novel_plot, filt_known_novel_plot, ncol=2))


##------------------------------------------------------------------------------
## Plots for the distribution of allele fraction of the variants before and
## after the second filtration.


histogram <- read.table('Tables/allele_fraction_BRIDGES.txt', sep = '\t')
filt_histogram <- read.table('Tables/filtered_allele_fraction_BRIDGES.txt', sep = '\t')

# Create the density plot for both tables.
density <- ggplot(data = histogram, aes(x=V3, fill=V2)) +
    stat_density() +
    facet_wrap(~V1, ncol = 1) +
    labs(title='Allele Fraction density', x="Allele fraction", y="Density") +
    scale_fill_discrete(name="Variant type") +
    theme_bw()

filt_density <- ggplot(data = filt_histogram, aes(x=V3, fill=V2)) +
    stat_density() +
    facet_wrap(~V1, ncol = 1) +
    labs(title='Allele Fraction density after filtration', x="Allele fraction",
         y="Density")+
    scale_fill_discrete(name="Variant type") +
    theme_bw()

# Join both plots and save it as a pdf.
arrangeGrob(density, filt_density)
ggsave('Plots/allele_fraction_density_BRIDGES.pdf', 
       arrangeGrob(density, filt_density))
