## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Mon Mar 14 11:03:05 2022
##
## GitHub: https://github.com/aboffelli/
##
## Description: Plots for quality control of the BRIDGES files.
##
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## ----------------------------------------------------------------------------- 

library(ggplot2)
library(ggrepel)
library(tidyverse)

setwd("~/Box/Notes/TestData/Bridges/ClinVar")

pie_chart <- function(file_table, plot_name) {
    colnames(file_table)[1] <- 'V1'
    x <- ggplot(data=file_table, aes(x='', y=Perc, 
                                     fill=paste(V1,V2))) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        geom_label_repel(data = file_table,
                         aes(y = pos, label = paste0(Perc, "%")),
                         size = 4.5, nudge_x = 0.6, show.legend = FALSE) +
        labs(title=plot_name)
    return(x)
}



known_novel <- read.table('novel_known_count.txt', sep='\t')

known_novel <- transform(known_novel, Perc = ave(V3, V1, FUN = function(x) round(x/sum(x), 2)*100))

known_novel <- known_novel %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

known_novel_plot <- pie_chart(known_novel, "Known vs novel percentage")
print(known_novel_plot)
ggsave('pie_chart.pdf', known_novel_plot)


################################################################################

histogram <- read.table('table.txt', sep = '\t')

density <- ggplot(data = histogram, aes(x=V3, fill=V2)) +
    stat_density() +
    facet_wrap(~V1, ncol = 1) +
    labs(title='Allele Fraction density', x="Allele fraction") +
    scale_fill_discrete(name="Variant type") +
    theme_bw()

print(density)
ggsave('density_212_files.pdf', density)
