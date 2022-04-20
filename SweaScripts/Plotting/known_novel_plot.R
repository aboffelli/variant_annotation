## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Mon Mar 14 11:03:05 2022
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description: Plot for quality control of the percentage of known and novel
##  variants in the SWEA files.
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

setwd("~/Box/Notes/Tables/SWEA/")

pie_chart <- function(file_table, plot_name) {
    colnames(file_table)[1] <- 'V1'
    x <- ggplot(data=file_table, aes(x='', y=Perc, 
                                     fill=V1)) +
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



known_novel <- read.table('novel_known_count_SWEA.txt', sep='\t') %>% 
    arrange(V1)


known_novel <- transform(known_novel, Perc = ave(V2, FUN = function(x) round(x/sum(x), 2)*100))

known_novel <- known_novel %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

known_novel_plot <- pie_chart(known_novel, "Known vs novel percentage - SWEA")
print(known_novel_plot)

ggsave('Plots/known_vs_novel_piechart_SWEA.pdf', known_novel_plot)
