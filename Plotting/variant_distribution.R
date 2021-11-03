library(ggplot2)
library(reshape2)
library(gridExtra)



setwd("~/Box/Notes/Scripts/Plotting")
known = read.table('known_variant_distribution.txt', sep='\t')
novel = read.table('novel_variant_distribution.txt', sep='\t')


print(ggplot(known, aes(x='', y=V2, fill=V1)) +
    geom_col() +
    geom_text(aes(label=V2, x=1.5), position = position_stack(vjust = 0.5), 
              size = 2, angle=55) +
    coord_polar(theta = 'y') +
    theme(legend.key.size  = unit(0.2, 'cm')) +
    guides(fill=guide_legend(ncol=2)) +
    theme_void())


print(ggplot(novel, aes(x='', y=V2, fill=V1)) +
    geom_col() +
    geom_text(aes(label=V2, x=1.48), position = position_stack(vjust = 0.5), 
              size = 2, angle=55) +
    coord_polar(theta = 'y') +
    theme(legend.key.size  = unit(0.2, 'cm')) +
    guides(fill=guide_legend(ncol=2)) +
    theme_void())
