## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Mon Nov 15 15:11:54 2021
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description:
##  Script to plot the distribution of variants by type and genes.
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## ----------------------------------------------------------------------------- 

library(ggplot2)
library(gridExtra)
options(scipen = 100)

pie_chart <- function(df, name) {
    ggplot(df, aes(x='', y=V2, fill=V1)) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        theme(legend.key.size  = unit(0.7, 'cm')) +
        guides(fill=guide_legend(ncol=1)) +
        labs(title=name) +
        scale_fill_discrete(name='Consequence')
        
}

bar_plot <- function(df, name, type='stack') {
    ggplot(df, aes(x=reorder(V1, -df$V3, sum), y=V3, fill=V2))+
        geom_bar(stat='identity', position=type, width=0.5) +
        theme_classic() + 
        labs(title=name) +
        scale_fill_discrete(name='Consequence') +
        theme(text=element_text(size=15),
            axis.text.x = element_text(size=7, angle=90), 
              legend.key.size = unit(0.3, 'cm')) +
        labs(x= 'Gene', y='Percentage')
}


setwd("~/Box/Notes/Tables/SWEA/VariantDistribution")

known = read.table('known_variant_distribution_SWEA.txt', sep='\t')
novel = read.table('novel_variant_distribution_SWEA.txt', sep='\t')
known$V1 <- paste(known$V1, ' (', known$V2, ')', sep='')
novel$V1 <- paste(novel$V1, ' (', novel$V2, ')', sep='')

known_pie <- pie_chart(known, "Known Variants")
# print(known_pie)

novel_pie <- pie_chart(novel, "Novel Variants")

to_save <- arrangeGrob(known_pie, novel_pie)
# print(novel_pie)
ggsave('Plots/pies_SWEA.pdf', to_save)
ggsave('Plots/pies_SWEA.png', to_save)

# Version without the introns
known_pie_no_intron <- pie_chart(known[-c(10:11),], 'Known variants')
novel_pie_no_intron <- pie_chart(novel[-c(10:11),], "Novel Variants")

to_save_no_intron <- to_save <- arrangeGrob(known_pie_no_intron,
                                            novel_pie_no_intron)
ggsave('Plots/pies_no_intron_SWEA.png', to_save_no_intron, width=17, height=18, unit='cm')
ggsave('Plots/pies_no_intron_SWEA.pdf', to_save_no_intron, width=17, height=18, unit='cm')

gene_known <- read.table('genes_known_variant_distribution_SWEA.txt', sep='\t')
gene_novel <- read.table('genes_novel_variant_distribution_SWEA.txt', sep='\t')
# no_intron <- FALSE

# # Without introns
gene_known_no_intron <- gene_known[!(gene_known$V2=='intron') & 
                                       !(gene_known$V2=='intergenic'),]
gene_novel_no_intron <- gene_novel[!(gene_novel$V2=='intron') & 
                                       !(gene_known$V2=='intergenic'),]
# no_intron <- TRUE

## -----------------------------------------------------------------------------
## Total tables
# x <- 1
# if (!no_intron) {
#     y <- 1008
# } else {
#     y <- 924
# }
# for (i in 0:2) {
#     if (x < y*2) {
#         gene_k_plot <- bar_plot(gene_known[x:(x + y - 1),], 
#                                 "Variants by gene (Known variants)")
#     }
#     else {
#         gene_k_plot <- bar_plot(gene_known[x:nrow(gene_known),], 
#                                 "Variants by gene (Known variants)")
#     }
#     nam <- paste('gene_k_plot', i, sep='')
#     assign(nam, gene_k_plot)
#     x = x + y
# }
# grid_save <- arrangeGrob(gene_k_plot0, gene_k_plot1, gene_k_plot2)
# if (!no_intron) {
#     ggsave('Plots/known_variants_distribution_by_gene_SWEA.pdf', grid_save, limitsize = FALSE)
# } else {
#     ggsave('Plots/known_variants_distribution_by_gene_no_intron_SWEA.pdf', grid_save, limitsize = FALSE)
# }
# 
# x <- 1
# for (i in 0:2) {
#     if (x < y*2) {
#         gene_n_plot <- bar_plot(gene_novel[x:(x + y - 1),], 
#                                 "Variants by gene (Novel variants)")
#     }
#     else {
#         gene_n_plot <- bar_plot(gene_novel[x:nrow(gene_novel),], "Variants by gene (Novel variants)")
#     }
#     nam <- paste('gene_n_plot', i, sep='')
#     assign(nam, gene_n_plot)
#     x = x + y
# }
# 
# grid_save <- arrangeGrob(gene_n_plot0, gene_n_plot1, gene_n_plot2)
# if (!no_intron) {
#     ggsave('Plots/novel_variants_distribution_by_gene_SWEA.pdf', grid_save, limitsize = FALSE)
# } else {
#     ggsave('Plots/novel_variants_distribution_by_gene_no_intron_SWEA.pdf', grid_save, limitsize = FALSE)
# }
# 

## -----------------------------------------------------------------------------
## Targeted genes tables
count_gene_k_plot <- bar_plot(gene_known, 'Variants by gene (Known variants)')
count_gene_n_plot <- bar_plot(gene_novel, 'Variants by gene (Novel variants)')
count_gene_k_plot_no_intron <- bar_plot(gene_known_no_intron, 'Variants by gene (Known variants) - No introns')
count_gene_n_plot_no_intron <- bar_plot(gene_novel_no_intron, 'Variants by gene (Novel variants) - No introns')

p_gene_k_plot <- bar_plot(gene_known, 'Variants by gene (Known variants)', 'fill')
p_gene_n_plot <- bar_plot(gene_novel, 'Variants by gene (Novel variants)', 'fill')
p_gene_k_plot_no_intron <- bar_plot(gene_known_no_intron, 'Variants by gene (Known variants) - No introns', 'fill')
p_gene_n_plot_no_intron <- bar_plot(gene_novel_no_intron, 'Variants by gene (Novel variants) - No introns', 'fill')

to_save <- arrangeGrob(count_gene_k_plot, count_gene_n_plot)
ggsave('Plots/count_variant_distribution_by_gene_SWEA.png', to_save, width=20, height=18, unit='cm')
ggsave('Plots/count_variant_distribution_by_gene_SWEA.pdf', to_save, width=20, height=18, unit='cm')
to_save <- arrangeGrob(count_gene_k_plot_no_intron, count_gene_n_plot_no_intron)
ggsave('Plots/count_variant_distribution_by_gene_no_intron_SWEA.png', to_save, width=20, height=18, unit='cm')
ggsave('Plots/count_variant_distribution_by_gene_no_intron_SWEA.pdf', to_save, width=20, height=18, unit='cm')

to_save <- arrangeGrob(p_gene_k_plot, p_gene_n_plot)
ggsave('Plots/percentage_variant_distribution_by_gene_SWEA.png', to_save, width=20, height=18, unit='cm')
ggsave('Plots/percentage_variant_distribution_by_gene_SWEA.pdf', to_save, width=20, height=18, unit='cm')
to_save <- arrangeGrob(p_gene_k_plot_no_intron, p_gene_n_plot_no_intron)
ggsave('Plots/percentage_variant_distribution_by_gene_no_intron_SWEA.png', to_save, width=20, height=18, unit='cm')
ggsave('Plots/percentage_variant_distribution_by_gene_no_intron_SWEA.pdf', to_save, width=20, height=18, unit='cm')