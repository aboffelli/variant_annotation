library(ggplot2)
library(gridExtra)
options(scipen = 100)

pie_chart <- function(df, name) {
    ggplot(df, aes(x='', y=V2, fill=V1)) +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        theme(legend.key.size  = unit(0.2, 'cm')) +
        guides(fill=guide_legend(ncol=1)) +
        labs(title=name) +
        scale_fill_discrete(name='Consequence') +
        theme_void()
}

bar_plot <- function(df, name) {
    ggplot(df, aes(x=reorder(V1, -df$V3, sum), y=V3, fill=V2))+
        geom_bar(stat='identity', position='stack', width=0.5) +
        theme_classic() + 
        labs(title=name) +
        scale_fill_discrete(name='Consequence') +
        theme(axis.text.x = element_text(size = 3, angle=90), 
              legend.key.size = unit(0.3, 'cm')) +
        labs(x= 'Gene', y='Count')
}


setwd("~/Box/Notes/Tables")
# known = read.table('known_variant_distribution.txt', sep='\t')
# novel = read.table('novel_variant_distribution.txt', sep='\t')
# known$V1 <- paste(known$V1, ' (', known$V2, ')', sep='')
# novel$V1 <- paste(novel$V1, ' (', novel$V2, ')', sep='')
# 
# known_pie <- pie_chart(known, "Known Variants")
# # print(known_pie)
# 
# novel_pie <- pie_chart(novel, "Novel Variants")
# 
# to_save <- arrangeGrob(known_pie, novel_pie)
# # print(novel_pie)
# ggsave('Plots/pies.pdf', to_save, limitsize=F)
# 
# # Version without the introns
# known_pie_no_intron <- pie_chart(known[-10,], 'Known variants')
# novel_pie_no_intron <- pie_chart(novel[-10,], "Novel Variants")
# 
# to_save_no_intron <- to_save <- arrangeGrob(known_pie_no_intron, 
#                                             novel_pie_no_intron)
# ggsave('Plots/pies_no_intron.pdf', to_save_no_intron, limitsize=F)

gene_known <- read.table('genes_known_variant_distribution.txt', sep='\t')
gene_novel <- read.table('genes_novel_variant_distribution.txt', sep='\t')
x <- aggregate(gene_known$V3, list(gene_known$V1), sum)
gene_known$V4[gene_known$V1==x$Group.1] <- x$x 
# no_intron <- FALSE

# # Without introns
gene_known <- gene_known[!gene_known$V2=='intron',]
gene_novel <- gene_novel[!gene_novel$V2=='intron',]
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
#     ggsave('Plots/known_variants_distribution_by_gene.pdf', grid_save, limitsize = FALSE)
# } else {
#     ggsave('Plots/known_variants_distribution_by_gene_no_intron.pdf', grid_save, limitsize = FALSE)
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
#     ggsave('Plots/novel_variants_distribution_by_gene.pdf', grid_save, limitsize = FALSE)
# } else {
#     ggsave('Plots/novel_variants_distribution_by_gene_no_intron.pdf', grid_save, limitsize = FALSE)
# }
# 

## -----------------------------------------------------------------------------
## Targeted genes tables
gene_k_plot <- bar_plot(gene_known, 'Variants by gene (Known variants)')
gene_n_plot <- bar_plot(gene_novel, 'Variants by gene (Novel variants)')

to_save <- arrangeGrob(gene_k_plot, gene_n_plot)
ggsave('Plots/count_variant_distribution_by_gene_no_intron.pdf', to_save,  limitsize=F)
