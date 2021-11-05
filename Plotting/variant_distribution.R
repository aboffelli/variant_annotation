library(ggplot2)
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
    ggplot(df, aes(x=gene_name, y=V2, fill=V1, group=-V2))+
        geom_bar(stat='identity', position='stack', width=0.5) +
        labs(title=name) +
        scale_fill_discrete(name='Consequence') +
        theme_classic()
}


#setwd("~/Box/Notes/Tables")
known = read.table('known_variant_distribution.txt', sep='\t')
novel = read.table('novel_variant_distribution.txt', sep='\t')
known$V1 <- paste(known$V1, ' (', known$V2, ')', sep='')
novel$V1 <- paste(novel$V1, ' (', novel$V2, ')', sep='')

test <- cbind(gene_name=c('geneA'), novel)
test2 <- cbind(gene_name=c('geneB'), novel)
test2$V2 <- test2$V2/20 * 5
all <- rbind(test, test2)

known_pie <- pie_chart(known, "Known Variants")
print(known_pie)
#ggsave('known.pdf', y)

novel_pie <- pie_chart(novel, "Novel Variants")
print(novel_pie)
#ggsave('novel.pdf', x)

test_by_gene_plot <- bar_plot(all, "Variants by gene")
print(test_by_gene_plot)
