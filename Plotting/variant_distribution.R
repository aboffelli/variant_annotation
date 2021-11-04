library(ggplot2)


setwd("~/Box/Notes/Tables")
known = read.table('known_variant_distribution.txt', sep='\t')
novel = read.table('novel_variant_distribution.txt', sep='\t')
known$V1 <- paste(known$V1, ' (', known$V2, ')', sep='')
novel$V1 <- paste(novel$V1, ' (', novel$V2, ')', sep='')

y <- ggplot(known, aes(x='', y=V2, fill=V1)) +
    geom_col(col='black', size=0.05) +
    coord_polar(theta = 'y') +
    theme(legend.key.size  = unit(0.2, 'cm'), legend. = 'Consequence') +
    guides(fill=guide_legend(ncol=1)) +
    labs(title="Known variants") +
    theme_void()
ggsave('known.pdf', y)


x <- ggplot(novel, aes(x='', y=V2, fill=V1)) +
    geom_col(col='black', size=0.05) +
    coord_polar(theta = 'y') +
    theme(legend.key.size  = unit(0.2, 'cm'), legend.title = 'Consequence') +
    guides(fill=guide_legend(ncol=1)) +
    labs(title="Novel variants") +
    theme_void()
ggsave('novel.pdf', x)


