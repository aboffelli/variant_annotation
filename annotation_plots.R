# Plots for the Annotation part
library(ggplot2)
library(reshape2)
library(gridExtra)

# Mac
setwd("~/Box/Notes/Tables")


fail_plot <- function(table, plot_name='') {
    p <- ggplot(data=melt(table), aes(x=variable, y=value)) +
        #geom_jitter(shape=21, size=1, aes(fill=variable), alpha=0.1) +
        geom_boxplot(aes(fill=variable), alpha=0.7) +
        theme_classic() +
        labs(x="Filter", y="Number of variants removed", title=plot_name) +
        scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet", "indianred1", "lightskyblue3"))
    return(p)
}


plot_box <- function(table, plot_name='') {
    p <- ggplot(melt(table), aes(x=variable, y=value)) +
        geom_boxplot(aes(fill=variable)) +
        scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet")) +
        labs(x="Group", y="Number of variants", title=plot_name) +
        theme_classic()
    return(p)
}


percentage <- function(table) {
    options(scipen = 100)
    melted <- melt(table)
    p <- aggregate(melted$value, list(melted$variable), sum)
    p <- cbind(p, Percentage=(p$x/sum(p$x))*100)
    colnames(p) <- c("Group", "Number", "Percentage")
    return(p)
}

double_percentage_plot <- function(p_table, plot_name='', cols=c('darkred', 'darkblue')) {
    p <- ggplot(p_table, aes(x=Group, y=Percentage, fill=Group)) +
        geom_bar(stat='identity', position = 'dodge') +
        scale_fill_manual(values=cols) +
        geom_text(aes(label=Percentage), position=position_dodge(width=0.9), vjust=-0.25) +
        theme_classic()
    return(p)
}
################################################################################
type_table <- read.table('type_comparison.txt', header = T, sep="\t")[,-1]
type_table_S4 <- read.table('Sor4/type_comparison.txt', header = T, sep="\t")[,-1]

type_comparison <- plot_box(type_table)
type_comparison_S4 <- plot_box(type_table_S4) 
# print(type_comparison)
# ggsave("Plots/type_comparison.pdf", plot=type_comparison)


a_table <- percentage(type_table)
a_table_S4 <- percentage(type_table_S4)
#write.table(a_table, file="percentage_novel_vs_known.txt", sep="\t", row.names = F, quote = F)


p_type_comparison <- ggplot(a_table, aes(x=Group, y=Percentage, fill=Group)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet")) +
    theme_classic()

p_type_comparison_S4 <- ggplot(a_table_S4, aes(x=Group, y=Percentage, fill=Group)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet")) +
    theme_classic()

#ggsave("Plots/p_type_table.pdf", plot=p_type_comparison)

grid.arrange(type_comparison, type_comparison_S4)
grid.arrange(p_type_comparison, p_type_comparison_S4)

known_on_sum <- double_percentage_plot(percentage(type_table[,1:2]), 'Known On', c("aquamarine", "coral"))
known_on_sum_S4 <- double_percentage_plot(percentage(type_table_S4[,1:2]), 'Known On', c("aquamarine", "coral"))
known_off_sum <- double_percentage_plot(percentage(type_table[,3:4]), 'Known Off', c("cadetblue1", "lightpink1"))
known_off_sum_S4 <- double_percentage_plot(percentage(type_table_S4[,3:4]), 'Known Off', c("cadetblue1", "lightpink1"))
novel_on_sum <- double_percentage_plot(percentage(type_table[,5:6]), "Novel On", c("darkolivegreen1", "lightgoldenrod1"))
novel_on_sum_S4 <- double_percentage_plot(percentage(type_table_S4[,5:6]), "Novel On", c("darkolivegreen1", "lightgoldenrod1"))
novel_off_sum <- double_percentage_plot(percentage(type_table[,7:8]), 'Novel Off', c("cornflowerblue", "blueviolet"))
novel_off_sum_S4 <- double_percentage_plot(percentage(type_table_S4[,7:8]), 'Novel Off', c("cornflowerblue", "blueviolet"))

a <- grid.arrange(known_on_sum, known_on_sum_S4)
b <- grid.arrange(known_off_sum, known_off_sum_S4)
c <- grid.arrange(novel_on_sum, novel_on_sum_S4)
d <- grid.arrange(novel_off_sum, novel_off_sum_S4)

grid.arrange(a, b, c, d)

################################################################################
# Difference of each filter
filters_known_on <- read.table('known_on_target_fail.txt', header = T, sep="\t")[,-1]

filters_known_off <- read.table('known_off_target_fail.txt', header = T, sep="\t")[,-1]

filters_novel_on <- read.table('novel_on_target_fail.txt', header = T, sep="\t")[,-1]

filters_novel_off <- read.table('novel_off_target_fail.txt', header = T, sep="\t")[,-1]

filters_k_on <- fail_plot(filters_known_on, "Known On Target")
# print(filters_k_on)


filters_k_off <- fail_plot(filters_known_off, "Known Off Target")
# print(filters_k_off)

filters_n_on <- fail_plot(filters_novel_on, "Novel On Target")
# print(filters_n_on)

filters_n_off <- fail_plot(filters_novel_off, "Novel Off Target")
# print(filters_n_off)

grid.arrange(filters_k_on, filters_k_off, filters_n_on, filters_n_off)


################################################################################
# Histogram for specific filters
sor_known_on <- read.table('KnownOnFail_sor3_values_table.txt', sep='\t')
sor_known_off <- read.table('KnownOffFail_sor3_values_table.txt', sep='\t')
sor_novel_on <- read.table('NovelOnFail_sor3_values_table.txt', sep='\t')
total_sor <- rbind(sor_known_on, sor_novel_on, sor_known_off)

n_bin = 30
total_sor_plot <- ggplot(data=total_sor, aes(x=V2, y=after_stat(density))) +
    geom_histogram(fill='gray', alpha=0.3, col='black', bins = n_bin) +
    geom_freqpoly(aes(col=V1)) +
    labs(x="SOR Value", y="Density", title = "SOR3 Filter") +
    theme_classic()


# ggplot(data=sor_known_on, aes(x=V2)) +
#     geom_histogram(fill='gray', alpha=0.3, col='black', bins = n_bin) +
#     geom_freqpoly(bins=n_bin, aes(col=V1)) +
#     labs(x="SOR Value", y="Count") +
#     theme_classic()

# off_and_novel <- rbind(sor_known_off, sor_novel_on)
# ggplot(data=off_and_novel, aes(x=V2)) +
#     geom_histogram(fill='gray', alpha=0.3, col='black', bins = n_bin) +
#     geom_freqpoly(bins=n_bin, aes(col=V1)) +
#     labs(x="SOR Value", y="Count") +
#     theme_classic()

################################################################################
# QD2S Filter
qds_known_on <- read.table('KnownOnFail_qd2s_values_table.txt', sep='\t')
qds_known_off <- read.table('KnownOffFail_qd2s_values_table.txt', sep='\t')
qds_novel_on <- read.table('NovelOnFail_qd2s_values_table.txt', sep='\t')
qds_novel_off <- read.table('NovelOffFail_qd2s_values_table.txt', sep='\t')
total_qds <- rbind(qds_known_on, qds_known_off, qds_novel_on)


total_qds_plot <- ggplot(data=total_qds, aes(x=V2, y=after_stat(density))) +
    geom_histogram(fill='gray', alpha=0.3, col='black', bins = n_bin) +
    geom_freqpoly(aes(col=V1)) +
    labs(x="QD Value", y="Density", title="QD2S Filter") +
    theme_classic()


################################################################################
# QD2I Filter

qdi_known_on <- read.table('KnownOnFail_qd2i_values_table.txt', sep='\t')
qdi_known_off <- read.table('KnownOffFail_qd2i_values_table.txt', sep='\t')
qdi_novel_on <- read.table('NovelOnFail_qd2i_values_table.txt', sep='\t')
qdi_novel_off <- read.table('NovelOffFail_qd2i_values_table.txt', sep='\t')
total_qdi <- rbind(qdi_known_on, qdi_known_off, qdi_novel_on)


total_qdi_plot <- ggplot(data=total_qdi, aes(x=V2, y=after_stat(density))) +
    geom_histogram(fill='gray', alpha=0.3, col='black', bins = n_bin) +
    geom_freqpoly(aes(col=V1)) +
    labs(x="QD Value", y="Density", title="QD2I Filter") +
    theme_classic()


grid.arrange(total_sor_plot, total_qdi_plot, total_qds_plot)
