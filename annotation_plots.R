# Plots for the Annotation part
library(ggplot2)
library(reshape2)

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
    melted <- melt(table)
    p <- aggregate(melted$value, list(melted$variable), sum)
    p <- cbind(p, Percentage=as.numeric(format((p$x/sum(p$x))*100, scientific=F)))
    colnames(p) <- c("Group", "Number", "Percentage")
    return(p)
}


################################################################################
type_table <- read.table('type_comparison.txt', header = T, sep="\t")[,-1]

type_comparison <- plot_box(type_table)
print(type_comparison)
ggsave("Plots/novel_vs_known.pdf", plot=novel_vs_known)

a_table <- percentage(type_table)
write.table(a_table, file="percentage_novel_vs_known.txt", sep="\t", row.names = F, quote = F)

p_type_table <- ggplot(a_table, aes(x=Group, y=Percentage, fill=Group)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet")) +
    theme_classic()
print(p_type_table)
ggsave("Plots/p_novel_vs_known.pdf", plot=p_novel_vs_known)


################################################################################
sor_known_on <- read.table('KnownOnFail_sor3_values_table.txt', sep='\t')
sor_known_off <- read.table('KnownOffFail_sor3_values_table.txt', sep='\t')
sor_novel_on <- read.table('NovelOnFail_sor3_values_table.txt', sep='\t')
total <- rbind(sor_known_on, sor_novel_on)

n_bin = 50
ggplot(data=total, aes(x=V2)) +
    geom_histogram(fill='gray', alpha=0.3, col='black', bins = n_bin) +
    geom_freqpoly(bins=n_bin, aes(col=V1)) +
    labs(x="SOR3 Value", y="Count") +
    theme_classic()
