# Plots for the Annotation part
library(ggplot2)
library(reshape2)

# Mac
setwd("~/Box/Notes/Tables")

fail_plot <- function(table) {
    p <- ggplot(data=melt(table), aes(x=variable, y=value)) +
        geom_jitter(shape=21, size=1, aes(fill=variable), alpha=0.1) +
        geom_boxplot(aes(fill=variable), alpha=0.7) +
        theme_classic() +
        labs(x="Filter", y="Number of variants removed") +
        scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet", "indianred1", "lightskyblue3"))
    return(p)
}


plot_box <- function(table) {
    p <- ggplot(melt(table), aes(x=variable, y=value)) +
        geom_boxplot(aes(fill=variable)) +
        scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
        xlab("Group") +
        ylab("Number of variants") +
        theme_classic()
    return(p)
}


percentage <- function(table) {
    melted <- melt(table)
    p <- aggregate(melted$value, list(melted$variable), sum)
    p <- cbind(p, Percentage=p$x/sum(p$x)*100)
    colnames(p) <- c("Group", "Number", "Percentage")
    return(p)
}

novel_vs_known_table <- read.table('novel_vs_known.txt', header = T, sep="\t")[,-1]
on_vs_off_table <- read.table('on_vs_off.txt', header = T, sep="\t")[,-1]
consequence_table <- read.table('consequence_list.txt', header=T, sep="\t")[,-1]


novel_vs_known <- plot_box(novel_vs_known_table)
print(novel_vs_known)
ggsave("Plots/novel_vs_known.pdf", plot=novel_vs_known)

on_vs_off <- plot_box(on_vs_off_table)
print(on_vs_off)
ggsave("Plots/on_vs_off.pdf", plot=on_vs_off)

# consq <- ggplot(melt(consequence_table), aes(x=variable, y=value)) +
#     geom_boxplot(aes(fill=variable)) +
#     theme_classic()
# print(consq)
# ggsave("Plots/consequence.pdf", plot=consq)



a_table <- percentage(novel_vs_known_table)
write.table(a_table, file="percentage_novel_vs_known.txt", sep="\t", row.names = F, quote = F)

p_novel_vs_known <- ggplot(a_table, aes(x=Group, y=Percentage, fill=Group)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
    theme_classic()
print(p_novel_vs_known)
ggsave("Plots/p_novel_vs_known.pdf", plot=p_novel_vs_known)



a_table2 <- percentage(on_vs_off_table)
write.table(a_table2, file="percentage_on_vs_off.txt", sep="\t", row.names = F, quote = F)

p_on_off <- ggplot(a_table2, aes(x=Group, y=Percentage, fill=Group)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
    theme_classic()
print(p_on_off)
ggsave("Plots/p_on_vs_off.pdf", plot=p_on_off)


on_fail <- read.table("on_target_fail.txt", header=T, sep="\t")[,-1]
off_fail <- read.table("off_target_fail.txt", header=T, sep="\t")[,-1]

known_fail <- read.table("known_fail.txt", header=T, sep="\t")[,-1]
novel_fail <- read.table("novel_fail.txt", header=T, sep="\t")[,-1]



on_fail_plot <- fail_plot(on_fail)
print(on_fail_plot)
off_fail_plot <- fail_plot(off_fail)
print(off_fail_plot)
