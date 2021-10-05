# Plots for the Annotation part
library(ggplot2)
library(reshape2)

# Mac
setwd("~/Box/Notes/Tables")


novel_table <- read.table('novel_vs_know_filter.txt', header = T, sep="\t")[,-1]
on_target_table <- read.table('number_var.txt', header = T, sep="\t")[,-1]
consequence_table <- read.table('consequence_list.txt', header=T, sep="\t")[,-1]


novel_vs_known <- ggplot(melt(novel_table), aes(x=variable, y=value)) +
    geom_boxplot(aes(fill=variable)) +
    scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
    xlab("Group") +
    ylab("Number of variants") +
    theme_classic()
print(novel_vs_known)
ggsave("novel_vs_known.pdf", plot=novel_vs_known)

on_vs_off <- ggplot(melt(on_target_table), aes(x=variable, y=value)) +
    geom_boxplot(aes(fill=variable)) +
    scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
    xlab("Group") +
    ylab("Number of variants") +
    theme_classic()
print(on_vs_off)
ggsave("on_vs_off.pdf", plot=on_vs_off)

# consq <- ggplot(melt(consequence_table), aes(x=variable, y=value)) +
#     geom_boxplot(aes(fill=variable)) +
#     theme_classic()
# print(consq)
# ggsave("consequence.pdf", plot=consq)


melted <- melt(novel_table)
a_table <- aggregate(melted$value, list(melted$variable), sum)
a_table <- cbind(a_table, Percentage=a_table$x/sum(a_table$x)*100)

p_novel_vs_known <- ggplot(a_table, aes(x=Group.1, y=Percentage, fill=Group.1)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
    xlab("Group") +
    theme_classic()
print(p_novel_vs_known)
ggsave("p_novel_vs_known.pdf", plot=p_novel_vs_known)


melted2 <- melt(on_target_table)
a_table2 <- aggregate(melted2$value, list(melted2$variable), sum)
a_table2 <- cbind(a_table2, Percentage=a_table2$x/sum(a_table2$x)*100)

p_on_off <- ggplot(a_table2, aes(x=Group.1, y=Percentage, fill=Group.1)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("cadetblue1","lightgoldenrod1", "coral","darkolivegreen1")) +
    xlab("Group") +
    theme_classic()
print(p_on_off)
ggsave("p_on_vs_off.pdf", plot=p_on_off)