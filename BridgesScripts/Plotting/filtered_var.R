# Create plots for the variant filtering
library(ggplot2)
library(reshape2)
library(plotly)


# Mac
setwd("~/Box/Notes")

# Windows
# setwd("C:\\Users\\Arthu\\Box\\Notes")

# Load tables
t_table <- as.data.frame(read.table('hf-ml_comparison.txt', header=T, sep="\t"))
t_table$Filtering_Type <- factor(t_table$Filtering_Type)

f_table <- as.data.frame(read.table('filters_tab.txt', header=T, sep="\t"))
f_table <- f_table[,-1]


# Stack bar for the percentage of variants that passed or failed in the filtering comparing the hard filtering and the machine learning approaches.
a_table <- aggregate(cbind(Pass=t_table$Pass, Fail=t_table$Fail), list(Filter=t_table$Filtering_Type), sum)

x <- (a_table$Pass*100)/(a_table$Pass + a_table$Fail)

p_table <- data.frame(Filter=a_table$Filter, Pass=round(x, 2), Fail=round(100-x, 2))


stackbar <- ggplot(data=melt(p_table), aes(y=value, x=Filter, fill=variable)) +
    theme_classic() +
    geom_bar(stat = "identity", width = 0.4) +
    scale_fill_manual(values= c("salmon3","darkred")) +
    ylab("Percentage (%)") +
    theme(
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.text=element_text(size=15))
    
print(stackbar)
ggsave("stackbar.png", plot=stackbar)
ggsave("stackbar.pdf", plot=stackbar)


# Box plot with the number of variables for each approach
total_box <- ggplot(data=melt(t_table), aes(x=Filtering_Type, y=value, fill=variable)) +
    geom_boxplot(alpha=0.7) +
    theme_classic() +
    scale_fill_manual(values= c("salmon3","darkred")) +
    labs(x="Filter", y="Number of variants") +
    theme(
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.text=element_text(size=15))

print(total_box)
ggsave("total_box.png", plot=total_box)
ggsave("total_box.pdf", plot=total_box)

# Plot for each filter in the Hard filtering approach

filter_plot <- ggplot(data=melt(f_table), aes(x=variable, y=value)) +
    geom_jitter(shape=21, size=1, aes(fill=variable), alpha=0.1) +
    geom_boxplot(aes(fill=variable), alpha=0.7) +
    theme_classic() +
    labs(x="Filter", y="Number of variants removed") +
    scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet", "indianred1", "lightskyblue3" ))

print(filter_plot)
ggsave("filter_plot.png", plot=filter_plot)
ggsave("filter_plot.pdf", plot=filter_plot)
    
no_jitter <- ggplot(data=melt(f_table), aes(x=variable, y=value)) +
    geom_boxplot(aes(fill=variable), alpha=0.7) +
    theme_classic() +
    labs(x="Filter", y="Number of variants removed") +
    scale_fill_manual(values=c("aquamarine", "coral", "cadetblue1", "lightpink1", "darkolivegreen1", "lightgoldenrod1", "cornflowerblue", "blueviolet", "indianred1", "lightskyblue3" ))

print(no_jitter)
ggsave("no_jitter.png", plot=no_jitter)
ggsave("no_jitter.pdf", plot=no_jitter)

# Interactive plots
# ggplotly(total_box)
# ggplotly(no_jitter)
