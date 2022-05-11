## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Thu Dec 16 14:59:07 2021
##
## GitHub: https://github.com/aboffelli/
##
## Description: Script to create plots for the ClinVar information in the
##  annotation.
##
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## -----------------------------------------------------------------------------

# TODO: Add comments
library(tidyverse)

pie_chart <- function(file_table, plot_name) {
    colnames(file_table)[1] <- 'V1'
    x <- ggplot(data=file_table, aes(x='', y=Perc, fill=paste0(V1, ' - ', Perc , '%'))) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        labs(title=plot_name)
    return(x)
}

setwd("C:/Users/Arthu/Box/Notes/Tables/SWEA/ClinvarTables")
setwd("~/Box/Notes/Tables/SWEA/ClinvarTables")

pathogenic <- read.table('pathogenic_count_SWEA.txt', sep = '\t')[,-1]
pathogenic <- as.data.frame(table(pathogenic))
pathogenic <- cbind(pathogenic, Perc=round(pathogenic$Freq/sum(pathogenic$Freq)*100, 4))


pathogenic_plot <- pie_chart(pathogenic, "Pathogenic Percentage")
ggsave('Plots/pathogenic_percentage_SWEA.pdf', pathogenic_plot, width=30, height=20, units='cm')
ggsave('Plots/pathogenic_percentage_SWEA.png', pathogenic_plot, width=30, height=20, units='cm')

################################################################################

clinical_type <- read.table('clinical_type_SWEA.txt', sep = '\t')[,-1]
clinical_type <- as.data.frame(table(clinical_type))
clinical_type <- cbind(clinical_type, Perc=round(clinical_type$Freq/sum(clinical_type$Freq)*100, 2))

clinical_type_plot <- pie_chart(clinical_type, 'Clinical Type Percentage')
ggsave('Plots/clinical_type_percentage_SWEA.pdf', clinical_type_plot, width=30, height=20, units='cm')
ggsave('Plots/clinical_type_percentage_SWEA.png', clinical_type_plot, width=30, height=20, units='cm')

################################################################################

most_common <- read_tsv('most_common_pathogenic_var_SWEA.txt', col_names=F) %>% 
    separate(X1, c(NA, NA, NA, "Gene"), "_") %>% 
    group_by(Gene) %>% 
    summarise(n = sum(X2))

most_common_plot <- ggplot(data=most_common, aes(x=reorder(Gene, -n), y=n))+
    geom_bar(stat='identity') +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x='Variant', y='Count', 
         title='Most common genes with pathogenic variants') +
    theme()
ggsave('Plots/most_common_genes_pathogenic_variant_SWEA.pdf', 
       most_common_plot, width=35, height=20, units='cm')
ggsave('Plots/most_common_genes_pathogenic_variant_SWEA.png', 
       most_common_plot, width=35, height=20, units='cm')

################################################################################

samples_pathogenic <- read.table('samples_pathogenic_SWEA.txt', sep='\t')
samples_pathogenic <- as.data.frame(table(samples_pathogenic$V1))

samples_pathogenic <- samples_pathogenic[order(samples_pathogenic$Freq),]
write.table(samples_pathogenic, file='number_of_pathogenic_var_by_sample_SWEA.txt', sep='\t', row.names=F, col.names=F, quote=F)
samples_pathogenic <- as.data.frame(table(samples_pathogenic$Freq))
samples_pathogenic$Var1 <- as.character(samples_pathogenic$Var1)
samples_pathogenic <- rbind(samples_pathogenic, c('0', 3043-sum(samples_pathogenic$Freq)))
samples_pathogenic$Var1 <- as.factor(samples_pathogenic$Var1)
samples_pathogenic$Freq <- as.integer(samples_pathogenic$Freq)
samples_pathogenic <- cbind(samples_pathogenic, Perc=round(samples_pathogenic$Freq/sum(samples_pathogenic$Freq)*100, 2))

samples_pathogenic_plot <- pie_chart(samples_pathogenic, "Pathogenic breast cancer variants by sample")
ggsave('Plots/pathogenic_variants_percentage_SWEA.pdf', samples_pathogenic_plot, width=30, height=20, units='cm')
ggsave('Plots/pathogenic_variants_percentage_SWEA.png', samples_pathogenic_plot, width=30, height=20, units='cm') 
    
    
