## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Thu Dec 16 14:59:07 2021
##
## GitHub: https://github.com/aboffelli/variant_annotation
##
## Description: Script to create plots for the ClinVar information in the
##  annotation for the BRIDGES data. Four plots are created
##
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## -----------------------------------------------------------------------------
library(ggplot2)
library(tidyverse, quietly = T)
library(ggrepel)

pie_chart <- function(file_table, plot_name) {
    colnames(file_table)[2] <- 'V2'
    x <- ggplot(data=file_table, aes(x='', y=Perc, 
                                     fill=V2)) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        labs(title=plot_name) +
        facet_wrap(~V1, ncol = 2) +
        geom_label_repel(data = file_table,
                         aes(y = pos, label = paste0(Perc, "%")),
                         size = 4.5, nudge_x = 0.6, show.legend = FALSE)
        
    return(x)
}

# setwd("C:/Users/Arthu/Box/Notes/Tables/ClinvarTables")
setwd("/Users/student/Box/Notes/TestData/Bridges/FilteredClinVarTables")

sample_type <- "Control"
fam_hist_num <- 3243
total_num <- 53306 - fam_hist_num

# sample_type <- "Samples"
# fam_hist_num <- 10927
# total_num <- 60239 - fam_hist_num

pathogenic <- read.table(paste0(sample_type, '/pathogenic_count.txt'), 
                         sep = '\t')[,-2]
pathogenic <- as.data.frame(table(pathogenic))
pathogenic <- transform(pathogenic, Perc = ave(Freq, V1, FUN = function(x) round(x/sum(x), 4)*100))

family_patho <- subset(pathogenic, V1 == "Family_hist") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_patho <- subset(pathogenic, V1 != "Family_hist") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
pathogenic <- bind_rows(family_patho, no_family_patho)

pathogenic_plot <- pie_chart(pathogenic, "Pathogenic Percentage")
print(pathogenic_plot)
ggsave(paste0(sample_type, '/pathogenic_percentage.pdf'), 
       pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/pathogenic_percentage.png'), 
       pathogenic_plot, width=30,height=20, units='cm')

################################################################################

clinical_type <- read.table(paste0(sample_type, '/clinical_type.txt'), 
                            sep = '\t')[,-2]
clinical_type <- as.data.frame(table(clinical_type))
clinical_type <- transform(clinical_type, Perc = ave(Freq, V1, FUN = function(x) round(x/sum(x), 4)*100))

family_clin <- subset(clinical_type, V1 == "Family_hist") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_clin <- subset(clinical_type, V1 != "Family_hist") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
clinical_type <- bind_rows(family_clin, no_family_clin)

clinical_type_plot <- pie_chart(clinical_type, 'Clinical Type Percentage')
print(clinical_type_plot)
ggsave(paste0(sample_type, '/clinical_type_percentage.pdf'), 
       clinical_type_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/clinical_type_percentage.png'), 
       clinical_type_plot, width=30, height=20, units='cm')

################################################################################

most_common <- read.table(paste0(
    sample_type,'/most_common_pathogenic_var.txt'), sep='\t')

most_common <- most_common %>% separate(V2, c(NA, NA, NA, "V2"), "_") 
most_common <- aggregate(most_common$V3, by=list(Gene=most_common$V2, V1=most_common$V1), FUN=sum)

most_common_plot <- ggplot(data=most_common, 
                           aes(x=reorder(Gene, -x), y=x))+
    geom_bar(stat='identity') +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x='Variant', y='Count', 
         title='Most common genes with pathogenic variants') +
    facet_wrap(~V1, ncol = 1) +
    theme()
print(most_common_plot)
ggsave(paste0(sample_type, '/most_common_genes_pathogenic_variant.pdf'), 
       most_common_plot, width=35, height=20, units='cm')
ggsave(paste0(sample_type, '/most_common_genes_pathogenic_variant.png'),
       most_common_plot, width=35, height=20, units='cm')

################################################################################

samples_pathogenic <- read.table(paste0(
    sample_type, '/samples_pathogenic.txt'), sep='\t')
samples_pathogenic <- as.data.frame(table(samples_pathogenic[,1:2]))

samples_pathogenic <- samples_pathogenic[order(samples_pathogenic$Freq),]
write.table(samples_pathogenic, 
            file=paste0(sample_type, '/number_of_pathogenic_var_by_sample.txt'),
            sep='\t', row.names=F, col.names=F, quote=F)
samples_pathogenic <- as.data.frame(table(samples_pathogenic$V1, samples_pathogenic$Freq))
samples_pathogenic$Freq[samples_pathogenic$Var1 == 'Family_hist' & samples_pathogenic$Var2=='0' ] <- fam_hist_num - sum(samples_pathogenic$Freq[samples_pathogenic$Var1=='Family_hist' & samples_pathogenic$Var2 != '0'])
samples_pathogenic$Freq[samples_pathogenic$Var1 != 'Family_hist' & samples_pathogenic$Var2=='0' ] <- fam_hist_num - sum(samples_pathogenic$Freq[samples_pathogenic$Var1 !='Family_hist' & samples_pathogenic$Var2 != '0'])

samples_pathogenic <- transform(samples_pathogenic, Perc = ave(Freq, Var1, FUN = function(x) round(x/sum(x), 4)*100))
colnames(samples_pathogenic)[1:2] <- c('V1', 'V2' )

family_samples <- subset(samples_pathogenic, V1 == "Family_hist") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_samples <- subset(samples_pathogenic, V1 != "Family_hist") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
samples_pathogenic <- bind_rows(family_samples, no_family_samples)

samples_pathogenic_plot <- pie_chart(
    samples_pathogenic, "Pathogenic breast cancer variants by sample")
print(samples_pathogenic_plot)
ggsave(paste0(sample_type, '/pathogenic_variants_percentage.pdf'), 
       samples_pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/pathogenic_variants_percentage.png'), 
       samples_pathogenic_plot, width=30, height=20, units='cm') 
    
    
