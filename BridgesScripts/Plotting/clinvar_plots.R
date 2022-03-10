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

pie_chart <- function(file_table, plot_name) {
    colnames(file_table)[1] <- 'V1'
    x <- ggplot(data=file_table, aes(x='', y=Perc, 
                                     fill=paste0(
                                         reorder(V1, sort(as.numeric(V1)) -1),
                                                 ' - ', Perc , '%'))) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        labs(title=plot_name)
    return(x)
}

# setwd("C:/Users/Arthu/Box/Notes/Tables/ClinvarTables")
setwd("/Users/student/Box/Notes/Tables/BRIDGES/ClinVarTables")

sample_type <- "Controls/BCFH"
total_num <- 3243
# total_num <- 53306

# sample_type <- "Cases"
# total_num <- 10927
# total_num <- 60239

pathogenic <- read.table(paste0(sample_type, '/pathogenic_count.txt'), 
                         sep = '\t')[,-1]
pathogenic <- as.data.frame(table(pathogenic))
pathogenic <- cbind(pathogenic, 
                    Perc=round(pathogenic$Freq/sum(pathogenic$Freq)*100, 3))


pathogenic_plot <- pie_chart(pathogenic, "Pathogenic Percentage")
ggsave(paste0(sample_type, '/pathogenic_percentage.pdf'), 
       pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/pathogenic_percentage.png'), 
       pathogenic_plot, width=30,height=20, units='cm')

################################################################################

clinical_type <- read.table(paste0(sample_type, '/clinical_type.txt'), 
                            sep = '\t')[,-1]
clinical_type <- as.data.frame(table(clinical_type))
clinical_type <- cbind(
    clinical_type, 
    Perc=round(clinical_type$Freq/sum(clinical_type$Freq)*100, 3))

clinical_type_plot <- pie_chart(clinical_type, 'Clinical Type Percentage')
ggsave(paste0(sample_type, '/clinical_type_percentage.pdf'), 
       clinical_type_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/clinical_type_percentage.png'), 
       clinical_type_plot, width=30, height=20, units='cm')

################################################################################

most_common <- read.table(paste0(
    sample_type,'/most_common_pathogenic_var.txt'), sep='\t')

most_common_plot <- ggplot(data=most_common[most_common$V2>5,], 
                           aes(x=reorder(V1, -V2), y=V2))+
    geom_bar(stat='identity') +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x='Variant', y='Count', 
         title='Most common pathogenic breast cancer variants') +
    theme()
ggsave(paste0(sample_type, '/most_common_pathogenic_variant.pdf'), 
       most_common_plot, width=35, height=20, units='cm')
ggsave(paste0(sample_type, '/most_common_pathogenic_variant.png'),
       most_common_plot, width=35, height=20, units='cm')

################################################################################

samples_pathogenic <- read.table(paste0(
    sample_type, '/samples_pathogenic.txt'), sep='\t')
samples_pathogenic <- as.data.frame(table(samples_pathogenic$V1))

samples_pathogenic <- samples_pathogenic[order(samples_pathogenic$Freq),]
write.table(samples_pathogenic, 
            file=paste0(sample_type, '/number_of_pathogenic_var_by_sample.txt'),
            sep='\t', row.names=F, col.names=F, quote=F)
samples_pathogenic <- as.data.frame(table(samples_pathogenic$Freq))
samples_pathogenic$Var1 <- as.character(samples_pathogenic$Var1)
samples_pathogenic <- rbind(samples_pathogenic, 
                            c('0', total_num-sum(samples_pathogenic$Freq)))
samples_pathogenic$Var1 <- factor(samples_pathogenic$Var1)
samples_pathogenic$Freq <- as.integer(samples_pathogenic$Freq)
samples_pathogenic <- cbind(
    samples_pathogenic, 
    Perc=round(samples_pathogenic$Freq/sum(samples_pathogenic$Freq)*100, 3))

samples_pathogenic_plot <- pie_chart(
    samples_pathogenic, "Pathogenic breast cancer variants by sample")
ggsave(paste0(sample_type, '/pathogenic_variants_percentage.pdf'), 
       samples_pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/pathogenic_variants_percentage.png'), 
       samples_pathogenic_plot, width=30, height=20, units='cm') 
    
    
