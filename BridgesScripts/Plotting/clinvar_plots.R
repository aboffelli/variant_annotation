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
library(tidyverse, quietly = T)
library(ggrepel)

# TODO: Add comments.
pie_chart <- function(file_table, plot_name) {
    x <- ggplot(data=file_table, aes(x='', y=Perc, 
                                     fill=Type)) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        labs(title=plot_name) +
        facet_wrap(~Hist, ncol = 2) +
        geom_label_repel(data = file_table,
                         aes(y = pos, label = paste0(Perc, "%")),
                         size = 4.5, nudge_x = 0.6, show.legend = FALSE,
                         max.overlaps = 30)
        
    return(x)
}

# setwd("C:/Users/Arthu/Box/Notes/Tables/ClinvarTables")
setwd("/Users/student/Box/Notes/Tables/BRIDGES/FilteredClinVarTables")

sample_type <- "Controls"
fam_hist_num <- 3698
total_num <- 53306 - fam_hist_num

sample_type <- "Cases"
fam_hist_num <- 11518
total_num <- 60239 - fam_hist_num

pathogenic <- read.table(paste0(sample_type, 
                                '/Tables/pathogenic_count_BRIDGES.txt'), 
                         sep = '\t')[,-2]
colnames(pathogenic)[2] <- 'Type'
pathogenic <- pathogenic %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    select(-V1)

pathogenic <- as.data.frame(table(pathogenic))
pathogenic <- transform(pathogenic, Perc = ave(Freq, Hist, FUN = function(x) round(x/sum(x), 4)*100))

family_patho <- subset(pathogenic, Hist == "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_patho <- subset(pathogenic, Hist != "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
pathogenic <- bind_rows(family_patho, no_family_patho)


pathogenic_plot <- pie_chart(pathogenic, paste("Pathogenic Percentage -",
                                               sample_type))

print(pathogenic_plot)

ggsave(paste0(sample_type, '/Plots/pathogenic_percentage_BRIDGES.pdf'), 
       pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/pathogenic_percentage_BRIDGES.png'), 
       pathogenic_plot, width=30,height=20, units='cm')

################################################################################

clinical_type <- read.table(paste0(sample_type, 
                                   '/Tables/clinical_type_BRIDGES.txt'), 
                            sep = '\t')[,-2]
colnames(clinical_type)[2] <- "Type"
clinical_type <- clinical_type %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    select(-V1)
clinical_type <- as.data.frame(table(clinical_type))
clinical_type <- transform(clinical_type, 
                           Perc = ave(Freq, Hist, 
                                      FUN = function(x) round(x/sum(x), 4)*100))

family_clin <- subset(clinical_type, Hist == "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_clin <- subset(clinical_type, Hist != "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
clinical_type <- bind_rows(family_clin, no_family_clin)

clinical_type_plot <- pie_chart(clinical_type, 
                                paste('Clinical Type Percentage -',
                                                     sample_type))
print(clinical_type_plot)
ggsave(paste0(sample_type, '/Plots/clinical_type_percentage_BRIDGES.pdf'), 
       clinical_type_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/clinical_type_percentage_BRIDGES.png'), 
       clinical_type_plot, width=30, height=20, units='cm')

################################################################################

most_common <- read.table(paste0(
    sample_type,'/Tables/most_common_pathogenic_var_BRIDGES.txt'), sep='\t')

most_common <- most_common %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    select(-V1)
most_common <- most_common %>% separate(V2, c(NA, NA, NA, "V2"), "_") 
most_common <- aggregate(most_common$V3, 
                         by=list(Gene=most_common$V2,
                                 Hist = most_common$Hist), FUN=sum)

most_common_plot <- ggplot(data=most_common, 
                           aes(x=reorder(Gene, -x), y=x))+
    geom_bar(stat='identity') +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x='Variant', y='Count', 
         title=paste('Most common genes with pathogenic variants -',
                     sample_type)) +
    facet_wrap(~Hist, ncol = 1) +
    theme()
print(most_common_plot)
ggsave(paste0(sample_type, '/Plots/most_common_genes_pathogenic_variant_BRIDGES.pdf'), 
       most_common_plot, width=35, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/most_common_genes_pathogenic_variant_BRIDGES.png'),
       most_common_plot, width=35, height=20, units='cm')

################################################################################

samples_pathogenic <- read.table(paste0(
    sample_type, '/Tables/samples_pathogenic_BRIDGES.txt'), sep='\t')

samples_pathogenic <- samples_pathogenic %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history"))

samples_with_patho <- samples_pathogenic %>% 
    select(V2, Hist) %>% 
    distinct()
    
samples_with_patho_count <- tibble(
    Type=c("Pathogenic family history",
           "Pathogenic no family history",
           "No pathogenic family history",
           "No pathogenic no family history"),
    
    Count=c(nrow(filter(samples_with_patho, Hist == "Family history")),
            nrow(filter(samples_with_patho, Hist != "Family history")),
            fam_hist_num - nrow(filter(samples_with_patho, 
                                       Hist == "Family history")),
            total_num - nrow(filter(samples_with_patho, 
                                    Hist != "Family history"))
            )
    )

samples_with_patho_count <- samples_with_patho_count %>% 
    transform(Perc = ave(Count,
                         FUN = function(x) round(x/sum(x),
                                                 4)*100)) %>% 
    arrange(Type) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
    


samples_with_patho_plot <- ggplot(data=samples_with_patho_count, 
                                  aes(x='', y=Perc, fill=Type)) +
    theme_void() +
    geom_col(col='black', size=0.05) +
    coord_polar(theta = 'y') +
    scale_fill_discrete(name='') +
    labs(title=paste("Percentage of samples with pathogenic variants -",
                     sample_type)) +
    geom_label_repel(data = samples_with_patho_count,
                     aes(y = pos, label = paste0(Perc, "%")),
                     size = 4.5, nudge_x = 0.6, show.legend = FALSE,
                     max.overlaps = 20)

print(samples_with_patho_plot)
ggsave(paste0(sample_type, '/Plots/percentage_samples_pathogenic_BRIDGES.pdf'), 
       samples_with_patho_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/percentage_samples_pathogenic_BRIDGES.png'), 
       samples_with_patho_plot, width=30, height=20, units='cm')
    

################################################################################

samples_pathogenic <- as.data.frame(table(samples_pathogenic[,c(1,2,9)])) %>% 
    filter(Freq > 0) %>% 
    select(-Freq)



# samples_pathogenic %>% as_tibble %>% select(V2, Freq) %>% group_by(V2) %>% nest %>%
#     mutate(any_pathologic = map_lgl(data, function(df){}))
# Purrr

# samples_pathogenic <- samples_pathogenic[order(samples_pathogenic$Freq),]
# write.table(samples_pathogenic,
#             file=paste0(sample_type, '/number_of_pathogenic_var_by_sample_BRIDGES.txt'),
#             sep='\t', row.names=F, col.names=F, quote=F)
samples_pathogenic <- as.data.frame(table(samples_pathogenic$V1, 
                                          samples_pathogenic$Hist)) %>% 
    filter(Freq > 0)

colnames(samples_pathogenic)[1:2] <- c('Type', 'Hist')
samples_pathogenic <- transform(samples_pathogenic,
                                Perc = ave(Freq,
                                           FUN = function(x) round(x/sum(x),
                                                                   4)*100))


family_samples <- samples_pathogenic %>%
    filter(Hist == "Family history") %>%
    transform(Perc = ave(Freq, FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

no_family_samples <- samples_pathogenic %>%
    filter(Hist == "No family history") %>% 
    mutate(Perc=Freq/total_num*100) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

samples_pathogenic <- bind_rows(family_samples, no_family_samples)

samples_pathogenic_plot <- pie_chart(
    family_samples, paste("Type of family history with patogenic variants -",
                              sample_type))
print(samples_pathogenic_plot)
ggsave(paste0(sample_type, '/Plots/type_history_percentage_BRIDGES.pdf'), 
       samples_pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/type_history_percentage_BRIDGES.png'), 
       samples_pathogenic_plot, width=30, height=20, units='cm') 
    
    
