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

pie_chart <- function(file_table, plot_name) {
    ## -------------------------------------------------------------------------
    ## Function to create a pie chart using ggplot2.
    ## -------------------------------------------------------------------------

    x <- ggplot(data=file_table, aes(x='', y=Perc, 
                                     fill=Type)) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        labs(title=plot_name) +
        facet_wrap(~Hist, ncol = 2) +
        
        # Label flags with the percentage.
        geom_label_repel(data = file_table,
                         aes(y = pos, label = paste0(Perc, "%")),
                         size = 4.5, nudge_x = 0.6, show.legend = FALSE,
                         max.overlaps = 30)
        
    return(x)
}

setwd("/Users/student/Box/Notes/Tables/BRIDGES/FilteredClinVarTables")

# Since the number of Control and Cases are different, and the tables will be
# in their respect directories, un-comment the group you are using before
#  continuing.

sample_type <- "Controls"
fam_hist_num <- 3698
total_num <- 53306 - fam_hist_num

# sample_type <- "Cases"
# fam_hist_num <- 11518
# total_num <- 60239 - fam_hist_num


## -----------------------------------------------------------------------------
## Pie chart of the percentage of pathogenic variants.

pathogenic <- read.table(paste0(sample_type, 
                                '/Tables/pathogenic_count_BRIDGES.txt'), 
                         sep = '\t')[,-2]

colnames(pathogenic)[2] <- 'Type'

# Add a column that differentiates only as Family history or No family history,
# since we are no interested in the type of family history right now.
pathogenic <- pathogenic %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    select(-V1)

# Count the unique combinations.
pathogenic <- as.data.frame(table(pathogenic))

# Add the percentage column. Based on the family history.
pathogenic <- transform(pathogenic, Perc = ave(Freq, Hist, FUN = function(x) round(x/sum(x), 4)*100))

# Since we calculated the percentage based on the Family history, we split the
# data frame to calculate the positions of the flags in the pie chart.
family_patho <- subset(pathogenic, Hist == "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_patho <- subset(pathogenic, Hist != "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

# Join the dataframes together again.
pathogenic <- bind_rows(family_patho, no_family_patho)


# Create the plot and save it as a pdf and as a png.
pathogenic_plot <- pie_chart(pathogenic, paste("Pathogenic Percentage -",
                                               sample_type))

print(pathogenic_plot)

ggsave(paste0(sample_type, '/Plots/pathogenic_percentage_BRIDGES.pdf'), 
       pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/pathogenic_percentage_BRIDGES.png'), 
       pathogenic_plot, width=30,height=20, units='cm')

## -----------------------------------------------------------------------------
## Pie chart for the percentage of types of clinical relevance in all variants.


clinical_type <- read.table(paste0(sample_type, 
                                   '/Tables/clinical_type_BRIDGES.txt'), 
                            sep = '\t')[,-2]
colnames(clinical_type)[2] <- "Type"

# Add a column that differentiates only as Family history or No family history,
# since we are no interested in the type of family history right now.
clinical_type <- clinical_type %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    select(-V1)

# Count the unique combinations.
clinical_type <- as.data.frame(table(clinical_type))

# Add the percentage column. Based on the family history.
clinical_type <- transform(clinical_type, 
                           Perc = ave(Freq, Hist, 
                                      FUN = function(x) round(x/sum(x), 4)*100))

# Since we calculated the percentage based on the Family history, we split the
# data frame to calculate the positions of the flags in the pie chart.
family_clin <- subset(clinical_type, Hist == "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
no_family_clin <- subset(clinical_type, Hist != "Family history") %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
# Join the dataframes together again.
clinical_type <- bind_rows(family_clin, no_family_clin)

# Create the plot and save it as a pdf and as a png.
clinical_type_plot <- pie_chart(clinical_type, 
                                paste('Clinical Type Percentage -',
                                                     sample_type))
print(clinical_type_plot)
ggsave(paste0(sample_type, '/Plots/clinical_type_percentage_BRIDGES.pdf'), 
       clinical_type_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/clinical_type_percentage_BRIDGES.png'), 
       clinical_type_plot, width=30, height=20, units='cm')

## -----------------------------------------------------------------------------
## Histogram showing which genes contatain more pathogenic variants.
 
most_common <- read.table(paste0(
    sample_type,'/Tables/most_common_pathogenic_var_BRIDGES.txt'), sep='\t')

# Add a column that differentiates only as Family history or No family history,
# since we are no interested in the type of family history right now.
most_common <- most_common %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    select(-V1)

# Since we have too many variants, we can isolate just the gene name to see 
# which genes have more pathogenic variants. 
most_common <- most_common %>% separate(V2, c(NA, NA, NA, "V2"), "_")

# Sum the counts for the same gene and family history.
most_common <- aggregate(most_common$V3, 
                         by=list(Gene=most_common$V2,
                                 Hist = most_common$Hist), FUN=sum)

# Create the bar plot, and save it as a pfd and as a png.
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

## -----------------------------------------------------------------------------
## Pie chart with the percentage of samples that contain at least one pathogenic
## variant.

samples_pathogenic <- read.table(paste0(
    sample_type, '/Tables/samples_pathogenic_BRIDGES.txt'), sep='\t')

# Add a column that differentiates only as Family history or No family history,
# since we are no interested in the type of family history right now.
samples_pathogenic <- samples_pathogenic %>% 
    mutate(Hist = if_else(V1 !="No_family_hist", "Family history", 
                          "No family history"))

# Get only the unique samples.
samples_with_patho <- samples_pathogenic %>% 
    select(V2, Hist) %>% 
    distinct()

# Create a new data frame with the counts of samples.    
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

# Add the percentage column and the flag position.
samples_with_patho_count <- samples_with_patho_count %>% 
    transform(Perc = ave(Count,
                         FUN = function(x) round(x/sum(x),
                                                 4)*100)) %>% 
    arrange(Type) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
    

# Since this plot is not divided into family history and non family history, 
# we can't you the pie chart function.
samples_with_patho_plot <- ggplot(data=samples_with_patho_count, 
                                  aes(x='', y=Perc, fill=Type)) +
    theme_void() +
    geom_col(col='black', size=0.05) +
    coord_polar(theta = 'y') +
    scale_fill_discrete(name='') +
    labs(title=paste("Percentage of samples with pathogenic variants -",
                     sample_type)) +
    # Add the flags.
    geom_label_repel(data = samples_with_patho_count,
                     aes(y = pos, label = paste0(Perc, "%")),
                     size = 4.5, nudge_x = 0.6, show.legend = FALSE,
                     max.overlaps = 20)

# Save the plot as a pdf and as a png.
print(samples_with_patho_plot)
ggsave(paste0(sample_type, '/Plots/percentage_samples_pathogenic_BRIDGES.pdf'), 
       samples_with_patho_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/percentage_samples_pathogenic_BRIDGES.png'), 
       samples_with_patho_plot, width=30, height=20, units='cm')
    

## -----------------------------------------------------------------------------
## Plot with the percentage of familial history type that contain pathogenic 
## variants. 

# Use the same table as before and count the unique combinations.
samples_pathogenic <- as.data.frame(table(samples_pathogenic[,c(1,2,9)])) %>% 
    filter(Freq > 0) %>% 
    select(-Freq)

# Get the frequency for each type of family history.
samples_pathogenic <- as.data.frame(table(samples_pathogenic$V1, 
                                          samples_pathogenic$Hist)) %>% 
    filter(Freq > 0)

colnames(samples_pathogenic)[1:2] <- c('Type', 'Hist')

# Add the percentage column
samples_pathogenic <- transform(samples_pathogenic,
                                Perc = ave(Freq,
                                           FUN = function(x) round(x/sum(x),
                                                                   4)*100))

# Since we calculated the percentage based on the Family history, we split the
# data frame to calculate the positions of the flags in the pie chart.
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

# Join the data frames together again.
samples_pathogenic <- bind_rows(family_samples, no_family_samples)


# Create the plot and save it as a pdf and as a png.
samples_pathogenic_plot <- pie_chart(
    family_samples, paste("Type of family history with patogenic variants -",
                              sample_type))
print(samples_pathogenic_plot)
ggsave(paste0(sample_type, '/Plots/type_history_percentage_BRIDGES.pdf'), 
       samples_pathogenic_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/type_history_percentage_BRIDGES.png'), 
       samples_pathogenic_plot, width=30, height=20, units='cm') 
    
    
