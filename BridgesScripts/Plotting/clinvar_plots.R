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

    x <- ggplot(data=file_table, aes(x=0, y=Perc, 
                                     fill=Type)) +
        theme_void() +
        geom_col(col='black', size=0.05) +
        coord_polar(theta = 'y') +
        scale_fill_discrete(name='') +
        labs(title=plot_name) +
        facet_wrap(~Sample_type + Hist + Hist50, ncol = 2,
                   strip.position = "bottom") +
        
        # Label flags with the percentage.
        geom_label_repel(data = file_table,
                         aes(y = pos, label = paste0(Perc, "%")),
                         size = 4.5, nudge_x = 1, show.legend = FALSE,
                         max.overlaps = 30)
        
    return(x)
}

setwd("/Users/student/Box/Notes/Tables/BRIDGES/FilteredClinVarTables")
# TODO: Add the controls in the same plot.

# Since the number of Control and Cases are different, and the tables will be
# in their respect directories, un-comment the group you are using before
#  continuing.

fam_hist_num <- 11518
total_num <- 60239 - fam_hist_num
c_total_num <- 53306

## -----------------------------------------------------------------------------
## Pie chart of the percentage of pathogenic variants.

case_pathogenic <- read_tsv('Cases/Tables/pathogenic_count_BRIDGES.txt', 
                         col_names = F) %>% 
    select(-X2) %>%
    
    # Add a column that differentiates only as Family history or 
    # No family history, since we are no interested in the type of family 
    # history right now.
    mutate(Hist = if_else(X1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    
    # Add a column with family history below 50.
    mutate(Hist50=if_else(str_detect(X1, "50BC"), "Family history 50", 
                          '')) %>% 
    select(-X1) %>% 
    
    # Count the frequency
    count(X3, Hist, Hist50) %>% 
    rename(Type = X3, Freq = n) %>% 
    
    # Add Cases column.
    add_column(Sample_type = "Cases", .before = "Type") %>% 
    
    # Add the percentage column. Based on the family history.
    transform(Perc = ave(Freq, Hist, Hist50,
              FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    
    # Calculate the position of the flags grouping by family history.
    group_by(Hist, Hist50) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
    

control_pathogenic <- read_tsv('Controls/Tables/pathogenic_count_BRIDGES.txt',
                               col_names = F) %>%
    select(X3) %>%
    
    # Count the frequency
    count(X3) %>% 
    rename(Type = X3, Freq = n) %>% 
    
    # Add the columns
    add_column(Sample_type = "Controls", .before = "Type") %>% 
    add_column(Hist = '', Hist50 = '', .after="Type") %>% 
    
    # Calculate the percentage.
    transform(Perc = ave(Freq, Hist, Hist50,
                         FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    
    # Calculate the position of the flags.
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))


# Join Cases and Controls tables.
pathogenic <- case_pathogenic %>% full_join(control_pathogenic) 
    
# Create the plot and save it as a pdf and as a png.
pathogenic_plot <- pie_chart(pathogenic, "Pathogenic Percentage")

print(pathogenic_plot)

ggsave('Plots/pathogenic_percentage_BRIDGES.pdf', 
       pathogenic_plot, width=30, height=20, units='cm')
ggsave('Plots/pathogenic_percentage_BRIDGES.png', 
       pathogenic_plot, width=30,height=20, units='cm')

## -----------------------------------------------------------------------------
## Pie chart for the percentage of types of clinical relevance in all variants.


cases_clinical_type <- read_tsv('Cases/Tables/clinical_type_BRIDGES.txt', 
                            col_names = F) %>% 
    select(-X2) %>%
    
    # Add a column that differentiates only as Family history or 
    # No family history, since we are no interested in the type of family 
    # history right now.
    mutate(Hist = if_else(X1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    
    # Add a column with family history below 50.
    mutate(Hist50=if_else(str_detect(X1, "50BC"), "Family history 50", 
                          '')) %>% 
    select(-X1) %>% 
    
    # Count the frequency
    count(X3, Hist, Hist50) %>% 
    rename(Type = X3, Freq = n) %>% 
    
    # Add Cases column.
    add_column(Sample_type = "Cases", .before = "Type") %>% 
    
    # Add the percentage column. Based on the family history.
    transform(Perc = ave(Freq, Hist, Hist50,
                         FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    
    # Calculate the position of the flags grouping by family history.
    group_by(Hist, Hist50) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

controls_clinical_type <- read_tsv('Controls/Tables/clinical_type_BRIDGES.txt', 
                                   col_names = F) %>% 
    select(X3) %>%
    
    # Count the frequency
    count(X3) %>% 
    rename(Type = X3, Freq = n) %>% 
    
    # Add the columns
    add_column(Sample_type = "Controls", .before = "Type") %>% 
    add_column(Hist = '', Hist50 = '', .after="Type") %>% 
    
    # Calculate the percentage.
    transform(Perc = ave(Freq, Hist, Hist50,
                         FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    
    # Calculate the position of the flags.
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

clinical_type <- cases_clinical_type %>% full_join(controls_clinical_type)

# Create the plot and save it as a pdf and as a png.
clinical_type_plot <- pie_chart(clinical_type, 
                               'Clinical Type Percentage')
print(clinical_type_plot)

ggsave(paste0(sample_type, '/Plots/clinical_type_percentage_BRIDGES.pdf'), 
       clinical_type_plot, width=30, height=20, units='cm')
ggsave(paste0(sample_type, '/Plots/clinical_type_percentage_BRIDGES.png'), 
       clinical_type_plot, width=30, height=20, units='cm')

## -----------------------------------------------------------------------------
## Histogram showing which genes contatain more pathogenic variants.
 
cases_most_common <- read_tsv('Cases/Tables/most_common_pathogenic_var_BRIDGES.txt',
                        col_names = F) %>%
    
    # Add a column that differentiates only as Family history or 
    # No family history, since we are no interested in the type of family 
    # history right now.
    mutate(Hist = if_else(X1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    
    # Add a column with family history below 50.
    mutate(Hist50=if_else(str_detect(X1, "50BC"), "Family history 50", 
                          '')) %>% 
    select(-X1) %>% 
    
    # Since we have too many variants, we can isolate just the gene name to see 
    # which genes have more pathogenic variants.
    separate(X2, c(NA, NA, NA, "Gene"), "_") %>%
    
    # Sum the counts for the same gene and family history.
    group_by(Gene, Hist, Hist50) %>% 
    summarise(x = sum(X3)) %>% 
    add_column(Sample_type = "Cases", .before = "Gene")

controls_most_commom <- read_tsv("Controls/Tables/most_common_pathogenic_var_BRIDGES.txt",
                                 col_names = F) %>% 
    select(-X1) %>% 
    separate(X2, c(NA, NA, NA, "Gene"), "_") %>%
    group_by(Gene) %>% 
    summarise(x = sum(X3)) %>% 
    # Add the columns
    add_column(Sample_type = "Controls", .before = "Gene") %>% 
    add_column(Hist = '', Hist50 = '', .after="Gene")

most_common <- cases_most_common %>% full_join(controls_most_commom)

# Create the bar plot, and save it as a pfd and as a png.
most_common_plot <- ggplot(data=most_common, 
                           aes(x=reorder(Gene, -x), y=x))+
    geom_bar(stat='identity') +
    theme_classic() +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x='Variant', y='Count', 
         title='Most common genes with pathogenic variants') +
    facet_wrap(~Sample_type + Hist + Hist50, ncol = 1) +
    theme()

print(most_common_plot)

ggsave('Plots/most_common_genes_pathogenic_variant_BRIDGES.pdf', 
       most_common_plot, width=35, height=20, units='cm')
ggsave('Plots/most_common_genes_pathogenic_variant_BRIDGES.png',
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
    
    
