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

setwd("C:/Users/Arthu/Box/Notes/Tables/BRIDGES/FilteredClinVarTables")
setwd("/Users/student/Box/Notes/Tables/BRIDGES/FilteredClinVarTables")


## Functions ----
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

## Pathogenic variants -------------------------------------------------------
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
pathogenic <- case_pathogenic %>% full_join(control_pathogenic) %>% 
    mutate(Type = case_when(
        str_detect(Type, "Pathogenic") ~ "Pathogenic/Likely pathogenic",
        TRUE ~ Type))
    
# Create the plot and save it as a pdf and as a png.
pathogenic_plot <- pie_chart(pathogenic, "Percentage of pathogenic variants")

print(pathogenic_plot)

ggsave('Plots/pathogenic_percentage_BRIDGES.pdf', 
       pathogenic_plot, width=30, height=20, units='cm')
ggsave('Plots/pathogenic_percentage_BRIDGES.png', 
       pathogenic_plot, width=30,height=20, units='cm')

## Clinical relevance ---------------------------------------------------------
## Pie chart for the percentage of types of clinical relevance in all variants.


cases_clinical_type <- read_tsv('Cases/Tables/clinical_type_BRIDGES.txt', 
                            col_names = F) %>% 
    select(-X2) %>%
    
    # Add a column that differentiates only as Family history or 
    # No family history, since we are no interested in the type of family 
    # history right now. Add another column with family history below 50.
    mutate(Hist = if_else(X1 !="No_family_hist", "Family history", 
                          "No family history"),
           Hist50=if_else(str_detect(X1, "50BC"), "Family history 50", ''),
           Type=case_when(
               str_detect(X3, 
                          "Benign|Likely_benign") ~ "Benign/Likely benign",
               str_detect(X3, 
                          "Pathogenic|Likely_pathogenic") ~ 
                   "Pathogenic/Likely pathogenic",
               str_detect(X3, "Uncertain") ~ "Uncertain significance",
               str_detect(X3, "Conflicting") ~ "Conflicting interpretations",
               str_detect(X3, "drug") ~"Drug response",
               TRUE ~ "Not provided/No interpretation")) %>%

    select(-X1, -X3) %>% 
    
    # Count the frequency
    count(Type, Hist, Hist50) %>% 
    
    # Add Cases column.
    add_column(Sample_type = "Cases", .before = "Type") %>% 
    
    # Add the percentage column. Based on the family history.
    transform(Perc = ave(n, Hist, Hist50,
                         FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    
    # Calculate the position of the flags grouping by family history.
    group_by(Hist, Hist50) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

controls_clinical_type <- read_tsv('Controls/Tables/clinical_type_BRIDGES.txt', 
                                   col_names = F) %>% 
    mutate(Type=case_when(
        str_detect(X3, 
                   "Benign|Likely_benign") ~ "Benign/Likely benign",
        str_detect(X3, 
                   "Pathogenic|Likely_pathogenic") ~ 
            "Pathogenic/Likely pathogenic",
        str_detect(X3, "Uncertain") ~ "Uncertain significance",
        str_detect(X3, "Conflicting") ~ "Conflicting interpretations",
        str_detect(X3, "drug") ~"Drug response",
        TRUE ~ "Not provided/No interpretation")) %>% 
    
    # Count the frequency
    count(Type) %>% 
    
    # Add the columns
    add_column(Sample_type = "Controls", .before = "Type") %>% 
    add_column(Hist = '', Hist50 = '', .after="Type") %>% 
    
    # Calculate the percentage.
    transform(Perc = ave(n, Hist, Hist50,
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

ggsave('Plots/v2_clinical_type_percentage_BRIDGES.pdf', 
       clinical_type_plot, width=30, height=20, units='cm')
ggsave('Plots/v2_clinical_type_percentage_BRIDGES.png', 
       clinical_type_plot, width=30, height=20, units='cm')

## Most common genes -------------------------------------------------------
## Histogram showing which genes containing more pathogenic variants.
 
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
    facet_wrap(~Sample_type + Hist + Hist50, ncol = 1,
               strip.position = "bottom") +
    theme()

print(most_common_plot)

ggsave('Plots/most_common_genes_pathogenic_variant_BRIDGES.pdf', 
       most_common_plot, width=35, height=20, units='cm')
ggsave('Plots/most_common_genes_pathogenic_variant_BRIDGES.png',
       most_common_plot, width=35, height=20, units='cm')

## Samples with pathogenic ---------------------------------------------------

# Number of samples in each type.
total_n_samples <- read_tsv("Cases/family_history_cases.txt", col_names = F) %>% 
    mutate(Type=case_when(
        str_detect(X2, "50BC") ~ "Family history BC before 50 yo",
        str_detect(X2, "FirstBC") ~ "First Degree BC",
        str_detect(X2, "SecondBC") ~ "Second degree BC",
        str_detect(X2, "50OC") ~ "Family history OC before 50 yo",
        str_detect(X2, "FirstOC") ~ "First degree OC",
        str_detect(X2, "SecondOC") ~ "Second degree OC")) %>% 
    select(Type) %>% 
    count(Type)

total_n <- 60239    
fam_hist_num <- 11518
no_fam_hist_num <- 60239 - fam_hist_num
c_total_num <- 53306


cases_samples_pathogenic <- read_tsv('Cases/Tables/samples_pathogenic_BRIDGES.txt',
                                 col_names = F) %>% 
    
    # Add a column that differentiates only as Family history or No family history,
    # since we are no interested in the type of family history right now.
    mutate(Hist = if_else(X1 !="No_family_hist", "Family history", 
                          "No family history")) %>% 
    
    # Add a column with family history below 50.
    mutate(Hist50=if_else(str_detect(X1, "50BC"), "Family history 50", ''))

controls_samples_pathogenic <- read_tsv(
    "Controls/Tables/samples_pathogenic_BRIDGES.txt", col_names = F)



# Get only the unique samples.
cases_samples_with_patho <- cases_samples_pathogenic %>% 
    select(X2, Hist, Hist50) %>% 
    distinct() %>% 
    count(Hist, Hist50)

controls_with_patho <- controls_samples_pathogenic %>% 
    select(X2) %>% 
    distinct() %>% 
    count()

## Pathogenic per sample ----
## Pie chart with the percentage of samples that contain at least one pathogenic
## variant.
# Create a new data frame with the counts of samples.    
cases_samples_with_patho_count <- tibble(
    Sample_type = "Cases",
    Type=c("Pathogenic with family history",
           "Pathogenic with family history before 50 yo",
           "Pathogenic without family history",
           "No pathogenic with family history",
           "No pathogenic without family history"),
    
    Count=c(cases_samples_with_patho$n,
            fam_hist_num - sum(filter(cases_samples_with_patho, 
                                       Hist == "Family history")$n),
            no_fam_hist_num - sum(filter(cases_samples_with_patho, 
                                    Hist != "Family history")$n))) %>% 
    
    # Add the percentage column and the flag position.
    transform(Perc = ave(Count,
                         FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    arrange(Type) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))
 
controls_with_patho_count <- tibble(
    Sample_type = "Controls",
    Type = c("Pathogenic", "No pathogenic"),
    Count= c(controls_with_patho$n,
             c_total_num - controls_with_patho$n)) %>% 
    
    # Add the percentage column and the flag position.
    transform(Perc = ave(Count,
                         FUN = function(x) round(x/sum(x), 4)*100)) %>% 
    arrange(Type) %>% 
    mutate(csum = rev(cumsum(rev(Perc))), 
           pos = Perc/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Perc/2, pos))

samples_with_patho_count <- cases_samples_with_patho_count %>% 
    full_join(controls_with_patho_count)
    
    

# Since this plot is not divided into family history, 
# we can't use the pie chart function.
samples_with_patho_plot <- ggplot(data=samples_with_patho_count, 
                                  aes(x='', y=Perc, fill=Type)) +
    theme_void() +
    geom_col(col='black', size=0.05) +
    coord_polar(theta = 'y') +
    scale_fill_discrete(name='') +
    labs(title="Percentage of samples with pathogenic variants") +
    facet_wrap(~Sample_type, strip.position = "bottom") +
    # Add the flags.
    geom_label_repel(data = samples_with_patho_count,
                     aes(y = pos, label = paste0(Perc, "%")),
                     size = 4.5, nudge_x = 0.6, show.legend = FALSE,
                     max.overlaps = 20)

# Save the plot as a pdf and as a png.
print(samples_with_patho_plot)
ggsave('Plots/percentage_samples_pathogenic_BRIDGES.pdf', 
       samples_with_patho_plot, width=30, height=20, units='cm')
ggsave('Plots/percentage_samples_pathogenic_BRIDGES.png', 
       samples_with_patho_plot, width=30, height=20, units='cm')
    

## Pathogenic per family history type --------------------------------------
## Plot with the percentage of familial history type that contain pathogenic 
## variants. 

# TODO: Change to barplot with yes or no pathogenic variants based on the total
# number of samples for each type of family history.

history_type <- cases_samples_pathogenic  %>%
    
    # Filter out the No family history
    filter(Hist == "Family history") %>%
    
    # Set the priority of Family history types.
    mutate(Type=case_when(
        str_detect(X1, "50BC") ~ "Family history BC before 50 yo",
        str_detect(X1, "FirstBC") ~ "First Degree BC",
        str_detect(X1, "SecondBC") ~ "Second degree BC",
        str_detect(X1, "50OC") ~ "Family history OC before 50 yo",
        str_detect(X1, "FirstOC") ~ "First degree OC",
        str_detect(X1, "SecondOC") ~ "Second degree OC")) %>% 

    # Get unique samples
    distinct(X2, Type) %>%
    
    # Count the frequency of family history type.
    count(Type) %>% 
    
    add_column(Pathogenic="Yes", .before = "Type")

history_type <- history_type %>% 
    add_row(Pathogenic = "No",
            Type = total_n_samples$Type,
            n = total_n_samples$n - history_type$n) %>% 
    add_row(Pathogenic = c("Yes", "No"),
            Type = "Controls",
            n = c(controls_with_patho$n,
                  c_total_num - controls_with_patho$n)) %>%
    group_by(Type) %>% 
    # Add the percentage and position of the flags.
    mutate(Perc = round(n/sum(n), 4)*100,
           label_y = cumsum(Perc)) %>% 
    arrange(Type)
    
# Create the plot and save it as a pdf and as a png.
samples_pathogenic_plot <- ggplot(data=history_type, aes(x=Type, 
                                                         y=Perc,
                                                         fill=Pathogenic)) +
    geom_bar(stat='identity', col='black', size=0.05) +
    scale_fill_discrete(name="Presence of pathogenic variants") +
    labs(title="Percentage of pathogenic variants per familial history type",
         x = "Type of family history") +
    geom_label(aes(y=label_y, label=paste("n =", n)),
               vjust = 1.1, colour = "black") +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    theme_classic()

print(samples_pathogenic_plot)
ggsave('Plots/type_history_percentage_BRIDGES.pdf', 
       samples_pathogenic_plot, width=30, height=20, units='cm')
ggsave('Plots/type_history_percentage_BRIDGES.png', 
       samples_pathogenic_plot, width=30, height=20, units='cm') 
    
    
