## -----------------------------------------------------------------------------
##
## Author: Arthur Boffelli Castro
##
## Date created: Fri Jan  7 11:45:03 2022
##
## GitHub: https://github.com/aboffelli/
##
## Description:
##
##  
## -----------------------------------------------------------------------------
## 
## Notes:
##  
##  
## ----------------------------------------------------------------------------- 

library(ggplot2)


setwd("~/Box/Notes/Tables/ClinvarTables")
samples_table <- read.table("non_patho_other_variants.txt", sep='\t')
samples_table_unique <- unique(samples_table[,-2])
samples_table_unique <- cbind(samples_table_unique, V4='yes')

plot <- ggplot(data=samples_table_unique, aes(x=V3, y=V1)) +
    geom_tile(fill='darkblue') +
    theme_bw()
    
    
ggsave('test_heat.pdf', plot, width=40, height = 60, units = 'cm')

abundance <- sort(table(samples_table[,2]), decreasing=T)
