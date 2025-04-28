library(dplyr)
library(tidyverse)
library(readr)

mdata.dermotype.conf <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
pwy_clean <- read.csv("./Data/3_pwy/Cleaned_Pathway_Abund.csv", check.names = F) 
names(pwy_clean)[1] <- "LibraryID"
pwy_clean <- pwy_clean %>% 
  column_to_rownames(var = "LibraryID")
pwy_ctr_sub <- read.csv("./Data/3_pwy/pwy_ctr_sub.csv", check.names = F, row.names = NULL ) %>% 
  select(-1)
pwy_asso_phenotype <-read.csv("./Data/3_pwy/pwy_asso_phenotype.csv")
source("./Script/Analysis/UtilityFunction_prep_data_for_pwy_species_contribution.R")

#==== pick a pathway =====
i = "Ac"

j = "HISDEG-PWY" 
j = "LACTOSECAT-PWY"

output <- pwy_contri_by_group(i, # skin site of interest
                              j	, # pick a pathway (shortened name, not the full name)
                              mdata.dermotype.conf, #mdata files 
                              pwy_ctr_sub,
                              pwy_clean, 
                              phenotype = mdata.dermotype.conf,
                              pwy_asso_phenotype,
                              n_taxa = 5
                              ) 


library(cowplot)
# Create the combined plot
combined_plot <- plot_grid(
  output[["vis_relab"]],  # First plot on the left
  output[["vis_average_species_contribution"]],  # Second plot
  ncol = 2,  # Arrange into 2 columns
  rel_widths = c(0.5, 1), 
  align = "h"# Adjust relative widths if needed
)

# Add title while keeping plots aligned
final_plot <- plot_grid(
  ggdraw() + draw_label(output[["pathway_fullname"]], size = 12, hjust = 0.5),  
  combined_plot,  
  ncol = 1,  # Stack title and plot vertically
  rel_heights = c(0.3, 1)  # Adjust title height relative to plots
)
final_plot

  

  