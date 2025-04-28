library(dplyr)
library(tidyverse)
library(readr)
load("/Users/chengchenli/Desktop/ASMP/Data/Metadata_dermotype_BatchABCD_with_StabilityScore.RData")
# source("/Users/chengchenli/ASMP-Network/script/0_Clean_pwy_abund.R")
load("/Users/chengchenli/Desktop/ASMP/Data/Cleaned_Pathway_Abund.RData")

#=========================== DA analysis =====================================#
library(Maaslin2)
dir.create("/Users/chengchenli/Desktop/ASMP/Output/3_ML/MaAsLin2_InterDermotype_pwy")

unique_sites = unique(mdata.dermotype.conf$Site)
for (i in 1: length(unique_sites)){
  PickedSite = unique_sites[i]
  if (PickedSite %in% c("Ll", "Ch")) next # skip because no dermotype discovered
  
  dir.create(paste("/Users/chengchenli/Desktop/ASMP/Output/3_ML/MaAsLin2_InterDermotype_pwy", 
                   PickedSite, sep = "/"))
  input_metadata = mdata.dermotype.conf |>
    filter(Site == PickedSite & LibraryID %in% rownames(pwy_clean)) |>
    data.frame() |>
    # left_join(select(tmp_unmap, LibraryID, unmap), join_by(LibraryID)) |>
    column_to_rownames(var = "LibraryID")
  input_data = pwy_clean[rownames(pwy_clean) %in% rownames(input_metadata), ]
    
  # ---- Iterate over each unique level of Dermotype_size for OvA logistic regression --- # 
  for (dermotype in unique(input_metadata$Dermotype_size)) {
    #---- Create a binary outcome: 1 for current dermotype, 0 for all others ---#
    input_metadata$binary_outcome <- ifelse(input_metadata$Dermotype_size == dermotype, dermotype, "rest")
    input_metadata$binary_outcome <- factor(input_metadata$binary_outcome, levels = c("rest", dermotype))
    fit_data <- Maaslin2(
      input_data,
      input_metadata,
      output = paste("/Users/chengchenli/Desktop/ASMP/Output/3_ML/MaAsLin2_InterDermotype_pwy", PickedSite, as.character(dermotype), sep = "/"),
      fixed_effects = c("binary_outcome", "Ethnicity", "Age", "Gender"),
      random_effects = c("SubjectID"),
      min_abundance = 0.001,
      min_prevalence = 0.1,
      max_significance = 0.05,
      plot_scatter = F, 
      reference = paste("binary_outcome,rest;Ethnicity,East Asian", sep = ""))
  }
} 
  
# #### merge all of significant results from MaAsLin2 #####
path <- "/Users/chengchenli/Desktop/ASMP/Output/3_ML/MaAsLin2_InterDermotype_pwy"
file_list <- list.files(path, pattern = "significant_results\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
data_list <- lapply(file_list, read.delim)
combined_data <- do.call(rbind, data_list)

DA_pwy <- combined_data |>
  filter(metadata %in% c("binary_outcome")) |>
  mutate(feature = str_replace_all(feature, "\\.", "-")) |>
  mutate(feature = str_replace_all(feature, "X(\\d+)", "\\1")) |>
  mutate(Site = str_replace(value, "\\-(\\d+)", ""))
n_distinct(DA_pwy$feature)

#----- number of unique differential pathways in each site-----#
tmp1 <- DA_pwy %>%
  group_by(Site) %>%
  summarise(n_distinct_features = n_distinct(feature))
tmp1
#----- number of unique pathways in each site-----#
all_pwy <- pwy_clean %>% 
  rownames_to_column("LibraryID") %>% 
  right_join(select(mdata.dermotype.conf,LibraryID, Site), join_by(LibraryID)) %>% 
  select(-LibraryID)

library(dplyr)
library(purrr)

all_pwy %>%
  group_split(Site) %>%
  map_df(~{
    df <- .x
    site_name <- unique(df$Site)
    df_paths <- df %>% select(where(is.numeric))
    df_paths_nonzero <- df_paths[, colSums(df_paths != 0, na.rm = TRUE) > 0]
    tibble(Site = site_name, n_unique_pathways = ncol(df_paths_nonzero))
  }) %>% 
  right_join(tmp1) %>% 
  mutate(perc_diff = n_distinct_features*100/n_unique_pathways) %>% 
  summary(perc_diff)


ggplot(DA_pwy, aes(x = value, y = abs(coef))) +
  geom_boxplot() +
  facet_grid(~Site, space = "free", scales = "free")
write.csv(DA_pwy, file = "/Users/chengchenli/Desktop/ASMP/Output/3_ML/DiffAbundPWY_Dermotype.csv",
          row.names = F)

# DA_pwy_unmap <- combined_data |>
#   filter(metadata %in% c("unmap")) |>
#   mutate(feature = gsub("\\.", "-", feature)) |>
#   mutate(feature = case_when( feature =="X1CMET2-PWY" ~ "1CMET2-PWY", 
#                               T ~ feature)) 
# n_distinct(DA_pwy_unmap$feature)
  
