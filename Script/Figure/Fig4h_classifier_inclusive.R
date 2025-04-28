# Load necessary libraries
library(dplyr)
library(tidyr)
library(openxlsx)
library(purrr)
library(stringr)

# List all subfolders in the "Inclusive" folder
subfolders <- list.dirs("./Data/4_ml/Inclusive_2025-0226", full.names = TRUE, recursive = FALSE)

process_subfolder_performance <- function(subfolder_path) {
  # Read the files from the subfolder
  performance <- read.csv(file.path(subfolder_path, "scores.csv"), check.names = F)
  performance_sub <- performance %>%  mutate(subfolder = basename(subfolder_path))
  return(performance_sub) 
}


# Function to process each subfolder
process_subfolder <- function(subfolder_path) {
  
  importance <- read.csv(file.path(subfolder_path, "importance.csv"), check.names = F)
  # Process the importance data
  importance_sub <- importance %>%
    summarise(across(everything(), ~ mean(abs(.) , na.rm = TRUE))) %>%  # Calculate column-wise averages
    pivot_longer(everything(), names_to = "feature", values_to = "importance") %>% 
    mutate(category = case_when(
      grepl("Age|Ethnicity|Gender", feature, ignore.case = F) ~ "demographic",
      grepl("Ac_|Ax_|Vf_|Fo_|Ps_|Sc_|Ub_", feature, ignore.case = F) ~ "dermotype",
      TRUE ~ "pathway"
    )) %>% 
    arrange(desc(abs(importance))) %>% 
    group_by(category) %>% 
    mutate(class_importance = mean(abs(importance))) %>%
    mutate(subfolder = basename(subfolder_path))  # Add subfolder name to each row

  return(importance_sub)
}

# Initialize the workbook
wb <- createWorkbook()

# Loop through subfolders to process them and add each to the workbook
for (subfolder_path in subfolders) {
  # Process the current subfolder's data
  importance_sub <- process_subfolder(subfolder_path)
  
  # Extract the subfolder name to use as the sheet name
  subfolder_name <- basename(subfolder_path)
  
  # Add a new sheet for each subfolder with its name
  addWorksheet(wb, subfolder_name)
  
  # Write the subfolder's processed data into the corresponding sheet
  writeData(wb, subfolder_name, importance_sub)
}

# Save the workbook with all sheets
#saveWorkbook(wb, "~/Desktop/ASMP/Output/3_ML/0_Phenotype_asso/Supplementary_File_inclusive_classifier.xlsx", overwrite = TRUE)

#---visualize performance -----#
all_performance_data <- list()
for (subfolder_path in subfolders) {
  performance_sub <- process_subfolder_performance(subfolder_path)
  all_performance_data[[subfolder_path]] <- performance_sub
}
combined_performance_data <- bind_rows(all_performance_data)
head(combined_performance_data)

ggplot(combined_performance_data, aes(x = subfolder, y = test_roc_auc_scorer)) +
  geom_boxplot() +
  ggpubr::theme_pubr()

combined_performance_data %>% 
  group_by(subfolder) %>% 
  summarise(across(where(is.numeric), ~ mean(. , na.rm = TRUE))) %>% 
  arrange(desc(test_roc_auc_scorer), desc(test_accuracy))
  

#----visualization importance----#
all_importance_data <- list()
for (subfolder_path in subfolders) {
  importance_sub <- process_subfolder(subfolder_path)
  all_importance_data[[subfolder_path]] <- importance_sub
}
combined_importance_data <- bind_rows(all_importance_data)
head(combined_importance_data)

library(dplyr)
library(ggsignif)
library(rstatix)

# Ensure the abs_importance column exists
combined_importance_data <- combined_importance_data %>%
  mutate(abs_importance = abs(importance)) %>% 
  mutate(category = factor(category, levels = c("dermotype",
                                                "pathway",
                                                
                                                "demographic"))) %>% 
  separate(feature, into = c("group1", "group2"), sep = "_") %>% 
  group_by(group1, subfolder, category) %>% 
  dplyr::summarize(abs_importance = mean(abs_importance)) %>% 
  ungroup() %>% 
  group_by(subfolder) %>% 
  # Calculate min and max importance per subfolder
  mutate(min_importance = min(abs_importance, na.rm = TRUE),
         max_importance = max(abs_importance, na.rm = TRUE)) %>% 
  # Apply min-max normalization
  mutate(scaled_importance = (abs_importance - min_importance) / (max_importance - min_importance)) %>% 
  ungroup() %>%
  select(-min_importance, -max_importance) %>% 
  # Sort by site, category, and scaled importance for better interpretation
  arrange(subfolder, category, desc(scaled_importance))

pairwise_tests <- combined_importance_data %>%
  group_by(subfolder) %>%
  pairwise_wilcox_test(
    formula = scaled_importance ~ category,
    alternative = "greater", 
    p.adjust.method = "BH"
  ) %>%
  filter(p.adj < 0.05) %>%  # Only keep significant comparisons
  rowwise() %>%
  mutate(
    # Calculate the y position as the max abs_importance among the two compared categories in the same subfolder plus an offset
    y.position = max(combined_importance_data$scaled_importance[
      combined_importance_data$subfolder == subfolder & 
        combined_importance_data$category %in% c(group1, group2)
    ]) + 0.2
  ) %>%
  ungroup()

# Create the boxplot with the Wilcoxon test annotations
p_inclusive <- ggplot( combined_importance_data,
       aes(x = category, y = scaled_importance, fill = category)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  geom_boxplot() +
  facet_grid(~subfolder, space = "free") +
  geom_signif(
    data = pairwise_tests,
    inherit.aes = FALSE,  # Disable inheriting the main plot's aesthetics
    aes(xmin = group1, xmax = group2, y_position = y.position, annotations = p.adj.signif),
    manual = TRUE
  ) +
  ggpubr::theme_pubr(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.direction = "vertical", 
        legend.position = "right") +
  labs(y = "Normalized feature importance\n(Ridge classifier)", x = "")
p_inclusive

