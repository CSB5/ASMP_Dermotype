library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(mltools)
library(data.table)
library(ggpubr)

# Load data
mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.dermotype <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
load("./Data/Taxon_Decontam_AbundFilt_Renorm.RData")

# Function to run the analysis with flexible Ax grouping
run_analysis <- function(merge_staph = FALSE) {
  
  df <- mdata.dermotype |> 
    mutate(SiteID = Site,
           SubjectID = paste(SubjectID, LR, sep = "")) |> 
    dplyr::select(Dermotype_size, SiteID, SubjectID) |> 
    pivot_wider(names_from = SiteID, values_from = Dermotype_size) |>
    dplyr::select(-contains(c('Ll', 'Ch'))) |>
    column_to_rownames("SubjectID")
  
  # Adjust Ax grouping if merge_staph = TRUE
  df <- df |> mutate(Ax = case_when(
    merge_staph & Ax %in% c("Ax-3", "Ax-5") ~ "Staphylococcus-\ndominated\n(Ax-3, Ax-5)",
    TRUE ~ Ax
  ))
  
  df <- as.data.frame(lapply(df, factor))  # Ensure factors are updated
  
  # One-hot encoding
  df_one_hot <- mltools::one_hot(data.table::data.table(df[, c("Ax", "Ac")])) 
  predictor_vars <- df_one_hot %>% select(contains("Ac")) %>% names()
  response_vars <- df_one_hot %>% select(contains("Ax")) %>% names()
  
  # Store results
  results_df <- data.frame()
  
  # Logistic regression loop
  for (response_var in response_vars) {
    for (predictor_var in predictor_vars) {
      model <- glm(as.formula(paste0("`", response_var, "` ~ `", predictor_var, "`")), 
                   data = df_one_hot, family = binomial(link = "logit"))
      coefficients <- coef(summary(model))[2, 1]
      p_value <- coef(summary(model))[2, 4]
      
      results_df <- rbind(results_df, data.frame(
        Predictor = predictor_var, 
        Response = response_var, 
        Beta_coefficient = coefficients,
        P_Value = p_value
      ))
    }
  }
  
  # Post-processing of results
  results_df <- results_df %>% 
    separate(Response, c("Response_Site", "Response"), "_") %>% 
    separate(Predictor, c("Predictor_Site", "Predictor"), "_") %>% 
    mutate(p_asterisk = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  # Generate heatmap plot
  p <- results_df %>% 
    mutate(Response = case_when(
      Response %in% c("Ax-3", "Ax-5") ~ paste0( Response,  "\nStaphylococcus-\ndominated"),
      TRUE ~ Response
    ),
    Predictor = case_when(
      Predictor %in% c("Ac-2") ~ paste0(Predictor,  "\nStaphylococcus-\ndominated"),
      TRUE ~ Predictor
    )) %>%
    ggplot(aes(x = Predictor, y = Response, fill = Beta_coefficient)) +
    geom_tile(color = "white") +
    facet_grid(Response_Site ~ Predictor_Site, space = "free", scale = "free") +
    scale_fill_gradient2(low = "blue", mid = "white", midpoint = 0, high = "red") +
    geom_text(aes(label = p_asterisk), color = "black", vjust = 0.5, hjust = 0.5) +
    labs(x = "", y = "") +
    ggpubr::theme_pubr(base_size = 7)
  
  return(list(results = results_df, plot = p))
}

# Run both analyses
analysis_separate <- run_analysis(merge_staph = FALSE)  # Ax-3 and Ax-5 separate
analysis_merged <- run_analysis(merge_staph = TRUE)     # Ax-3 and Ax-5 grouped

# Access results and plots
analysis_separate$plot  # Plot for separate Ax-3 and Ax-5
analysis_merged$plot    # Plot for merged Ax-3 and Ax-5

#===========check S.hominis abundance correlation ========#
m_Staph <- mat %>% 
  select(contains("Staphy") & !contains("phage") ) %>% 
  #select(contains("Cuti")& contains("acnes") & !contains("phage") ) %>% 
  mutate(Staphylococcus =  rowSums(across(everything()))) %>% 
  rownames_to_column(var = "LibraryID") %>% 
  left_join(select(mdata.dermotype, LibraryID, SubjectID, Site, LR)) %>% 
  filter(Site %in% c("Ac", "Ax")) %>% 
  mutate(SubjectID_LR = paste(SubjectID, LR, sep = "")) 


plot_bug_abund_cor <- function(bug){
  m1 <- m_Staph %>% 
    select(any_of(bug), SubjectID_LR, Site) %>% 
    pivot_wider(names_from = Site, values_from = .data[[bug]])
  
  # Remove rows with missing values
  m1_clean <- na.omit(m1)
  
  m1_clean_non_zero <- m1 #%>% filter(Ax > 0.001 & Ac > 0.001)
  # Calculate the Spearman correlation coefficient
  cor_test <- cor.test(log10(m1_clean_non_zero$Ax), log10(m1_clean_non_zero$Ac), method = "spearman", exact = F)
  spearman_pval <- cor_test$p.value
  spearman_cor <- cor_test$estimate
  # Create the plot with geom_smooth and annotate the correlation coefficient
  p <- ggplot(m1_clean, aes(x = Ac, y = Ax)) +
    geom_jitter(color = "grey30") +
    geom_smooth(method = "lm", se = T, data = m1) + #[m1$Ax > 0.001 & m1$Ac > 0.001, ]) +
    scale_x_log10(labels = scales::percent) +
    scale_y_log10(labels = scales::percent) +
    # geom_vline(xintercept = 0.001, color = "red", linetype = 2) +
    # geom_hline(yintercept = 0.001, color = "red", linetype = 2) +
    ggtitle(paste("Spearman r =", round(spearman_cor, 3),
                  "\np-value =", signif(spearman_pval, 2))) +
    xlab(paste0("[Ac] ", str_replace_all(bug, "s__|_", " "))) +
    ylab(paste0("[Ax] ", str_replace_all(bug, "s__|_", " "))) +
    ggpubr::theme_pubr(base_size = 7)
  return(p)
}
p_sh <- plot_bug_abund_cor(bug = "s__Staphylococcus_hominis")
p_epi <- plot_bug_abund_cor(bug = "s__Staphylococcus_epidermidis")
p_all_staph <- plot_bug_abund_cor(bug = "Staphylococcus")

library(cowplot)
p_set <- plot_grid(p_all_staph, p_sh, p_epi, ncol = 3)
p_set



