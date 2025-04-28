library(reshape2)
library(stringr)
library(ggpubr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(poibin)
library(tibble)


d_all <- read.csv("./Data/2_basic/Ancillary_species_Site_Dermotype_minimal_3.csv")
mdata.dermotype <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
load("./Data/Taxon_Decontam_AbundFilt_Renorm.RData")
mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")

mat_auxillary <- mat %>% select(all_of(d_all$Taxa_og)) 

merged_data <- merge(mdata.dermotype[, c("LibraryID", "SubjectID", "Site")], 
                     mat_auxillary, 
                     by.x = "LibraryID", by.y = 0) %>% 
  mutate(SubjectID = str_pad(SubjectID , 5, pad = "0"))
  
# Step 1: Convert abundance data to presence/absence
df_presence <- merged_data %>%
  pivot_longer(cols = starts_with("s__"), names_to = "Taxon", values_to = "Abundance") %>%
  mutate(Present = if_else(Abundance > 0, 1, 0, missing = 0)) %>%  # Binary presence
  select(-Abundance) %>%  # Remove raw abundance
  pivot_wider(names_from = "Taxon", values_from = "Present")  # Reshape to wide format

# Step 2: Calculate prevalence by site (baseline probability for each taxon-site pair)
prevalence_by_site <- df_presence %>%
  pivot_longer(cols = starts_with("s__"), names_to = "Taxon", values_to = "Present") %>%
  group_by(Site, Taxon) %>%
  summarise(Prevalence = mean(Present, na.rm = TRUE), .groups = "drop")

# Step 3: Compute expected probabilities for each taxon (list of site-specific prevalences)
expected_probabilities <- prevalence_by_site %>%
  group_by(Taxon) %>%
  summarise(p_vec = list(Prevalence), .groups = "drop") %>%
  { setNames(.$p_vec, .$Taxon) }

# Step 4: Count positive sites for each subject-taxon pair
subject_counts <- df_presence %>%
  pivot_longer(cols = starts_with("s__"), names_to = "Taxon", values_to = "Present") %>%
  group_by(SubjectID, Taxon) %>%
  summarise(Count = sum(Present), .groups = "drop")  # Total positive sites per subject-taxon

# Step 5: Compute p-values for each subject-taxon pair
subject_enrichment <- subject_counts %>%
  left_join(expected_probabilities %>% enframe(name = "Taxon", value = "p_vec"), by = "Taxon") %>%
  rowwise() %>%
  mutate(
    PValue = if_else(
      all(!is.na(p_vec)) && length(p_vec) > 0,  # Check for valid p_vec
      1 - ppoibin(Count - 1, p = p_vec),  # Compute p-value
      NA_real_  # Return NA if p_vec is invalid
    )
  ) %>%
  ungroup()

# Step 6: Apply Benjamini-Hochberg correction globally
subject_enrichment <- subject_enrichment %>%
  mutate(AdjPValue = p.adjust(PValue, method = "BH"))

# Step 7: Filter significant subjects after correction
significant_subjects <- subject_enrichment %>%
  filter(AdjPValue < 0.05) %>% 
  mutate(Enriched = "T")
n_distinct(significant_subjects$SubjectID)
n_distinct(significant_subjects$Taxon)

significant_subjects_multi_ancillary_taxa <- 
  significant_subjects %>% 
  group_by(SubjectID) %>% 
  summarize(n = n())
  
# Step 8: Visualize results (optional)
library(ggplot2)
library(dplyr)
library(stringr)

# Step 1: Identify taxa that have at least one significant subject
taxa_sig <- significant_subjects %>% 
  pull(Taxon) %>%
  unique()

# Step 2: Create Taxon_clean and reorder based on significant subjects
# Clean the Taxon names (remove s__ and replace _ with space)
subject_counts_flagged <- subject_counts%>% 
  left_join(significant_subjects) %>% 
  mutate(
    Taxon_clean = gsub("s__", "", Taxon),
    Taxon_clean = gsub("_", " ", Taxon_clean)
  )

# Step 3: Count the number of significant subjects per taxon (after cleaning)
taxa_order <- subject_counts_flagged %>%
  filter(Taxon %in% taxa_sig) %>%
  mutate(
    Taxon_clean = gsub("s__", "", Taxon),
    Taxon_clean = gsub("_", " ", Taxon_clean)
  ) %>%
  group_by(Taxon_clean) %>%
  summarise(n_signif = n_distinct(SubjectID)) %>%
  arrange(desc(n_signif)) %>%
  pull(Taxon_clean)

# Step 4: Reorder the plot_data based on the new Taxon_clean levels
plot_data <- subject_counts_flagged %>%
  filter(Taxon %in% taxa_sig) %>%
  mutate(
    Taxon_clean = factor(Taxon_clean, levels = taxa_order)  # Reorder based on the number of significant subjects
  )

# Step 5: Create the faceted violin plot with wrapped facet labels
g_enrichment <- ggplot(plot_data, aes(x = "", y = Count)) +
  # Full violin for the overall distribution (light blue)
  geom_violin(fill = "lightblue", color = NA) +
  geom_jitter(
    data = filter(plot_data, Enriched == "T"),
    aes(y = Count),
    color = "red", width = 0.15, size = 1
  ) +
  # Facet by the cleaned Taxon name, with 2 rows
  facet_wrap(~ Taxon_clean, 
             scales = "free_y", 
             nrow = 2, 
             labeller = labeller(Taxon_clean = function(x) str_wrap(x, width = 10))) +
  labs(
    title = "",
    x = NULL,
    y = "Number of skin site", 
  ) +
  ggpubr::theme_pubr(base_size = 7) +
  theme(
    # Italic facet labels
    strip.text = element_text(face = "italic"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(breaks = c(0,3,6,9)) # Set Y axis to integers
g_enrichment



#=========== cross-site correlation ===================#
within_subject_prev <- merged_data %>%
  dplyr::group_by(SubjectID_Site = paste(SubjectID,Site, sep = "_"), .drop = F) %>%
  dplyr::summarise(across(contains(c("s__")), \(x) sum(x, na.rm = TRUE))) %>% 
  mutate(SubjectID = str_extract(SubjectID_Site, ".*(?=_)")) %>%
  dplyr::group_by(SubjectID) %>%
  dplyr::summarise(across(contains(c("s__")), ~ sum(. > 0, na.rm = TRUE))) 
melted_prev <- reshape2::melt(within_subject_prev, id.vars = "SubjectID") %>%
  filter(value > 0) %>%
  dplyr::rename(prev = value)

# Intra-individual: the non-zero median RA of these virus
within_subject_prev <- merged_data %>%
  dplyr::group_by(SubjectID_Site = paste(SubjectID,Site, sep = "_"), .drop = F) %>%
  #select(-LibraryID) %>%
  dplyr::summarise(across(contains(c("s__")), \(x) sum(x, na.rm = TRUE))) %>% 
  mutate(SubjectID = str_extract(SubjectID_Site, ".*(?=_)")) %>%
  dplyr::group_by(SubjectID) %>%
  dplyr::summarise(across(contains(c("s__")), ~ sum(. > 0, na.rm = TRUE))) 
melted_prev <- reshape2::melt(within_subject_prev, id.vars = "SubjectID") %>%
  filter(value > 0) %>%
  dplyr::rename(prev = value)
non_zero_median <- function(x) {
  median(x[x != 0], na.rm = T)
}
within_subject_RA <- merged_data %>%
  group_by(SubjectID) %>%
  dplyr::summarize(across(all_of(unique(melted_prev$variable)), non_zero_median, .names = "{.col}"))
melted_RA <- reshape2::melt(within_subject_RA, id.vars = "SubjectID") %>% 
  filter(!is.na(value))%>%
  dplyr::rename(RA = value)

# the number of individuals harboring these virus at >= 1 skin sites 
count_positive_obs <- within_subject_prev %>%
  dplyr::summarize(across(-c(SubjectID),  ~ sum(. > 0, na.rm = TRUE)))

#merge and subset to the interested taxa group  
melted_data <- melted_prev %>%
  left_join(melted_RA, join_by("SubjectID", "variable"))%>%
  left_join(., d_all[, c("Taxa", "Taxa_og", "Taxa_group")], join_by( variable == Taxa_og),
            relationship = "many-to-many") 

df = as.data.frame(t(count_positive_obs)) %>% 
  rownames_to_column(var = "variable") %>% 
  filter(variable %in% as.character(melted_data$variable) == T) %>%
  right_join(melted_data, join_by("variable")) %>%
  arrange(-V1, -prev, -RA) %>%
  distinct() 

# visualization #
p_fisher <- list()
df_gt10_carrier <- df %>% arrange(Taxa_group) %>% filter(V1 >10) 
list_microbe <- unique(df_gt10_carrier$variable)
#list_microbe <- c(unique(bloom_example$variable))
for (i in c(1:length(list_microbe))){
  bloom_taxa = list_microbe[i]
  tmp <- mat %>% dplyr::select(all_of(bloom_taxa)) %>%
    dplyr::rename("bloom_taxa" = bloom_taxa) %>%
    rownames_to_column(var = "LibraryID") %>%
    left_join(mdata.dermotype[, c("LibraryID", "SubjectID", "Site")], join_by("LibraryID")) %>%
    group_by(SubjectID, Site) %>%
    dplyr::summarise(bloom_taxa = mean(bloom_taxa)) %>%
    mutate(bloom_taxa = ifelse(bloom_taxa > 0, 1, 0)) %>%
    pivot_wider(names_from  = Site,
                values_from = bloom_taxa) %>%
    column_to_rownames(var = "SubjectID") %>%
    filter(rowSums(.) > 0)
  
  pairwise_fisher <- function(col1, col2) {
    tbl <- table(col1, col2)
    
    # Check if the table has at least 2 rows and 2 columns
    if (nrow(tbl) < 2 || ncol(tbl) < 2) {
      return(NULL)  # Skip the test
    }
    test <- fisher.test(tbl)
    return(data.frame(odds_ratio = test$estimate,
                      p_value_single = test$p.value))
  }
  
  # Calculate pairwise fisher test for every combination of columns
  results <- combn(names(tmp), 2, function(cols) {
    res <- pairwise_fisher(tmp[[cols[1]]], tmp[[cols[2]]])
    
    if (!is.null(res)) {  # Only add results if the test was performed
      res$col1 <- cols[1]
      res$col2 <- cols[2]
      return(res)
    } else {
      return(NULL)  # Skip if res is NULL
    }
  }, simplify = FALSE)
  
  # Remove NULL entries from results
  results <- Filter(Negate(is.null), results)
  
  # Combine results into a single data frame
  df_results <- bind_rows(results)
  df_results <- df_results %>% 
    mutate(p_value = p.adjust(p_value_single, method = "BH")) %>% 
    filter(p_value < 0.05 & odds_ratio != Inf) %>%
    arrange(p_value, -odds_ratio) %>% 
    select(-p_value_single)
  
  df_results_2 <- df_results |>
    dplyr::rename(c("col2" = "col1", "col1" = "col2")) |>
    dplyr::select(col1, col2, p_value, odds_ratio)
  
  results_sorted <- rbind(df_results, df_results_2)
  unique_sites <- unique(c(results_sorted$col1, results_sorted$col2))
  all_combinations <- expand.grid(col1 = unique_sites, col2 = unique_sites)
  full_data <- merge(all_combinations, results_sorted, by = c("col1", "col2"), all.x = TRUE)
  # Order factors for plotting
  desired_order <- c(  
    "Ax", "Ac","Ll", "Vf", "Ub", "Fo","Ch", "Ps", "Sc")
  full_data$col1 <- factor(full_data$col1, levels = desired_order)
  full_data$col2 <- factor(full_data$col2, levels = desired_order)
  full_data2 <- full_data[which(as.numeric(full_data$col1) < as.numeric(full_data$col2)), ]
  
  if (nrow(full_data2) > 0) {
    p = ggplot(full_data2, aes(x = col1, y = col2)) +
      geom_tile(aes(fill = log10(p_value)), color = "white") +
      #geom_text(aes(label = round(odds_ratio,1)), na.rm = TRUE) +
      scale_fill_gradient(low = "red", high = "yellow", 
                          na.value = "white", 
                          name = "Fisher's\ntest\nlog10(p-value)",
                          limits = c(-5, -1)) +
      scale_x_discrete(limits = levels(full_data2$col1)) + 
      scale_y_discrete(limits = levels(full_data2$col2)) +
      theme_pubr(base_size = 7, border = T, x.text.angle = 0) + 
      theme(
        axis.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.4),
        legend.box.background = element_rect(color = "black")
      ) +
      annotate(geom = "text", 
               label = str_replace_all(bloom_taxa, "s__|_", " "), 
               fontface = "italic", 
               x = -Inf, y = Inf, 
               hjust = 0, vjust = 2, 
               size = 3
      )
    p 
    p_fisher[[bloom_taxa]] <- p 
  }
  print(i)
  # ggsave(filename = paste("../Output/4_common/", bloom_taxa, "Site-Site_correlation.png"),
  #        plot = p, height = 7, width = 7)
}
plot_grid(p_fisher[[1]], p_fisher[[2]], p_fisher[[3]],  p_fisher[[5]],p_fisher[[6]], nrow  = 1)


#=========== association with host factor ===================#
results_df_1 <- read.csv("~/ASMP_Dermotype/Data/2_basic/Ancillary_taxa_carriage_host_factor.csv")
results_df_2 <- results_df_1 %>%
  mutate(p_asterisk = case_when(
    PValue <= 0.01 ~ "**",
    PValue <= 0.05 ~ "*",
    PValue <= 0.1 ~ "#",
    TRUE ~ ""
  )) %>% 
  mutate(
    Factor_Category = gsub("Skin|_skin", "", Column_Category_2), 
    Factor_Category = factor(Factor_Category, 
                             levels = c("Demographics", "Physio",
                                        "Behavior", "Discomfort",
                                        "Disease", "GeneralHealth")),
    Data_type = factor(Data_type, levels = c("Character", "Numeric"))) 


results_df_3 <- results_df_2 %>%
  filter(Factor_Category != "GeneralHealth") %>%
  filter(Formatted_Predictor != "Gender (Female)") %>% 
  filter(!str_detect(Formatted_Predictor, "not described above")) %>% 
  mutate(Formatted_Predictor = ifelse(Formatted_Predictor == "PH forearm", "pH forearm", Formatted_Predictor)) %>% 
  filter(PValue < 0.1) %>%
  arrange(Factor_Category, 
          Column_Category, 
          factor(Column_Subcategory, levels = c(setdiff(unique(Column_Subcategory), "General"), "General")), 
          desc(Data_type))

results_df_3$Formatted_Predictor<- factor(results_df_3$Formatted_Predictor, 
                                          levels = unique(results_df_3$Formatted_Predictor))


results_df_3 <- results_df_3 %>% 
  mutate(Taxa = str_replace_all(Taxa, "s__|_", " "))

p_carrier_phenotype <- results_df_3 %>%
  ggplot(aes(y = Taxa, x = Formatted_Predictor)) +
  geom_tile(aes(fill = EffectSize), color = "white") +
  geom_text(aes(label = p_asterisk), vjust = -0.3, hjust = 0.5, color = "black", size = 2) +
  geom_text(aes(label = ifelse(abs(EffectSize) < 0.2 & EffectSize < 0 , "-", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "grey30") +
  geom_text(aes(label = ifelse(abs(EffectSize) < 0.2 & EffectSize > 0 , "+", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "grey30") +
  geom_text(aes(label = ifelse(abs(EffectSize) >= 0.2 & EffectSize < 0 , "-", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "white") +
  geom_text(aes(label = ifelse(abs(EffectSize) >= 0.2 & EffectSize > 0 , "+", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "white") +
  scale_fill_gradient2(low = "#013220", mid = "grey90", high = "maroon", midpoint = 0, guide = "colorbar") +
  labs(title = "", x = "", y = "", fill = "Beta-coefficient") +
  facet_grid(cols = vars(Factor_Category),
             rows = vars(Taxa_category),
             space = "free", 
             scale = "free", 
             labeller = label_wrap_gen(width = 5)) + 
  theme_pubr(border = T, base_size = 7, x.text.angle = 40) +
  theme(axis.text.y = element_text(face = "italic"), 
        legend.position = "right", 
        strip.text.y = element_text(angle = 0),
        strip.text = element_text(size = rel(1.2)))
p_carrier_phenotype
