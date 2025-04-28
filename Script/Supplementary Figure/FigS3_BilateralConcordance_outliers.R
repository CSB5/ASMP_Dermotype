library(ggpubr)
library(vegan)
library(reshape2)
library(tidyverse)
library(tibble)
library(dplyr)
library(svglite)
library(tidyr)

# load("/Users/chengchenli/Desktop/ASMP/Data/Tax_BatchABCD_AbundFilt_PrevFilt.RData")
mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.ex <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
load("./Data/2_basic/Pearson.RData")

# Calculate the upper and lower fences for each group
boxplot_stats <- tmp2 %>%
  mutate(Var1_SubjectID = str_pad(Var1_SubjectID , 5, pad = "0"),
    Var2_SubjectID = str_pad(Var2_SubjectID , 5, pad = "0")) %>% 
  filter(type == "intra-individual") %>%
  group_by(Var1_Site) %>% 
  dplyr::summarize(lower_fence = quantile(value, 0.25) - 1.5 * IQR(value)) %>%
  #dplyr::summarize(lower_fence = 0.7) %>%
  data.frame()

# Add outlier information to the original data frame
tmp3 <- tmp2 %>%
  filter(type == "intra-individual") %>%
  distinct(Var1_Site, Var1_SubjectID, value) %>%
  dplyr::left_join(boxplot_stats, by = c("Var1_Site")) %>%
  mutate(outlier = case_when(value < lower_fence ~ "Discordant",
                             T ~ "Concordant")) %>% 
  mutate(Var1_SubjectID = str_pad(Var1_SubjectID , 5, pad = "0"))

#============= produce a summary table ==================#
#install.packages("pacman")
pacman::p_load(
  rio,          # File import
  here,         # File locator
  skimr,        # get overview of data
  tidyverse,    # data management + ggplot2 graphics 
  gtsummary,    # summary statistics and tests
  rstatix,      # summary statistics and statistical tests
  janitor,      # adding totals and percents to tables
  scales,       # easily convert proportions to percents  
  flextable     # converting tables to pretty images
)

g_table <- tmp3 %>%
  tabyl(Var1_Site, outlier) %>%
  adorn_percentages(denominator = "row") %>%
  adorn_pct_formatting() %>%
  adorn_ns(position = "front") %>%
  rename(Site = Var1_Site) %>%  # Rename first column
  flextable() %>%
  autofit() %>%
  align(align = "center", part = "body") %>%
  align(align = "center", part = "header") %>%
  footnote(
    i = 1, j = 3,  # Target the third column header
    value = as_paragraph("Discordant: Pearson correlation coefficient < 25th percentile - 1.5*IQR"),
    ref_symbols = "†",
    part = "header"
  ) %>%
  gen_grob(fit = "auto", just = "center")

plot(g_table)


tmp4 <- mdata.ex %>% 
  select(-c(LibraryID, Dermotype_size, Stability)) %>% 
  mutate(SubjectID = str_pad(SubjectID , 5, pad = "0")) %>% 
  distinct() %>% 
  dplyr::right_join(tmp3, by = join_by(SubjectID == Var1_SubjectID, 
                                       Site == Var1_Site)) %>%
  mutate(Site = factor(Site, levels = c("Ax", "Ll", "Vf", "Ac", "Ub",  
                                        "Fo", "Ch","Ps", "Sc"))) %>% 
  select(-value, -lower_fence)

#============enrichment of discordance within individual============#
m <- tmp3 %>% 
  mutate(discordance = case_when(outlier == "Discordant" ~ 1, 
                                 outlier == "Concordant" ~ 0,
                                 TRUE ~ NA_real_)) %>% 
  select(Var1_Site, Var1_SubjectID, discordance) %>% 
  pivot_wider(names_from = "Var1_Site", values_from = "discordance") %>% 
  column_to_rownames(var = "Var1_SubjectID")
m$n_discordance <- rowSums(m, na.rm = TRUE)
discordance_prevalence <- tmp3 %>%
  mutate(discordance = case_when(outlier == "Discordant" ~ 1, 
                                 outlier == "Concordant" ~ 0, 
                                 TRUE ~ NA_real_)) %>%  # Ensure the discordance column exists here
  group_by(Var1_Site) %>%
  summarise(Prevalence = mean(discordance, na.rm = TRUE), .groups = "drop") %>% 
  ungroup() %>%
  rename(Site = Var1_Site)

total_non_na <- m %>%
  select(-n_discordance) %>%  # Exclude the n_discordance column
  summarise_all(~sum(!is.na(.))) %>%  # Count non-NA values in each column
  sum()  # Sum up the non-NA counts for all columns
total_non_na
average_discordance_prevalence <- sum(m$n_discordance)/total_non_na
p_vec <- setNames(rep(average_discordance_prevalence, 9), discordance_prevalence$Site)

subject_df <- m %>% rownames_to_column(var = "SubjectID")

subject_enrichment <- subject_df %>%
  rowwise() %>%
  mutate(
    # ppoibin: cumulative probability of observing <= (n-1) discordances;
    # subtract from 1 gives the p-value for observing >= n discordances.
    p_value = 1 - poibin::ppoibin(n_discordance - 1, p = p_vec)
  ) %>%
  ungroup() %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Filter significant subjects after correction
significant_subjects <- subject_enrichment %>%
  filter(adj_p_value  < 0.05)

# Load necessary library
library(ggplot2)
library(ggrepel)

p <- ggplot(subject_enrichment, aes(x = n_discordance)) +
  # Histogram for the population distribution
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  # Add density curve for smoother visualization
  geom_density(aes(y = after_stat(count)), color = "darkblue", linewidth = 1) +
  scale_x_continuous(breaks = seq(0,9, by = 1 )) + 
  # Highlight significant subjects
  geom_point(data = subset(subject_enrichment, adj_p_value < 0.05), 
             aes(x = n_discordance, y = 0), 
             color = "red", size = 1, shape = 16) +
  # Add labels for significant subjects with ggrepel
  geom_text_repel(data = subset(subject_enrichment, adj_p_value < 0.05), 
                  aes(x = n_discordance, y = 0, label = SubjectID), 
                  color = "red", size = 3, box.padding = 1, max.overlaps = 10) +  # Box padding and max overlaps
  # Customize theme and labels
  labs(title = "",
       x = "Number of discordant sites within individuals",
       y = "Number of individuals") +
  ggpubr::theme_pubr(base_size = 7)
p


#===========check bilateral discordance distribution cross-sites ======#
############ Fisher test ############ 
df = tmp3 %>% 
  select(Var1_Site, Var1_SubjectID, outlier) %>% 
  pivot_wider(names_from = "Var1_Site", values_from = "outlier") %>% 
  column_to_rownames(var = "Var1_SubjectID")
cols <- colnames(df)
combinations <- combn(cols, 2)
test_pair <- function(index) {
  col1 <- combinations[1, index]
  col2 <- combinations[2, index]
  table_data <- table(df[, col1], df[, col2])
  
  # Check for empty tables
  if (any(dim(table_data) == 0)) {
    return(data.frame(S1 = col1, S2 = col2, 
                      P_Value = NA, OR = NA))
  }
  
  test_res <- tryCatch({
    fisher.test(table_data, alternative = "two.sided")
  }, error = function(e) NULL)
  
  # If test_res is NULL due to an error in fisher.test
  if (is.null(test_res)) {
    p_val <- NA
    odds <- NA
  } else {
    p_val <- test_res$p.value
    odds <- ifelse(is.null(test_res$estimate), 999, test_res$estimate)
  }
  
  return(data.frame(S1 = col1, S2 = col2, 
                    p_val = p_val, OR = odds))
}
# Apply the function to each combination
results <- do.call(rbind, lapply(1:ncol(combinations), test_pair))
results$P_Value <- p.adjust(results$p_val, method = "BH")

results_1 <- results %>% 
  mutate(P_Value = case_when(P_Value >= 0.05 ~ NA, 
                             T ~ P_Value)) %>% 
  mutate(OR = case_when(P_Value < 0.05 ~ OR, 
                        T ~ NA)) %>% 
  select(-p_val)

results_2 <- results_1 |>
  rename(c("S2" = "S1", "S1" = "S2")) |>
  dplyr::select(S1, S2, P_Value, OR)

results_sorted <- rbind(results_1, results_2)
results_sorted$log_P_Value <- log10(results_sorted$P_Value)
unique_sites <- unique(c(results_sorted$S1, results_sorted$S2))
all_combinations <- expand.grid(S1 = unique_sites, S2 = unique_sites)
full_data <- merge(all_combinations, results_sorted, by = c("S1", "S2"), all.x = TRUE)

desired_order <- c("Ps", "Sc", "Fo", "Ch", "Ub", "Vf", "Ll", "Ac","Ax")

full_data$S1 <- factor(full_data$S1, levels = desired_order)
full_data$S2 <- factor(full_data$S2, levels = desired_order)

p <- 
  ggplot(full_data, aes(x = S1, y = S2)) +
  geom_tile(aes(fill = log_P_Value), color = "white") +
  geom_text(aes(label = ifelse( OR != 999, round(OR,1), "")), 
            size = 3.5, na.rm = TRUE) +
  scale_fill_gradient(low = "red", high = "yellow", 
                      na.value = "transparent", 
                      name = "Bilateral\ndiscordance\ncross-site\ncorrelation\n \nlog10\n(Fisher's\np-value)") +
  scale_y_discrete(limits = levels(full_data$S1)) + 
  #theme_minimal() + 
  ggpubr::theme_pubr(margin = T, border = T, legend = "right") +
  theme(
    panel.background = element_blank(), 
    axis.ticks = element_blank(), 
    axis.title = element_blank(),
    legend.box.background = element_rect(colour = "black")) 

print(p)

#================Discordance vs host factors=============#
results_df_1 <- read.csv("./Data/2_basic/BilateralDiscordance_host_factor_assoc.csv")
## visualization ##
library(ggplot2)
# Convert PValue into asterisks for significance
results_df_2 <- results_df_1 %>%
  mutate(p_asterisk = case_when(
    PValue <= 0.01 ~ "***",
    PValue <= 0.05 ~ "**",
    PValue <= 0.1 ~ "*",
    TRUE ~ ""
  )) 

# Create the plot
results_df_3 <- results_df_2 %>%
  mutate(
    Factor_Category = gsub("Skin|_skin", "", Column_Category_2), 
    Factor_Category = factor(Factor_Category, 
                             levels = c("Demographic", "Physio",
                                        "Behavior", "Subclinical",
                                        "Disease", "GeneralHealth")),
    Data_type = factor(Data_type, levels = c("Character", "Numeric"))) %>%  
  mutate(Site = factor(Site, levels = c("Ax", "Ll", "Ac", "Vf", "Ub", "Fo", "Ch", "Ps", "Sc")))

#============== rearrange the order of phenotypes for vis ===========# 
results_df_3 <- results_df_3 %>%
  filter(grepl("not described above", Formatted_Predictor) == F) %>% 
  mutate(Formatted_Predictor = ifelse(Formatted_Predictor == "PH forearm", "pH forearm", Formatted_Predictor)) %>% 
  filter(p_Value < 0.05) %>%
  arrange(Factor_Category, 
          Column_Category, 
          factor(Column_Subcategory, levels = c(setdiff(unique(Column_Subcategory), "General"), "General")), 
          desc(Data_type))

# Reorder the NewTerm factor levels based on the arranged dataframe
results_df_3$Formatted_Predictor<- factor(results_df_3$Formatted_Predictor, 
                                          levels = unique(results_df_3$Formatted_Predictor))

p_outlier_phenotype_sub <- results_df_3 %>%
  filter(p_asterisk != "") %>% 
  ggplot(aes(x = Site, y = Formatted_Predictor)) +
  geom_tile(aes(fill = EffectSize), color = "white") +
  geom_text(aes(label = p_asterisk), vjust = 0.2, hjust = 0.5, size = 4, color = "black") +
  geom_text(aes(label = ifelse(abs(EffectSize) < 0.2 & EffectSize < 0 , "-", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "grey30") +
  geom_text(aes(label = ifelse(abs(EffectSize) < 0.2 & EffectSize > 0 , "+", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "grey30") +
  geom_text(aes(label = ifelse(abs(EffectSize) >= 0.2 & EffectSize < 0 , "-", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "white") +
  geom_text(aes(label = ifelse(abs(EffectSize) >= 0.2 & EffectSize > 0 , "+", "")), vjust = 0.8, hjust = 0.5, size = 4, colour = "white") +
  scale_fill_gradient2(low = "#013220", mid = "grey90", high = "maroon", midpoint = 0, guide = "colorbar") +
  labs(title = "", x = "", y = "", fill = "Beta\ncoefficient") +
  facet_grid(rows = vars(Factor_Category),
             labeller = label_wrap_gen(),
             cols = vars(Site),
             space = "free", 
             scale = "free") + 
  theme_pubr(border = T) + 
  theme(axis.text.x.bottom =  element_blank(), 
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.direction = "vertical",
        legend.position = "right",
        #legend.position.inside = c(-1.5, 0.5), 
        legend.background = element_rect(linewidth = 0.2, color = "black", fill = "transparent")
  )
p_outlier_phenotype_sub
