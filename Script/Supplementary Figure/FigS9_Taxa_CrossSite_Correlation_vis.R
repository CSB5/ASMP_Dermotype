library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

tmp <- read.csv("~/ASMP_Dermotype/Data/2_basic/Taxa_cross_site_correlation_all.csv")
tmp_out <- read.csv("~/ASMP_Dermotype/Data/2_basic/Taxa_cross_site_correlation.csv")

#-------visualize spearman r distirbution of each taxa---#
# Compute median values for each Taxa1
medians <- tmp_out %>%
  group_by(Taxa1) %>%
  summarise(med = median(Spearman_r, na.rm = TRUE),
            ave = mean(Spearman_r, na.rm = T), 
            max = max(Spearman_r, na.rm = T))

p_overall_distribution <- 
  ggplot(tmp_out, aes(x = reorder(Taxa1, Spearman_r, median, decreasing = T), 
                      y = Spearman_r)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = ifelse(Site1 == "Ax" | Site2 == "Ax", "red", "black")),
              alpha = 0.5) +
  scale_color_identity() +
  labs(x = "", y = "Cross-site correlation of \ncore species abundance (Spearman rho)") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  # Replace underscores with spaces
  ggpubr::theme_pubr(x.text.angle = 45, base_size = 10) +
  theme(axis.text.x = element_text(face = "italic")) +
  
  # Add median labels
  geom_text(data = medians, aes(x = Taxa1, y = med, label = round(med, 3)), 
            vjust = -0.5, size = 5, color = "blue")

p_overall_distribution


#-------visualize spearman r distirbution of each taxa---#
#----but exclude axilla site given its distinctness from other skin sites ---#
medians_no_Ax <- tmp_out %>%
  #filter(Site1 == "Ax" | Site2 == "Ax") %>% 
  filter(Site1 != "Ax" & Site2 != "Ax") %>% 
  filter(Taxa1 == "Malassezia_restricta") %>% 
  group_by(Taxa1) %>%
  summarize(med = median(Spearman_r, na.rm = TRUE),
         mean = mean(Spearman_r, na.rm = TRUE),
         max =  max (Spearman_r, na.rm = TRUE))


  tmp_out %>%
  filter(Site1 != "Ax" & Site2 != "Ax") %>% 
  ggplot(., aes(x = reorder(Taxa1, Spearman_r, median, decreasing = T), 
                      y = Spearman_r)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = ifelse(Site1 == "Ax" | Site2 == "Ax", "red", "black")),
              alpha = 0.5) +
  scale_color_identity() +
  labs(x = "", y = "Cross-site correlation of \ncore species abundance (Spearman rho)") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +  # Replace underscores with spaces
  ggpubr::theme_pubr(x.text.angle = 45, base_size = 10) +
  theme(axis.text.x = element_text(face = "italic")) +
  
  # Add median labels
  geom_text(data = medians_no_Ax, aes(x = Taxa1, y = med, label = round(med, 3)), 
            vjust = -0.5, size = 5, color = "blue") -> p_overall_distribution_no_Ax

p_overall_distribution_no_Ax

median_total <- medians %>% 
  rename( med_with_Ax = med ) %>% 
  left_join(medians_no_Ax)
view(median_total)

cowplot::plot_grid(
  p_overall_distribution,  
  p_overall_distribution_no_Ax, ncol = 1)

#----------visualize for each taxa ----------------------#
# Ensure symmetry by creating a unique key for combinations
tmp_filled <- tmp %>%
  rowwise() %>%
  mutate(
    combo_key = paste(sort(c(paste(Site1, Taxa1), paste(Site2, Taxa2))), collapse = " - ")
  ) %>%
  ungroup() %>%
  group_by(combo_key) %>% # Group by the symmetrical combination key
  mutate(
    p_value = ifelse(is.na(p_value), first(na.omit(p_value)), p_value),
    adjust_p = ifelse(is.na(adjust_p), first(na.omit(adjust_p)), adjust_p),
    Annotation = ifelse(is.na(Annotation), first(na.omit(Annotation)), Annotation)
  ) %>%
  ungroup() %>%
  select(-combo_key) # Remove the temporary combo_key column


# Function to plot a heatmap for each Taxa1-Taxa2 combination
plot_heatmap <- function(data, taxa1, taxa2) {
  heatmap_data <- data %>%
    filter(Taxa1 == taxa1 & Taxa2 == taxa2) %>%
    select(Site1, Site2, Spearman_r, Annotation)
  
  # Create the heatmap
  p <- ggplot(heatmap_data, aes(x = Site2, y = Site1, fill = Spearman_r)) +
    geom_tile(color = "white") +  # Heatmap tiles
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", na.value = "black",
      midpoint = 0, limits = c(-1, 1), name = "Spearman r"
    ) +
    geom_text(aes(label = Annotation), color = "black") +  # Add annotations
    labs(
      x = gsub("_", " ", taxa1), y = gsub("_", " ", taxa2)
    ) +
    ggpubr::theme_pubr(x.text.angle = 45, base_size = 10, legend = "right") +
    theme(legend.key.width = unit(0.2, "cm"), 
          axis.title = element_text(face = "italic"))
  
  return(p)
}

# Generate heatmaps for each Taxa1 and Taxa2 combination
unique_combinations <- tmp_filled %>%
  distinct(Taxa1, Taxa2) %>%
  arrange(Taxa1, Taxa2)

heatmaps <- lapply(1:nrow(unique_combinations), function(i) {
  taxa1 <- unique_combinations$Taxa1[i]
  taxa2 <- unique_combinations$Taxa2[i]
  plot_heatmap(tmp_filled, taxa1, taxa2)
})

for (i in seq_along(heatmaps)) {
  print(heatmaps[[i]])
}




