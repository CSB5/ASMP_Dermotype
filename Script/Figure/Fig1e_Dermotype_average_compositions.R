library(dplyr)
library(ggplot2)
library(scales)
library(tidyverse)
library(svglite)

mat <- read.csv("./Data/Taxonomic abundance table.csv")
mat = mat %>% 
  column_to_rownames(var = "X")
mdata.dermotype.conf <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")


colors =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000","#E6AB02",
             "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#F46D44", 
             "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", 
             "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6fabd0", 
             "#1B9E77", "#D95F02", "#7570B4", "#FB8072",  "#A6761D", "#CCCCCC")

tmp = mat |>
      rownames_to_column(var = "LibraryID") |>
      left_join(mdata.dermotype.conf[, c("Dermotype_size", "LibraryID", "Stability")], join_by(LibraryID)) %>%
      # filter(Stability >= 0.9 | grepl("Ch|Ll", Dermotype_size)) %>% 
      dplyr::group_by(Dermotype_size) %>% 
      dplyr::mutate(count = n()) %>% 
      ungroup() %>% 
      mutate(Dermotype_size = paste0(Dermotype_size, " (n= ",count, ")")) %>% 
      select(-count, -Stability)

# ======== average composition ========== #
tmp_long <- tmp %>% 
  gather(key = "Taxa", value = "Value", -LibraryID, -Dermotype_size)

# Calculate the median for each taxon within each Dermotype_size
median_values <- tmp_long %>% 
  dplyr::group_by(Dermotype_size, Taxa) %>% 
  dplyr::summarize(Median_Value = median(Value, na.rm = TRUE))

# For each Dermotype_size, get the top 15 taxa by median value
top_taxa <- median_values %>% 
  filter(Median_Value > 0 ) %>%
  dplyr::group_by(Dermotype_size) %>%
  arrange(Dermotype_size, desc(Median_Value)) %>% 
  slice_head(n = 15)

taxa_summary <- top_taxa %>%
  dplyr::group_by(Taxa) %>%
  dplyr::summarise(
    num_appearances = n(),
    mean_median_value = mean(Median_Value, na.rm = TRUE)
  ) %>% 
  arrange(desc(num_appearances), desc(mean_median_value))

# Calculate the mean relative abundance for each taxon within each Dermotype_size
mean_abundance <- tmp_long %>% 
  dplyr::group_by(Dermotype_size, Taxa) %>% 
  dplyr::summarise(Mean_Abundance = mean(Value, na.rm = TRUE))

# Filter to get only those taxa from the top_taxa
filtered_mean_abundance <- mean_abundance %>%
  semi_join(top_taxa, by = c("Taxa"))

# Row bind the result to a new dataframe
result <- filtered_mean_abundance %>% 
  dplyr:: select(Taxa, value = Mean_Abundance, Dermotype_size)

# Calculate the sum of value for each Dermotype_size category
sum_values <- result %>% 
  dplyr::group_by(Dermotype_size) %>% 
  dplyr::summarise(Total_Abundance = sum(value))

# Calculate the value for "Others"
sum_values <- sum_values %>% 
  mutate(value = 1 - Total_Abundance) %>% 
  dplyr::select(-Total_Abundance)

# Create a dataframe for "Others"
others_df <- sum_values %>% 
  mutate(Taxa = "Others")

# 4. Row bind "Others" dataframe to the original result
final_result <- bind_rows(result, others_df)
final_result$Taxa <- str_replace_all(final_result$Taxa, c('s__' = '', '_' = ' '))
final_result$Site <- sub("\\-.*", "", final_result$Dermotype_size)
final_result$Site <- sub("\\(.*", "", final_result$Site)
final_result$Site <- sub("\\ ", "", final_result$Site)
final_result$Site <- factor(final_result$Site, 
                            levels = c("Ax", "Ac", "Ll", "Vf", "Ub", "Fo", "Ch", "Ps", "Sc"))
p_stackplot = ggplot(final_result, aes(x = Dermotype_size, fill = Taxa)) + 
  geom_bar(aes(weight = value), position = position_fill(reverse = T)) +
  facet_grid(.~ Site, scales = "free", space = "free") +
  scale_fill_manual(values = colors[c(1:n_distinct(final_result$Taxa)-1, 30)]) +
  ylab("Relative Abundance") + 
  theme_pubr( legend = "right",
              x.text.angle = 45, 
              margin = F, 
              border = T, 
              base_size = 7) +  # )
  theme(
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        axis.title.x=element_blank(),   
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(face = "bold", angle = 0),       # Bold facet labels
        strip.background = element_rect( fill = "white"),
        ) + # Remove x-axis ticks
  guides(fill = guide_legend(keyheight = 0.3, reverse = F, position = "bottom", ncol = 3)) + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0))


 p_stackplot$data$Taxa = factor(p_stackplot$data$Taxa, ordered = TRUE, 
                               levels = c(str_replace_all(taxa_summary$Taxa, c('s__' = '', '_' = ' ')),
                                          "Others"))

 
 print(p_stackplot)
