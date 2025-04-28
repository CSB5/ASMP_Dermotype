library(ggpubr)
library(vegan)
library(reshape2)
library(tidyverse)
library(tibble)
library(dplyr)
library(svglite)
library(tidyr)
library(rio)
library(here)
library(skimr)
library(tidyverse)
library(gtsummary)
library(rstatix)
library(janitor)
library(scales)
library(flextable)
library(googlesheets4)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(broom)

mat <- read.csv("./Data/Taxonomic abundance table.csv")
mat = mat %>% 
  column_to_rownames(var = "X")
Library_to_keep <- rownames(mat)

#===== load the questionnaire data + basic dermographics ==== #
mdata.ex <- read.csv("~/Desktop/ASMP/Data/SkinHealth_Questionnaire_perSample.csv", stringsAsFactors = F)
mdata.ex <- mdata.ex %>%
  mutate(SubjectID = str_pad(SubjectID, 5, pad = "0")) %>% 
  filter(LibraryID %in% Library_to_keep) %>% 
  rowwise() %>% 
  mutate(Itchy_average_1_to_10 = 
           mean(c_across(contains("Itch") & ends_with("1_to_10")), na.rm = TRUE)
  ) %>% 
  ungroup() %>% data.frame()

rownames(mdata.ex) <- mdata.ex$LibraryID

#===== load the LR Pearson correlation coefficient ==== #
load("./Data/2_basic/Pearson.RData")

# Calculate the upper and lower fences for each group
boxplot_stats <- tmp2 %>% 
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
  mutate(Var1_SubjectID = str_pad(Var1_SubjectID, 5, pad = "0")) 

tmp4 <- mdata.ex %>% 
  select(-c(LibraryID, Dermotype_size, Stability)) %>% 
  distinct() %>% 
  dplyr::right_join(tmp3, by = join_by(SubjectID == Var1_SubjectID, 
                                Site == Var1_Site)) %>%
  mutate(Site = factor(Site, levels = c("Ax", "Ll", "Vf", "Ac", "Ub",  
                                        "Fo", "Ch","Ps", "Sc"))) %>% 
  select(-value, -lower_fence)


#========visualize specific associations =========#
names(tmp4)
#------------1. categorical variables ----------------#
plot_site_var <- function(site, var, var_level, var2 = NULL, var2_level = NULL) {

  data <- tmp4 %>%
  filter(Site == site) 

  if (!is.null(var2)) {
    data <- data %>%
      filter(.data[[var2]] == var2_level)
  }

proportion_data <- data %>% 
  group_by(outlier) %>%
  summarise(
    count_yes = sum(.data[[var]] == var_level),
    total = n(),
    proportion = count_yes / total,
  )

# Perform Fisher's exact test
test_data <- data %>% 
  group_by(outlier) %>%
  summarise(count_yes = sum(.data[[var]] == var_level), 
            count_no = sum(.data[[var]] != var_level))

# Prepare data for Fisher's test
contingency_table <- matrix(c(test_data$count_yes, test_data$count_no), nrow = 2)

# Perform the Fisher's exact test
fisher_test <- fisher.test(contingency_table)
p_value <- fisher_test$p.value

# Get the maximum proportion for correct y positioning
max_proportion <- max(proportion_data$proportion)

proportion_data$Site = site
proportion_data$var = var
proportion_data$pvalue = fisher_test$p.value

# Create the bar plot with error bars and p-value annotation
p1 <- proportion_data  %>%
  mutate(site = site) %>% 
  data.frame() %>% 
  mutate(var = str_replace_all(var, "_", " ")) %>% 
  mutate(label = paste0(var, " (", var_level, ")")) %>% 
  mutate(outlier = paste0(outlier, "\n(n=", total, ")")) %>% 
  ggplot(aes(x = outlier, y = proportion, fill = outlier)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("blue", "red")) +
  facet_grid(cols = vars(site)) +
  ylab(paste("Individuals with", tolower(var_level), tolower(str_replace_all(var, "_", " ")), "(%)")) +
  xlab("") +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), limit = c(0, max_proportion*1.1)) +
  ggpubr::theme_pubr(legend = "none", border = T, base_size = 7, x.text.angle = 30)+ 
  annotate("text", x = 1.2, y = max_proportion,
           label = paste0("Fisher's test\np = ", round(p_value, 3))) 
return(p1)
}
plot_site_var("Ps", "Hair_covering_frequency", "Daily") -> p1
#------------2. continous variables ----------------#
plot_site_con_var <- function(site, var,  var2 = NULL, var2_level = NULL) {
  data <- tmp4 %>%
    filter(Site == site)
  
  if (!is.null(var2)) {
    data <- data %>%
      filter(.data[[var2]] == var2_level)
  }
  max_var <- max(data[[var]])
  
  p_box <- data %>%
    group_by(outlier) %>% 
    mutate(total = n()) %>% 
    ungroup() %>% 
    mutate(outlier = paste0(outlier, "\n(n=", total, ")")) %>% 
    ggplot(aes(x = outlier, y = .data[[var]])) +
    geom_boxplot(outliers = F) +
    geom_jitter(alpha = 0.7, aes(color = outlier), show.legend = NA, size = 0.5) +
    scale_color_manual(values = c("red", "blue")) +
    scale_y_continuous(limits = c(0, max_var*1.1)) +
    facet_grid(cols = vars(Site)) + 
    ylab(str_replace_all(var, "_", " "))+
    xlab("") +
    stat_compare_means(method = "wilcox", label.x.npc = 0.4, label.y.npc = 0.95,
                       label.sep = "\n") +
    theme_pubr(legend = "none", border = T, base_size = 7, x.text.angle = 30)
  
  return(p_box)
}

plot_site_con_var("Fo", "Irritation_itching_1_to_10"
                  ) -> p4
#-----------putting everything together ----------#
library(cowplot)
p_set <- plot_grid(p1, 
                   #p2, 
                   #p3, 
                   p4, #p5, p6, 
                   nrow = 1)
p_set 

