library(dplyr)
library(tidyverse)
library(tibble)
library(readxl)
library(lmerTest)
library(ggpubr)
library(mltools)
library(data.table)
library(googlesheets4)
library(effectsize)
library(stats)

df_sig_one_hot <- read.csv(file = "./Data/4_host_factor/sig_demorgraphic_assoc.csv")
# vis #
p_demographic <- 
  df_sig_one_hot %>% 
  # --- aesthetic  purpose ----#
  mutate(Formatted_Predictor  = 
           factor(Formatted_Predictor, 
                  levels = c("Age", 
                             "Gender (Female)",
                             "Ethnicity (South Asian)", 
                             "Ethnicity (East Asian)", 
                             "Ethnicity (South East Asian)"
                             ))) %>% 
                        
  # --- mann-whitney u test ----#
  ggplot(., aes(x = Dermotype, y = Formatted_Predictor )) +
    geom_tile(
      data = subset(df_sig_one_hot, Test == "Mann-Whitney U"),
      aes(fill = EffectSize), color = "white"
    ) +
  scale_fill_gradient2(
    low = "#013220", mid = "white", high = "maroon", midpoint = 0, 
    na.value = "black", name = "Mann-Whitney U test\nrank biserial correlation", 
    limits = c(max(df_sig_one_hot$EffectSize[df_sig_one_hot$Test == "Mann-Whitney U"])*(-1.1),
      max(df_sig_one_hot$EffectSize[df_sig_one_hot$Test == "Mann-Whitney U"])*1.1), 
    guide = guide_colorbar(barwidth = 0.5, barheight = 5) 
  ) +
  ggnewscale::new_scale_fill() +
  # --- fisher's test ----#
  geom_tile(
    data = subset(df_sig_one_hot , Test == "Fisher's Exact"),
    aes(fill = EffectSize), color = "white"
  ) +
  scale_fill_gradient2(
    "low" = "purple", mid = "white", "high" = "orange", midpoint = 1,
    na.value = "black", name = "Fisher's Test\nodds ratio", 
    guide = guide_colorbar(barwidth = 0.5, barheight = 5) 
  ) +
  geom_text(aes(label = p_asterisk), color = "black", size = 4, hjust = 0.5, vjust = 0.5) +
  labs(title = "", x = "", y = "") +
  facet_grid(cols = vars(Site), scales = "free", space = "free") +
  ggpubr::theme_pubr(base_size = 12, legend = "right", x.text.angle = 45) +
  theme(
    legend.position = "right",         # Place legends on the right
    legend.box = "horizontal",         # Arrange legends side by side
    legend.spacing.x = unit(0.5, "cm") # Add spacing between legends
  )
p_demographic

 

