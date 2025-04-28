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
library(ggplot2)

df_sig_vis <- read.csv("./Data/4_host_factor/sig_skin_phenotype_assoc.csv")
p <- df_sig_vis %>%
  ggplot(aes(x = Dermotype, y = Formatted_Predictor)) +
  geom_tile(aes(fill = EffectSize), color = "white") +
  geom_text(aes(label = p_asterisk), vjust = -0.1, hjust = 0.5, color = "black", size.unit = "pt", size = 8) +
  geom_text(aes(label = ifelse(abs(EffectSize) < 0.5 & EffectSize < 0 , "-", "")), vjust = 0.7, hjust = 0.5, colour = "grey30", size.unit = "pt", size = 5) +
  geom_text(aes(label = ifelse(abs(EffectSize) < 0.5& EffectSize > 0 , "+", "")), vjust = 0.7, hjust = 0.5, colour = "grey30", size.unit = "pt", size = 5) +
  geom_text(aes(label = ifelse(abs(EffectSize) >= 0.5 & EffectSize < 0 , "-", "")), vjust = 0.7, hjust = 0.5, colour = "white", size.unit = "pt", size = 5) +
  geom_text(aes(label = ifelse(abs(EffectSize) >= 0.5& EffectSize > 0 , "+", "")), vjust = 0.7, hjust = 0.5, colour = "white", size.unit = "pt", size = 5) +
  scale_fill_gradient2(low = "#013220", mid = "grey90", high = "maroon", midpoint = 0, guide = "colorbar") +
  labs(title = "", x = "", y = "", fill = "Beta \ncoefficient") +
  ggh4x::facet_grid2(
    rows = vars(Factor_Category),
    labeller = label_wrap_gen(),
    cols = vars(Site),
    space = "free", 
    scale = "free") + 
  theme_pubr(border = T, base_size = 12, x.text.angle = 45) + 
  theme(
    legend.position = "right", 
    legend.direction = "vertical"
    #legend.position.inside = c(-1.2, 0.9)
  )
p


