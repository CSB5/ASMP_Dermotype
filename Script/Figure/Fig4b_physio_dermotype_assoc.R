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
library(patchwork)

tmp <- read.csv("./Data/4_host_factor/physio_measurement.csv")
df_sig <- read.csv("./Data/4_host_factor/sig_physio_variation.csv")
# ----approach 1: we do pair-wise comparison between dermotypes for sites with global significance
physio_list <- c("TEWL_forearm",
                 "Surface_hydration_forearm",
                 "PH_forearm",
                 "TEWL_cheek",
                 "Sebum_cheek")
#============ other sites  ==========#
p_physio <- list()
list_site <- c("Ax", "Ac", "Vf",  "Ub", "Fo", "Ps","Sc")
color_dermotype <- list(c('#66A61E', '#FF0000', '#5E4FA2',  '#F46D44' ,'#4288BD'),
                        c("green", "#728c69"),
                        c("blue",  "#73c2fb"), 
                        c("purple",  "#C64B8C"), 
                        c("pink",  "#FF007F"),
                        c("#F46D44", "#00008B"),
                        c("#F46D44", "#00008B"))

# Get only the significant Site-Phenotype pairs
sig_pairs <- df_sig %>% select(Site, Phenotype)
p_physio <- list()
for (i in seq_along(list_site)) {
  Pickedsite <- list_site[i]
  sig_physio_for_site <- sig_pairs %>%
    filter(Site == Pickedsite) %>%
    pull(Phenotype)
  
  # Skip sites that have no significant phenotypes
  if (length(sig_physio_for_site) == 0) next 
    site_data <- tmp %>%
      filter(Site == Pickedsite) 
  for (j in seq_along(sig_physio_for_site)) {
    Pickedphysio <- sig_physio_for_site[j]
    
    y_label <- str_replace_all(Pickedphysio, "_", " ")
    
    p_physio[[Pickedphysio]][[Pickedsite]] <- site_data %>%
      ggplot(aes(x= Dermotype_size, y= .data[[Pickedphysio]], fill= Dermotype_size)) + 
      stat_boxplot(geom = "errorbar", width = 0.5) + 
      geom_boxplot(outliers = FALSE) + 
      geom_jitter(size = 0.5, width = 0.2, height = 0, alpha = 0.3) +
      scale_fill_manual(values= color_dermotype[[i]]) + 
      geom_pwc(method = 'wilcox_test',
               p.adjust.method = "BH",
               label = "p.signif", hide.ns = TRUE, 
               tip.length = 0.01, label.size = 8, 
               vjust = 1, color = "black") +
      ggpubr::theme_pubr(x.text.angle = 30, base_size = 10) + 
      labs(y= y_label) +
      theme(
            legend.position="none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = rel(1.2)))
  }
}

all_plots <- flatten(p_physio)
final_plot <- wrap_plots(all_plots, nrow = 1) +
  plot_layout(widths = c(2, rep(1, length(all_plots)-1)))
final_plot  

# Combine plots by Pickedphysio and facet by Site
combined_plots <- list()
for (physio in names(p_physio)) {
  plots_for_physio <- p_physio[[physio]]
  
  # Calculate the number of unique Dermotype_size levels for each Site
  dermotype_counts <- sapply(plots_for_physio, function(plot) {
    length(unique(plot$data$Dermotype_size))
  })
  
  # Use the counts to set relative widths
  relative_widths <- dermotype_counts / sum(dermotype_counts)
  
  # Combine plots with proportional widths
  combined_plot <- wrap_plots(plots_for_physio, nrow = 1, widths = relative_widths) +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'none')
  
  combined_plots[[physio]] <- combined_plot
}

# Combine all combined_plots into one final plot
final_plot <- wrap_plots(combined_plots, nrow  = 1, 
                         widths = c(9,1.5,4,1.5))
final_plot


