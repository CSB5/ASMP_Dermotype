library(dplyr)
library(tidyverse)
library(readr)
local_path <- "/Users/chengchenli/Desktop/ASMP/"
load(paste0(local_path, "Data/Cleaned_Pathway_Abund.RData"))
load(paste0(local_path,"Data/Tax_BatchABCD_AbundFilt.RData"))
DA_taxa <- read.csv(paste0(local_path, "Output/3_ML/DiffAbundTaxa_Dermotype.csv"))
phenotype <- read.csv(paste0(local_path,"Data/SkinHealth_Questionnaire_perSample.csv"))
source("~/ASMP_network_pwy/script/5_DA_pathway_lineage_reformating.R")
load( "./pwy_asso_phenotype.RData")

#==========  relab (TSS after removing unmapped + unintergrated + highly associated w/ kitome) ===== 
load("./Data/Cleaned_Pathway_Abund.RData")
pwy_clean_site <- pwy_clean %>% 
  rownames_to_column("LibraryID") %>% 
  left_join(mdata.dermotype.conf) %>% 
  filter(Site == "Ac") 
pwy_clean_site %>% 
  #filter(`LACTOSECAT-PWY` >0) %>% 
ggplot(., aes(x = `ARGDEG-PWY` , y = `HISDEG-PWY`)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Dermotype_size) +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(
    method = "spearman",
    color = "black",
    #label.x = df$min_abundance,  # Use group-specific min (original scale)
    #label.y = df$max_value,      # Use group-specific max (original scale)
    hjust = -0.1,                # Adjust horizontal position
    vjust = 1.5                  # Adjust vertical position
  ) +
  theme_pubr()


cor.test(pwy_clean_site$`LACTOSECAT-PWY`, pwy_clean_site$`HISDEG-PWY`)
#===== check pathway contribution  ===
pwy_ctr <- read_tsv(paste0(local_path,"Output/3_ML/Pathway_from_Indrik/combined_pwy_relab.tsv"))
names(pwy_ctr) <- gsub(".humann2_Abundance", "", names(pwy_ctr))
names(pwy_ctr)[1] <- "pwy"
pwy_ctr_sub <- pwy_ctr %>%
  data.frame() %>%
  select(pwy, any_of(mdata.dermotype.conf$LibraryID)) %>%
  #filter(!str_detect(pwy, "\\|")) %>%
  filter(!str_detect(pwy, "UNMAP")) %>% 
  filter(!str_detect(pwy, "UNINTEGRATED")) %>%
  separate(pwy, c("ID", "pwy_taxa"), ":" ) %>%
  #separate(pwy_taxa, c("pwy", "taxa"), "\\|") %>% 
  # filter(is.na(taxa) == T) %>%
  mutate_at(vars(matches("WMS")), function(x) replace_na(x, 0))

source("~/ASMP_network_pwy/script/UtilityFunction_prep_data_for_pwy_species_contribution.R")

#==== pick a pathway =====
i = "Ac"

j = "HISDEG-PWY" 
j = "LACTOSECAT-PWY"
j = "PWY-5028"
j = "PWY-5030"

output <- pwy_contri_by_group(i, # skin site of interest
                              j	, # pick a pathway (shortened name, not the full name)
                              mdata.dermotype.conf, #mdata files 
                              pwy_ctr_sub,
                              pwy_clean, 
                              phenotype,
                              pwy_asso_phenotype,
                              n_taxa = 5
                              ) 


library(cowplot)
# Create the combined plot
combined_plot <- plot_grid(
  output[["vis_relab"]],  # First plot on the left
  output[["vis_average_species_contribution"]],  # Second plot
  ncol = 2,  # Arrange into 2 columns
  rel_widths = c(0.5, 1), 
  align = "h"# Adjust relative widths if needed
)

# Add title while keeping plots aligned
final_plot <- plot_grid(
  ggdraw() + draw_label(output[["pathway_fullname"]], size = 12, hjust = 0.5),  
  combined_plot,  
  ncol = 1,  # Stack title and plot vertically
  rel_heights = c(0.3, 1)  # Adjust title height relative to plots
)
final_plot
ggsave(final_plot, 
       filename = paste0("pwy_", i, "_", j, ".svg" ),
       width = 7, 
       height = 2)


#==== check pathway-phenotype associations ==========
df <- pwy_clean_site %>%
  left_join(select(phenotype, -SubjectID), join_by(LibraryID)) %>% 
  rowwise() %>% 
  mutate(Irritation_severity = 
           mean(c_across(starts_with("Irritation") & ends_with("1_to_10") ), na.rm = TRUE)) %>% 
  group_by(SubjectID) %>% 
  mutate(ave_HISDEG_PWY = mean(`HISDEG-PWY`))
            
df1 <- df %>% 
  distinct(SubjectID, ave_HISDEG_PWY, Sensitive_skin_overall) 
table(df1$Sensitive_skin_overall)

df1 %>% 
  ggplot(., aes(x = Sensitive_skin_overall, 
                 y = ave_HISDEG_PWY)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.6) + 
    scale_y_log10(labels = scales::percent_format(suffix = "")) +
    theme_pubr() +
  geom_pwc(method = "wilcox_test", 
           method.args = list(#alternative = "greater", 
                             p.adjust.method = "BH"), 
                             label = "p.adj.signif", 
                             color = "red",
                             vjust = 1) +
  labs(x = "Sensitive skin", y = "HISDEG-PWY abundance (%)")

# df %>% 
#   distinct(SubjectID, ave_HISDEG_PWY, Irritation_severity) %>% 
#   #filter(Irritation_severity >0) %>% 
#   ggplot(., aes(x = Irritation_severity, 
#                 y = ave_HISDEG_PWY)) +
#   geom_point() +
#   scale_y_log10() +
#   theme_pubr()

#======check histidine degradation pwy and galactose degradation pwy======#
pwy_clean_site %>% 
  ggplot(., aes(y = `HISDEG-PWY`, x = `LACTOSECAT-PWY`)) +
  geom_point(alpha = 0.6) + 
  scale_x_log10(label = scales::percent_format(suffix = "")) +
  scale_y_log10(label = scales::percent_format(suffix = "")) +
  facet_wrap(~ Dermotype_size) +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(
    method = "spearman",
    color = "black",
    #label.x = df$min_abundance,  # Use group-specific min (original scale)
    #label.y = df$max_value,      # Use group-specific max (original scale)
    hjust = -0.1,                # Adjust horizontal position
    vjust = 1.5                  # Adjust vertical position
  ) +
  ggpubr::theme_pubr(base_size = 12, border = T) 
  

  