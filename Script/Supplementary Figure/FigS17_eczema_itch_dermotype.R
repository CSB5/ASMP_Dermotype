library(dplyr)
library(tidyverse)
library(tibble)
library(readxl)
library(lmerTest)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(broom)

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.ex <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
Library_to_keep <- rownames(mat)

#---questionnaire data is accessible upon approval from HELIOS cohort, hence not included here----#
tmp <- read.csv("./SkinHealth_Questionnaire_perSample.csv", stringsAsFactors = FALSE)
tmp <- tmp %>%
  mutate(SubjectID = str_pad(SubjectID, 5, pad = "0")) %>% 
  filter(LibraryID %in% Library_to_keep) %>% 
  rowwise() %>% 
  mutate(Irritation_severity = 
           mean(c_across(starts_with("Irritation") & ends_with("1_to_10") ), na.rm = TRUE), 
         
    Irritation_non_itch_average_1_to_10 = 
           mean(c_across(!contains("Itch") & starts_with("Irritation") & ends_with("1_to_10") ), na.rm = TRUE), 
         Itchy_average_1_to_10 = 
           mean(c_across(contains("Itch") & ends_with("1_to_10")), na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(eczema_hx_AND_current_condition = 
           case_when(Currently_have_a_skin_condition == "No" & 
                       Ever_have_had_eczema == "No" ~ "No", 
                     T ~ "Yes")) 



j <- "Ever_have_had_eczema" 
#==========severity measurement ========#
site_data <-  tmp %>%
  filter(Site == "Ac") %>%
  select( -LibraryID)%>%
  filter(Stability >= 0.9) %>% 
  select(-Stability) %>%
  group_by(SubjectID, Site) %>%
  mutate(num_obs = n()) %>%
  filter( num_obs == 1 |
            (num_obs > 1 & duplicated(paste(SubjectID, Site, Dermotype_size)))
  ) %>%
  ungroup() %>% 
  group_by(Ever_have_had_eczema) %>%
  mutate(Ever_have_had_eczema_n = n()) %>% 
  group_by(Ever_have_had_eczema, Dermotype_size) %>% 
  mutate(Dermotype_size_n = n()) %>% 
  ungroup() %>% 
  mutate(facet_label = case_when(Ever_have_had_eczema == "No" ~ paste0("No eczema history\n(n=", Ever_have_had_eczema_n, ")"),
                                 Ever_have_had_eczema == "Yes" ~ paste0("Have had eczema\n(n=", Ever_have_had_eczema_n, ")"), 
                                 T ~ "Other")) %>% 
  mutate(Dermotype_size_label = paste0(Dermotype_size, "\n(n=", Dermotype_size_n, ")"))

p_continous_1 <- site_data %>% 
  ggplot(., aes(x = Dermotype_size, y = Irritation_severity,fill = Dermotype_size))+ 
    geom_boxplot(outliers = F) +
    geom_jitter(aes(alpha = 0.7)) +
    scale_fill_manual(values = c("blue", "red"))+
    facet_grid(~ facet_label) +
    ggpubr::geom_pwc(method = "wilcox.test",
                     method.args = list(alternative = "less"), 
                     label = "p.adj.signif") +
    theme_pubr(legend = "none") +
    labs(x = "", y = "Irritation severity score (1 to 10)") 


p_continous_2 <- site_data %>%
  ggplot(., aes(x = Dermotype_size, y = Itchy_average_1_to_10,fill = Dermotype_size))+ 
  geom_boxplot(outliers = F) +
  geom_jitter(aes(alpha = 0.7)) +
  scale_fill_manual(values = c("blue", "red"))+
  facet_grid(~ facet_label) +
  ggpubr::geom_pwc(method = "wilcox.test",
                   method.args = list(alternative = "less"), 
                   label = "p.adj.format") +
  theme_pubr(legend = "none") +
  labs(x = "", y = "Itch severity score (1 to 10)") 


p_set_3 <- gridExtra::grid.arrange(p_continous_1, p_continous_2, nrow = 1)
p_set_3



