library(dplyr)
library(tidyverse)
library(readr)

df <- read.csv("./Data/3_assoc_diff/Hubs_using_weighted_degree.csv")
View(df)

oral <- "oralis|parvula|matruchotii|dentocariosa|sanguinis|musila|massilensis"
skin_commensal <- "oslensis|tuberculo|pseudogeni|hominis|acnes"

hubs <- df %>% 
  mutate(Site = gsub("[-0-9]*", "", Dermotype)) %>% 
  mutate(Site = factor(Site, levels = c("Ax", "Ac", "Vf", "Ub", "Fo", "Sc"))) %>% 
  mutate(taxa_group = 
           case_when(str_detect(Hub, "oralis|parvula|matruchotii|dentocariosa|sanguinis|musila|massilensis") ~ "2_oral",
                     str_detect(Hub, "osloensis|tuberculo|pseudogeni|hominis|epiderm|acnes") ~ "3_skin",
                     T ~ "1_other")) %>% 
  group_by(Hub) %>% 
  mutate(n_hub = n()) %>% 
  ungroup() %>% 
  arrange(taxa_group, n_hub, desc(Site))
order <- hubs %>% distinct(Hub) %>% unlist()

p_hub_species <- hubs %>% mutate(Hub = factor(Hub, levels = order)) %>% 
  ggplot(., aes(x = Dermotype, y = Hub)) +
  geom_tile(fill = "#007191", color = "white") +
  facet_wrap(~Site, scales   = "free_x", nrow = 1) +
  ggpubr::theme_pubr(base_size = 12, x.text.angle = 45, border = T) +
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(face = "italic"))
p_hub_species 