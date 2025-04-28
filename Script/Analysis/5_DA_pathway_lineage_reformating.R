
setwd("~/ASMP_network_pwy/")
# source('./script/utils_pathway.R')
library(tidyverse)
library(ggh4x)
library(ggpubr)
library(dplyr)

color.list <- #colorblind friendly palette
  c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    #"black", 
    "gold1",
    "skyblue2",
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    #"gray70", 
    "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )


source("~/ASMP_network_pwy/script/1_pwy_translator_ID_to_pathwayname.R")

#------------------- pathway lineage information -------------------#
pwy_class <- read_tsv("./Data/map_metacyc-pwy_lineage.tsv", col_names = c("pwd", "class")) %>%
  group_by(pwd) %>%
  mutate(count = n()) %>%
  filter(!(count > 1 & str_detect(class, 
  "Super-Pathways|Metabolic-Clusters|Photosynthesis|Electron|Alcohol|Deoxyribonucleotide-Biosynthesis|Other-|Cell-Structure|GALACTOSE|SECONDARY|CO2-Fixation|Polymer-Degradation"))) %>%
  select(-count) %>%
  ungroup() %>%
  group_by(pwd) %>%
  mutate(count = n()) %>%
  arrange(-count, pwd) %>%
  select(-count)
#------------------- Load metadata & community prevalence data -------------------#
load("~/Desktop/ASMP/Data/Metadata_dermotype_BatchABCD_with_StabilityScore.RData")
load("~/Desktop/ASMP/Data/pathway_community_prevalence_per_dermotype.RData")
load("~/Desktop/ASMP/Data/pathway_community_ave_abund_per_dermotype.RData")
ggplot(pwy_ave_abund, aes(x = Dermotype_size, y = ave_abund)) +
  geom_boxplot()

pwy_ave_abund_ranked <- pwy_ave_abund %>%
  group_by(Dermotype_size) %>%
  mutate(rank = rank(-ave_abund),  # Rank variables by ave_abund from high to low
         percentile = percent_rank(ave_abund)) 

#------------------- Load maaslin2 results -------------------#
All_res <- read.csv("~/Desktop/ASMP/Output/3_ML/DiffAbundPWY_Dermotype.csv")
#----------------reformat pwy class, add scripts to highlight skin related pathways:
#------ for example: moisturizing factors, vitamin, Heme -------# 
# Function to convert a string to lowercase but capitalize the first letter
library(stringr)
convert_to_proper_case <- function(string) {
  string <- str_to_lower(string) %>%
    str_to_title() %>%
    str_replace("Biosynthesis", "Syn") %>%
    str_replace("Biosyn", "Syn") %>%
    str_replace("Degradation", "Deg") %>%
    str_replace("Other-Energy", "Other Energy-Metabolism") %>%
    str_replace("-And-", "&") %>%
    str_replace(" And", "&") %>%
    str_replace("Pyr-Nuc", "Pyrimidine") %>%
    str_replace("Pur-Nuc", "Purine")
  return(string)
}
Reformat_pwy_class <- function(site_res){
  site_res <- site_res %>%
    rename(., c("ID" = "feature", 
                "Dermotype_size" = "value")) %>%
    mutate(Site = gsub("-.*", "", Dermotype_size)) %>%
    left_join(., pwy_class, c("ID" = "pwd"), relationship = "many-to-many") %>%
    filter(!str_detect(class, "Photosynthesis")) %>%
    # mutate(sim = stringdist( Superclass2, class)) %>%
    group_by(Dermotype_size, ID) %>%
    mutate(count = n()) %>%
    arrange(-count, Dermotype_size, ID) %>%
    #-------handle pathways with more than 1 class in the metaCyc database------#
    mutate(random_index = ifelse(count > 1, sample(row_number(), 1), 1)) %>%
    slice_head(n = 1) %>%
    select(!c(count, random_index)) %>%
    mutate(
      Superclass1 = sapply(strsplit(class, "\\|"), function(x) x[1]), 
      Superclass2 = sapply(strsplit(class, "\\|"), function(x) x[2]), 
      Superclass3 = sapply(strsplit(class, "\\|"), function(x) x[3]),
      Superclass4 = sapply(strsplit(class, "\\|"), function(x) x[4])) %>%
    #------create labels for plotting --------#
    mutate(label_1 = case_when(
      str_detect(Superclass2, "Cofactor") & str_detect(Superclass3, "Vitamin") ~ Superclass3,
      str_detect(Superclass2, "Cofactor") & !str_detect(Superclass3, "Vitamin") ~ Superclass2,
      str_detect(Superclass2, "Amino-Acid")  ~ Superclass2,
      is.na(Superclass2) == T & str_detect(Superclass1, "Superpathway") ~ Superclass1,
      is.na(Superclass2) == T & ! str_detect(Superclass1, "Superpathway") ~ paste("Other", Superclass1), 
      T ~ Superclass2)) %>%
    
    mutate(label_2 = case_when(
      str_detect(Superclass2, "Cofactor") & str_detect(Superclass3, "Vitamin") ~ Superclass4,
      str_detect(Superclass2, "Cofactor") & !str_detect(Superclass3, "Vitamin") ~ Superclass3,
      str_detect(Superclass2, "Amino-Acid") & str_detect(Superclass3, "IND-AMINO-ACID") ~ Superclass4,
      # str_detect(Superclass2, "Amino-Acid") & !str_detect(Superclass3, " IND-AMINO-ACID") ~ paste("Other", Superclass2),
      str_detect(Superclass2, "Noncarbon") ~ Superclass3,
      ID == "PWY0-781" ~ "aspartate superpathway", 
      ID == "SULFATE-CYS-PWY" ~ "sulfate assimilation", 
      ID == "MET-SAM-PWY" ~ "S-adenosyl-L-methionine biosynthesis",
      ID == "PRPP-PWY" ~ "histidine, purine, and pyrimidine biosynthesis",
      is.na(Superclass3) == T ~ " ", 
      T ~ Superclass3)) %>%
    mutate(label_1 = convert_to_proper_case(label_1),
           label_2 = convert_to_proper_case(label_2)) %>%
    left_join(pwy_translator, join_by(ID))
  return(site_res)
}

All_res_1 <- Reformat_pwy_class(All_res)
All_res_1 <- All_res_1 %>%
  left_join(pwy_prevalence, join_by(Site, Dermotype_size, ID == variable))  %>%
  left_join(pwy_ave_abund_ranked, join_by(Dermotype_size, ID == variable))  %>%
  mutate(Site = factor(Site, levels = c(
    "Ax", "Ac", "Vf", "Ub", "Fo", "Ps", "Sc"
  ))) 

