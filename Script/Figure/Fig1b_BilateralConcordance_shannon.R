library(ggpubr)
library(vegan)
library(reshape2)
library(tidyverse)
library(tibble)
library(dplyr)
library(svglite)

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.ex <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

#========================== Pearson Corr ==========================#
getPearson <- function (fn) {
    PearsonCorr = cor(t(fn)) %>% as.matrix()
    subject_Pearson_pairwise <- subset(reshape2::melt(PearsonCorr), Var1 != Var2)
    
    subject_Pearson_pairwise <- mdata.ex %>%
      subset(select = c(Site, SubjectID,  LibraryID, SiteID)) %>%
      dplyr::rename(Var1_Site = Site, Var1_SubjectID = SubjectID,
                    Var1 = LibraryID, Var1_SiteID = SiteID) %>%
      merge(subject_Pearson_pairwise, by = "Var1")
    
    subject_Pearson_pairwise <- mdata.ex %>%
      subset(select = c(Site, SubjectID,  LibraryID, SiteID)) %>%
      dplyr::rename(Var2_Site = Site, Var2_SubjectID = SubjectID,
                    Var2 = LibraryID, Var2_SiteID = SiteID) %>%
      merge(subject_Pearson_pairwise, by = "Var2")
    
    within_subject_Pearson_pairwise <- subject_Pearson_pairwise %>%
      filter(Var2_SubjectID == Var1_SubjectID) %>%
      filter(Var2_Site == Var1_Site) %>%
      mutate (type = "intra-individual")
    
    between_subject_Pearson_pairwise <- subject_Pearson_pairwise %>%
      filter(Var2_Site == Var1_Site) %>%
      filter(Var2_SubjectID != Var1_SubjectID) %>%
      mutate (type = "inter-individual")
  
    tmp = rbind(within_subject_Pearson_pairwise, between_subject_Pearson_pairwise)
    return(tmp)
}

tmp2 <- getPearson(mat)
save(tmp2, file = "./Data/2_basic/Pearson.RData")

#========== boxplot Figure 1b full profile ==========================# 
load("./Data/2_basic/Pearson.RData")
# Calculate the upper and lower fences for each group
boxplot_stats <- tmp2 %>% 
  group_by(Var1_Site, type) %>% 
  dplyr::summarize(lower_fence = quantile(value, 0.25) - 1.5 * IQR(value),
            upper_fence = quantile(value, 0.75) + 1.5 * IQR(value))

# Add outlier informatio the original data frame
tmp3 <- left_join(tmp2, boxplot_stats, by = c("Var1_Site", "type"))
tmp3$outlier <- ifelse(tmp3$value < tmp3$lower_fence | tmp3$value > tmp3$upper_fence, TRUE, FALSE)
tmp3$outlier[tmp3$type == "inter-individual"] = FALSE


# Create a temporary dataset with only "intra-individual" type
tmp_intra_individual <- tmp2%>% 
  filter(type == "intra-individual") %>% 
  distinct(Var1_Site, Var1_SubjectID, value) 
# Reorder Var1_Site based on the median value
ordered_sites <- with(tmp_intra_individual, reorder(Var1_Site, value, FUN = median, na.rm = TRUE))
# Extract the levels of the ordered Var1_Site
order <- levels(ordered_sites)

p1a = 
  tmp2 %>%   
  mutate(Var1_Site = factor(Var1_Site, levels = order)) %>% 
  ggplot(., aes(x= Var1_Site, 
                y = value, fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c( "#7FB0FF", "#E53975")) +
  # geom_point(data = subset(tmp3, outlier == TRUE), aes(x= Var1_Site, y = value, fill = type), shape = 20, size = 0.5) + 
  labs(x = paste("Site"),
       y = paste("Pearson correlation (", "\u03C1", ")", sep = ""))  +
  stat_compare_means(
    aes(group = type),
    label.y = 1.05,  # Adjust position of comparison statistics
    label = "p.signif",  # Display p-values with formating
    hide.ns = TRUE,  # Hide non-significant comparisons
    method = "wilcox.test",
    method.args = list(alternative = "less"), 
    symnum.args=list(
      cutpoints = c(0, 0.001, 0.005, 0.05, 1), 
      symbols = c( "\u2217","<0.005" ,"<0.05", "ns")),
    family='mono') + 
  ggpubr::theme_pubr( border = T, legend = "top", margin = F, base_size = 7) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm")
    )

print(p1a)


############## calculate shanon index ################
pal.1 <- c("#ffe119","orange",  "purple","maroon", "pink", 
           "blue","cyan", "green", "#DC0000FF")
site_colors <- setNames(pal.1, c('Sc', "Ps", "Fo", "Ch", "Ub", "Vf", "Ll", "Ac", "Ax"))

shannon = data.frame(diversity(fn, index = "shannon")) 
names(shannon) = "shannon"
alpha_meta = shannon %>% 
  rownames_to_column(var = "LibraryID") %>%
  left_join(mdata.ex, join_by("LibraryID")) %>% 
  mutate(Site = factor(Site, levels=order)) 

p1a.1 = ggplot(alpha_meta, 
               aes(x = Site, 
                   y = shannon)) + 
  geom_boxplot(aes(fill = Site), outlier.size = 0.5) + 
  scale_fill_manual(values =  site_colors) + 
  labs (y = "Shannon Index") +
  ggpubr::theme_pubr(border = T, legend = "none", margin = F, base_size = 7) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p1a_all = p1a + patchwork::inset_element(p1a.1, left = 0.5, right = 0.97, bottom = 0.01, top = 0.51)
p1a_all 
