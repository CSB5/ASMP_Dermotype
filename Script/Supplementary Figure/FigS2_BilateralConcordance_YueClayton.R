library(ggpubr)
library(vegan)
library(reshape2)
library(tidyverse)
library(tibble)
library(dplyr)
library(svglite)

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
fn = mat
mdata.ex <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
#========================== Yue-Clayton Theta Index ==========================#
# Yue-Clayton Theta formula: https://mothur.org/wiki/thetayc/ #
getYC <- function (fn) {
  yc.dist <- vegan::designdist(fn,
                             method="J/(A+B-J)",
                             terms = c( "quadratic"),
                             abcd = FALSE)

  m <- yc.dist %>% as.matrix() 
  yc_index <- reshape2::melt(m)[reshape2::melt(upper.tri(m))$value,]
            
  subject_YC_pairwise <- mdata.ex %>%
    subset(select = c(Site, SubjectID, LibraryID, SiteID)) %>%
    dplyr::rename(Var1_Site = Site, Var1_SubjectID = SubjectID,
                  Var1 = LibraryID, Var1_SiteID = SiteID) %>%
    merge(yc_index, by = "Var1")

  subject_YC_pairwise <- mdata.ex %>%
    subset(select = c(Site, SubjectID, LibraryID, SiteID)) %>%
    dplyr::rename(Var2_Site = Site, Var2_SubjectID = SubjectID,
                  Var2 = LibraryID, Var2_SiteID = SiteID) %>%
    merge(subject_YC_pairwise, by = "Var2")
  
  within_subject_YC_pairwise <- subject_YC_pairwise %>%
    filter(Var2_SubjectID == Var1_SubjectID) %>%
    filter(Var2_Site == Var1_Site) %>%
    mutate (type = "intra-individual")
  
  between_subject_YC_pairwise <- subject_YC_pairwise %>%
    filter(Var2_Site == Var1_Site) %>%
    filter(Var2_SubjectID != Var1_SubjectID) %>%
    mutate (type = "inter-individual")
  
  tmp2 = rbind(within_subject_YC_pairwise, between_subject_YC_pairwise)
  return(tmp2)
}

#========= get YC theta for the full profile ==========#
tmp2_w_Cacnes <- getYC(mat)

#========= get YC theta with removal of C.acnes ==========#
fn <- mat %>% select(-s__Cutibacterium_acnes)
fn <- fn[]/rowSums(fn[])
tmp2_wo_Cacnes <- getYC (fn)

#========= Merge YC thetas  ==========#
tmp2_w_Cacnes$group <- "Full profile"
# Create a temporary dataset with only "intra-individual" type
tmp_intra_individual <- tmp2_w_Cacnes %>% 
  filter(type == "intra-individual") %>% 
  distinct(Var1_Site, Var1_SubjectID, value) 
# Reorder Var1_Site based on the median value
ordered_sites <- with(tmp_intra_individual, reorder(Var1_Site, value, FUN = median, na.rm = TRUE))
# Extract the levels of the ordered Var1_Site
order <- levels(ordered_sites)

p1a <-
  # Apply the ordered levels to the original dataset
  tmp2_w_Cacnes %>%
  mutate(Var1_Site = factor(Var1_Site, levels = order)) %>% 
  ggplot(., aes(x= Var1_Site, 
                       y = value, fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c( "#7FB0FF", "#E53975")) +
  # geom_point(data = subset(tmp3, outlier == TRUE), aes(x= Var1_Site, y = value, fill = type), shape = 20, size = 0.5) + 
  labs(x = paste("Site"),
       y = paste("Yue-Clayton theta"))  +
  facet_grid(cols = vars(group), scales = "free", space = "free") + 
  stat_compare_means(
    aes(group = type),
    label.y = 1.05,  # Adjust position of comparison statistics
    label = "p.signif",  # Display p-values with forfnting
    hide.ns = TRUE,  # Hide non-significant comparisons
    method = "wilcox.test",
    method.args = list(alternative = "less"), 
    symnum.args=list(
      cutpoints = c(0, 0.001, 0.005, 0.05, 1), 
      symbols = c( "\u2217","<0.005" ,"<0.05", "ns")),
    family='mono') + 
  theme(panel.grid.major = 
          element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth =  1),
        axis.text = element_text(size = rel(1.1)),  # Increase axis font size
        axis.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = rel(1.1)),  # Increase legend font size
        legend.title = element_blank(),  # Increase legend title font size
        legend.background = element_blank(), 
        legend.position = "top", 
        legend.direction = "horizontal")

print(p1a)
#==========two-way hierarhical clustering of sites and subjecs ==========#
# Convert the data to a wide format
heatmap_data <- tmp_intra_individual %>%
  pivot_wider(names_from = Var1_Site, values_from = value)

# Convert to matrix, removing the SubjectID column
heatmap_matrix <- as.matrix(heatmap_data[,-1])

# Set rownames to SubjectID
rownames(heatmap_matrix) <- heatmap_data$Var1_SubjectID

# Perform hierarchical clustering
row_clustering <- hclust(dist(heatmap_matrix), "ward.D2")  # Clustering by rows (subjects)
col_clustering <- hclust(dist(t(heatmap_matrix)), "ward.D2") # Clustering by columns (sites)
# Install and load ComplexHeatmap
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  install.packages("ComplexHeatmap")
}
library(ComplexHeatmap)

# Plot the heatmap
Heatmap(heatmap_matrix, 
        cluster_rows = row_clustering, 
        cluster_columns = col_clustering,
        show_row_names = FALSE,  # Do not label rows
        show_column_names = TRUE)  # Label columns


#======= boxplot FigureS full profile + C.acnes removed  ================# 
getPearson <- function (fn) {
  PearsonCorr = cor(t(fn))
  subject_Pearson_pairwise <- subset(melt(as.matrix(PearsonCorr)), Var1 != Var2)
  
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
fn <- mat %>% select(-s__Cutibacterium_acnes)
fn <- fn[]/rowSums(fn[])
tmp2_wo_Cacnes <- getPearson(fn)

tmp2_wo_Cacnes$group <- "After removing C. acnes"
tmp2_wo_Cacnes <- tmp2_wo_Cacnes |>
  mutate(Var1_Site = factor(Var1_Site, levels = order))
p = ggplot(tmp2_wo_Cacnes, 
           aes(x=Var1_Site,
               y = value, fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c( "#7FB0FF", "#E53975")) +
  # geom_point(data = subset(tmp3, outlier == TRUE), aes(x= Var1_Site, y = value, fill = type), shape = 20, size = 0.5) + 
  labs(x = "Site",
       y = paste("Pearson correlation (", "\u03C1", ")", sep = ""))   +
  facet_grid(cols = vars(group), scales = "free", space = "free") + 
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
  theme(panel.grid.major = 
          element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth =  1),
        strip.text = element_text(face = "italic", size = rel(1.1)), 
        axis.text = element_text(size = rel(1.1)),  # Increase axis font size
        axis.title = element_text(size = rel(1.1)),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = rel(1.1)),  # Increase legend font size
        legend.title = element_blank(),  # Increase legend title font size
        legend.background = element_blank(), 
        legend.position = "top", 
        legend.direction = "horizontal")
print(p)

library(patchwork)
p_Set <- p1a | p
p_Set

