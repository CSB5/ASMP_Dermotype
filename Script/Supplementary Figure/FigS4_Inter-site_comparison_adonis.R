#==========pairwise adonis -> inter-site differences ========#
library(dplyr)
library(vegan)
library(pairwiseAdonis)
library(tidyverse)

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
fn = mat 
mdata.ex <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

mdata.ex <- mdata.ex %>% data.frame()
rownames(mdata.ex) = mdata.ex$LibraryID

mdata.ex_sorted <- mdata.ex[match(rownames(fn), mdata.ex$LibraryID), ]
identical(mdata.ex_sorted$LibraryID, rownames(fn))

BetaDiv <- vegdist(fn, method = "bray")
set.seed(12345)
result <- pairwise.adonis2(fn ~ Site, 
                           data = mdata.ex_sorted, 
                           method = "bray")

df <- data.frame(matrix(NA, nrow = 36, ncol = 3))
names(df)= c("comparison", "R2", "P")
for (i in seq_along(result)[-1]){
  df$comparison[i-1] <- names(result)[i]
  df$R2[i-1] <- result[[i]]$R2[1]
  df$P[i-1] <- result[[i]]$`Pr(>F)`[1]
}

df1 <- df %>% separate(comparison, c("S1", "S2"), "_vs_") 
df2 <- df1 %>% dplyr::rename("S2" = "S1", "S1" = "S2")
df12 <- rbind(df1, df2)
# Unique values
unique_sites <- unique(c(df1$S1, df1$S2))
# Create a dataframe of all combinations
all_combinations <- expand.grid(S1 = unique_sites, S2 = unique_sites)
# Left join with results_sorted
full_data <- all_combinations %>% 
  left_join(df12, by = c("S1", "S2")) %>%
  mutate(R2 = case_when(S1 == S2 ~ 0,
                        T ~ R2),
         # p = case_when(S1 == S2 ~ 1,
         #               T ~ P),
         p_sig = case_when(P <= 0.001 ~ "***",
                           P > 0.001 & P <= 0.01 ~ "**",
                           P > 0.01 & P <= 0.05 ~ "*",
                           T ~ ""))

full_data <- full_data %>% 
  mutate(padj = p.adjust(p, method = "BH")) %>% 
  mutate(p_sig = case_when(padj <= 0.001 ~ "***",
                           padj > 0.001 & padj <= 0.01 ~ "**",
                           padj > 0.01 & padj <= 0.05 ~ "*",
                    T ~ ""))
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist(cormat)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#====== get R2 value =============#
full_data_2 <- full_data %>% 
  select(-P, -p_sig) %>%
  rename(var1 = S1, var2 = S2, cor = R2) %>%
  rstatix::cor_spread(value = "cor") %>%
  column_to_rownames(var = "rowname")
#full_data_3 <- reorder_cormat(full_data_2)
full_data_2[full_data_2 == 0] <- NA
#====== get R2 value =============#
full_data_3 <- full_data %>% 
  mutate(label = case_when(
    R2 == 0 ~  NA,
    R2 <= 0.01 ~  paste(p_sig, "\n", round(R2, 3)),
    T ~  paste(p_sig, "\n", round(R2, 2)))) %>%
  select(-R2, -P, -p_sig) %>%
  rename(var1 = S1, var2 = S2, cor = label) %>%
  rstatix::cor_spread(value = "cor") %>%
  column_to_rownames(var = "rowname")

#==================== use pheatmap ==========================#
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#

library(pheatmap)
library(svglite)
 ComplexHeatmap::pheatmap(as.matrix(full_data_2),
                         display_numbers = as.matrix(full_data_3), 
         clustering_distance_cols = function(m) as.dist(m), 
         clustering_distance_rows = function(m) as.dist(m), 
         color = colorRampPalette(c('white','red'))(100), 
         fontsize_number = 8,
         number_format = "%.3f",
         number_color = "black",
         show_row_den = F, 
         column_title_rot = 0,
         legend = T,
         angle_col = "0",

         heatmap_legend_param = list(title = "Pairwise \nadonis R2", at = c(0,0.1, 0.2,0.3,0.4)))


