library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(tidyverse)
library(svglite)
library(ggpubr)

#========= refine the post-prevalence filtering relative abundance table ========= #
mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.dermotype.conf <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

colors =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", 
             "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", 
             "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", 
             "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6fabd0", 
             "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")

tmp = mat |>
  rownames_to_column(var = "LibraryID") |>
  left_join(mdata.dermotype.conf[, c("Dermotype_size", "LibraryID")], join_by(LibraryID))

#======= calculate shannon index =======#
library(ggpubr)
df_Shannon <- vegan::diversity(mat, "shannon") %>% 
  data.frame() %>%
  rownames_to_column(var = "LibraryID") %>% 
  dplyr::rename(c("Shannon"= ".")) %>%
  left_join(., mdata.dermotype.conf, by = "LibraryID") 
p_shannon <- list()
#============ Axilla ============##

df <- df_Shannon %>% filter(Site == "Ax")

# Conduct pairwise Wilcoxon tests with BH correction
test_results <- pairwise.wilcox.test(df$Shannon, df$Dermotype_size, 
                                     #, alternative = "greater"
                                     p.adjust.method = "BH")
p_values_df <- as.data.frame(as.table(test_results$p.value)) %>%
  rename(group1 = Var1, group2 = Var2, p.adj = Freq)

# Filter significant comparisons (p < 0.05)
significant_comparisons <- p_values_df %>% filter(p.adj < 0.05)

# Create a mapping of p-values to significance levels
significant_comparisons <- significant_comparisons %>%
  mutate(p.signif = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))
significant_comparisons$y.position <- seq(max(df$Shannon) * 1.1, length.out = nrow(significant_comparisons), by = 0.1)

# Plot with BH-corrected p-values
p_shannon[["Ax"]] <- df %>%
  ggplot(aes(x= Dermotype_size, y=Shannon)) + 
  #stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill= Dermotype_size)) + 
  geom_jitter(size = 0.5, width = 0.2, height = 0, alpha = 0.5) +
  labs(y='Shannon Index')+ 
  scale_fill_manual(values=c('#66A61E', '#FF0000', '#5E4FA2',  '#F46D44' ,'#4288BD')) + 
  theme_pubr(base_size = 10, x.text.angle = 45) + 
  theme(
        legend.position="none",
        axis.title.x = element_blank()) +
  stat_pvalue_manual(significant_comparisons, 
                     label = "p.signif", 
                     step.increase = 0.05)

print(p_shannon[["Ax"]])


#============ other sites  ==========#
list_site <- c("Ac", "Vf",  "Ub", "Fo", "Ps","Sc")
color_dermotype <- list(c("green", "#728c69"),
                        c("blue",  "#73c2fb"), 
                        c("purple",  "#C64B8C"), 
                        c("pink",  "#FF007F"),
                        c("#F46D44", "#00008B"),
                        c("#F46D44", "#00008B"))

for (i in seq_along(list_site)){
  Pickedsite <- list_site[i]
  df <- df_Shannon %>% filter(Site == Pickedsite)
  library(ggpubr)
  p_shannon[[Pickedsite]] <- df %>%
    ggplot(aes(x= Dermotype_size, y= Shannon)) + 
    geom_boxplot(outlier.shape = NA, aes(fill= Dermotype_size)) + 
    geom_jitter(size = 0.5, width = 0.2, height = 0, alpha = 0.5) +
    #scale_y_log10() + 
    labs(y='')+ 
    scale_fill_manual(values= color_dermotype[[i]]) + 
    stat_compare_means(comparisons = list(c(sort(unique(df$Dermotype_size)))), 
                       method = 'wilcox.test', 
                       label = "p.signif", 
                       hide.ns = T, 
                       method.args = list( #alternative = "greater",
                                           p.adjust.methods = "BH")) +
    theme_pubr(base_size = 10, x.text.angle = 45) + 
    theme(
          legend.position="none",
          axis.title.x = element_blank())   
  print(p_shannon[[Pickedsite]])

}

p_shannon_all <- gridExtra::grid.arrange(grobs = p_shannon, nrow = 1, 
                                        widths = c(4,2,2,2,2,2,2))
plot(p_shannon_all)


#==========Differential abundant taxa ======#
# ==== keep the long taxa name here ================== #
DA_taxa <- read.csv("./Data/1_inter_dermotype/DiffAbundTaxa_Dermotype_TSS.csv")
DA_taxa$ShortTaxa <- gsub("_", " ", DA_taxa$Taxa)
DA_taxa_1 <- DA_taxa %>%
  #filter(Taxa %in% filtered_data$Taxa) %>%
  mutate(
    # ShortTaxa = factor(ShortTaxa, ordered = F, 
    # levels = desired_order),
    Site = factor(Site, levels = 
                    c(
                      "Ax",
                      "Ll", 
                      "Ac",
                      "Vf", 
                      "Ub", 
                      "Fo",  
                      "Ch", 
                      "Ps",
                      "Sc")))

qval_to_asterisk <- function(qval) {
  if (is.na(qval)) {
    return("")
  } else if (qval < 0.01) {
    return("***")
  } else if (qval < 0.05) {
    return("**")
  } else if (qval < 0.1) {
    return("*")
  } else {
    return("")
  }
}

DA_taxa_1 <- DA_taxa_1 %>%
  mutate(asterisks = purrr::map_chr(qval, qval_to_asterisk))

tmp <- DA_taxa %>% 
  dplyr::group_by(ShortTaxa) %>%
  dplyr::summarise(N = n()) %>%
  arrange(N)

sorter <- tmp$ShortTaxa

p_log2fc <- ggplot(DA_taxa_1, 
                   aes(x = value, y = factor(ShortTaxa, levels = c(tmp$ShortTaxa[tmp$ShortTaxa %in% sorter == F], sorter)))) +
  geom_tile(aes(fill = coef)) + 
  geom_text(aes(label = asterisks), vjust = 0.5, hjust = 0.5, size = 2.5) + 
  theme_classic() + 
  scale_fill_gradient2(
    name = "Log2FC",
    low = "blue", 
    high = "red", 
    mid = "white", 
    midpoint = 0, 
    limits = c(min(DA_taxa_1$coef, na.rm = TRUE), 
               max(DA_taxa_1$coef, na.rm = TRUE)),
    guide = guide_colourbar(barwidth = 0.7)
  ) + labs(x = "", y = "") +
  facet_grid(
    cols = vars(Site), 
    scale = "free", space = "free") + 
  theme_classic() +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(position = "top") +
  facet_grid(cols = vars(Site), space = "free", scales = "free") +
  theme(
    axis.text.x = element_text(angle = 90, size = rel(1.3), face = "bold"),
    axis.text.y = element_text(size = rel(1.3), face = "italic"),
    axis.ticks = element_blank(),
    strip.text.x.top = element_text(size = rel(1.3)),
    strip.text.y.right = element_blank(),
    legend.position.inside = c(0.9, 0.3))
p_log2fc
