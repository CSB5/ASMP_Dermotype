library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(tidyverse)
library(svglite)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggtext)
library(ggh4x)
library(ggforce)
library(glue)

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.dermotype.conf <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")
pwy_clean <- read.csv("./Data/3_pwy/Cleaned_Pathway_Abund.csv", check.names = F) 
names(pwy_clean)[1] <- "LibraryID"
pwy_clean <- pwy_clean %>% 
  column_to_rownames(var = "LibraryID")

colors =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", 
             "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", 
             "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", 
             "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6fabd0", 
             "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")


Taxa_Shannon_median <- vegan::diversity(mat, "shannon") %>% 
  data.frame() %>%
  dplyr::rename("Shannon" = ".") %>%
  merge(mdata.dermotype.conf , by.x=0, by.y= "LibraryID") %>% 
  dplyr::group_by(Dermotype_size) %>% 
  dplyr::summarize(Median_tax_shannon = median(Shannon))

#======= calculate pathway shannon index =======#
tmp = pwy_clean |>
  rownames_to_column(var = "LibraryID") |>
  left_join(mdata.dermotype.conf [, c("Dermotype_size", "LibraryID")], join_by(LibraryID))

library(ggpubr)
df_Shannon <- vegan::diversity(pwy_clean, "shannon") %>% 
  data.frame() %>%
  dplyr::rename("Shannon" = ".")%>%
  merge(mdata.dermotype.conf , by.x=0, by.y= "LibraryID") %>% 
  left_join(Taxa_Shannon_median, join_by(Dermotype_size))

p_shannon <- list()
#============ Axilla ============##
df <- df_Shannon %>% filter(Site == "Ax")
# Conduct pairwise Wilcoxon tests
test_results <- pairwise.wilcox.test(df$Shannon, df$Dermotype_size, 
                                     p.adjust.method = "none", # Adjust this if needed
                                     paired = FALSE)
p_values_df <- as.data.frame(as.table(test_results$p.value))
# Filter comparisons based on p-value threshold
significant_comparisons <- p_values_df %>% 
  filter(Freq < 0.05) %>% 
  dplyr::select(Var1, Var2, Freq)
print(significant_comparisons)

p_shannon[["Ax"]] <- df %>%
  ggplot(aes(x= Dermotype_size, y=Shannon, fill= Dermotype_size)) + 
  stat_boxplot(geom = "errorbar", width = 0.5, aes()) + 
  geom_boxplot(outlier.size = 0.5) + 
  geom_jitter(size = 0.5, width = 0.2, height = 0, alpha = 0.5) +
  #scale_y_log10() + 
  labs(y='Shannon Index (pathway)')+ 
  scale_fill_manual(values=c('#66A61E', '#FF0000', '#5E4FA2',  '#F46D44' ,'#4288BD')) + 
  stat_compare_means(comparisons = list(#c("Ax-2", "Ax-1"),
                                        #c("Ax-3", "Ax-1"), 
    c("Ax-4", "Ax-3"),
    c("Ax-4", "Ax-5"),
    c("Ax-3", "Ax-5"),
    c("Ax-2", "Ax-5"),
    c("Ax-2", "Ax-4"),
    c("Ax-4", "Ax-1"), 
    c("Ax-5", "Ax-1")), 
  label = "p.signif", hide.ns = T,  step.increase = 0.05, vjust = 1) +
  geom_point(aes(x = Dermotype_size, y = Median_tax_shannon),
             shape = 23, size = 3, fill = "cyan", color = "black") +
  theme_classic() + 
  theme(axis.text.x = element_text(size=rel(1.2), angle = 90),
        axis.text.y = element_text(size=rel(1.2)),
        legend.position="none",
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=rel(1.2)),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = "black"))    
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
    ggplot(aes(x= Dermotype_size, y= Shannon, fill= Dermotype_size)) + 
    stat_boxplot(geom = "errorbar", width = 0.5, aes()) + 
    geom_boxplot(outlier.size = 0.5) + 
    geom_jitter(size = 0.5, width = 0.2, height = 0, alpha = 0.5) +
    #scale_y_log10() + 
    labs(y='')+ 
    scale_fill_manual(values= color_dermotype[[i]]) + 
    stat_compare_means(comparisons = list(c(unique(df$Dermotype_size))), 
                       #method = 'wilcox.test',
                       label = "p.signif", hide.ns = T) +
    geom_point(aes(x = Dermotype_size, y = Median_tax_shannon),
               shape = 23, size = 3, fill = "cyan", color = "black") +
    theme_classic() + 
    theme(axis.text.x = element_text(size=rel(1.2), angle = 90),
          axis.text.y = element_text(size=rel(1.2)),
          legend.position="none",
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=rel(1.2)),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = "black"))   
  print(p_shannon[[Pickedsite]])

}

p_shannon_all <- gridExtra::grid.arrange(grobs = p_shannon, nrow = 1, 
                                        widths = c(4,2,2,2,2,2,2))
plot(p_shannon_all)


#=========DA pathway and their subclass ========#
All_res_1 <- read.csv("~/ASMP_Dermotype/Data/3_pwy/DA_pwy_dermotypes_All_res_1.csv", row.names = NULL)
pwy_translator <- read.csv("~/ASMP_Dermotype/Data/3_pwy/pwy_translator_ID.csv", row.names = NULL) %>% 
  select(all_of(c("ID", "pwy")))
df <- All_res_1|> 
  mutate(Dermotype_size = factor(Dermotype_size, levels = c(
    "Ax-1", "Ax-2", "Ax-3", "Ax-4", "Ax-5", 
    "Ac-1", "Ac-2", 
    "Vf-1", "Vf-2", 
    "Ub-1", "Ub-2", 
    "Fo-1", "Fo-2", 
    "Ps-1", "Ps-2", 
    "Sc-1", "Sc-2"
  ))) |>
  #filter(! (str_detect(Dermotype_size, "-1") & Site != "Ax")) |>
  select(ID, Superclass1, coef, Subclass, Dermotype_size) |>
  spread(key = Dermotype_size, value = coef) |>
  left_join(pwy_translator, join_by(ID)) |>
  arrange(Superclass1, Subclass) |>
  column_to_rownames("pwy") |>
  select(-ID, -Superclass1, -Subclass) |> 
  mutate_all(~ifelse(is.na(.), 0, .)) 

max_abs_value <- max(abs(df), na.rm = TRUE)
breaks <- seq(-max_abs_value, max_abs_value, length.out = 101)
breaks <- c(breaks[1:50], 0, breaks[52:101])  # Insert 0 at the middle of breaks
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
annot_df <-   All_res_1|>
  ungroup() |>
  distinct(ID, Superclass1, Subclass) |>
  left_join(pwy_translator, join_by(ID)) |>
  arrange( Superclass1, Subclass) |>
  select(-ID) |>
  column_to_rownames("pwy")

n <- n_distinct(annot_df$Subclass)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
annot_color <- sample(col_vector, n)
names(annot_color) <- unique(annot_df$Subclass)
annot_color <- list(Subclass = annot_color)


pheatmap::pheatmap(t(df), 
                   #filename = paste0("~/Desktop/ASMP/Output/5_PWY/Heatmap_log2FC_full", ".pdf"),
                   breaks = breaks, 
                   color = colors,
                   annotation_col = annot_df,
                   annotation_colors  = annot_color,
                   cluster_cols = F, 
                   cluster_rows = F,
                   gaps_row = c(5, 7, 9, 10, 11, 13, 15), 
                   height = 12, 
                   width = 27)

# subset tpheatmap()# subset the heatmap --to emphasize one dermotype of each site 
to_plot_list <- c(    "Ax-1", "Ax-2", "Ax-3", "Ax-4", "Ax-5", 
                      "Ac-2", 
                      "Vf-2", 
                      "Ub-2",  
                      "Fo-2", 
                      "Ps-2",  
                      "Sc-2")
df1 <- df |> select(any_of(to_plot_list))
pheatmap::pheatmap(df1, 
                   #filename = paste0("~/Desktop/ASMP/Output/5_PWY/Heatmap_log2FC_subset", ".pdf"),
                   breaks = breaks, 
                   color = colors,
                   annotation_row = annot_df,
                   annotation_colors  = annot_color,
                   show_rownames = F,
                   cluster_cols = F, 
                   cluster_rows = F,
                   gaps_col  = c(5, 6), 
                   height = 7, 
                   width = 5.5)

#========selected DA pwys ==========#
log2FC_threshold <- 2
prev_threshold <- 0.95
percentile_threshold <- 0 
n_top <- 10
top_pwy <- All_res_1 %>%
  filter(Type == "conserved") |>
  filter(prev > prev_threshold) |>
  filter(!(Site != "Ax" & str_detect(Dermotype_size, "-1"))) |>
  dplyr::select(ID, pwy, 
                coef, Site, Dermotype_size,
                Superclass1, Superclass2, Superclass3, Superclass4,
                label_1, label_2, Subclass,
                prev, HighPrev, Type) %>%
  group_by(Dermotype_size) %>%
  arrange(Dermotype_size, -abs(coef)) %>%
  slice_max(n= n_top, abs(coef)) |>
  ungroup() 

tmp <- All_res_1 |>
  filter(prev > prev_threshold) %>%
  filter(!(Site != "Ax" & str_detect(Dermotype_size, "-1"))) |>
  filter(ID %in% top_pwy$ID)
n_distinct(tmp$ID)
top_pwy <- tmp

color.list <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A")
color_mapping <- data.frame(
  Superclass1 = unique(top_pwy$Superclass1),
  color = color.list[1:length(unique(top_pwy$Superclass1))]
)
top_pwy_sort <- top_pwy %>%
  ungroup() %>%
  select(Superclass1, label_1)%>%
  distinct(Superclass1, label_1) %>%
  arrange(Superclass1, label_1) 

input <- top_pwy %>%
  mutate(pwy = str_replace_all(pwy, "N-acetylglucosamine", "GlcNAc")) |>
  mutate(pwy = str_replace_all(pwy, "N-acetylmannosamine", "ManNAc")) |>
  mutate(pwy = str_replace_all(pwy, "N-acetylneuraminate", "NeuAc")) |>
  mutate(pwy = str_replace_all(pwy, "superpathway of", "SP")) |>
  mutate(pwy = str_replace_all(pwy, " \\s*\\([^\\)]+\\)", "")) %>%
  left_join(color_mapping, join_by(Superclass1)) %>%
  mutate(text_color = glue("<span style='color:{color}'>{pwy}</span>")) %>%
  mutate(Site = factor(Site, levels = c(
    "Ax", "Ac", "Vf", "Ub", "Fo", "Ps", "Sc"
  ))) %>%
  mutate(label_1 = factor(label_1, levels = top_pwy_sort$label_1)) 
n_distinct(input$ID)

heatmap <- ggplot(input, aes(x = text_color, y = Dermotype_size)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", guide = "colorbar") + 
  geom_point(data = . %>% filter(Type == "conserved"), aes(x = text_color, y = Dermotype_size), 
             shape = 21, size = 4, fill = "white", color = "black") +
  facet_nested(Site ~ Subclass,
               labeller = label_wrap_gen(width = 5, multi_line = TRUE),
               space = "free",
               scales = "free",
               switch = "y"
  ) +
  labs(x = '', y = '', fill = 'log2FC') +
  theme_pubr(base_size = 12, x.text.angle = 45, border = T) +
  theme(
    strip.text.y.left  = element_text(angle = 0),
    strip.placement = "inside",
    strip.text.x.top = element_text(angle = 0, size = rel(0.8)),
    axis.text.y = element_markdown(),  # <-- This ensures text_color renders as HTML
    axis.text.x = element_markdown(),  # Apply if text_color is on x-axis
    legend.position = "left",
    legend.direction = "vertical",
    legend.margin = margin(t = 5, b = 5)
  )  
heatmap
