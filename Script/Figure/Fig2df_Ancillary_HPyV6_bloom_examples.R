library(dplyr)
library(scales)
library(tibble)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggpubr)
library(tidyverse)
library(tidyr)

d_all <- read.csv("./Data/2_basic/Ancillary_species_Site_Dermotype_minimal_3.csv")
mdata.dermotype <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv") %>% 
  mutate(SubjectID = str_pad(SubjectID, 5, pad = "0"))

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")

mat_auxillary <- mat %>% 
  select(all_of(d_all$Taxa_og)) 

merged_data <- merge(mdata.dermotype[, c("LibraryID", "SubjectID", "Site")], 
                     mat_auxillary, 
                     by.x = "LibraryID", by.y = 0)

# =========== Intra-individual: the number of skin sites that harboring these species =========== #
within_subject_prev <- merged_data %>%
  dplyr::group_by(SubjectID_Site = paste(SubjectID,Site, sep = "_"), .drop = F) %>%
  dplyr::summarise(across(contains(c("s__")), \(x) sum(x, na.rm = TRUE))) %>% 
  mutate(SubjectID = str_extract(SubjectID_Site, ".*(?=_)")) %>%
  dplyr::group_by(SubjectID) %>%
  dplyr::summarise(across(contains(c("s__")), ~ sum(. > 0, na.rm = TRUE))) 
melted_prev <- reshape2::melt(within_subject_prev, id.vars = "SubjectID") %>%
  filter(value > 0) %>%
  dplyr::rename(prev = value)

# =========== Intra-individual: the non-zero median RA of these virus =========== # 
# Function to calculate non-zero median
non_zero_median <- function(x) {
  median(x[x != 0], na.rm = T)
}
within_subject_RA <- merged_data %>%
  group_by(SubjectID) %>%
  dplyr::summarize(across(all_of(unique(melted_prev$variable)), non_zero_median, .names = "{.col}"))
melted_RA <- reshape2::melt(within_subject_RA, id.vars = "SubjectID") %>% 
  filter(!is.na(value))%>%
  dplyr::rename(RA = value)

# ===========  the number of individuals harboring these virus at >= 1 skin sites =========== #  
count_positive_obs <- within_subject_prev %>%
  dplyr::summarize(across(-c(SubjectID),  ~ sum(. > 0, na.rm = TRUE)))

# ===========  merge and subset to the interested taxa group =========== # 
melted_data <- melted_prev %>%
  left_join(melted_RA, join_by("SubjectID", "variable"))%>%
  left_join(., d_all[, c("Taxa", "Taxa_og", "Taxa_group")], join_by( variable == Taxa_og),
            relationship = "many-to-many") 

df = as.data.frame(t(count_positive_obs)) %>% 
  rownames_to_column(var = "variable") %>% 
  filter(variable %in% as.character(melted_data$variable) == T) %>%
  right_join(melted_data, join_by("variable")) %>%
  arrange(-V1, -prev, -RA) %>%
  distinct() 

#================plotting bloom examples ================#
library(svglite)
bloom_example <- df %>%
  group_by(Taxa_group) %>%
  arrange(desc(prev), desc(RA), desc(V1)) %>%
  slice_head(n = 10)

view(bloom_example)

example <- c("01274", "00708")
bloom_example <- bloom_example %>% 
  filter(SubjectID %in% example ) %>% 
  filter(Taxa_group != "Fungi")

PickedID = bloom_example$SubjectID
ID_list = mdata.dermotype$LibraryID[mdata.dermotype$SubjectID %in% PickedID]
  
input <- mat  %>%
  filter(rownames(mat) %in% ID_list) %>% 
  merge(mdata.dermotype[, c("LibraryID", "SubjectID", "Site")], 
        by.x = 0, by.y = "LibraryID") %>%
  group_by(SubjectID, Site) %>%
  dplyr::summarise(across(all_of(contains("s__")), mean)) %>%
  ungroup() %>%
  mutate(SubjectID_Site = paste(SubjectID, Site, sep = "_")) %>%
  column_to_rownames(var = "SubjectID_Site") %>%
  select(-SubjectID, -Site)
  
n_taxa <- 6
who = names(sort(colMeans(input), decreasing = TRUE))[1:n_taxa]
f = input[,names(input) %in% who]
f = f[order(colMeans(f), decreasing = T)]
f$Others = 1-rowSums(f)
f$Others = ifelse(f$Others < 0, 0, f$Others)
who = c(who, "Others")
f = data.frame(f)
f = f[order(f[,1], decreasing = T ), ]
f = data.frame(t(f), check.names = FALSE)
head(f)
f$Taxa <- row.names(f)
m = melt(f)
col.30 =   c("#4288BD", "#E7298A",  "#FCCDE5",  "#BC80BD", "#FEE08B","#CCEBC5", "#D54E4F","#444444", 
             "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", "#00008B",  "#9E0142",
             "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", "#80B1D4", "#B4DE69",
              "#6fabd0", "#66C2A5", "#FDAE61",
             "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")

m$Taxa <- str_replace_all(m$Taxa, c('s__' = '', "_" = " ")) 
m <- m %>% tidyr::separate(variable, c("SubjectID", "Site"), "_") %>%
     left_join(bloom_example[,c("SubjectID", "Taxa_group")], join_by(SubjectID)) %>% 
     mutate(Site = factor(Site, levels = c("Ax", "Ac", "Ll", "Vf", "Ub", "Fo", "Ch", "Ps", "Sc")),
            Taxa_group = factor(Taxa_group, levels = c("Eukaryotic virus", "Bacteria", "Bacteriophage")),
            SubjectID = factor(SubjectID, levels = example))

p_tmp <- ggplot(m, aes(Site, fill = Taxa)) + 
  geom_bar(aes(weight = value), position = position_fill(reverse = TRUE)) +  
  facet_grid( rows= vars(SubjectID),
              scales = "free", 
              space = "free"
  ) +
  scale_fill_manual(values = col.30[c(1:n_taxa, 30)]) + 
  ylab("Relative abundance (%)") +  
  guides(fill = guide_legend(ncol = 3, reverse=F, keyheight = 0.7, keywidth = 0.7)) +
  scale_y_continuous(labels = percent_format(suffix = ""), limits = c(0,1), expand = c(0, 0)) + 
  theme_pubr(base_size = 7, legend = "bottom") + 
  theme(axis.text.x = element_text(),
        axis.title.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank(), 
        strip.text = element_text(size = rel(1)))
p_tmp$data$Taxa = factor(p_tmp$data$Taxa, ordered = TRUE, 
                   levels = c(str_replace_all(who, c('s__' = '', '_' = ' '))))

p_tmp

###### site-site correlation in terms of harboring outlier viruses ##
p_fisher <- list()
df_gt10_carrier <- df %>% arrange(Taxa_group) %>% filter(V1 >10) 
list_microbe <- unique(df_gt10_carrier$variable)
#list_microbe <- c(unique(bloom_example$variable))
for (i in c(1:length(list_microbe))){
bloom_taxa = list_microbe[i]
tmp <- mat %>% dplyr::select(all_of(bloom_taxa)) %>%
  dplyr::rename("bloom_taxa" = bloom_taxa) %>%
  rownames_to_column(var = "LibraryID") %>%
  left_join(mdata.dermotype[, c("LibraryID", "SubjectID", "Site")], join_by("LibraryID")) %>%
  group_by(SubjectID, Site) %>%
  dplyr::summarise(bloom_taxa = mean(bloom_taxa)) %>%
  mutate(bloom_taxa = ifelse(bloom_taxa > 0, 1, 0)) %>%
  pivot_wider(names_from  = Site,
              values_from = bloom_taxa) %>%
  column_to_rownames(var = "SubjectID") %>%
  filter(rowSums(.) > 0)

pairwise_fisher <- function(col1, col2) {
  tbl <- table(col1, col2)
  
  # Check if the table has at least 2 rows and 2 columns
  if (nrow(tbl) < 2 || ncol(tbl) < 2) {
    return(NULL)  # Skip the test
  }
  test <- fisher.test(tbl)
  return(data.frame(odds_ratio = test$estimate,
                    p_value_single = test$p.value))
}

# Calculate pairwise fisher test for every combination of columns
results <- combn(names(tmp), 2, function(cols) {
  res <- pairwise_fisher(tmp[[cols[1]]], tmp[[cols[2]]])
  
  if (!is.null(res)) {  # Only add results if the test was performed
    res$col1 <- cols[1]
    res$col2 <- cols[2]
    return(res)
  } else {
    return(NULL)  # Skip if res is NULL
  }
}, simplify = FALSE)

# Remove NULL entries from results
results <- Filter(Negate(is.null), results)

# Combine results into a single data frame
df_results <- bind_rows(results)
df_results <- df_results %>% 
  mutate(p_value = p.adjust(p_value_single, method = "BH")) %>% 
  filter(p_value < 0.05 & odds_ratio != Inf) %>%
  arrange(p_value, -odds_ratio) %>% 
  select(-p_value_single)

df_results_2 <- df_results |>
  dplyr::rename(c("col2" = "col1", "col1" = "col2")) |>
  dplyr::select(col1, col2, p_value, odds_ratio)

results_sorted <- rbind(df_results, df_results_2)
unique_sites <- unique(c(results_sorted$col1, results_sorted$col2))
all_combinations <- expand.grid(col1 = unique_sites, col2 = unique_sites)
full_data <- merge(all_combinations, results_sorted, by = c("col1", "col2"), all.x = TRUE)
desired_order <- c(  
  "Ax", "Ac","Ll", "Vf", "Ub", "Fo","Ch", "Ps", "Sc")
full_data$col1 <- factor(full_data$col1, levels = desired_order)
full_data$col2 <- factor(full_data$col2, levels = desired_order)
full_data2 <- full_data[which(as.numeric(full_data$col1) < as.numeric(full_data$col2)), ]

if (nrow(full_data2) > 0) {
  p = ggplot(full_data2, aes(x = col1, y = col2)) +
    geom_tile(aes(fill = log10(p_value)), color = "white") +
    #geom_text(aes(label = round(odds_ratio,1)), na.rm = TRUE) +
    scale_fill_gradient(low = "red", high = "yellow", 
                        na.value = "white", 
                        name = "Fisher's\ntest\nlog10(p-value)",
                        limits = c(-5, -1)) +
    scale_x_discrete(limits = levels(full_data2$col1)) + 
    scale_y_discrete(limits = levels(full_data2$col2)) +
    theme_pubr(base_size = 7, border = T, x.text.angle = 0) + 
    theme(
          axis.title = element_blank(),
          legend.direction = "vertical",
          legend.position = "inside",
          legend.position.inside = c(0.8, 0.4),
          legend.box.background = element_rect(color = "black")
          ) +
    annotate(geom = "text", 
             label = str_replace_all(bloom_taxa, "s__|_", " "), 
             fontface = "italic", 
              x = -Inf, y = Inf, 
              hjust = 0, vjust = 2)
  p 
p_fisher[[bloom_taxa]] <- p 
}
print(i)
}
p_fisher[[1]]

