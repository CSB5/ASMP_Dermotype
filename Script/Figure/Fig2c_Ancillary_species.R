library(dplyr)
library(scales)
library(tibble)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ====================== read in data ==========================
mdata.dermotype.conf <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv") %>% 
  mutate(SubjectID = str_pad(SubjectID, 5, pad = "0"))
  
load("./Data/Taxon_Decontam_AbundFilt_Renorm.RData")

load("./Data/2_basic/Outliers_LowPrev_Dermotype_percentile.RData")
d_all$Taxa = ifelse(d_all$Taxa == "Acinetobacter pittii calcoaceticus nosocomialis",
                  "Acinetobacter pittii/calcoaceticus/nosocomialis", d_all$Taxa)
d_all_Dermotype <- d_all


load("./Data/2_basic/Outliers_LowPrev_Site_percentile.RData")
d_all$Taxa = ifelse(d_all$Taxa == "Acinetobacter pittii calcoaceticus nosocomialis",
                    "Acinetobacter pittii/calcoaceticus/nosocomialis", d_all$Taxa)
d_all_Site <- d_all %>% 
  dplyr::rename(Dermotype = Site) 

d_all <- rbind(d_all_Site, d_all_Dermotype)

minimal_n <- 3
d_all <- d_all %>% 
  mutate(Site = gsub("-.*", "", Dermotype)) %>% 
  filter(n_non_zero >= minimal_n)
write.csv(d_all, paste0("./Data/2_basic/Ancillary_species_Site_Dermotype_minimal_", 
  minimal_n, ".csv"))


mat_ancillary <- mat %>% select(all_of(unique(d_all$Taxa_og)))
merged_data <- merge(mdata.dermotype.conf[, c("LibraryID", "SubjectID", "Site" )], 
                     mat_ancillary, 
                     by.x = "LibraryID", by.y = 0)

# =========== Intra-individual: the number of skin sites that harboring these species =========== #
within_subject_prev <- merged_data %>%
  dplyr::group_by(SubjectID_Site = paste(SubjectID,Site, sep = "_"), .drop = F) %>%
  select(-LibraryID) %>%
  summarise(across(contains(c("s__")), \(x) sum(x, na.rm = TRUE))) %>% 
  mutate(SubjectID = str_extract(SubjectID_Site, ".*(?=_)")) %>%
  dplyr::group_by(SubjectID) %>%
  summarise(across(contains(c("s__")), ~ sum(. > 0, na.rm = TRUE))) 

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

# ===========  the number of indiviudals harboring these virus at >= 1 skin sites =========== #  
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

# ===========  draw figures =========== # 
p1 = ggplot(df, aes(x= reorder(Taxa, V1), y = V1)) +
  geom_bar(stat = "identity", position = "identity", fill = "grey") +
  labs(x = "", y = paste0("Number of \ncarriers")) +
  scale_x_discrete(expand = expansion(add = .2)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(rows = vars(Taxa_group), space = "free", scales = "free", switch = "y")+ 
  theme_pubr(base_size = 7, margin = F, border = T, x.text.angle = 90) +
  theme(
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(),
    strip.text.y = element_blank()
    ) +
  coord_flip()+
  guides(fill = "none") 

p1


#-------- draw the within-individual prevalence --------#
# Create discrete abundance categories in your data
df <- df %>%
  mutate(RA_category = case_when(
    RA < 0.1 ~ "<10%",
    RA >= 0.1 & RA < 0.3 ~ "10-30%",
    RA >= 0.3 ~ ">=30%"
  ), 
  RA_category = factor(RA_category, levels = c(">=30%", "10-30%", "<10%")))

# Create the plot
label_vector <- setNames(df$V1, df$Taxa)
p2 <- ggplot(df, aes(x = reorder(Taxa, V1), 
                     y = prev)) +
  geom_boxplot(color = "grey", outliers = F) +
  geom_jitter(aes(color = RA_category), alpha = 0.7, size = 0.1) +
  scale_color_manual(
    values = c("<10%" = "#F2E5BF", "10-30%" = "#257180", ">=30%" = "#FF6500"),
    name = "Intra-individual\nmedian\nabundance (%)"
  ) +
  labs(x = "", y = paste0("Number of sites", "\n", "within an individual")) +
  scale_y_continuous(breaks = c(1,3,5,7,9), expand = c(0,0.5)) +
  facet_grid(rows = vars(Taxa_group), space = "free", scales = "free") + 
  theme_pubr(legend = "right", base_size = 5, border = TRUE, margin = FALSE, x.text.angle = 90) +
  coord_flip() +
  theme(
    axis.text.y.right = element_text(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 0)
  ) +
  scale_x_discrete(labels = label_vector)
p2

# ========= the non-0 median RA of short-listed bugs, regardless of sites ========= #

df3 <- d_all %>%
       filter(Taxa_og %in% df$variable)

Desired_order <- df %>% arrange(desc(V1), desc(Taxa)) %>% distinct(Taxa)
df3$Taxa <- factor(df3$Taxa, levels = Desired_order$Taxa)
df3$Site= factor(df3$Site, levels =
                c(
                  "Ax",
                  "Ac",
                  "Ll",
                  "Vf",
                  "Ub",
                  "Fo",
                  "Ch",
                  "Ps",
                  "Sc"))



df3 <- df3 %>%
  mutate(
    Abundance_category = case_when(
      median_abund < 0.1 ~ "<10%",
      median_abund >= 0.1 & median_abund < 0.5 ~ "10-50%",
      median_abund >= 0.5 ~ ">50%",
      TRUE ~ NA_character_
    ),
  Abundance_category = factor(Abundance_category, levels = c(">50%", "10-50%", "<10%")))


# Set color palette for categories
category_colors <- c(">50%" = "#C80000",
                     "10-50%" =  "#FFA07A", 
                     "<10%" = "#E0BBE4"
                     )  

# Create the heatmap with categories
p3 <- ggplot(df3, aes(x = Dermotype, y = Taxa)) + 
  geom_tile(aes(fill = Abundance_category), color = "white") + 
  scale_fill_manual(values = category_colors, 
                    na.value = "white", 
                    name = "Within-niche\nmedian\nrelative\nabundance (%)") +
  scale_y_discrete(limits = rev) +
  facet_grid(cols = vars(Site), 
             rows = vars(Taxa_group), 
             scales = "free", space = 'free',
             switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_pubr(base_size = 5, margin = F, border = T, legend = "left", x.text.angle = 90 ) +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"), 
    strip.text.x.top = element_blank(), 
    strip.text.y.left = element_blank()
  )

p3
p3_without_legend <- p3 + 
  theme(legend.position = "none") +
  xlab ("Median abundance among carriers (%)")

# =========== assemble figures =================== #
p_combine_1 = p3 + p2 + plot_layout(ncol = 2,widths = c(4, 1.3))
p_combine_1 
