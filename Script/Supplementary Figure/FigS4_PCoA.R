library(dplyr)
library(tibble)
library(vegan)
library(reshape2)
library(tidyr)
library(ggpubr)

mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
fn = mat 
mdata.ex <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

##########calculate average of each subject######################
mat.subj <-  fn %>% 
  rownames_to_column("LibraryID") 

# Convert mat to long format
mat_long <- mat.subj %>% pivot_longer(cols = -LibraryID, names_to = "taxa", values_to = "abundance")

# Join mat_long and mdata.ex
joined_df <- inner_join(mat_long, mdata.ex, by = "LibraryID")

# Calculate average
averaged_df <- joined_df %>%
  group_by(SubjectID, Site, taxa) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

# Pivot the data to wide format to get the averaged relative abundance table
averaged_mat <- averaged_df %>%
  pivot_wider(names_from = taxa, values_from = abundance) %>%
  mutate(NewID = paste(SubjectID, Site, sep = "-")) %>% 
  column_to_rownames(var = "NewID") %>% 
  select(-SubjectID) %>% 
  select(-Site) %>% 
  mutate(across(everything(), as.numeric, .names = "{.col}")) %>% 
  data.frame()


BetaDiv <- vegdist(averaged_mat, method = "bray")

#perform multi-dimensional scaling on the distance matrix using this function
cmds <- cmdscale(BetaDiv, k=5, eig =TRUE, x.ret = TRUE)
#calculate the amount of variation each axis in the MDS plot accounts for, using the eigen values
eigen <- round(cmds$eig / sum(cmds$eig) * 100, 1) 
#Format the data for making MDS plots in ggplot
mds.values  <- cmds$points

mdata.ex.2 = mdata.ex %>% mutate(NewID = paste(SubjectID, Site, sep = "-")) %>%
  distinct(NewID, Site)

mds.data <- data.frame(NewID=rownames(mds.values), 
                       X=mds.values[,1], Y=mds.values[,2], Z=mds.values[,3]) %>%
  merge(mdata.ex.2, by = "NewID" ) %>% 
  column_to_rownames(var = "NewID")

################# 2D_PCoA ################# 
pal.1 <- c("#ffe119","orange", "purple", "maroon", "pink",
            "#DC0000FF", "green", "blue","cyan")
mds.data.1 = mds.data %>% mutate(Site = factor(Site, levels = c('Sc', "Ps",  "Fo", "Ch","Ub", "Ax","Ac", "Vf", "Ll" ))) %>%
  mutate(Site_fullname = case_when(Site == "Sc" ~ "(Sc) Scalp",
                                   Site == "Ps" ~ "(Ps) Parietal scalp",
                                   Site == "Fo" ~ "(Fo) Forehead",
                                   Site == "Ch" ~ "(Ch) Cheek", 
                                   Site == "Ub" ~ "(Ub) Upper back", 
                                   Site == "Vf" ~ "(Vf) Volar forearm", 
                                   Site == "Ll" ~ "(Ll) Lower leg", 
                                   Site == "Ac" ~ "(Ac) Antecubital fossa", 
                                   Site == "Ax" ~ "(Ax) Axilla", T ~ "NotApplicable"), 
         Site_fullname = factor(Site_fullname, levels = c("(Sc) Scalp",
                                                          "(Ps) Parietal scalp",
                                                          "(Fo) Forehead",
                                                       "(Ch) Cheek", 
                                                          "(Ub) Upper back", 
                                                          "(Ax) Axilla", 
                                                          "(Ac) Antecubital fossa", 
                                                          "(Vf) Volar forearm", 
                                                          "(Ll) Lower leg", 
                                                          T ~ "NotApplicable")))


site_colors <- setNames(pal.1, levels(mds.data.1$Site_fullname)[1:9])

p0 <- ggplot(mds.data.1 ,aes(x = X, y = Y)) + 
  geom_point( aes(color = Site_fullname), size = 1) +
  scale_color_manual(values=site_colors) + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
      legend.key.size = unit(2, "line"),  
      legend.position = "left",
      legend.text = element_text(size = rel(1)),  # Increase legend font size
      legend.title = element_blank(),  # Increase legend title font size
      legend.key= element_blank(), 
      legend.box.spacing = unit(0, "pt")) +  
  guides(color = guide_legend(override.aes = list(size=3)))
p_legend <- as_ggplot(get_legend(p0))

main_data <- mds.data.1 
library(ggside)
p1 = ggplot(main_data ,aes(x = X, y = Y)) + 
  geom_point( aes(color = Site_fullname), size = 1) +
  #scale_shape_manual(values=c(7, 15, 17, 8, 25, 3, 23, 20, 20))+
  stat_ellipse(aes(color = Site_fullname), level = 0.5) + 
  scale_color_manual(values=site_colors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth =  1),
        axis.text = element_text(size = rel(1)),  # Increase axis font size
        axis.title = element_text(size = rel (1.2)),
        legend.key.size = unit(2, "line"),  
        legend.position = "none",
        legend.text = element_text(size = rel(1)),  # Increase legend font size
        legend.title = element_blank(),  # Increase legend title font size
        legend.key=element_rect(fill="white"), 
        legend.box.spacing = unit(0, "pt")) + 
  labs(x = paste("PCoA 1 (", round(eigen[1], digits =2), "%)", sep = ''),
       y = paste("PCoA 2 (", round(eigen[2], digits =2), "%)", sep = '')) + 
  guides(color = guide_legend(override.aes = list(size=2)))

p1

####=========== PCs =========== #####
# mds.data.2 <- mds.data.1 %>% filter(!Site %in% c("Ax","Ac", "Vf", "Ll","Ub"))
mds.data.2 <- mds.data.1
ggplot(mds.data.2, aes(y = reorder(Site, X), x = X, fill = Site)) +
  ggridges::geom_density_ridges2(scale = 1.5, quantile_lines = T, quantiles = 0.5, alpha = 0.9) +
  scale_fill_manual(values=pal.1) + 
  theme_classic() +
  xlab("PCoA 1") + ylab("Density") +
  theme(legend.position = "none",
        axis.text = element_text(size=rel(1.2)),  
        axis.title = element_text(size=rel(1.2)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth =  1)) -> p_PC1
p_PC1

ggplot(mds.data.2, aes(y = reorder(Site, X), x = Y, fill = Site)) +
  ggridges::geom_density_ridges2(scale = 1.5, 
                                 quantile_lines = T, 
                                 quantiles = 0.5, 
                                 alpha = 0.9
  ) +
  scale_fill_manual(values=pal.1) + 
  theme_classic() +
  xlab("PCoA 2") + ylab("Density") +
  theme(legend.position = "none",
        axis.text.x = element_text(size=rel(1.2)),  
        axis.title.x = element_text(size=rel(1.2)),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth =  1)) -> p_PC2
p_PC2  

p_subs = gridExtra::arrangeGrob(p_PC1, p_PC2, ncol = 2, widths = c(1.1, 0.9))
plot(p_subs)
p_all = gridExtra::arrangeGrob(p_legend, p1, p_subs, ncol = 3, 
                               widths = c(1.4,3.5,2.8))
plot(p_all)

