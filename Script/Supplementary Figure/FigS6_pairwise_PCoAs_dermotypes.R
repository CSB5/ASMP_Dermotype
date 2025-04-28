######## Data Prep ###################
mat <- read.csv("./Data/Taxonomic abundance table.csv") %>% 
  column_to_rownames(var = "X")
mdata.dermotype.conf <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

load("./Data/Taxon_Decontam_AbundFilt_Renorm.RData")
colors =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", 
             "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", 
             "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", 
             "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6fabd0", 
             "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")

plot_pcoa_pairwise_with_legend <- function(PickedSite, my_colors, stability_cutoff = 0.9) {
  # Prep data
  fn <- tdata_per_site[[PickedSite]]
  data.dist <- vegan::vegdist(fn, "bray")  
  cmds <- cmdscale(data.dist, k = 5, eig = TRUE, x.ret = TRUE)
  eigen <- round(cmds$eig / sum(cmds$eig) * 100, 1)
  
  # Extract top 3 dimensions and merge metadata
  input <- as.data.frame(cmds$points[,1:3]) 
  names(input) <- c("PCoA1", "PCoA2", "PCoA3")
  input <- input %>% 
    tibble::rownames_to_column("LibraryID") %>% 
    dplyr::left_join(mdata.dermotype.conf, by = "LibraryID") %>% 
    dplyr::filter(Stability >= stability_cutoff)
  
  # Ensure factor and map colors
  input$Dermotype_size <- factor(input$Dermotype_size)
  if (is.null(names(my_colors))) {
    if (length(my_colors) < length(levels(input$Dermotype_size))) {
      stop("Not enough colors provided for the number of dermotypes.")
    }
    names(my_colors) <- levels(input$Dermotype_size)
  }
  
  # ggpairs plot
  p <- GGally::ggpairs(
    input[, c("PCoA1", "PCoA2", "PCoA3", "Dermotype_size")],
    columns = 1:3,
    mapping = aes(color = Dermotype_size, fill = Dermotype_size, alpha = 0.5),
    upper = list(continuous = "blank"),
    lower = list(continuous = "points"),
    diag = list(continuous = "densityDiag")
  )
  
  # Apply theme and color to all subplots
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol) {
      p[i, j] <- p[i, j] + 
        scale_color_manual(values = my_colors) +
        scale_fill_manual(values = my_colors) +
        theme_pubr(border = TRUE, base_size = 7)
    }
  }
  
  # Create dummy for legend
  legend_plot <- ggplot(input, aes(x = PCoA1, y = PCoA2, color = Dermotype_size, fill = Dermotype_size)) +
    geom_point() +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    theme_minimal()
  
  # Extract and insert legend
  legend <- cowplot::get_legend(legend_plot + 
                                  theme(legend.position = "right", 
                                        legend.title = element_blank(), 
                                        legend.text = element_text(size = rel(1.2))))
  p[1, 3] <- cowplot::ggdraw(legend)
  
  return(p)
}


custom_colors <- c("lightblue", "red", "orange", "blue", "purple")
plot_pcoa_pairwise_with_legend(PickedSite = "Ac", my_colors = c("green", "#728c69")) -> p_ac
plot_pcoa_pairwise_with_legend(PickedSite = "Vf", my_colors = c( "#73c2fb", "red")) -> p_vf
plot_pcoa_pairwise_with_legend(PickedSite = "Fo", my_colors = c("#FFED6F","purple"))-> p_fo
plot_pcoa_pairwise_with_legend(PickedSite = "Ub", my_colors = c("#00A087FF", "#CF4E9CFF"))-> p_ub
plot_pcoa_pairwise_with_legend(PickedSite = "Sc", my_colors = c("orange","blue"))-> p_sc
plot_pcoa_pairwise_with_legend(PickedSite = "Ps", my_colors = c("orange","blue"))-> p_ps

library(GGally)
library(gridExtra)

# Convert ggmatrix to grobs
g1 <- ggmatrix_gtable(p_ac)
g2 <- ggmatrix_gtable(p_vf)
g3 <- ggmatrix_gtable(p_fo)
g4 <- ggmatrix_gtable(p_ub)
g5 <- ggmatrix_gtable(p_sc)
g6 <- ggmatrix_gtable(p_ps)

# Arrange them
grid.arrange(g1, g2, g3, g4, g5, g6, nrow = 2)
