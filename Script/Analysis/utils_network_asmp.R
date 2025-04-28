# Script to plot networks using igraph
# date: 20231128

library(tidyverse)
library(igraph)
library(ggh4x)

color.list <- #colorblind friendly palette
  c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2",
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )

colors =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44",
              "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072",
              "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142",
              "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6FABD0",
              "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")
oral <- c('S. mitis oralis pneumoniae', 'S. salivarius', 'R. dentocariosa',
          'V. parvula', 'C. matruchotii', 'R. aeria', 'C. durum',
          'A. viscosus', 'R. mucilaginosa', 'P. melaninogenica',
          'S. sanguinis', "N. flavescens", "A. massiliensis")
corynebact <- c('C. tuberculostearicum', 'C. pseudogenitalium',
                'C. kroppenstedtii', 'C. accolens', 'C. jeikeium', 'C. durum', 'C. matruchotii')
staph <- c('S. hominis', 'S. epidermidis', 'S. haemolyticus', 'S. saprophyticus',
           'S. caprae capitis')
cuti <- c('C. acnes', 'C. avidum', 'C. granulosum')
malass <- c('M. restricta', 'M. furfur', 'M. globosa', 'M. japonica', 'M. dermatis', 'M. slooffiae')
strep <- c('S. mitis oralis pneumoniae', 'S. salivarius', 'S. sanguinis')
rothia <- c('R. dentocariosa','R. aeria','R. mucilaginosa') 
microco <- c('M. luteus')
morax <- c('M. osloensis')

# function to generate nodes and edges from association matrix
generate_nodes_edges <- function(assoMat, Cor) {
  assoMat <- assoMat %>% dplyr::select(all_of(c("Taxa1", "Taxa2", paste0(Cor), "Padj")))
  nodes <- data.frame(id = unique(c(assoMat$Taxa1, assoMat$Taxa2)),
                           label = unique(c(assoMat$Taxa1, assoMat$Taxa2))) %>% 
    unique()
  # color the same genus
  p_filt <- 0.1
  nodes <- nodes %>% mutate(color = 
                              case_when(label %in% corynebact ~ "deeppink1",
                                        label %in% staph ~ "#A6761D",
                                        label %in% cuti ~ "gold1",
                                        label %in% malass ~ "green4",
                                        label %in% strep ~ "steelblue4",
                                        label %in% rothia ~ "#CAB2D6",
                                        label %in% microco ~ "#FB8072",
                                        label %in% morax ~ "#FF0000",
                                        T ~ "grey"))
                            
  edges <- assoMat %>% 
    rename(from = 'Taxa1', to = 'Taxa2', cor = Cor) %>% 
    mutate(value = abs(cor)) %>% 
    # add weight: if significant higher weight
    mutate(weight = case_when( Padj < p_filt ~ 4,
                              TRUE ~ 1))
  return(list(nodes = nodes,
              edges = edges))
  
}

# function to generate graph and modify graph properties
modif_graph <- function(attr) {
  g = graph_from_data_frame(d = attr$edges, 
                            vertices = attr$nodes,
                            directed = F)
  
  E(g)$lty = 5
  E(g)$lty[E(g)$value < 0.5 & E(g)$Padj < p_filt] = 2
  E(g)$lty[E(g)$value >= 0.5 & E(g)$Padj < p_filt] = 1
  E(g)$lty[E(g)$value < 0.3] = 0
  
  E(g)$Weight = 10
  
  E(g)$color = "black"
  E(g)$color[E(g)$cor >= 0.3] = "#6A3D9A" #"dodgerblue2"
  E(g)$color[E(g)$cor <= - 0.3] = "#FF7F00" #E31A1C"
  # E(g)$color[E(g)$cor >= 0.3 & E(g)$Significance == 'N'] = "#7570B4"
  # E(g)$color[E(g)$cor <= - 0.3 & E(g)$Significance == 'N'] = "#FF9F90"
  # E(g)$color[E(g)$cor >= 0.3 & E(g)$Significance == 'Y'] = "#6A3D9A" #"dodgerblue2"
  # E(g)$color[E(g)$cor <= - 0.3 & E(g)$Significance == 'Y'] = "#FF7F00" #E31A1C"
  
  sub.g = delete_edges(g, which(E(g)$value < 0.3))
  return(sub.g)
}

# helper function to standardize node size, nodeSizeSpread = 4, cexNodes = 1
getsize <- function(x, nodeSizeSpread, cexNodes) {
  (x - min(x)) / (max(x) - min(x)) * nodeSizeSpread + cexNodes
}

# helper function to rescale node size based on degrees
rescale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}

# function to get node size based on relabund
set_nodesize <- function(normCounts, subgraph) {
  normCounts <- data.frame(normCounts)
  normCounts_sum <- colSums(normCounts) %>% data.frame()
  
  # Clean up rownames by replacing all non-first dots with spaces
  rownames(normCounts_sum) <- sapply(strsplit(rownames(normCounts_sum), "\\."), function(x) {
    if (length(x) > 1) {
      paste(x[1], paste(x[-1], collapse = " "), sep = ".")
    } else {
      x
    }
  })
  # Step 1: Replace the first period (.) with a space, leaving the rest unchanged
  rownames(normCounts_sum) <- gsub("^([^\\.]+\\.)", "\\1 ", rownames(normCounts_sum))
  # filter only edges from subgraph
  normCounts_sum <- normCounts_sum[rownames(normCounts_sum) %in% vertex_attr(subgraph)$name,]
  nodeSize <- getsize(normCounts_sum, nodeSizeSpread = 4, cexNodes = 8)
  return(nodeSize)
}

# https://stackoverflow.com/questions/46139709/plot-two-igraph-networks-using-the-same-coordinates-and-same-placement-in-the-pl
# 
# adjust node size
# https://stackoverflow.com/questions/12058556/adjusting-the-node-size-in-igraph-using-a-matrix
# 
# edge thickness
# https://stackoverflow.com/questions/69765986/change-edge-thickness-according-to-a-column-of-a-dataframe-containing-relations
# 
# arrange using tkplot
# https://igraph.org/r/doc/tkplot.html

# label position
# https://stackoverflow.com/questions/28542966/igraph-positioning-labels-and-removing-empty-space-in-grid-layout
# 
# #Get the coordinates of the Nodes
# set.seed(1)
# Coords <- layout_with_fr(sub.g.ac.1) %>% 
#   as_tibble() %>%
#   bind_cols(data_frame(names = names(V(sub.g.ac.1))))
# 
# #get the coordinates of the remaining Nodes
# NetCoords <- data_frame(names = names(V(sub.g.ac.2))) %>%
#   left_join(Coords, by= "names")
# coords_rescaled <- sapply(Coords[-3], function(x) -1+((x-min(x))*2)/diff(range(x)))
# rownames(coords_rescaled) <- Coords$names
