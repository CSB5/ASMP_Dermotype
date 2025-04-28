######## Data Prep ###################
#load("./Metadata_dermotype_BatchABCD.RData")
load("./Data/Taxon_Decontam_AbundFilt_Renorm.RData")
colors =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", 
             "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", 
             "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", 
             "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6fabd0", 
             "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")

PickedSite =  "Ax"
SUB_mat = tdata_per_site[[PickedSite]]
#SUB_mdata.ex = mdata.ex[mdata.ex$LibraryID %in% rownames(SUB_mat) ,]

## PCoA ########################
fn = SUB_mat
library(MASS)
library(plotly)
library(plot3D)
library(vegan)
library(dplyr)
library(rgl)
library(ks)

fn = SUB_mat
data.dist <- vegdist(fn, "bray")  
SUB_cmds <- cmdscale(data.dist, k = 5, eig =TRUE, x.ret = TRUE)
SUB_eigen <- round(SUB_cmds$eig / sum(SUB_cmds$eig) * 100, 1)
#Format the data for making MDS plots in ggplot
SUB_mds.values  <- SUB_cmds$points
SUB_mds.data <- data.frame(LibraryID=rownames(SUB_mds.values), 
                           X= SUB_mds.values[,1], 
                           Y= SUB_mds.values[,2], 
                           Z = SUB_mds.values[,3], 
                           W = SUB_mds.values[,4], 
                           Q = SUB_mds.values[,5]) 


input = SUB_cmds$points[,1:2]
H = Hscv(input, deriv.order = 1)
fhat <- kde(x = input) 
# plot(fhat, display = "persp", phi = 25, theta = 50,
#      xlab = "PCoA 1", ylab = "PCoA 2", zlab = "PCoA 3", 
#      col.fun = viridis::viridis, thin = 3)
kms <- kms(x = input, 
           keep.path = T,
           merge = F
)
summary(kms)
tmp = data.frame(kms$mode)
estimate <- kde(x = input[,1:2], 
                H = Hscv(input, deriv.order = 1)[1:2, 1:2],
                eval.points = kms$mode[, 1:2]) 
tmp$z_values_modes <- estimate$estimate %>% slice_head(n=5)

########### 3D plot ########### 
breaks <- pretty(fhat$estimate, n = 20)
colors <- viridis::viridis(length(breaks)-1)
# Match density values to the color gradient
col_matrix <- matrix(colors[cut(fhat$estimate, breaks = breaks,
                                include.lowest = TRUE, labels = FALSE)],
                     nrow = length(fhat$eval.points[[1]]))
library(rgl)
bg3d(color = "white", fogtype = "none")
 material3d(specular="black")
# Plot without axes
persp3d( x = fhat$eval.points[[1]], 
         y = fhat$eval.points[[2]],
        z = as.matrix(fhat$estimate),
        col = col_matrix,
        shininess = 100,
        aspect = c(1.5, 1.5, 0.8),
        theta = 60, phi = 40,
        xlab = "", ylab = " ", zlab = " ",
        bty = "f") # Remove labels
# Add the mode points
points3d(tmp$X1-0.03, tmp$X2+0.03, tmp$z_values_modes+0.03,
         col = "red", pch = 19, size = 15)
widget <- rglwidget()
rgl.snapshot("./Output/Axilla.png", fmt = "png")


