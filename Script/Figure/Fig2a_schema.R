library(dplyr)
library(tibble)
library(ggplot2)
library(svglite)
library(matrixStats)
library(scales)
library(ggExtra)
library(ggpubr)

load("./Data/Taxon_Decontam_AbundFilt_Renorm.RData")

#========= threshold to use ========#
HighPrev_th <- 0.5
LowPrev_th <- 0.1
trendline_th <- 0.05

#========= fit the prevalence-abundance trajectory ======#
PickedGroup <- "Ll"
SUB_mat = tdata_per_site[[PickedGroup]]
t.mat = data.frame(t(SUB_mat)) 
nsample = ncol(t.mat)
d2 = t.mat %>% 
  mutate (variance = rowVars(as.matrix(.)),
          prev = apply(t.mat[,1:nsample], 1, function(x) {sum(x > 0)/nsample}), 
          median_abund = apply(t.mat[, 1:nsample], 1, function(x) median(x[x > 0],na.rm = T)),
          mean_abund =  apply(t.mat[, 1:nsample], 1, function(x) mean(x[x > 0]))) %>% 
  rownames_to_column(var = "Taxa") %>% 
  dplyr::select(Taxa, prev, median_abund, mean_abund, variance) %>% 
  filter(prev > 0) %>%
  mutate(Taxa_group = case_when(grepl("virus", Taxa) == T ~ "Eukaryotic virus",
                                grepl("phage", Taxa) == T ~ "Bacteriophage",
                                grepl("Malassezia", Taxa) == T ~ "Fungi",
                                grepl("Aspergilla", Taxa) == T ~ "Fungi",
                                grepl("Candida_",Taxa) == T ~ "Fungi",
                                T ~ "Bacteria")) 

# Calculating the 95th percentile of median_abund
percentile_95 <- quantile(d2$median_abund, 0.95)

# Linear regression with prev >= 0.05
reg_data <- subset(d2, prev >= 0.05)
fit <- lm(log(median_abund) ~ log(prev), data = reg_data)

# Prediction for the entire dataset
d2$fit <- exp(predict(fit, newdata = data.frame(prev = d2$prev)))

# 2 standard deviations from the regression line
std_dev <- sd(resid(fit))
d2$lower_bound <- exp(log(d2$fit) - 2 * std_dev)
d2$upper_bound <- exp(log(d2$fit) + 2 * std_dev)

d2$group <- ifelse(d2$prev >= 0.5, "core species",
                   ifelse(d2$prev < 0.1 & d2$median_abund > percentile_95, "ancillary species",
                          "others"))
# Create the base plot
p1 <- 
  ggplot(d2, aes(x = prev, y = median_abund)) +
  geom_jitter(color = "grey30", alpha = 0.7, size = 0.5) +
  #scale_color_brewer(palette = "Dark2") +
  scale_x_log10(labels = percent_format(suffix = ""),
                breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  scale_y_log10(labels = percent_format(suffix = "")) +
  theme_pubr(base_size = 7, border = T, margin = F) +
  geom_line(data = subset(d2, prev >= trendline_th), aes(y = fit), color = '#6666FF') +
  geom_line(data = subset(d2, prev < trendline_th), aes(y = fit), color = '#6666FF', linetype = "dashed") +
  annotate(geom = "rect",
           xmin = HighPrev_th, xmax = 1, 
           ymin = 0, ymax = 1,
           colour = "white", fill = "#2C77BF", alpha = 0.4) +
  labs(x ="Species prevalence (%)", y = "Non-zero median abundance (%)") +
  annotate("text", x = (0.5 +1)/2, y = 0.3, label = "core \nspecies", color = "#2C77BF", size = 3)
ggsave("/Users/chengchenli/Desktop/ASMP/Output/4_common/Fig2b_schema_core_only.svg",
       p1, 
       width = 4, height = 4,
       bg = "transparent")


p2 <- p1 +
  annotate(geom = "rect",
           xmin = 3/400, xmax = LowPrev_th, 
           ymin = percentile_95, ymax = 1,
           colour = "white", fill = "#D95F02", alpha = 0.4)+
  annotate("text", x = (3/400 + 0.1)/4, y = 0.3 , label = "ancillary \nspecies", color = "#D95F02", size = 3) +

  geom_vline(xintercept = LowPrev_th, linetype = "dotted", color = "#D95F02")  +

  geom_vline(xintercept = 3/400, linetype = "dotted", color = "#D95F02") +
  annotate("text", x = 3/400, y = 0.001, label = "3 samples", color = "#D95F02", size = 3) +
  
  geom_hline(yintercept = percentile_95, linetype = "dotted", color = "#D95F02") +
  annotate("text", x = 2/400, y = percentile_95 + 0.01, label = "95th percentile", color = "#D95F02", size = 3)  +
  theme_pubr(base_size = 7, border = T, margin = F)

p_schema <- ggExtra::ggMarginal(p2, type = "density",
                    fill = "#7570B3",
                    color = "#7570B3",
                    alpha = 0.5,
                    size = 5)


p_schema
