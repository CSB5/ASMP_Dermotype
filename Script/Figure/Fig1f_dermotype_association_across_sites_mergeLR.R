library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggplot2)

mdata.dermotype <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

######### Site-Site similarity ###################
df <- mdata.dermotype |> 
  mutate(SiteID = paste(Site)) |>
  mutate(SubjectID = paste0(SubjectID, LR)) |> 
  dplyr::select(Dermotype_size, SiteID, SubjectID) |> 
  pivot_wider(names_from = SiteID, values_from = Dermotype_size) |>
  dplyr::select(-contains(c('Ll', 'Ch'))) |>
  column_to_rownames("SubjectID")

############ Fisher test ############ 
df = data.frame(df) 
cols <- colnames(df)
combinations <- combn(cols, 2)
test_pair <- function(index) {
  col1 <- combinations[1, index]
  col2 <- combinations[2, index]
  table_data <- table(df[, col1], df[, col2])
  
  if (any(dim(table_data) == 0)) {
    return(data.frame(S1 = col1, S2 = col2, 
                      P_Value = NA, OR = NA))
  }
  
  # If SAX is detected in column names
  if (str_detect(col1, "SAX") | str_detect(col2, "SAX")) {
    test_res <- tryCatch({
      fisher.test(table_data, B = 2000, simulate.p.value = TRUE)
    }, error = function(e) NULL)
  } else {
    test_res <- tryCatch({
      fisher.test(table_data, alternative = "two.sided")
    }, error = function(e) NULL)
  }
  
  # If test_res is NULL due to an error in fisher.test
  if (is.null(test_res)) {
    p_val <- NA
    odds <- NA
  } else {
    p_val <- test_res$p.value
    odds <- ifelse(is.null(test_res$estimate), 999, test_res$estimate)
  }
  
  return(data.frame(S1 = col1, S2 = col2, 
                    p_val = p_val, OR = odds))
}
results <- do.call(rbind, lapply(1:ncol(combinations), test_pair))
results$P_Value <- p.adjust(results$p_val, method = "BH")

results_1 <- results |> 
  filter(P_Value < 0.05) |>
  select(-p_val)

results_2 <- results_1 |>
  rename(c("S2" = "S1", "S1" = "S2")) |>
  dplyr::select(S1, S2, P_Value, OR)

results_sorted <- rbind(results_1, results_2)
results_sorted$log_P_Value <- log10(results_sorted$P_Value)
unique_sites <- unique(c(results_sorted$S1, results_sorted$S2))
all_combinations <- expand.grid(S1 = unique_sites, S2 = unique_sites)
full_data <- merge(all_combinations, results_sorted, by = c("S1", "S2"), all.x = TRUE)

desired_order <- c("Ps", "Sc", "Fo", "Vf", "Ub", "Ac","Ax")

full_data$S1 <- factor(full_data$S1, levels = desired_order)
full_data$S2 <- factor(full_data$S2, levels = desired_order)

#========vis===========#
p <- 
  ggplot(full_data, aes(x = S1, y = S2)) +
  geom_tile(aes(fill = log_P_Value), color = "white") +
  geom_text(aes(label = ifelse( OR != 999, round(OR,1), "")), 
            size = 7, size.unit = "pt", na.rm = TRUE) +
  scale_fill_gradient(low = "red", high = "yellow", 
                      na.value = "transparent", 
                      name = "Dermotype cross-site correlation\n log10(Fisher's p-value)",
                      guide = "colourbar") +
  scale_y_discrete(limits = levels(full_data$S1)) + 
  ggpubr::theme_pubr(border = T, legend = "bottom", base_size = 7)+
  theme(
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.key.height =  rel(0.3),
        legend.key.width = rel(0.8))

print(p)