library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
Stackplot.generic.pwy <- function(input, # stratified by species 
                                  input_pwy_total, # community-level, not stratified by species 
                                  n_taxa, metadata_input, 
                                  V = "variable", unique_ID, 
                                  group, 
                                  title_add) {

  # input = input / as.matrix(input_pwy_total)
  input = input %>% data.frame() %>%  mutate_all(~ifelse(is.na(.), 0, .))
  n_taxa = min(ncol(input), n_taxa)
  # t_input = t(input) |> data.frame(check.names = F)
  # t_input$mean = apply(t_input, 1, function(x) mean(x))
  # n_mean_gt0 = nrow(t_input[t_input$mean > 0, ])
  # n_taxa = min(n_mean_gt0, n_taxa)
  
  #Creating bar plots for visualizations
  who = names(sort(apply(input, 2, function(x) mean(x)), decreasing = TRUE))[1:n_taxa]
  # who = names(sort(colMeans(input), decreasing = TRUE))[1:n_taxa]
  f = input[,names(input) %in% who]
  f$Other = rowSums(input[, !(names(input) %in% who)])
  #Sort f by most abundant taxa
  f = f[order(apply(f, 2, function(x) mean(x)), decreasing = T)]
  # adding in other taxa as a column for display
  f$Total = input_pwy_total[match(rownames(f), rownames(input_pwy_total)),1]
  f$Unstratified = f$Total - rowSums(f[, names(f) != "Total"])
  f$Unstratified = ifelse(f$Unstratified < 0, 0, f$Unstratified)
  f$Total <- NULL
  who = c(who,  "Other", "Unstratified")
  # now order samples (rows) by the most abundant taxon/Sort the x-axis by most abundant
  f = f[order(f[,1], decreasing = T ), ]
  #transpose f to add in the sample names for the melt
  f = data.frame(t(f), check.names = FALSE)
  # head(f)
  # add Taxon column prior to melting
  f$Taxa <- row.names(f)
  # melt
  m = melt(f)
  #make the colors for 29 taxa + 1 "Other" (30 total), making the "Other" gray.
  col.30 =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", "#E7298A", "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", "#80B1D4", "#B4DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D54E4F", "#FDAE61", "#FEE08B", "#444444", "#66C2A5", "#6fabd0", "#1B9E77", "#D95F02", "#7570B4", "#E6AB02", "#A6761D", "#CCCCCC")
  
  col = scale_fill_manual(values = c(col.30[c(1:n_taxa, 30)], "grey20"))
  
  metadata_input = data.frame(metadata_input)

  m[[group]] = metadata_input[ match(m$variable, metadata_input[[unique_ID]]), group]
  p = ggplot(m, aes(m[[V]], y = value)) + 
    geom_bar(aes(fill=Taxa), stat = "identity", position = position_stack(reverse = TRUE)) +  
    facet_grid(.~ m[[group]], scales = "free", space ="free") +
    col + xlab("Sample") + ylab("") + 
    theme_pubr(legend = "left") + 
    theme(axis.text.x = element_blank(),  
          #axis.title = element_text(size = rel(1.1)),
          legend.text = element_text(face = "italic", size = rel(0.9))) + 
    guides(fill = guide_legend(ncol = 1, 
                               title = paste("Top", n_taxa, "contributing taxa"),
                               reverse=TRUE, 
                               keyheight = 0.8)) +
    scale_y_continuous(
      #trans = ,
      #labels = percent_format(), 
                      expand=c(0,0) ) +
    labs(title=paste(title_add)) + 
    #theme(plot.title = element_text(size = rel(1.2))) 
  p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = who)
  plot(p)
  return(p)
  }
