library(reshape2)
library(dplyr)
library(tidyverse)
library(readr)
pwy_contri_by_group <- function(i, # skin site of interest
                                j, # pathway ID of interest 
                                mdata.dermotype.conf, #mdata files 
                                pwy_ctr_sub, # NON-TSS pathway relab file contains both a) total community-level, & b) species-stratified
                                pwy_clean, #TSS normalized pathway relab file just contain total community-level
                                phenotype,
                                pwy_asso_phenotype,
                                n_taxa){

  output <- list()
  library_list <- mdata.dermotype.conf$LibraryID[mdata.dermotype.conf$Site == i]
  pwy_list <- j
  top_pwy_total_ctr <- pwy_ctr_sub %>%
    select(any_of(c("ID", "pwy_taxa",library_list))) %>%
    filter(ID %in% j) %>%
    mutate_all(~ifelse(is.na(.), 0, .)) %>%
    mutate(pwy_taxa = str_replace(pwy_taxa, "Enhydrobacter_aerosaccus", "Moraxella_oslensis"),
           pwy_taxa= str_replace(pwy_taxa, "Propionibacterium_acnes", "Cutibacterium_acnes"))
  
  #----1. to plot relab difference -------#
  tss_relab <- pwy_clean %>% 
    select(any_of(j)) %>% 
    rownames_to_column("LibraryID") %>% 
    filter(LibraryID %in% library_list) %>% 
    left_join(phenotype, join_by(LibraryID)) 
  output[["tss_relab"]] <- tss_relab
  
  pathway_fullname <- paste(j,  pwy_ctr_sub$pwy_taxa[pwy_ctr_sub$ID == j][1], sep = "\n")
  title_add <- pathway_fullname
  output[["pathway_fullname"]] <- pathway_fullname
  
  # library(rlang)
  p_relab <- ggplot(tss_relab, aes(x = Dermotype_size, y = get(j))) +
    geom_boxplot(aes(fill = Dermotype_size)#, outlier.shape = NA
                 ) +
    #geom_jitter(color = "black", alpha = 0.1) +
    scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    scale_fill_manual(values = c("blue", "red", "orange", "green", "purple")) +
    ggpubr::theme_pubr(base_size = 10) + 
    labs(x = "", y = "Pathway abundance (%)") + # paste0(j, "\n", top_pwy_total_ctr$pwy_taxa[1])) +
    ggpubr::stat_compare_means(method = "wilcox", 
                               method.args = list(alternative = "two.sided"), 
                               label = "p.signif") +
    theme(legend.position = "none") +
    coord_flip() 
  output[["vis_relab"]] <- p_relab
   # 
  #----2. to plot contribution by each species-------#
  #-----2.0 read in kitome species -----#
  kitome <- readxl::read_excel("./Data/Kitome.xlsx")
  metadata_input <- mdata.dermotype.conf %>% filter(LibraryID %in% library_list)
  unique_ID <- "LibraryID"
  group <- "Dermotype_size"
  title_add <- pathway_fullname
  n_taxa <- n_taxa
  #----2.1 calculate % contributed by each species among stratifiable -------#
  tmp <- top_pwy_total_ctr |>
    mutate_all(~ifelse(is.na(.), 0, .)) |>
    filter(str_detect(pwy_taxa, "\\|")) |>
    mutate(pwy_taxa = sapply(strsplit(pwy_taxa, "\\|"), function(x) x[2])) |>
    mutate(pwy_taxa = str_extract(pwy_taxa, "s__.*")) |>
    # get ride of kitome taxa #
    filter(!pwy_taxa %in% kitome$Taxa) %>% 
    mutate(pwy_taxa = gsub( "s__", "", pwy_taxa)) |>
    column_to_rownames(var = "pwy_taxa") |>
    select(-ID) |>
    t() |>
    data.frame(check.names = F) |>
    rownames_to_column("LibraryID") %>%
    rowwise() %>%
    # mutate(
    #   sum_numeric = sum(c_across(where(is.numeric)), na.rm = TRUE)
    #   #Unstratified = ifelse(sum_numeric == 0, 0, max(0, 1 - sum_numeric))
    # ) %>%
    # ungroup() %>%
    # select(-sum_numeric) %>%
    # rowwise() %>%
    mutate(
      sum_numeric = sum(c_across(where(is.numeric)), na.rm = TRUE),
      across(where(is.numeric), ~ ifelse(sum_numeric == 0, ., ./sum_numeric), .names = "{col}")
    ) %>%
    filter(sum_numeric > 0) %>% 
    select(-sum_numeric) %>% 
    ungroup()
  output[["species_contri_percent"]] <- tmp
  
  #----2.2 rescale ---- #
  tmp_rescale <- tmp %>%
    # rowwise() %>%
    # mutate(
    #   sum_numeric = sum(c_across(where(is.numeric)), na.rm = TRUE),
    #   Unstratified = ifelse(sum_numeric == 0, 0, max(0, 1 - sum_numeric))
    # ) %>%
    # ungroup() %>%
    # select(-sum_numeric) %>%
    # rowwise() %>%
    # mutate(
    #   sum_numeric = sum(c_across(where(is.numeric)), na.rm = TRUE),
    #   across(where(is.numeric), ~ ifelse(sum_numeric == 0, ., ./sum_numeric), .names = "{col}")
    # ) %>%
    # ungroup() %>%
    left_join(tss_relab[, c(j, "LibraryID")], by = "LibraryID") %>%
    rowwise() %>%
    # mutate(
    #   Unstratified = ifelse(sum_numeric == 0 & .data[[j]] != 0, 1, Unstratified)) %>%
    mutate(
      across(where(is.numeric), ~ .*.data[[j]], .names = "{col}")
    ) %>%
    ungroup() %>%
    select(-any_of(c("sum_numeric", j)))
  
  output[["species_contri_rescale"]] <- tmp_rescale
  
  
  #----2.3 mean species contribution by groups---- #
  ave_contrib <- tmp_rescale %>% 
    left_join(select_at(mdata.dermotype.conf, c("LibraryID", "Dermotype_size"))) %>% 
    column_to_rownames("LibraryID") %>% 
    reshape2::melt() %>% 
    #filter(str_detect(variable, "mitis")) %>% 
    group_by(Dermotype_size, variable) %>% 
    #   mutate(mean = mean(value)) %>% 
    #   filter(mean > 0.005) %>% 
    #   ungroup %>% 
    #   ggplot(., aes(x = Dermotype_size, y = value*100)) +
    #   geom_boxplot() +
    #   scale_y_log10() +
    #   facet_wrap(~variable, nrow = 3) +
    #   ggpubr::stat_compare_means(method = "wilcox", hide.ns = F, label = "p.format")
    #group_by(Dermotype_size, variable) %>% 
    # ave_contrib 
    summarise(mean = mean(value, na.rm = T)) %>%
    mutate(mean = round(mean*100, 2)) %>%
    arrange(desc(mean), variable) %>%
    filter(mean > 0.1) %>%
    ggplot(., aes(fill= variable, y=mean, x = Dermotype_size)) +
    geom_bar(position="dodge", stat="identity") +
    facet_wrap(~variable)
  output[["mean_spec_contri_percent_vis"]] <- ave_contrib 
  
  
  library(tidyverse)
  
  # #------ 2.3.1 Calculate mean differences
  # mean_diff <- tmp %>% 
  #   select(-sum_numeric) %>% 
  #   left_join(select(mdata.dermotype.conf, LibraryID, Dermotype_size), by = "LibraryID") %>% 
  #   pivot_longer(-c(LibraryID, Dermotype_size), names_to = "variable", values_to = "value") %>% 
  #   group_by(variable, Dermotype_size) %>% 
  #   summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>% 
  #   pivot_wider(names_from = Dermotype_size, values_from = mean_value, names_prefix = "Dermotype_") %>% 
  #   mutate(diff = `Dermotype_Ac-1` - `Dermotype_Ac-2`) # Adjust category names as needed
  # 
  # # Step 2: Prepare data for stacked bar plot
  # plot_data <- mean_diff %>% 
  #   mutate(direction = ifelse(diff > 0, "Positive", "Negative")) %>% 
  #   mutate(abs_diff = abs(diff)) %>% 
  #   mutate(Site = i)
  # 
  # # Step 3: Plot
  # ggplot(plot_data, aes(x = Site, y = diff, fill = variable)) +
  #   geom_bar(stat = "identity", position = "stack") +
  #   coord_flip() +
  #   geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  #   labs(
  #     x = "Variable",
  #     y = "Difference in Mean",
  #     fill = "Direction",
  #     title = "Differences in Mean Between Dermotype Sizes"
  #   ) +
  #   ggpubr::theme_pubr(legend = "none")
  
  
  #---2.4 vis  % contributed by each species  -----#
  input = #tmp_rescale 
    tmp %>% 
    # filter(sum_numeric >0) %>% 
    # select(-sum_numeric) %>% 
    data.frame() %>%  
    mutate_all(~ifelse(is.na(.), 0, .)) %>% 
    column_to_rownames("LibraryID")
  n_taxa = min(ncol(input), n_taxa)
  
  #Creating bar plots for visualizations
  who = names(sort(apply(input, 2, function(x) mean(x)), decreasing = TRUE))[1:n_taxa]
  # who = names(sort(colMeans(input), decreasing = TRUE))[1:n_taxa]
  f = input[,names(input) %in% who]
  f$Other = rowSums(input[, !(names(input) %in% who)])
  #Sort f by most abundant taxa
  f = f[order(apply(f, 2, function(x) mean(x)), decreasing = T)]
  # adding in other taxa as a column for display
  who = c(who, "Other")
  # now order samples (rows) by the most abundant taxon/Sort the x-axis by most abundant
  f = f[order(f[,1], decreasing = T ), ]
  #transpose f to add in the sample names for the melt
  f = data.frame(t(f), check.names = FALSE)
  # head(f)
  # add Taxon column prior to melting
  f$Taxa <- row.names(f)
  # melt
  m = reshape2::melt(f)
  #make the colors for 29 taxa + 1 "Other" (30 total), making the "Other" gray.
  col.30 =   c("#4288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D44", "#E7298A", 
               "#00008B", "#8DD4C7", "#FFFFB4", "#BEBADA", "#FB8072", "#80B1D4", "#B4DE69", 
               "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D54E4F", "#FDAE61", "#FEE08B", 
               "#444444", "#66C2A5", "#6fabd0", "#1B9E77", "#D95F02", "#7570B4", "#E6AB02",
               "#A6761D", "#CCCCCC")
  col = scale_fill_manual(values = c(col.30[1:n_taxa], "grey90"))
  metadata_input = data.frame(metadata_input)
  V = "variable"
  m[[group]] = metadata_input[ match(m$variable, metadata_input[[unique_ID]]), group]
  p = ggplot(m, aes(m[[V]], y = value)) + 
    geom_bar(aes(fill=Taxa), stat = "identity", position = position_stack(reverse = TRUE)) +  
    facet_grid(.~ m[[group]], scales = "free", space ="free") +
    col + xlab("Sample") + ylab("") + 
    scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    theme_pubr(legend = "right", base_size = 10) + 
    theme(axis.text.x = element_blank(),  
          #axis.title = element_text(size = rel(1.1)),
          legend.text = element_text(face = "italic")) + 
    guides(fill = guide_legend(ncol = 1, 
                               title = paste("Top", n_taxa, "contributing taxa"),
                               reverse=TRUE, 
                               keyheight = 0.55)) +
    scale_y_continuous(
      #trans = ,
      #labels = percent_format(), 
      expand=c(0,0) ) +
    labs(title=paste(title_add)) + 
    theme(#plot.title = element_text(size = rel(1.2)),
          panel.spacing = unit(0, "lines")) 
  p$data$Taxa = factor(p$data$Taxa, ordered = TRUE, levels = who)
  plot(p)
  output[["vis_species_contri_rescale"]] <- p
  
  
  #---2.4.1 vis average % contributed by each species  -----#
  m_summary <- m %>%
    group_by(!!sym(group), Taxa) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  p <- ggplot(m_summary, aes(x = !!sym(group), y = value)) + 
    geom_bar(aes(fill = Taxa), stat = "identity", position = position_stack(reverse = TRUE)) +  
    col + 
    xlab("") + 
    ylab("Species contribution (%)") + 
    theme_pubr(legend = "right", base_size = 10) + 
    theme(
      axis.text.y = element_blank(),  
      legend.text = element_text(face = "italic"),  # smaller legend text
      #legend.title = element_text(size = rel(0.9)),                # smaller legend title
      legend.key.size = unit(0.3, "cm"),                           # smaller legend keys
      plot.margin = margin(5, 5, 5, 5, "pt"),                      # reduce plot margins if needed
      panel.spacing = unit(0, "lines")
    ) + 
    guides(fill = guide_legend(
      ncol = 1, 
      title = paste("Top", n_taxa, "contributing taxa"),
      reverse = TRUE, 
      keyheight = 0.4                                           # reduce key height
    )) +
    scale_y_continuous(
      expand = c(0, 0), labels = scales::percent_format(suffix = "")
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    coord_flip()
  
  # Adjust Taxa levels for ordered display
  p$data$Taxa <- factor(p$data$Taxa, ordered = TRUE, levels = who) 
  plot(p)
  output[["vis_average_species_contribution"]] <- p 
  

  #----3. plot pathway and phenotype associations----#
  df <- pwy_asso_phenotype |>
    filter(Site %in% i & ID %in% j) |>
    arrange(desc(abs(coef)) ) |>
    mutate(metadata_value = 
             ifelse(metadata == value, metadata, 
                    paste0(metadata, "\n(", value, ")"))) |>
    mutate(p.format = case_when(
      qval < 0.05 & qval >= 0.01  ~ "*",
      qval < 0.01 & qval >= 0.001~ "**", 
      qval < 0.001 ~ "***", 
      T ~ ""))
  p_pwy_pheno <- ggplot(df, aes(x = pwy, y = metadata_value)) +
    geom_tile(aes(fill = coef)) +
    scale_fill_gradient2(low = "#013220", high = "maroon", mid = "white", 
                         midpoint = 0,
                         guide=guide_colorbar(barwidth=.5, barheight = 3)) +
    geom_text(aes(label = p.format)) +
    theme_pubr(legend = "right", base_size = 10) +
    labs(x= "Pwy", y = "") +
    theme(axis.text.x = element_blank(),
          #axis.text.y = element_text(size = rel(1.1))
          )
  output[["pwy_pheno_association"]] <-   p_pwy_pheno
  g_2 <- gridExtra::arrangeGrob(output[["vis_relab"]],
                                output[["vis_species_contri_rescale"]], 
                                widths = c(1,4))
  g_3 <- 
    gridExtra::arrangeGrob(
      g_2,
      p_pwy_pheno,
      widths  = c(4,1))
  
  output[["vis_all"]] <- g_3

  return(output)
}
