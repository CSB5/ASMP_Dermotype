library(dplyr)
library(ggplot2)
library(tidyverse)
library(scales)

mdata.dermotype <- read.csv("./Data/Metadata_dermotype_BatchABCD_with_StabilityScore.csv")

mdata.dermotype = mdata.dermotype %>% 
  add_count(Dermotype_size) %>%
  mutate(Dermotype_n = paste(Dermotype_size, ' (n=', n, ')', sep = "")) 
tmp = mdata.dermotype %>% 
  distinct(Dermotype_size)


# load previously-shortlisted High Prevalence Taxa at site- & dermotype- level
load("./Data/2_basic/HighPrev_site.RData")
d_all_site = d_all %>% dplyr::rename(Category = Site)
rm(d_all)

load("./Data/2_basic/HighPrev_Dermotype.RData")

d_all_dermotype <- d_all %>% 
  left_join(tmp, 
            join_by(Dermotype == Dermotype_size)) %>% 
  dplyr::rename(Category = Dermotype)
rm(d_all)

d_Category = rbind(d_all_site, d_all_dermotype)
d_all = d_Category %>% 
  mutate(Genus = sub(" .*", "", Taxa))

d_all = d_Category %>% 
  mutate(Genus = sub(" .*", "", Taxa)) %>%
  mutate(Site = sub("\\-.*", "", Category)) %>%
  filter(grepl("Ch-|Ll-", Category) == F) 

#========to sort Taxa in an order that makes sense + easier to spot pattern =========#
#--------1. based on how many times a taxon categorized as a core -------------------#
d_all.2 = d_all
d_all = d_all %>% 
  group_by(Taxa) %>% 
  summarise(N = n(), where = paste(Site, collapse = "&"), 
            mean_abund = median(median_abund), mean_prev = median(prev)) %>% 
  arrange(-N, where, -mean_prev, -mean_abund) 

#--------2. seperate Malassezia from bacteria--------------------------------------#
order_by_count <-  unlist(d_all$Taxa[!d_all$Taxa %in% 
                                           c("Malassezia globosa",
                                             "Malassezia restricta")])

#--------3. clustering taxa based on their presence/absence in core species--------#
#------Not necessary, just to better summarize trends ---------------------#
library(tidyr)
wide_tab <- d_all.2 %>% 
  select(Taxa, Category, median_abund) %>%
  filter(!grepl("Malassezia", Taxa)) %>%
  spread(., key = Category, value = median_abund) %>%
  column_to_rownames(var = "Taxa") %>%
  mutate_all(., ~as.integer(!is.na(.))) %>%
  mutate_all(., ~replace_na(.,0)) 
library(ggdendro)
library(dendextend)
hc <- hclust(dist(wide_tab), "average")
dend <- hc %>% as.dendrogram() %>%
  #rev.dendrogram() %>% 
  rotate(., c( 2:4, 5:6,1, 13:15, 11:12,10,9,8,7)) %>%
  dendro_data(., type = "rectangle")
g_dend <- ggplot() +
  geom_segment(data = segment(dend),
               aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  geom_text(data = label(dend),
            aes(x = x, y = y, label = label, hjust = 0),
            size = 3
  ) +
  coord_flip() +
  scale_y_reverse(expand = c(1, 0)) +
  theme_dendro()

#-------------4. extract the correct order ----------------#
taxa_order_plot <- c(rev(dend$labels$label),
                     "Malassezia globosa",
                     "Malassezia restricta")

#------------5. resort the original data and plot it ----------------#
d_all.2 = d_all.2 %>% 
  mutate(Taxa = factor(d_all.2$Taxa, levels = taxa_order_plot),
         SiteType.1 = case_when(
           Site %in% c("Ll", "Vf", "Ch", "Fo", "Ub") ~ "Dry", 
           Site %in% c("Ac", "Ax", "Sc", "Ps") ~ "Moist", 
           T~"Others"),
         SiteType.2 = case_when(
           Site %in% c("Ac", "Ax", "Ll", "Vf") ~ "Non-Sebaceous", 
           Site %in% c("Sc", "Ps", "Ch", "Fo", "Ub") ~ "Sebaceous", 
           T~"Others"),
         Site = factor(Site, levels = 
                         c(
                           "Ax",
                           "Ac",
                           "Ll", 
                           "Vf",
                           "Ub", 
                           "Fo",  
                           "Ch", 
                           "Ps",
                           "Sc")),
         highlight_group = case_when(median_abund < 0.005 ~ "ExtremelyLow",
                                     T ~ "Others"),
         prev_group = case_when(prev < 0.75 ~ "Prevalent", 
                                #prev < 0.9 ~ "Shell", 
                                T ~ "Core"))
        
library(khroma)
library(scales)
p = ggplot(d_all.2, 
           aes(x = Category, y = Taxa, fill = median_abund*100
               )) + 
  geom_point(aes(
    size = prev, 
    #color = highlight_group,
    stroke = ifelse(highlight_group == "ExtremelyLow", 1.5, 0.5)),
    shape= 21)+
  scale_size_continuous(name = "Prevalence\n(%)", 
                        breaks = c(0.5,0.75, 1),
                        labels = c(50, 75, 100),
                        range = c(1, 6)) +
  #scale_color_manual(values = c("purple","black"))+
  #scale_shape_manual(values=c("ExtremelyLow" = 1, "Others" = 21)) + 
  scale_fill_gradient2(low = "white", 
                       #high = "#08306B", mid = "#9ECAE1", 
                       high = "#00008B", mid = "skyblue", 
                       midpoint = 0.3,
                       name = "Non-zero\nmedian\nabundance\n(%)",    # Name for the color legend
                       #labels = scales::percent_format(scale = 100, accuracy = 0.1, suffix = ""),
                       trans = "log10",
                       limits = c(min(d_all.2$median_abund*100), 100), 
                       # breaks = c(0.005, 0.01, 0.1, 0.5, 1), 
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', function(x) round(10^x,1)),
                       guide = guide_colourbar(nbin = 100, barwidth = 1, 
                                               barheight = 8, title.position = "top")) +
  scale_x_discrete(position = "top") +   
  scale_y_discrete(limits=rev) +
  facet_grid(cols = vars(Site), 
             rows = vars(Taxa_group), 
             scales = "free", space='free') +
  facet_grid(cols = vars(Site), 
             rows = vars(Taxa_group),
             scales = "free", space='free') +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(reverse = TRUE)) +
  theme_pubr(legend = "right") +
  theme(axis.text.x = element_text(angle = 50,
                                   hjust = 0, vjust = 0,  face = "bold"),
        axis.text.y = element_text(face = "italic", 
                                   colour = "black"),
        legend.position.inside = c(0.5, 0.9), 
        legend.direction = "vertical",
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.line.y.left = element_blank(), 
        strip.text.y.right = element_blank(),
        panel.grid.major.x = element_line(color = "grey20", linewidth = 0.25)
        ) 
           
p
