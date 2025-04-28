n_top <- 2
top_degradation <- All_res_1 %>%
  # filter(Site == PickedSite) %>% 
  #filter(prev > prev_threshold) |>
  filter(!(Site != "Ax" & str_detect(Dermotype_size, "-1"))) |>
  dplyr::select(ID, pwy, 
                coef, Site, Dermotype_size,
                Superclass1, Superclass2, Superclass3, Superclass4,
                label_1, label_2, Subclass,
                prev, HighPrev, Type) %>% 
  filter(Superclass1 == "Degradation") %>% 
  filter(Type != "all_low") %>% 
  group_by(Dermotype_size) %>%
  arrange(Dermotype_size, -abs(coef)) %>%
  slice_max(n= n_top, abs(coef)) |>
  ungroup() 

input <- All_res_1 %>%
  filter(ID %in% top_degradation$ID) %>% 
  filter(!(Site != "Ax" & str_detect(Dermotype_size, "-1"))) |>
  mutate(is_shared = ifelse(Type == "conserved", TRUE, FALSE)) |>
  mutate(pwy = str_replace_all(pwy, "N-acetylglucosamine", "GlcNAc")) |>
  mutate(pwy = str_replace_all(pwy, "N-acetylmannosamine", "ManNAc")) |>
  mutate(pwy = str_replace_all(pwy, "N-acetylneuraminate", "NeuAc")) |>
  # mutate(pwy = str_replace_all(pwy, "superpathway of", "SP")) |>
  mutate(pwy = str_replace_all(pwy, " \\s*\\([^\\)]+\\)", "")) %>%
  # left_join(color_mapping, join_by(Superclass1)) %>%
  # mutate(text_color = glue("<span style='color:{color}'>{pwy}</span>")) %>%
  mutate(Site = factor(Site, levels = c(
    "Ax", "Ac", "Vf", "Ub", "Fo", "Ps", "Sc"
  ))) %>% 
  arrange(Site, desc(Subclass))

  

heatmap <-  ggplot(input, aes(x = Dermotype_size, y = pwy)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", guide = "colorbar", 
                       name = "Beta coefficient") + 
  geom_point(data = . %>% filter(is_shared), aes(x = Dermotype_size, y = pwy), 
             shape = 21, size = 3, fill = "white", color = "black") +
  facet_nested( Subclass~Site,
                # Superclass4 ~Site + Dermotype_size,
                labeller = label_wrap_gen(width= 5, multi_line = T),
                space = "free",
                scales = "free",
                switch = "y"
  ) +
  labs(x = '', y = '', fill = 'log2FC') +
  theme_pubr(base_size = 10, x.text.angle = 30) +
  theme(strip.text.y.left  = element_blank(),
  strip.placement = "inside",
  strip.text.x.top = element_text(angle = 0),
  panel.spacing.y = unit(0.2, "lines"),
  panel.spacing.x = unit(0, "lines"),
  # legend.position = "top",
  # legend.direction = "horizontal",
  # legend.box = "horizontal",
  legend.position = "right",
  legend.direction = "vertical"
  #legend.margin = margin(t = 5, b = 5),
  #legend.position.inside = c(0, 0,9)
  # legend.text = element_text(size = rel(1.4)),
  # legend.title = element_text(size = rel(1.4))
  ) 

  