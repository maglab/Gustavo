library(dplyr)
library(tidyr)
library(ggplot2)


GenArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::select(Gene,ARC_Meaning) %>% unique() %>% dplyr::rename(Arc=ARC_Meaning)

PltFrm = GenArcFrm %>% pull(Gene) %>% table() %>% as.data.frame()#%>% sort() %>% rev()
colnames(PltFrm) = c("Gene", "Plt")

PltArcFrm = merge(GenArcFrm,PltFrm)

PltArcFrm = PltArcFrm %>% rename(Gen = Gene)

# Asegurar tipos y niveles
PltArcFrm <- PltArcFrm %>%
  mutate(
    Plt = factor(Plt, levels = 1:6),
    Arc = factor(Arc)
  )

# Conteos por ARC × Plt
head(PltArcFrm)

plot_df <- PltArcFrm %>%
  dplyr::count(Arc, Plt, name = "n") %>%
  tidyr::complete(Arc, Plt, fill = list(n = 0))

#plot_df <- PltArcFrm %>%
#  count(Arc, Plt, name = "n") %>%
#  complete(Arc, Plt, fill = list(n = 0))

# Etiquetas (a., b., c., ...)
arc_levels <- levels(plot_df$Arc)
strip_labs <- setNames(
  paste0(letters[seq_along(arc_levels)], ". ", arc_levels),
  arc_levels
)

# Totales por mitad
counts_long <- plot_df %>%
  mutate(plt_num = as.numeric(Plt),
         side = if_else(plt_num <= 3, "Low", "High")) %>%
  group_by(Arc, side) %>%
  summarise(total = sum(n), .groups = "drop")

# Posiciones verticales altas
arc_max <- plot_df %>%
  group_by(Arc) %>%
  summarise(m = max(n), .groups = "drop")

labels_df <- counts_long %>%
  left_join(arc_max, by = "Arc") %>%
  mutate(
    x = if_else(side == "Low", 2, 5),
    y_label = m * 1.55,       # más cerca del tope
    y_count = m * 1.43,       # justo debajo
    label = if_else(side == "Low", "Low ARC-Pleiotropy", "High ARC-Pleiotropy"),
    countlab = paste0("n = ", total)
  )

# Gráfica final
ggplot(plot_df, aes(x = as.numeric(Plt), y = n)) +
  geom_col(fill = "#a6cee3", color = "grey30", width = 0.85) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey40") +
  geom_text(data = labels_df, inherit.aes = FALSE,
            aes(x = x, y = y_label, label = label),
            size = 2.8, vjust = -0.5, fontface = "bold") +
  geom_text(data = labels_df, inherit.aes = FALSE,
            aes(x = x, y = y_count, label = countlab),
            size = 2.6, vjust = 1) +
  facet_wrap(~ Arc, ncol = 2, scales = "free_y",
             labeller = labeller(Arc = strip_labs)) +
  scale_x_continuous(breaks = 1:6, labels = 1:6,
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    x = "ARC-Pleiotropy level",
    y = "Number of genes",
    title = "Distribution of ARC-Pleiotropy levels per ARC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "#f2f2f2", colour = NA),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(12, 12, 12, 12),
    # 👇 aumenta la distancia vertical entre filas:
    panel.spacing.y = unit(1.3, "lines")  # puedes probar 1.5 o 2 para más separación
  )


# Guardar en SVG tamaño carta (8.5 x 11 pulgadas)
ggsave(
  filename = "Last_Figures/Supp_ARC_Pleiotropy_Distribution.svg",
  width = 8.5,     # pulgadas
  height =7,     # pulgadas
  units = "in",
  dpi = 300,
  device = "svg"
)
