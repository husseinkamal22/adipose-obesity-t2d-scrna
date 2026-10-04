# ============================================================
# 06_final_figures.R
# Final figures for adipose scRNA-seq project
# ============================================================

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

results_dir <- here::here("results")

figures_dir <- file.path(
  results_dir,
  "figures"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load broad cell-type composition
# ------------------------------------------------------------

composition <- read.csv(
  file.path(
    results_dir,
    "broad_cell_type_composition_by_sample.csv"
  )
)

cat("\nFigure preparation started successfully.\n")

# ------------------------------------------------------------
# 3. Figure 1 — Broad cell-type composition
# ------------------------------------------------------------

p1 <- ggplot(
  composition,
  aes(
    x = sample_id,
    y = percentage,
    fill = broad_cell_type
  )
) +
  geom_col(width = 0.75) +
  facet_wrap(~ tissue, scales = "free_x") +
  labs(
    title = "Broad cell-type composition across adipose samples",
    x = NULL,
    y = "Cell proportion (%)",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p1)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure1_cell_type_composition.png"
  ),
  plot = p1,
  width = 10,
  height = 6,
  dpi = 300
)

cat("\nFigure 1 saved successfully.\n")

# ------------------------------------------------------------
# 4. Figure 2 — Inflammation and related gene-program scores
# ------------------------------------------------------------

inflam <- read.csv(
  file.path(
    results_dir,
    "inflammation_cell_level_scores.csv"
  )
)

inflam_summary <- inflam %>%
  group_by(
    sample_name,
    disease,
    tissue,
    broad_cell_type
  ) %>%
  summarise(
    InflammatoryScore = mean(
      InflammatoryScore,
      na.rm = TRUE
    ),
    ChemokineScore = mean(
      ChemokineScore,
      na.rm = TRUE
    ),
    ECMScore = mean(
      ECMScore,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Figure 2A — Inflammatory score
# ------------------------------------------------------------

p2a <- ggplot(
  inflam_summary,
  aes(
    x = broad_cell_type,
    y = InflammatoryScore,
    fill = disease
  )
) +
  geom_col(
    position = "dodge"
  ) +
  facet_wrap(
    ~ tissue,
    scales = "free_x"
  ) +
  labs(
    title = "Inflammatory gene-program scores",
    x = NULL,
    y = "Mean score",
    fill = "Disease state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p2a)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure2A_inflammatory_score.png"
  ),
  plot = p2a,
  width = 11,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Figure 2B — Chemokine score
# ------------------------------------------------------------

p2b <- ggplot(
  inflam_summary,
  aes(
    x = broad_cell_type,
    y = ChemokineScore,
    fill = disease
  )
) +
  geom_col(
    position = "dodge"
  ) +
  facet_wrap(
    ~ tissue,
    scales = "free_x"
  ) +
  labs(
    title = "Chemokine-related gene-program scores",
    x = NULL,
    y = "Mean score",
    fill = "Disease state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p2b)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure2B_chemokine_score.png"
  ),
  plot = p2b,
  width = 11,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Figure 2C — ECM/stromal score
# ------------------------------------------------------------

p2c <- ggplot(
  inflam_summary,
  aes(
    x = broad_cell_type,
    y = ECMScore,
    fill = disease
  )
) +
  geom_col(
    position = "dodge"
  ) +
  facet_wrap(
    ~ tissue,
    scales = "free_x"
  ) +
  labs(
    title = "ECM/stromal gene-program scores",
    x = NULL,
    y = "Mean score",
    fill = "Disease state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p2c)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure2C_ECM_score.png"
  ),
  plot = p2c,
  width = 11,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Save summary table
# ------------------------------------------------------------

write.csv(
  inflam_summary,
  file.path(
    results_dir,
    "inflammation_score_summary_by_broad_cell_type.csv"
  ),
  row.names = FALSE
)

cat("\nFigure 2A, 2B, and 2C saved successfully.\n")
# ------------------------------------------------------------
# 5. Load CellChat focused pathways
# ------------------------------------------------------------

cellchat_focus <- read.csv(
  file.path(
    results_dir,
    "cellchat_focused_pathways.csv"
  )
)

cat("\nCellChat focused pathways loaded successfully.\n")

# ------------------------------------------------------------
# 6. Figure 3 — Focused CellChat pathways
# ------------------------------------------------------------

p3 <- ggplot(
  cellchat_focus,
  aes(
    x = pathway,
    y = communication_probability,
    fill = sample
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "Selected inferred cell-cell communication pathways",
    x = "CellChat pathway",
    y = "Communication probability",
    fill = "Sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p3)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure3_focused_CellChat_pathways.png"
  ),
  plot = p3,
  width = 12,
  height = 7,
  dpi = 300
)

cat("\nFigure 3 saved successfully.\n")

# ------------------------------------------------------------
# 7. Load CXCL top sender-receiver interactions
# ------------------------------------------------------------

cxcl_top10 <- read.csv(
  file.path(
    results_dir,
    "cellchat_CXCL_top10_sender_receiver.csv"
  )
)

cxcl_top10 <- cxcl_top10 %>%
  mutate(
    interaction = paste(
      source,
      "→",
      target
    )
  )

# ------------------------------------------------------------
# 8. Figure 4 — CXCL sender-receiver patterns
# ------------------------------------------------------------

p4 <- ggplot(
  cxcl_top10,
  aes(
    x = reorder(
      interaction,
      total_prob
    ),
    y = total_prob,
    fill = sample
  )
) +
  geom_col() +
  facet_wrap(
    ~ sample,
    scales = "free_y"
  ) +
  coord_flip() +
  labs(
    title = "Top inferred CXCL sender-receiver interactions",
    x = "Sender → receiver",
    y = "Total communication probability",
    fill = "Sample"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "none"
  )

print(p4)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure4_CXCL_sender_receiver.png"
  ),
  plot = p4,
  width = 12,
  height = 9,
  dpi = 300
)

cat("\nFigure 4 saved successfully.\n")