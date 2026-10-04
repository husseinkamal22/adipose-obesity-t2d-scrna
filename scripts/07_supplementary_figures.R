# ============================================================
# 07_supplementary_figures.R
# Supplementary figures: QC, embeddings, annotation validation
# Project: Adipose obesity / metabolic disease scRNA-seq (GSE278526)
#
# READ-ONLY PLOTTING SCRIPT
# Flow: readRDS -> extract existing information -> plot -> ggsave
#
# FORBIDDEN IN THIS SCRIPT (do not add):
#   NormalizeData, FindVariableFeatures, ScaleData, RunPCA,
#   FindNeighbors, FindClusters, RunUMAP, FindAllMarkers, SingleR,
#   PercentageFeatureSet, saveRDS, merge/integration, any RDS edit
#
# Skipped on purpose (decision recorded in the figure plan):
#   FigureS2_QC_percent_mt.png  (percent.mt saved only in BRI-1456)
#
# Samples (tissue per official NCBI GEO GSE278526 sample records):
#   BRI-1456 = Healthy Omentum
#   BRI-1457 = Unhealthy Omentum
#   BRI-1458 = Healthy SQ
#   BRI-1459 = Unhealthy SQ
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)


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
# 2. Load annotated objects (read only, never written back)
# ------------------------------------------------------------

rds_files <- file.path(
  results_dir,
  sprintf("seurat%d_annotated.rds", 1:4)
)

stopifnot(all(file.exists(rds_files)))

seurat_list <- lapply(rds_files, readRDS)

required_meta <- c(
  "sample_id", "disease", "tissue",
  "nFeature_RNA", "nCount_RNA",
  "doublet_score", "doublet_class",
  "seurat_clusters", "SingleR_label",
  "cell_type", "broad_cell_type"
)

for (i in seq_along(seurat_list)) {
  s <- seurat_list[[i]]

  stopifnot(all(required_meta %in% colnames(s@meta.data)))
  stopifnot(all(c("pca", "umap") %in% names(s@reductions)))

  cat(
    "\nLoaded seurat", i, ":", ncol(s), "cells |",
    s@meta.data$sample_id[1], "\n"
  )
}

cat("\nAll four annotated objects loaded (read only).\n")


# ------------------------------------------------------------
# 3. Shared helpers: panel labels and ordered factors
# ------------------------------------------------------------

# Panel labels are built from GEO, NOT from the tissue column in
# the annotated RDS objects.
#
# Reason: data/GSE278526_cell_barcodes_metadata.tsv.gz carries a
# tissue column that is swapped relative to the official GEO
# sample records, and that value is what sits inside the saved
# annotated RDS objects (RDS files are read-only here and are
# never rewritten).
#
# Authoritative source: NCBI GEO GSE278526 family SOFT file
# (!Sample_characteristics_ch1, GSM8548235-GSM8548238):
#   GSM8548235 / BRI-1456 = Omental  / Healthy
#   GSM8548236 / BRI-1457 = Omental  / Unhealthy
#   GSM8548237 / BRI-1458 = Subcutaneous / Healthy
#   GSM8548238 / BRI-1459 = Subcutaneous / Unhealthy
geo_tissue <- c(
  "BRI-1456" = "Omentum",
  "BRI-1457" = "Omentum",
  "BRI-1458" = "SQ",
  "BRI-1459" = "SQ"
)

make_label <- function(meta_data) {
  sample_id <- meta_data$sample_id[1]

  paste0(
    sample_id,
    " (",
    meta_data$disease[1],
    " ",
    geo_tissue[[sample_id]],
    ")"
  )
}

sample_order <- vapply(
  seurat_list,
  function(x) x@meta.data$sample_id[1],
  character(1)
)

stopifnot(all(sample_order %in% names(geo_tissue)))

sample_labels <- vapply(
  seurat_list,
  function(x) make_label(x@meta.data),
  character(1)
)

names(seurat_list) <- sample_labels

sample_levels <- sample_labels[order(sample_order)]

cat("\nPanel labels (from GEO GSE278526 sample records):\n")
cat(paste0("  ", sample_levels, collapse = "\n"), "\n")

broad_levels <- c(
  "Stromal_fibroblast",
  "Myeloid_macrophage",
  "Endothelial",
  "T_cells",
  "NK_cells",
  "Dendritic_cells",
  "B_cells",
  "Epithelial_mesothelial",
  "Vascular_smooth_muscle",
  "Unresolved"
)


# ------------------------------------------------------------
# 4. Extract metadata table (single pass, reused below)
# ------------------------------------------------------------

meta_df <- bind_rows(
  lapply(names(seurat_list), function(nm) {
    md <- seurat_list[[nm]]@meta.data

    data.frame(
      sample_label = nm,
      nFeature_RNA = md$nFeature_RNA,
      nCount_RNA = md$nCount_RNA,
      doublet_score = md$doublet_score,
      doublet_class = as.character(md$doublet_class),
      seurat_clusters = as.character(md$seurat_clusters),
      SingleR_label = as.character(md$SingleR_label),
      cell_type = as.character(md$cell_type),
      broad_cell_type = as.character(md$broad_cell_type),
      stringsAsFactors = FALSE
    )
  })
)

meta_df$sample_label <- factor(
  meta_df$sample_label,
  levels = sample_levels
)

meta_df$broad_cell_type <- factor(
  meta_df$broad_cell_type,
  levels = broad_levels
)

cat("\nMetadata extracted:", nrow(meta_df), "cells.\n")


# ------------------------------------------------------------
# 5. FigureS1 — QC metrics (nFeature_RNA, nCount_RNA)
#     Post-filter values, from saved metadata only.
# ------------------------------------------------------------

qc_long <- meta_df %>%
  select(sample_label, nFeature_RNA, nCount_RNA) %>%
  pivot_longer(
    cols = c(nFeature_RNA, nCount_RNA),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      nFeature_RNA = "Detected genes per cell (nFeature_RNA)",
      nCount_RNA = "UMI counts per cell (nCount_RNA)"
    )
  )

p_s1 <- ggplot(
  qc_long,
  aes(
    x = sample_label,
    y = value,
    fill = sample_label
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    alpha = 0.85,
    linewidth = 0.2
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.9,
    fill = "white",
    linewidth = 0.2
  ) +
  scale_y_log10() +
  facet_wrap(~ metric, scales = "free_y") +
  guides(fill = "none") +
  labs(
    title = "Quality metrics per sample",
    subtitle = "Cells retained after upstream QC filtering",
    x = NULL,
    y = "Log10 value",
    caption = paste(
      "Source: seurat1-4_annotated.rds metadata.",
      "QC filtering was applied upstream in 01_primary_*.R.",
      "Distributions are post-filter."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "grey92"),
    plot.caption = element_text(hjust = 0)
  )

print(p_s1)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS1_QC_metrics.png"
  ),
  plot = p_s1,
  width = 10,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 6. FigureS3 — scDblFinder doublet scores (retained singlets)
#     doublet_class is "singlet" for every retained cell.
# ------------------------------------------------------------

p_s3 <- ggplot(
  meta_df,
  aes(
    x = sample_label,
    y = doublet_score,
    fill = sample_label
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    alpha = 0.85,
    linewidth = 0.2
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.9,
    fill = "white",
    linewidth = 0.2
  ) +
  guides(fill = "none") +
  labs(
    title = "scDblFinder doublet scores per sample",
    subtitle = "Only singlets are retained in the annotated objects",
    x = NULL,
    y = "Doublet score",
    caption = paste(
      "Source: seurat1-4_annotated.rds metadata (doublet_score, doublet_class).",
      "Doublet removal was applied upstream in 01_primary_*.R."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.caption = element_text(hjust = 0)
  )

print(p_s3)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS3_doublet_score.png"
  ),
  plot = p_s3,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# 7. Extract saved UMAP coordinates (one embedding per sample;
#    no integration, no re-embedding)
# ------------------------------------------------------------

umap_df <- bind_rows(
  lapply(names(seurat_list), function(nm) {
    s <- seurat_list[[nm]]

    cells <- colnames(s)
    emb <- Embeddings(s, "umap")[cells, , drop = FALSE]
    md <- s@meta.data[cells, , drop = FALSE]

    data.frame(
      sample_label = nm,
      umap_1 = emb[, 1],
      umap_2 = emb[, 2],
      seurat_clusters = as.character(md$seurat_clusters),
      broad_cell_type = as.character(md$broad_cell_type),
      stringsAsFactors = FALSE
    )
  })
)

umap_df$sample_label <- factor(
  umap_df$sample_label,
  levels = sample_levels
)

umap_df$broad_cell_type <- factor(
  umap_df$broad_cell_type,
  levels = broad_levels
)

cat("\nUMAP coordinates extracted:", nrow(umap_df), "cells.\n")


# ------------------------------------------------------------
# 8. FigureS4 — UMAP by clusters (panel per sample)
# ------------------------------------------------------------

cluster_pos <- umap_df %>%
  group_by(sample_label, seurat_clusters) %>%
  summarise(
    umap_1 = mean(umap_1),
    umap_2 = mean(umap_2),
    .groups = "drop"
  )

p_s4 <- ggplot(
  umap_df,
  aes(
    x = umap_1,
    y = umap_2,
    color = seurat_clusters
  )
) +
  geom_point(size = 0.3, alpha = 0.75) +
  geom_text_repel(
    data = cluster_pos,
    aes(x = umap_1, y = umap_2, label = seurat_clusters),
    size = 3.2,
    box.padding = 0.6,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~ sample_label, scales = "free") +
  labs(
    title = "UMAP by clustering resolution 0.5",
    subtitle = "One saved embedding per sample (no integrated UMAP)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cluster",
    caption = paste(
      "Source: reductions$umap and seurat_clusters stored in",
      "seurat1-4_annotated.rds."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    aspect.ratio = 1,
    plot.caption = element_text(hjust = 0)
  )

print(p_s4)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS4_UMAP_clusters.png"
  ),
  plot = p_s4,
  width = 13,
  height = 9,
  dpi = 300
)


# ------------------------------------------------------------
# 9. FigureS5 — UMAP by final broad cell type (panel per sample)
# ------------------------------------------------------------

p_s5 <- ggplot(
  umap_df,
  aes(
    x = umap_1,
    y = umap_2,
    color = broad_cell_type
  )
) +
  geom_point(size = 0.3, alpha = 0.75) +
  facet_wrap(~ sample_label, scales = "free") +
  labs(
    title = "UMAP by final broad cell type",
    subtitle = "One saved embedding per sample (no integrated UMAP)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Broad cell type",
    caption = paste(
      "Source: reductions$umap and broad_cell_type stored in",
      "seurat1-4_annotated.rds."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    aspect.ratio = 1,
    plot.caption = element_text(hjust = 0)
  )

print(p_s5)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS5_UMAP_broad_cell_type.png"
  ),
  plot = p_s5,
  width = 13,
  height = 9,
  dpi = 300
)


# ------------------------------------------------------------
# 10. FigureS6 — PCA variance explained (elbow)
#      Stdev comes from the saved pca reduction.
# ------------------------------------------------------------

pca_df <- bind_rows(
  lapply(names(seurat_list), function(nm) {
    stdev <- seurat_list[[nm]][["pca"]]@stdev

    data.frame(
      sample_label = nm,
      pc = seq_along(stdev),
      variance_pct = (stdev^2) / sum(stdev^2) * 100,
      stringsAsFactors = FALSE
    )
  })
)

pca_df$sample_label <- factor(
  pca_df$sample_label,
  levels = sample_levels
)

p_s6 <- ggplot(
  pca_df,
  aes(
    x = pc,
    y = variance_pct
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    color = "red"
  ) +
  facet_wrap(~ sample_label) +
  labs(
    title = "PCA variance explained per sample",
    subtitle = "Saved pca reduction; dashed line = dims 1:10 used upstream",
    x = "Principal component",
    y = "Variance explained (%)",
    caption = paste(
      "Source: reductions$pca stored in seurat1-4_annotated.rds.",
      "Neighbors/UMAP were computed on dims 1:10 in 01_primary_*.R."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    plot.caption = element_text(hjust = 0)
  )

print(p_s6)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS6_PCA_elbow.png"
  ),
  plot = p_s6,
  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 11. FigureS7 — Marker validation dot plot
#      DotPlot is used only to read the saved normalized data;
#      scale = FALSE keeps raw values comparable across samples.
# ------------------------------------------------------------

marker_panel <- c(
  "CD3D", "CD3E",
  "NKG7", "GNLY",
  "FCN1", "C1QB",
  "COL1A1", "COL3A1",
  "MYH11",
  "PTPRB", "ERG",
  "CLEC9A",
  "LILRA4", "CLEC4C",
  "MS4A1", "CD79A",
  "MSLN", "KRT19"
)

dot_df <- bind_rows(
  lapply(names(seurat_list), function(nm) {
    s <- seurat_list[[nm]]

    stopifnot(all(marker_panel %in% rownames(s[["RNA"]])))

    # scale = FALSE keeps raw mean(expm1) values, so the colour
    # scale stays comparable across the four samples.
    dp <- DotPlot(
      s,
      features = marker_panel,
      group.by = "broad_cell_type",
      scale = FALSE
    )$data

    data.frame(
      sample_label = nm,
      broad_cell_type = as.character(dp$id),
      feature = as.character(dp$features.plot),
      avg_exp = dp$avg.exp,
      pct_exp = dp$pct.exp,
      stringsAsFactors = FALSE
    )
  })
)

dot_df$sample_label <- factor(
  dot_df$sample_label,
  levels = sample_levels
)

dot_df$broad_cell_type <- factor(
  dot_df$broad_cell_type,
  levels = broad_levels
)

dot_df$feature <- factor(
  dot_df$feature,
  levels = marker_panel
)

# Display transform only: log compresses the raw mean(expm1)
# range so low and high expressing markers stay comparable.
dot_df$avg_log <- log1p(dot_df$avg_exp)

p_s7 <- ggplot(
  dot_df,
  aes(
    x = feature,
    y = broad_cell_type,
    color = avg_log,
    size = pct_exp
  )
) +
  geom_point() +
  facet_wrap(~ sample_label) +
  scale_color_gradient(
    low = "grey85",
    high = "darkblue"
  ) +
  scale_size(range = c(0.6, 4.5)) +
  labs(
    title = "Marker validation across broad cell types",
    subtitle = "Average expression and percentage of expressing cells",
    x = NULL,
    y = "Broad cell type",
    color = "Mean expression\n(log1p)",
    size = "% expressing",
    caption = paste(
      "Read from the saved RNA data slot of seurat1-4_annotated.rds;",
      "no normalization or marker testing was recomputed."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    strip.background = element_rect(fill = "grey92"),
    plot.caption = element_text(hjust = 0)
  )

print(p_s7)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS7_marker_validation.png"
  ),
  plot = p_s7,
  width = 15,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 12. FigureS8 — SingleR labels vs manual cluster annotation
#      Table only: stored SingleR_label vs stored cell_type.
#      SingleR is NOT re-run.
# ------------------------------------------------------------

singler_df <- bind_rows(
  lapply(names(seurat_list), function(nm) {
    tab <- as.data.frame(
      table(
        singleR = as.character(
          seurat_list[[nm]]@meta.data$SingleR_label
        ),
        manual = as.character(
          seurat_list[[nm]]@meta.data$cell_type
        )
      ),
      stringsAsFactors = FALSE
    )

    colnames(tab)[3] <- "count"

    tab$sample_label <- nm

    tab
  })
)

singler_df <- singler_df %>%
  filter(count > 0) %>%
  group_by(sample_label, manual) %>%
  mutate(
    proportion = count / sum(count)
  ) %>%
  ungroup()

singler_df$sample_label <- factor(
  singler_df$sample_label,
  levels = sample_levels
)

p_s8 <- ggplot(
  singler_df,
  aes(
    x = singleR,
    y = manual,
    fill = proportion
  )
) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ sample_label) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(0, 1)
  ) +
  labs(
    title = "SingleR labels vs manual cluster annotation",
    subtitle = "Row-normalized: share of each manual cell type per SingleR label",
    x = "SingleR label (Human Primary Cell Atlas)",
    y = "Manual cell type annotation",
    fill = "Proportion",
    caption = paste(
      "Table of stored SingleR_label and cell_type columns in",
      "seurat1-4_annotated.rds. SingleR was not re-run."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    strip.background = element_rect(fill = "grey92"),
    plot.caption = element_text(hjust = 0)
  )

print(p_s8)

ggsave(
  filename = file.path(
    figures_dir,
    "FigureS8_SingleR_vs_manual.png"
  ),
  plot = p_s8,
  width = 16,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 13. Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("Supplementary figures saved:\n")
cat("- FigureS1_QC_metrics.png\n")
cat("- FigureS3_doublet_score.png\n")
cat("- FigureS4_UMAP_clusters.png\n")
cat("- FigureS5_UMAP_broad_cell_type.png\n")
cat("- FigureS6_PCA_elbow.png\n")
cat("- FigureS7_marker_validation.png\n")
cat("- FigureS8_SingleR_vs_manual.png\n")
cat("(FigureS2 percent.mt intentionally skipped)\n")
cat("============================================\n")
