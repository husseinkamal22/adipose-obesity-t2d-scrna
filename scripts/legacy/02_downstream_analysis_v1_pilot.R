# ============================================================
# 02 — DOWNSTREAM ANALYSIS
# ============================================================

# RESEARCH QUESTION
#
# Which adipose cell populations and intercellular signaling
# programs are associated with the inflammatory state of
# unhealthy obesity/metabolic disease?
#
# Study design:
#
# BRI-1456 = Healthy Omentum
# BRI-1457 = Unhealthy Omentum
# BRI-1458 = Healthy SQ
# BRI-1459 = Unhealthy SQ
#
# This script starts from the already processed Seurat objects.
# Primary processing has already been completed and saved.
# ============================================================


# ============================================================
# 01 — LOAD PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)


# ============================================================
# 02 — LOAD THE FOUR PROCESSED SEURAT OBJECTS
# ============================================================

seurat1 <- readRDS(
  "../results/seurat1_primary_processed.rds"
)

seurat2 <- readRDS(
  "../results/seurat2_primary_processed.rds"
)

seurat3 <- readRDS(
  "../results/seurat3_primary_processed.rds"
)

seurat4 <- readRDS(
  "../results/seurat4_primary_processed.rds"
)


# ============================================================
# 03 — DEFINE SAMPLE INFORMATION
# ============================================================

seurat1$disease <- "Healthy"
seurat1$tissue <- "Omentum"
seurat1$sample_name <- "BRI-1456"

seurat2$disease <- "Unhealthy"
seurat2$tissue <- "Omentum"
seurat2$sample_name <- "BRI-1457"

seurat3$disease <- "Healthy"
seurat3$tissue <- "SQ"
seurat3$sample_name <- "BRI-1458"

seurat4$disease <- "Unhealthy"
seurat4$tissue <- "SQ"
seurat4$sample_name <- "BRI-1459"


# ============================================================
# 04 — CHECK SAMPLE SIZES
# ============================================================

cat(
  "BRI-1456 | Healthy | Omentum:",
  ncol(seurat1),
  "cells\n"
)

cat(
  "BRI-1457 | Unhealthy | Omentum:",
  ncol(seurat2),
  "cells\n"
)

cat(
  "BRI-1458 | Healthy | SQ:",
  ncol(seurat3),
  "cells\n"
)

cat(
  "BRI-1459 | Unhealthy | SQ:",
  ncol(seurat4),
  "cells\n"
)


# ============================================================
# 05 — CHECK AVAILABLE METADATA
# ============================================================

cat("\n--- BRI-1456 metadata ---\n")
print(colnames(seurat1@meta.data))

cat("\n--- BRI-1457 metadata ---\n")
print(colnames(seurat2@meta.data))

cat("\n--- BRI-1458 metadata ---\n")
print(colnames(seurat3@meta.data))

cat("\n--- BRI-1459 metadata ---\n")
print(colnames(seurat4@meta.data))


# ============================================================
# 06 — CHECK CELL TYPE ANNOTATIONS
# ============================================================

cat("\n--- BRI-1456 cell_type ---\n")

if ("cell_type" %in% colnames(seurat1@meta.data)) {
  print(table(seurat1$cell_type))
} else {
  cat("No cell_type column found.\n")
}


cat("\n--- BRI-1457 cell_type ---\n")

if ("cell_type" %in% colnames(seurat2@meta.data)) {
  print(table(seurat2$cell_type))
} else {
  cat("No cell_type column found.\n")
}


cat("\n--- BRI-1458 cell_type ---\n")

if ("cell_type" %in% colnames(seurat3@meta.data)) {
  print(table(seurat3$cell_type))
} else {
  cat("No cell_type column found.\n")
}


cat("\n--- BRI-1459 cell_type ---\n")

if ("cell_type" %in% colnames(seurat4@meta.data)) {
  print(table(seurat4$cell_type))
} else {
  cat("No cell_type column found.\n")
}


# ============================================================
# 07 — CHECK SINGLER ANNOTATIONS
# ============================================================

cat("\n--- BRI-1456 SingleR ---\n")

if ("SingleR_label" %in% colnames(seurat1@meta.data)) {
  print(table(seurat1$SingleR_label))
} else {
  cat("No SingleR_label column found.\n")
}


cat("\n--- BRI-1457 SingleR ---\n")

if ("SingleR_label" %in% colnames(seurat2@meta.data)) {
  print(table(seurat2$SingleR_label))
} else {
  cat("No SingleR_label column found.\n")
}


cat("\n--- BRI-1458 SingleR ---\n")

if ("SingleR_label" %in% colnames(seurat3@meta.data)) {
  print(table(seurat3$SingleR_label))
} else {
  cat("No SingleR_label column found.\n")
}


cat("\n--- BRI-1459 SingleR ---\n")

if ("SingleR_label" %in% colnames(seurat4@meta.data)) {
  print(table(seurat4$SingleR_label))
} else {
  cat("No SingleR_label column found.\n")
}


# ============================================================
# END OF INITIAL DOWNSTREAM SETUP
# ============================================================
#
# At this stage we have:
#
# 1. Loaded all four processed samples.
# 2. Added disease and tissue information.
# 3. Checked the available metadata.
# 4. Checked the existing cell-type annotations.
#
# We will NOT perform comparative analysis yet.
#
# The next step is to harmonize cell-type annotations across
# all four samples before comparing Healthy vs Unhealthy.
# ============================================================

# ============================================================
# 08 — BRI-1457: CLUSTER AND MARKER INSPECTION
# ============================================================

# BRI-1457 = Unhealthy Omentum
#
# Goal:
# Inspect the existing clusters and their marker genes
# before assigning unified cell-type annotations.
# ============================================================


# ============================================================
# 08 — BRI-1457: CLUSTER AND MARKER INSPECTION
# ============================================================

# BRI-1457 = Unhealthy Omentum
#
# Goal:
# Inspect the existing clusters and their marker genes
# before assigning unified cell-type annotations.
# ============================================================


# ============================================================
# 08.1 — CHECK CLUSTER SIZES
# ============================================================

cat("\n--- BRI-1457 cluster sizes ---\n")

cluster_sizes_2 <- table(
  seurat2$seurat_clusters
)

print(cluster_sizes_2)


# ============================================================
# 08.2 — CHECK NUMBER OF CLUSTERS
# ============================================================

cat(
  "\nNumber of clusters in BRI-1457:",
  length(unique(seurat2$seurat_clusters)),
  "\n"
)


# ============================================================
# 08.3 — FIND MARKERS FOR BRI-1457
# ============================================================

markers_2 <- FindAllMarkers(
  seurat2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


# ============================================================
# 08.4 — TOP 10 MARKERS PER CLUSTER
# ============================================================

top10_markers_2 <- markers_2 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()


# ============================================================
# 08.5 — PRINT TOP MARKERS
# ============================================================

print(top10_markers_2)


# ============================================================
# 08.6 — SAVE BRI-1457 MARKERS
# ============================================================

write.csv(
  markers_2,
  "../results/BRI-1457_all_markers.csv",
  row.names = FALSE
)

write.csv(
  top10_markers_2,
  "../results/BRI-1457_top10_markers.csv",
  row.names = FALSE
)


# ============================================================
# 08.7 — VISUALIZE TOP MARKERS
# ============================================================

DoHeatmap(
  seurat2,
  features = unique(top10_markers_2$gene)
) +
  NoLegend()
# ============================================================
# 09 — ANNOTATE BRI-1458
# ============================================================

cluster_annotations_3 <- c(
  `0`  = "Unresolved",
  `1`  = "Stromal",
  `2`  = "T_cells",
  `3`  = "Fibroblast_stromal",
  `4`  = "NK_cells",
  `5`  = "T_cells",
  `6`  = "Inflammatory_monocyte_like",
  `7`  = "Epithelial_mesothelial_like",
  `8`  = "Myeloid_monocyte_macrophage_like",
  `9`  = "Endothelial",
  `10` = "Vascular_smooth_muscle",
  `11` = "B_cells",
  `12` = "Stromal_mesenchymal",
  `13` = "Macrophage_myeloid",
  `14` = "Epithelial_mesothelial_like",
  `15` = "pDC_like"
)

seurat3$cell_type <- unname(
  cluster_annotations_3[as.character(seurat3$seurat_clusters)]
)

cat("\n--- BRI-1458 final cluster annotation ---\n")

print(
  data.frame(
    Cluster = names(cluster_annotations_3),
    Cell_Type = unname(cluster_annotations_3),
    Cells = as.numeric(
      table(
        factor(
          seurat3$seurat_clusters,
          levels = names(cluster_annotations_3)
        )
      )
    )
  )
)

cat("\n--- BRI-1458 cell-type summary ---\n")

print(
  sort(
    table(seurat3$cell_type),
    decreasing = TRUE
  )
)
# ============================================================
# 10 — HARMONIZE CELL-TYPE ANNOTATIONS
# ============================================================

# ------------------------------------------------------------
# BRI-1456 — Healthy Omentum
# ------------------------------------------------------------

cluster_annotations_1 <- c(
  `0`  = "Stromal",
  `1`  = "Stromal_mesenchymal",
  `2`  = "T_cells",
  `3`  = "NK_cells",
  `4`  = "Unresolved",
  `5`  = "Myeloid_monocyte_macrophage_like",
  `6`  = "Fibroblast_stromal",
  `7`  = "Vascular_smooth_muscle",
  `8`  = "Endothelial",
  `9`  = "DC_like",
  `10` = "pDC_like"
)

seurat1$cell_type <- unname(
  cluster_annotations_1[as.character(seurat1$seurat_clusters)]
)


# ------------------------------------------------------------
# BRI-1457 — Unhealthy Omentum
# ------------------------------------------------------------

cluster_annotations_2 <- c(
  `0`  = "Stromal",
  `1`  = "Fibroblast_stromal",
  `2`  = "Unresolved",
  `3`  = "Endothelial",
  `4`  = "T_cells",
  `5`  = "NK_cells",
  `6`  = "DC_like",
  `7`  = "Stromal_mesenchymal",
  `8`  = "Vascular_smooth_muscle",
  `9`  = "Inflammatory_monocyte_like",
  `10` = "Myeloid_monocyte_macrophage_like",
  `11` = "Endothelial_ACKR1",
  `12` = "pDC_like"
)

seurat2$cell_type <- unname(
  cluster_annotations_2[as.character(seurat2$seurat_clusters)]
)


# ------------------------------------------------------------
# BRI-1458 — Healthy SQ
# ------------------------------------------------------------

# Already assigned above.
# We keep the same harmonized vocabulary.


# ------------------------------------------------------------
# BRI-1459 — Unhealthy SQ
# ------------------------------------------------------------

cluster_annotations_4 <- c(
  `0`  = "Unresolved",
  `1`  = "Stromal",
  `2`  = "Epithelial_mesothelial_like",
  `3`  = "T_cells",
  `4`  = "Myeloid_monocyte_macrophage_like",
  `5`  = "T_cells",
  `6`  = "Cytotoxic_T_cells",
  `7`  = "Stromal_mesenchymal",
  `8`  = "NK_cells",
  `9`  = "Fibroblast_stromal",
  `10` = "Endothelial",
  `11` = "B_cells",
  `12` = "Epithelial_mesothelial_like",
  `13` = "Vascular_smooth_muscle"
)

seurat4$cell_type <- unname(
  cluster_annotations_4[as.character(seurat4$seurat_clusters)]
)


# ============================================================
# 11 — CHECK HARMONIZED CELL TYPES
# ============================================================

cat("\n--- BRI-1456 ---\n")
print(sort(table(seurat1$cell_type), decreasing = TRUE))

cat("\n--- BRI-1457 ---\n")
print(sort(table(seurat2$cell_type), decreasing = TRUE))

cat("\n--- BRI-1458 ---\n")
print(sort(table(seurat3$cell_type), decreasing = TRUE))

cat("\n--- BRI-1459 ---\n")
print(sort(table(seurat4$cell_type), decreasing = TRUE))


# ============================================================
# 12 — SAVE ANNOTATED OBJECTS
# ============================================================

saveRDS(
  seurat1,
  "../results/seurat1_annotated.rds"
)

saveRDS(
  seurat2,
  "../results/seurat2_annotated.rds"
)

saveRDS(
  seurat3,
  "../results/seurat3_annotated.rds"
)

saveRDS(
  seurat4,
  "../results/seurat4_annotated.rds"
)

cat("\nAnnotated Seurat objects saved successfully.\n")

# ============================================================
# 13 — CELL-TYPE COMPOSITION ACROSS SAMPLES
# ============================================================

# Combine metadata from the four annotated samples

meta1 <- seurat1@meta.data
meta2 <- seurat2@meta.data
meta3 <- seurat3@meta.data
meta4 <- seurat4@meta.data

meta1$sample_name <- "BRI-1456"
meta1$disease <- "Healthy"
meta1$tissue <- "Omentum"

meta2$sample_name <- "BRI-1457"
meta2$disease <- "Unhealthy"
meta2$tissue <- "Omentum"

meta3$sample_name <- "BRI-1458"
meta3$disease <- "Healthy"
meta3$tissue <- "SQ"

meta4$sample_name <- "BRI-1459"
meta4$disease <- "Unhealthy"
meta4$tissue <- "SQ"


# ------------------------------------------------------------
# Combine the four metadata tables
# ------------------------------------------------------------

all_meta <- bind_rows(
  meta1,
  meta2,
  meta3,
  meta4
)


# ------------------------------------------------------------
# Cell counts
# ------------------------------------------------------------

cell_counts <- all_meta %>%
  count(
    sample_name,
    disease,
    tissue,
    cell_type,
    name = "n_cells"
  )


cat("\n--- CELL COUNTS ---\n")

print(cell_counts)


# ------------------------------------------------------------
# Cell percentages within each sample
# ------------------------------------------------------------

cell_composition <- cell_counts %>%
  group_by(sample_name) %>%
  mutate(
    total_cells = sum(n_cells),
    percentage = 100 * n_cells / total_cells
  ) %>%
  ungroup()


cat("\n--- CELL COMPOSITION (%) ---\n")

print(
  cell_composition %>%
    arrange(sample_name, desc(percentage))
)


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

write.csv(
  cell_counts,
  "../results/cell_type_counts_by_sample.csv",
  row.names = FALSE
)

write.csv(
  cell_composition,
  "../results/cell_type_composition_by_sample.csv",
  row.names = FALSE
)


cat(
  "\nCell-type composition tables saved successfully.\n"
)
# ============================================================
# 14 — PRINT COMPLETE CELL COMPOSITION TABLE
# ============================================================

cat("\n--- COMPLETE CELL COMPOSITION TABLE ---\n")

print(
  cell_composition %>%
    arrange(
      tissue,
      disease,
      desc(percentage)
    ),
  n = Inf
)
# ============================================================
# 16 — INFLAMMATION-RELATED CELL-STATE SCORING
# ============================================================
#
# Research question:
# Which adipose cell populations and intercellular signaling
# programs are associated with the inflammatory state of
# unhealthy obesity/metabolic disease?
#
# Purpose:
# Estimate three biological programs at single-cell level:
# 1. Inflammatory response
# 2. Chemokine signaling
# 3. ECM / stromal remodeling
#
# Important:
# - These are descriptive / hypothesis-generating scores.
# - They are NOT statistical tests.
# - They do NOT prove causality.
# - We will compare samples descriptively before any formal
#   inferential analysis.
# ============================================================


cat("\n============================================================\n")
cat("16 — INFLAMMATION-RELATED CELL-STATE SCORING\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 16.1 — Define gene programs
# ------------------------------------------------------------

cat("\n--- DEFINING BIOLOGICAL GENE PROGRAMS ---\n")

inflammatory_genes <- c(
  "IL1B",
  "IL6",
  "TNF",
  "NFKBIA",
  "NFKBIZ",
  "TNFAIP3",
  "CXCL8",
  "CCL2",
  "CCL3",
  "CCL4",
  "CCL5",
  "S100A8",
  "S100A9",
  "S100A10",
  "S100A12",
  "ICAM1",
  "VCAM1",
  "PTGS2",
  "NAMPT",
  "IRF1"
)


chemokine_genes <- c(
  "CCL2",
  "CCL3",
  "CCL4",
  "CCL5",
  "CCL7",
  "CCL8",
  "CXCL1",
  "CXCL2",
  "CXCL3",
  "CXCL8",
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "CXCL12",
  "CXCL13",
  "CXCL16"
)


ecm_stromal_genes <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL5A2",
  "COL6A1",
  "COL6A2",
  "COL6A3",
  "FN1",
  "DCN",
  "LUM",
  "SPARC",
  "SPARCL1",
  "FBN1",
  "FBN2",
  "POSTN",
  "VCAN",
  "MMP2",
  "MMP9",
  "TIMP1",
  "TIMP2"
)


# ------------------------------------------------------------
# 16.2 — Function to keep genes actually present in object
# ------------------------------------------------------------

get_present_genes <- function(seurat_object, gene_set) {
  
  present <- intersect(
    gene_set,
    rownames(seurat_object)
  )
  
  return(present)
}


# ------------------------------------------------------------
# 16.3 — Score one Seurat object
# ------------------------------------------------------------

score_inflammation_programs <- function(seurat_object, sample_id) {
  
  cat("\n--- Processing:", sample_id, "---\n")
  
  # Find genes available in this dataset
  inflammatory_present <- get_present_genes(
    seurat_object,
    inflammatory_genes
  )
  
  chemokine_present <- get_present_genes(
    seurat_object,
    chemokine_genes
  )
  
  ecm_present <- get_present_genes(
    seurat_object,
    ecm_stromal_genes
  )
  
  
  # Print number of genes available
  cat(
    "Inflammatory genes present:",
    length(inflammatory_present),
    "\n"
  )
  
  cat(
    "Chemokine genes present:",
    length(chemokine_present),
    "\n"
  )
  
  cat(
    "ECM/stromal genes present:",
    length(ecm_present),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Add module scores
  # ----------------------------------------------------------
  
  if (length(inflammatory_present) >= 3) {
    
    seurat_object <- AddModuleScore(
      object = seurat_object,
      features = list(inflammatory_present),
      name = "InflammatoryScore"
    )
    
    # Rename generated column
    colnames(seurat_object@meta.data)[
      ncol(seurat_object@meta.data)
    ] <- "InflammatoryScore"
    
  } else {
    
    seurat_object$InflammatoryScore <- NA_real_
    
  }
  
  
  if (length(chemokine_present) >= 3) {
    
    seurat_object <- AddModuleScore(
      object = seurat_object,
      features = list(chemokine_present),
      name = "ChemokineScore"
    )
    
    colnames(seurat_object@meta.data)[
      ncol(seurat_object@meta.data)
    ] <- "ChemokineScore"
    
  } else {
    
    seurat_object$ChemokineScore <- NA_real_
    
  }
  
  
  if (length(ecm_present) >= 3) {
    
    seurat_object <- AddModuleScore(
      object = seurat_object,
      features = list(ecm_present),
      name = "ECMScore"
    )
    
    colnames(seurat_object@meta.data)[
      ncol(seurat_object@meta.data)
    ] <- "ECMScore"
    
  } else {
    
    seurat_object$ECMScore <- NA_real_
    
  }
  
  
  return(seurat_object)
}


# ------------------------------------------------------------
# 16.4 — Apply scoring to the four samples
# ------------------------------------------------------------

cat("\n--- APPLYING SCORES TO ALL FOUR SAMPLES ---\n")


seurat1 <- score_inflammation_programs(
  seurat1,
  "BRI-1456"
)


seurat2 <- score_inflammation_programs(
  seurat2,
  "BRI-1457"
)


seurat3 <- score_inflammation_programs(
  seurat3,
  "BRI-1458"
)


seurat4 <- score_inflammation_programs(
  seurat4,
  "BRI-1459"
)


# ------------------------------------------------------------
# 16.5 — Build combined metadata table
# ------------------------------------------------------------

cat("\n--- BUILDING COMBINED SCORE TABLE ---\n")


score_meta1 <- seurat1@meta.data %>%
  mutate(
    sample_name = "BRI-1456",
    disease = "Healthy",
    tissue = "Omentum"
  )


score_meta2 <- seurat2@meta.data %>%
  mutate(
    sample_name = "BRI-1457",
    disease = "Unhealthy",
    tissue = "Omentum"
  )


score_meta3 <- seurat3@meta.data %>%
  mutate(
    sample_name = "BRI-1458",
    disease = "Healthy",
    tissue = "SQ"
  )


score_meta4 <- seurat4@meta.data %>%
  mutate(
    sample_name = "BRI-1459",
    disease = "Unhealthy",
    tissue = "SQ"
  )


all_score_meta <- bind_rows(
  score_meta1,
  score_meta2,
  score_meta3,
  score_meta4
)


# ------------------------------------------------------------
# 16.6 — Keep relevant variables
# ------------------------------------------------------------

score_table <- all_score_meta %>%
  select(
    sample_name,
    disease,
    tissue,
    cell_type,
    broad_cell_type,
    InflammatoryScore,
    ChemokineScore,
    ECMScore
  )


# ------------------------------------------------------------
# 16.7 — Cell-level score summary by broad cell type
# ------------------------------------------------------------

cat("\n--- SCORE SUMMARY BY SAMPLE AND BROAD CELL TYPE ---\n")


score_summary <- score_table %>%
  group_by(
    sample_name,
    disease,
    tissue,
    broad_cell_type
  ) %>%
  summarise(
    n_cells = n(),
    
    inflammatory_mean =
      mean(InflammatoryScore, na.rm = TRUE),
    
    inflammatory_median =
      median(InflammatoryScore, na.rm = TRUE),
    
    chemokine_mean =
      mean(ChemokineScore, na.rm = TRUE),
    
    chemokine_median =
      median(ChemokineScore, na.rm = TRUE),
    
    ECM_mean =
      mean(ECMScore, na.rm = TRUE),
    
    ECM_median =
      median(ECMScore, na.rm = TRUE),
    
    .groups = "drop"
  )


print(
  score_summary,
  n = Inf
)


# ------------------------------------------------------------
# 16.8 — Focus on biologically relevant populations
# ------------------------------------------------------------

cat("\n--- FOCUSED INFLAMMATION-RELATED CELL POPULATIONS ---\n")


focused_score_summary <- score_summary %>%
  filter(
    broad_cell_type %in% c(
      "Myeloid_macrophage",
      "Stromal_fibroblast",
      "Dendritic_cells",
      "Endothelial",
      "T_cells",
      "NK_cells"
    )
  ) %>%
  arrange(
    tissue,
    disease,
    broad_cell_type
  )


print(
  focused_score_summary,
  n = Inf
)


# ------------------------------------------------------------
# 16.9 — Save cell-level scores
# ------------------------------------------------------------

write.csv(
  score_table,
  "../results/inflammation_cell_level_scores.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 16.10 — Save summary tables
# ------------------------------------------------------------

write.csv(
  score_summary,
  "../results/inflammation_score_summary_by_cell_type.csv",
  row.names = FALSE
)


write.csv(
  focused_score_summary,
  "../results/inflammation_score_summary_focused.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 16.11 — Save updated Seurat objects
# ------------------------------------------------------------

saveRDS(
  seurat1,
  "../results/seurat1_inflammation_scored.rds"
)


saveRDS(
  seurat2,
  "../results/seurat2_inflammation_scored.rds"
)


saveRDS(
  seurat3,
  "../results/seurat3_inflammation_scored.rds"
)


saveRDS(
  seurat4,
  "../results/seurat4_inflammation_scored.rds"
)


# ------------------------------------------------------------
# 16.12 — Final message
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "Inflammation-related scoring completed successfully.\n",
  "\nSaved files:\n",
  "1. inflammation_cell_level_scores.csv\n",
  "2. inflammation_score_summary_by_cell_type.csv\n",
  "3. inflammation_score_summary_focused.csv\n",
  "4. seurat1_inflammation_scored.rds\n",
  "5. seurat2_inflammation_scored.rds\n",
  "6. seurat3_inflammation_scored.rds\n",
  "7. seurat4_inflammation_scored.rds\n",
  "============================================================\n"
)