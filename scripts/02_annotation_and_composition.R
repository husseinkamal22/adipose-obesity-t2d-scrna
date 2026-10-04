# ============================================================
# 02 — ANNOTATION AND CELL COMPOSITION
# GSE278526 adipose tissue scRNA-seq
#
# Samples (tissue per official NCBI GEO GSE278526 sample records):
# BRI-1456 = Healthy Omentum
# BRI-1457 = Unhealthy Omentum
# BRI-1458 = Healthy SQ
# BRI-1459 = Unhealthy SQ
# ============================================================


# ============================================================
# 01 — LOAD PACKAGES
# ============================================================

library(Seurat)
library(dplyr)


# ============================================================
# 02 — DEFINE PROJECT PATHS
# ============================================================

results_dir <- here::here("results")


# ============================================================
# 03 — LOAD PRIMARY PROCESSED OBJECTS
# ============================================================

seurat1 <- readRDS(
  file.path(results_dir, "seurat1_primary_processed.rds")
)

seurat2 <- readRDS(
  file.path(results_dir, "seurat2_primary_processed.rds")
)

seurat3 <- readRDS(
  file.path(results_dir, "seurat3_primary_processed.rds")
)

seurat4 <- readRDS(
  file.path(results_dir, "seurat4_primary_processed.rds")
)

cat("All four primary processed Seurat objects loaded.\n")


# ============================================================
# 04 — ADD SAMPLE INFORMATION
# ============================================================

seurat1$sample_id <- "BRI-1456"
seurat1$disease <- "Healthy"
seurat1$tissue <- "Omentum"

seurat2$sample_id <- "BRI-1457"
seurat2$disease <- "Unhealthy"
seurat2$tissue <- "Omentum"

seurat3$sample_id <- "BRI-1458"
seurat3$disease <- "Healthy"
seurat3$tissue <- "SQ"

seurat4$sample_id <- "BRI-1459"
seurat4$disease <- "Unhealthy"
seurat4$tissue <- "SQ"


# ============================================================
# 05 — DEFINE FINAL CLUSTER ANNOTATIONS
# ============================================================

cluster_annotations_1 <- c(
  `0` = "Stromal",
  `1` = "Stromal_mesenchymal",
  `2` = "T_cells",
  `3` = "NK_cells",
  `4` = "Unresolved",
  `5` = "Myeloid_monocyte_macrophage_like",
  `6` = "Fibroblast_stromal",
  `7` = "Vascular_smooth_muscle",
  `8` = "Endothelial",
  `9` = "DC_like",
  `10` = "pDC_like"
)

cluster_annotations_2 <- c(
  `0` = "Stromal",
  `1` = "Fibroblast_stromal",
  `2` = "Unresolved",
  `3` = "Endothelial",
  `4` = "T_cells",
  `5` = "NK_cells",
  `6` = "DC_like",
  `7` = "Stromal_mesenchymal",
  `8` = "Vascular_smooth_muscle",
  `9` = "Inflammatory_monocyte_like",
  `10` = "Myeloid_monocyte_macrophage_like",
  `11` = "Endothelial_ACKR1",
  `12` = "pDC_like"
)

cluster_annotations_3 <- c(
  `0` = "Unresolved",
  `1` = "Stromal",
  `2` = "T_cells",
  `3` = "Fibroblast_stromal",
  `4` = "NK_cells",
  `5` = "T_cells",
  `6` = "Inflammatory_monocyte_like",
  `7` = "Epithelial_mesothelial_like",
  `8` = "Myeloid_monocyte_macrophage_like",
  `9` = "Endothelial",
  `10` = "Vascular_smooth_muscle",
  `11` = "B_cells",
  `12` = "Stromal_mesenchymal",
  `13` = "Macrophage_myeloid",
  `14` = "Epithelial_mesothelial_like",
  `15` = "pDC_like"
)

cluster_annotations_4 <- c(
  `0` = "Unresolved",
  `1` = "Stromal",
  `2` = "Epithelial_mesothelial_like",
  `3` = "T_cells",
  `4` = "Myeloid_monocyte_macrophage_like",
  `5` = "T_cells",
  `6` = "Cytotoxic_T_cells",
  `7` = "Stromal_mesenchymal",
  `8` = "NK_cells",
  `9` = "Fibroblast_stromal",
  `10` = "Endothelial",
  `11` = "B_cells",
  `12` = "Epithelial_mesothelial_like",
  `13` = "Vascular_smooth_muscle"
)


# ============================================================
# 06 — APPLY CLUSTER ANNOTATIONS
# ============================================================

apply_cluster_annotation <- function(seurat_object, annotations) {
  
  clusters <- as.character(Idents(seurat_object))
  
  cell_types <- unname(
    annotations[clusters]
  )
  
  names(cell_types) <- colnames(seurat_object)
  
  seurat_object$cell_type <- cell_types
  
  return(seurat_object)
}

seurat1 <- apply_cluster_annotation(
  seurat1,
  cluster_annotations_1
)

seurat2 <- apply_cluster_annotation(
  seurat2,
  cluster_annotations_2
)

seurat3 <- apply_cluster_annotation(
  seurat3,
  cluster_annotations_3
)

seurat4 <- apply_cluster_annotation(
  seurat4,
  cluster_annotations_4
)


# ============================================================
# 07 — DEFINE BROAD CELL TYPES
# ============================================================

define_broad_cell_type <- function(cell_type) {
  
  case_when(
    cell_type %in% c(
      "Stromal",
      "Stromal_mesenchymal",
      "Fibroblast_stromal"
    ) ~ "Stromal_fibroblast",
    
    cell_type %in% c(
      "Myeloid_monocyte_macrophage_like",
      "Macrophage_myeloid",
      "Inflammatory_monocyte_like"
    ) ~ "Myeloid_macrophage",
    
    cell_type %in% c(
      "Endothelial",
      "Endothelial_ACKR1"
    ) ~ "Endothelial",
    
    cell_type %in% c(
      "T_cells",
      "Cytotoxic_T_cells"
    ) ~ "T_cells",
    
    cell_type == "NK_cells" ~ "NK_cells",
    
    cell_type %in% c(
      "DC_like",
      "pDC_like"
    ) ~ "Dendritic_cells",
    
    cell_type == "B_cells" ~ "B_cells",
    
    cell_type == "Epithelial_mesothelial_like" ~
      "Epithelial_mesothelial",
    
    cell_type == "Vascular_smooth_muscle" ~
      "Vascular_smooth_muscle",
    
    cell_type == "Unresolved" ~ "Unresolved",
    
    TRUE ~ "Other"
  )
}


seurat1$broad_cell_type <-
  define_broad_cell_type(seurat1$cell_type)

seurat2$broad_cell_type <-
  define_broad_cell_type(seurat2$cell_type)

seurat3$broad_cell_type <-
  define_broad_cell_type(seurat3$cell_type)

seurat4$broad_cell_type <-
  define_broad_cell_type(seurat4$cell_type)


# ============================================================
# 08 — CHECK ANNOTATIONS
# ============================================================

cat("\n============================================\n")
cat("BRI-1456 cell types\n")
cat("============================================\n")
print(table(seurat1$cell_type))

cat("\n============================================\n")
cat("BRI-1457 cell types\n")
cat("============================================\n")
print(table(seurat2$cell_type))

cat("\n============================================\n")
cat("BRI-1458 cell types\n")
cat("============================================\n")
print(table(seurat3$cell_type))

cat("\n============================================\n")
cat("BRI-1459 cell types\n")
cat("============================================\n")
print(table(seurat4$cell_type))


# ============================================================
# 09 — BUILD CELL COMPOSITION TABLE
# ============================================================

build_composition_table <- function(seurat_object) {
  
  data.frame(
    sample_id = seurat_object$sample_id,
    disease = seurat_object$disease,
    tissue = seurat_object$tissue,
    broad_cell_type = seurat_object$broad_cell_type
  ) %>%
    count(
      sample_id,
      disease,
      tissue,
      broad_cell_type,
      name = "n_cells"
    ) %>%
    group_by(sample_id) %>%
    mutate(
      percentage = 100 * n_cells / sum(n_cells)
    ) %>%
    ungroup()
}


composition <- bind_rows(
  build_composition_table(seurat1),
  build_composition_table(seurat2),
  build_composition_table(seurat3),
  build_composition_table(seurat4)
)


# ============================================================
# 10 — SAVE CELL COMPOSITION TABLES
# ============================================================

write.csv(
  composition,
  file.path(
    results_dir,
    "broad_cell_type_composition_by_sample.csv"
  ),
  row.names = FALSE
)

counts <- composition %>%
  select(
    sample_id,
    disease,
    tissue,
    broad_cell_type,
    n_cells
  )

write.csv(
  counts,
  file.path(
    results_dir,
    "broad_cell_type_counts_by_sample.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11 — SAVE ANNOTATED SEURAT OBJECTS
# ============================================================

saveRDS(
  seurat1,
  file.path(results_dir, "seurat1_annotated.rds")
)

saveRDS(
  seurat2,
  file.path(results_dir, "seurat2_annotated.rds")
)

saveRDS(
  seurat3,
  file.path(results_dir, "seurat3_annotated.rds")
)

saveRDS(
  seurat4,
  file.path(results_dir, "seurat4_annotated.rds")
)


# ============================================================
# 12 — FINAL MESSAGE
# ============================================================

cat("\n============================================\n")
cat("Annotation and composition workflow completed.\n")
cat("============================================\n")

cat("\nSaved files:\n")
cat("- seurat1_annotated.rds\n")
cat("- seurat2_annotated.rds\n")
cat("- seurat3_annotated.rds\n")
cat("- seurat4_annotated.rds\n")
cat("- broad_cell_type_counts_by_sample.csv\n")
cat("- broad_cell_type_composition_by_sample.csv\n")