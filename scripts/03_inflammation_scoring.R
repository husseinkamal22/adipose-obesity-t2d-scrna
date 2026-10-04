# ============================================================
# 03 — INFLAMMATION-RELATED CELL-STATE SCORING
# ============================================================

library(Seurat)
library(dplyr)


# ============================================================
# 01 — LOAD ANNOTATED SEURAT OBJECTS
# ============================================================

seurat1 <- readRDS(
  here::here("results", "seurat1_annotated.rds")
)

seurat2 <- readRDS(
  here::here("results", "seurat2_annotated.rds")
)

seurat3 <- readRDS(
  here::here("results", "seurat3_annotated.rds")
)

seurat4 <- readRDS(
  here::here("results", "seurat4_annotated.rds")
)


# ============================================================
# 02 — ADD SAMPLE INFORMATION
# ============================================================

seurat1$sample_name <- "BRI-1456"
seurat1$disease <- "Healthy"
seurat1$tissue <- "SQ"

seurat2$sample_name <- "BRI-1457"
seurat2$disease <- "Unhealthy"
seurat2$tissue <- "SQ"

seurat3$sample_name <- "BRI-1458"
seurat3$disease <- "Healthy"
seurat3$tissue <- "Omentum"

seurat4$sample_name <- "BRI-1459"
seurat4$disease <- "Unhealthy"
seurat4$tissue <- "Omentum"


# ============================================================
# 03 — DEFINE BROAD CELL TYPES
# ============================================================

make_broad_cell_type <- function(seurat_object) {
  
  seurat_object$broad_cell_type <- case_when(
    
    seurat_object$cell_type %in% c(
      "Stromal",
      "Stromal_mesenchymal",
      "Fibroblast_stromal"
    ) ~ "Stromal_fibroblast",
    
    seurat_object$cell_type %in% c(
      "Myeloid_monocyte_macrophage_like",
      "Macrophage_myeloid",
      "Inflammatory_monocyte_like"
    ) ~ "Myeloid_macrophage",
    
    seurat_object$cell_type %in% c(
      "Endothelial",
      "Endothelial_ACKR1"
    ) ~ "Endothelial",
    
    seurat_object$cell_type %in% c(
      "T_cells",
      "Cytotoxic_T_cells"
    ) ~ "T_cells",
    
    seurat_object$cell_type == "NK_cells" ~
      "NK_cells",
    
    seurat_object$cell_type %in% c(
      "DC_like",
      "pDC_like"
    ) ~ "Dendritic_cells",
    
    seurat_object$cell_type == "B_cells" ~
      "B_cells",
    
    seurat_object$cell_type ==
      "Epithelial_mesothelial_like" ~
      "Epithelial_mesothelial",
    
    seurat_object$cell_type ==
      "Vascular_smooth_muscle" ~
      "Vascular_smooth_muscle",
    
    seurat_object$cell_type ==
      "Unresolved" ~
      "Unresolved",
    
    TRUE ~ "Other"
  )
  
  return(seurat_object)
}


seurat1 <- make_broad_cell_type(seurat1)
seurat2 <- make_broad_cell_type(seurat2)
seurat3 <- make_broad_cell_type(seurat3)
seurat4 <- make_broad_cell_type(seurat4)


# ============================================================
# 04 — DEFINE BIOLOGICAL GENE PROGRAMS
# ============================================================

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


# ============================================================
# 05 — FUNCTION TO IDENTIFY AVAILABLE GENES
# ============================================================

get_present_genes <- function(seurat_object, gene_set) {
  
  intersect(
    gene_set,
    rownames(seurat_object)
  )
}


# ============================================================
# 06 — FUNCTION TO CALCULATE THREE SCORES
# ============================================================

score_programs <- function(
    seurat_object,
    sample_id
) {
  
  cat("\n--------------------------------------------\n")
  cat("Processing:", sample_id, "\n")
  cat("--------------------------------------------\n")
  
  
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
  
  
  cat(
    "Inflammatory genes:",
    length(inflammatory_present),
    "\n"
  )
  
  cat(
    "Chemokine genes:",
    length(chemokine_present),
    "\n"
  )
  
  cat(
    "ECM/stromal genes:",
    length(ecm_present),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Inflammatory score
  # ----------------------------------------------------------
  
  seurat_object <- AddModuleScore(
    object = seurat_object,
    features = list(inflammatory_present),
    name = "InflammatoryScore"
  )
  
  
  colnames(seurat_object@meta.data)[
    ncol(seurat_object@meta.data)
  ] <- "InflammatoryScore"
  
  
  # ----------------------------------------------------------
  # Chemokine score
  # ----------------------------------------------------------
  
  seurat_object <- AddModuleScore(
    object = seurat_object,
    features = list(chemokine_present),
    name = "ChemokineScore"
  )
  
  
  colnames(seurat_object@meta.data)[
    ncol(seurat_object@meta.data)
  ] <- "ChemokineScore"
  
  
  # ----------------------------------------------------------
  # ECM score
  # ----------------------------------------------------------
  
  seurat_object <- AddModuleScore(
    object = seurat_object,
    features = list(ecm_present),
    name = "ECMScore"
  )
  
  
  colnames(seurat_object@meta.data)[
    ncol(seurat_object@meta.data)
  ] <- "ECMScore"
  
  
  return(seurat_object)
}


# ============================================================
# 07 — SCORE ALL FOUR SAMPLES
# ============================================================

seurat1 <- score_programs(
  seurat1,
  "BRI-1456"
)

seurat2 <- score_programs(
  seurat2,
  "BRI-1457"
)

seurat3 <- score_programs(
  seurat3,
  "BRI-1458"
)

seurat4 <- score_programs(
  seurat4,
  "BRI-1459"
)


# ============================================================
# 08 — COMBINE METADATA
# ============================================================

score_meta1 <- seurat1@meta.data
score_meta2 <- seurat2@meta.data
score_meta3 <- seurat3@meta.data
score_meta4 <- seurat4@meta.data


all_score_meta <- bind_rows(
  score_meta1,
  score_meta2,
  score_meta3,
  score_meta4
)


# ============================================================
# 09 — CREATE SCORE TABLE
# ============================================================

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


# ============================================================
# 10 — SUMMARIZE SCORES BY CELL TYPE
# ============================================================

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
      mean(
        InflammatoryScore,
        na.rm = TRUE
      ),
    
    inflammatory_median =
      median(
        InflammatoryScore,
        na.rm = TRUE
      ),
    
    chemokine_mean =
      mean(
        ChemokineScore,
        na.rm = TRUE
      ),
    
    chemokine_median =
      median(
        ChemokineScore,
        na.rm = TRUE
      ),
    
    ECM_mean =
      mean(
        ECMScore,
        na.rm = TRUE
      ),
    
    ECM_median =
      median(
        ECMScore,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ============================================================
# 11 — FOCUS ON RELEVANT CELL POPULATIONS
# ============================================================

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


# ============================================================
# 12 — PRINT RESULTS
# ============================================================

cat("\n============================================================\n")
cat("SCORE SUMMARY BY SAMPLE AND BROAD CELL TYPE\n")
cat("============================================================\n")


print(
  score_summary,
  n = Inf
)


cat("\n============================================================\n")
cat("FOCUSED INFLAMMATION-RELATED CELL POPULATIONS\n")
cat("============================================================\n")


print(
  focused_score_summary,
  n = Inf
)


# ============================================================
# 13 — SAVE RESULTS
# ============================================================

write.csv(
  score_table,
  here::here("results", "inflammation_cell_level_scores.csv"),
  row.names = FALSE
)


write.csv(
  score_summary,
  here::here("results", "inflammation_score_summary_by_cell_type.csv"),
  row.names = FALSE
)


write.csv(
  focused_score_summary,
  here::here("results", "inflammation_score_summary_focused.csv"),
  row.names = FALSE
)


# ============================================================
# 14 — SAVE SCORED SEURAT OBJECTS
# ============================================================

saveRDS(
  seurat1,
  here::here("results", "seurat1_inflammation_scored.rds")
)

saveRDS(
  seurat2,
  here::here("results", "seurat2_inflammation_scored.rds")
)

saveRDS(
  seurat3,
  here::here("results", "seurat3_inflammation_scored.rds")
)

saveRDS(
  seurat4,
  here::here("results", "seurat4_inflammation_scored.rds")
)


# ============================================================
# 15 — FINISHED
# ============================================================

cat("\n============================================================\n")
cat("INFLAMMATION SCORING COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nResults saved in results/\n")