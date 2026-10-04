# ============================================================
# 04_cellchat_analysis.R
# CellChat analysis of adipose tissue scRNA-seq
#
# Samples:
# BRI-1456 = Healthy Omentum
# BRI-1457 = Unhealthy Omentum
# BRI-1458 = Healthy SQ
# BRI-1459 = Unhealthy SQ
#
# Tissue labels follow the official NCBI GEO GSE278526 sample
# records. The deposited per-cell metadata file
# (GSE278526_cell_barcodes_metadata.tsv.gz) carries the opposite
# assignment; that file is not tracked and is not modified here.
#
# Purpose:
# Infer potential cell-cell communication patterns
# using CellChat.
#
# Important:
# Results are descriptive/inferential predictions based on
# ligand-receptor expression. They are NOT causal evidence.
#
# NOTE — DO NOT EXECUTE DURING REPOSITORY CLEANUP.
# This script rebuilds and overwrites results/cellchat_BRI-*.rds.
# Those four objects already exist from the original run and are
# the inputs to 05_cellchat_results.R. They were deliberately not
# regenerated while preparing this repository, because doing so
# would mean re-running part of the analysis. This file documents
# the workflow that produced them; it is not to be run as part of
# the cleanup.
# ============================================================

library(Seurat)
library(CellChat)
library(dplyr)

# ------------------------------------------------------------
# 1. Load inflammation-scored Seurat objects
# ------------------------------------------------------------

seurat1 <- readRDS(here::here("results", "seurat1_inflammation_scored.rds"))
seurat2 <- readRDS(here::here("results", "seurat2_inflammation_scored.rds"))
seurat3 <- readRDS(here::here("results", "seurat3_inflammation_scored.rds"))
seurat4 <- readRDS(here::here("results", "seurat4_inflammation_scored.rds"))

# ------------------------------------------------------------
# 2. Remove unresolved cells
# ------------------------------------------------------------

seurat1_cc <- subset(
  seurat1,
  subset = broad_cell_type != "Unresolved"
)

seurat2_cc <- subset(
  seurat2,
  subset = broad_cell_type != "Unresolved"
)

seurat3_cc <- subset(
  seurat3,
  subset = broad_cell_type != "Unresolved"
)

seurat4_cc <- subset(
  seurat4,
  subset = broad_cell_type != "Unresolved"
)

# ------------------------------------------------------------
# 3. Create CellChat object
# ------------------------------------------------------------

create_cellchat_object <- function(seurat_object, sample_id) {
  
  seurat_object$samples <- sample_id
  
  cellchat <- createCellChat(
    object = seurat_object,
    group.by = "broad_cell_type"
  )
  
  cellchat@DB <- CellChatDB.human
  
  return(cellchat)
}

# ------------------------------------------------------------
# 4. Create CellChat objects
# ------------------------------------------------------------

cellchat1 <- create_cellchat_object(
  seurat1_cc,
  "BRI-1456"
)

cellchat2 <- create_cellchat_object(
  seurat2_cc,
  "BRI-1457"
)

cellchat3 <- create_cellchat_object(
  seurat3_cc,
  "BRI-1458"
)

cellchat4 <- create_cellchat_object(
  seurat4_cc,
  "BRI-1459"
)

# ------------------------------------------------------------
# 5. CellChat analysis function
# ------------------------------------------------------------

run_cellchat <- function(cellchat) {
  
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# ------------------------------------------------------------
# 6. Run CellChat
# ------------------------------------------------------------

cellchat1 <- run_cellchat(cellchat1)

saveRDS(
  cellchat1,
  here::here("results", "cellchat_BRI-1456.rds")
)

cellchat2 <- run_cellchat(cellchat2)

saveRDS(
  cellchat2,
  here::here("results", "cellchat_BRI-1457.rds")
)

cellchat3 <- run_cellchat(cellchat3)

saveRDS(
  cellchat3,
  here::here("results", "cellchat_BRI-1458.rds")
)

cellchat4 <- run_cellchat(cellchat4)

saveRDS(
  cellchat4,
  here::here("results", "cellchat_BRI-1459.rds")
)

cat("\nAll four CellChat analyses completed and saved.\n")
