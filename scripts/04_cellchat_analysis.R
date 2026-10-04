library(Seurat)
library(CellChat)
library(dplyr)

seurat1 <- readRDS(here::here("results", "seurat1_inflammation_scored.rds"))
seurat2 <- readRDS(here::here("results", "seurat2_inflammation_scored.rds"))
seurat3 <- readRDS(here::here("results", "seurat3_inflammation_scored.rds"))
seurat4 <- readRDS(here::here("results", "seurat4_inflammation_scored.rds"))

seurat1_cc <- subset(seurat1, subset = broad_cell_type != "Unresolved")
seurat2_cc <- subset(seurat2, subset = broad_cell_type != "Unresolved")
seurat3_cc <- subset(seurat3, subset = broad_cell_type != "Unresolved")
seurat4_cc <- subset(seurat4, subset = broad_cell_type != "Unresolved")

create_cellchat_object <- function(seurat_object, sample_id) {
  
  seurat_object$samples <- sample_id
  
  cellchat <- createCellChat(
    object = seurat_object,
    group.by = "broad_cell_type"
  )
  
  cellchat@DB <- CellChatDB.human
  
  return(cellchat)
}

cellchat1 <- create_cellchat_object(seurat1_cc, "BRI-1456")
cellchat2 <- create_cellchat_object(seurat2_cc, "BRI-1457")
cellchat3 <- create_cellchat_object(seurat3_cc, "BRI-1458")
cellchat4 <- create_cellchat_object(seurat4_cc, "BRI-1459")

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
cellchat1 <- run_cellchat(cellchat1)