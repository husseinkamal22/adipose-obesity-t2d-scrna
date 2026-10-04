# ============================================================
# 02 — DOWNSTREAM ANALYSIS
# ============================================================

# RESEARCH QUESTION
#
# Which adipose cell populations and intercellular signaling
# programs are associated with the inflammatory state of
# unhealthy obesity/metabolic disease?
#
# Current analysis:
# BRI-1456 — cell-type characterization and communication
# pipeline development.
#
# Future comparative analysis:
# Healthy vs Unhealthy obesity/metabolic disease.
# ============================================================


library(Seurat)
library(dplyr)
library(ggplot2)

seurat1 <- readRDS(
  "../results/seurat1_primary_processed.rds"
)

seurat1
# ============================================================
# 02 — CHECK ANNOTATED CELL TYPES
# ============================================================

table(seurat1$cell_type)
# ============================================================
# 03 — CELL TYPE COMPOSITION
# ============================================================

cell_composition <- as.data.frame(
  prop.table(table(seurat1$cell_type)) * 100
)

colnames(cell_composition) <- c(
  "cell_type",
  "percentage"
)

cell_composition
# ============================================================
# 04 — CELL TYPE COMPOSITION PLOT
# ============================================================

ggplot(
  cell_composition,
  aes(
    x = reorder(cell_type, percentage),
    y = percentage
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Cell Type Composition",
    x = "Cell Type",
    y = "Percentage of Cells"
  ) +
  theme_minimal()

ggsave(
  "../results/figure_cell_type_composition.png",
  width = 8,
  height = 6,
  dpi = 300
)
# ============================================================
# 05 — DEFINE BROAD CELL COMPARTMENTS
# ============================================================

seurat1$cell_compartment <- case_when(
  seurat1$cell_type %in% c(
    "Stromal",
    "Stromal_mesenchymal",
    "Fibroblast_stromal"
  ) ~ "Stromal",
  
  seurat1$cell_type %in% c(
    "Myeloid_monocyte_macrophage_like",
    "cDC1_like",
    "pDC_like"
  ) ~ "Myeloid",
  
  seurat1$cell_type %in% c(
    "T_cells",
    "NK_cells"
  ) ~ "Lymphoid",
  
  seurat1$cell_type == "Endothelial" ~ "Endothelial",
  
  seurat1$cell_type == "Vascular_smooth_muscle" ~ "Vascular_smooth_muscle",
  
  seurat1$cell_type == "Low_complexity_unresolved" ~ "Unresolved"
)
# ============================================================
# 06 — BROAD CELL COMPARTMENT COMPOSITION
# ============================================================

compartment_composition <- as.data.frame(
  prop.table(table(seurat1$cell_compartment)) * 100
)

colnames(compartment_composition) <- c(
  "cell_compartment",
  "percentage"
)

compartment_composition
table(seurat1$cell_compartment)
# ============================================================
# 07 — SAVE BROAD COMPARTMENT SUMMARY
# ============================================================

write.csv(
  compartment_composition,
  file = "../results/broad_compartment_composition.csv",
  row.names = FALSE
)
# ============================================================
# CELLCHAT — STEP 1: LOAD PACKAGES
# ============================================================

library(CellChat)
library(patchwork)
# ============================================================
# CELLCHAT — STEP 2: PREPARE CELL TYPES
# ============================================================

seurat1_cellchat <- subset(
  seurat1,
  subset = cell_type != "Low_complexity_unresolved"
)

table(seurat1_cellchat$cell_type)
# ============================================================
# CELLCHAT — STEP 3: CREATE CELLCHAT OBJECT
# ============================================================

cellchat <- createCellChat(
  object = seurat1_cellchat,
  group.by = "cell_type"
)

cellchat
# ============================================================
# CELLCHAT — STEP 4: SET HUMAN DATABASE
# ============================================================

cellchat@DB <- CellChatDB.human
# ============================================================
# CELLCHAT — STEP 5: PREPARE SIGNALING DATA
# ============================================================

cellchat <- subsetData(cellchat)
# ============================================================
# CELLCHAT — STEP 6: IDENTIFY OVEREXPRESSED GENES
# ============================================================

cellchat <- identifyOverExpressedGenes(cellchat)
# ============================================================
# CELLCHAT — STEP 7: IDENTIFY OVEREXPRESSED INTERACTIONS
# ============================================================

cellchat <- identifyOverExpressedInteractions(cellchat)
# ============================================================
# CELLCHAT — STEP 8: COMPUTE COMMUNICATION PROBABILITY
# ============================================================

cellchat <- computeCommunProb(cellchat)
# ============================================================
# CELLCHAT — STEP 9: FILTER COMMUNICATION
# ============================================================

cellchat <- filterCommunication(
  cellchat,
  min.cells = 10
)
# ============================================================
# CELLCHAT — STEP 10: COMPUTE SIGNALING PATHWAYS
# ============================================================

cellchat <- computeCommunProbPathway(cellchat)
# ============================================================
# CELLCHAT — STEP 11: AGGREGATE COMMUNICATION NETWORK
# ============================================================

cellchat <- aggregateNet(cellchat)
# ============================================================
# CELLCHAT — STEP 12: OVERALL COMMUNICATION NETWORK
# ============================================================

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE
)
# ============================================================
# CELLCHAT — STEP 13: EXTRACT COMMUNICATION RESULTS
# ============================================================

lr_results <- subsetCommunication(cellchat)

head(lr_results)

# ============================================================
# CELLCHAT — STEP 14: NUMBER OF COMMUNICATIONS
# ============================================================

netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE
)