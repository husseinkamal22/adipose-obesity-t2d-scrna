# ============================================================
# 01 — PRIMARY PROCESSING
# BRI-1458
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)
library(SingleR)
library(celldex)

# ============================================================
# 01 — READ 10X DATA
# ============================================================

sample3_path <- here::here("data", "GSM8548237_BRI-1458")

counts3 <- Read10X(
  data.dir = sample3_path
)

seurat3 <- CreateSeuratObject(
  counts = counts3,
  project = "GSM8548237_BRI-1458"
)

# ============================================================
# 02 — QC
# ============================================================

percent.mt <- PercentageFeatureSet(
  seurat3,
  pattern = "^MT-"
)

seurat3 <- subset(
  seurat3,
  subset = nFeature_RNA > 100 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

# ============================================================
# 03 — DOUBLETS
# ============================================================

sce_d3 <- as.SingleCellExperiment(seurat3)

set.seed(100)

sce_d3 <- scDblFinder(sce_d3)

seurat3$doublet_score <- colData(sce_d3)$scDblFinder.score

seurat3$doublet_class <- colData(sce_d3)$scDblFinder.class

seurat3 <- subset(
  seurat3,
  subset = doublet_class == "singlet"
)

# ============================================================
# 04 — NORMALIZATION
# ============================================================

seurat3 <- NormalizeData(seurat3)

# ============================================================
# 05 — VARIABLE FEATURES
# ============================================================

seurat3 <- FindVariableFeatures(
  seurat3,
  nfeatures = 2000
)

# ============================================================
# 06 — SCALING
# ============================================================

seurat3 <- ScaleData(
  seurat3,
  features = VariableFeatures(seurat3)
)

# ============================================================
# 07 — PCA
# ============================================================

seurat3 <- RunPCA(
  seurat3,
  features = VariableFeatures(seurat3)
)

# ============================================================
# 08 — NEIGHBORS
# ============================================================

seurat3 <- FindNeighbors(
  seurat3,
  dims = 1:10
)

# ============================================================
# 09 — CLUSTERING
# ============================================================

seurat3 <- FindClusters(
  seurat3,
  resolution = 0.5
)

# ============================================================
# 10 — UMAP
# ============================================================

seurat3 <- RunUMAP(
  seurat3,
  dims = 1:10
)

# ============================================================
# 11 — MARKERS
# ============================================================

markers3 <- FindAllMarkers(
  seurat3,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers3,
  file = here::here("results", "cluster_markers_BRI1458.csv"),
  row.names = FALSE
)

# ============================================================
# 12 — TOP 10 MARKERS
# ============================================================

top_markers3 <- markers3 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

write.csv(
  top_markers3,
  file = here::here("results", "top_markers_per_cluster_BRI1458.csv"),
  row.names = FALSE
)

# ============================================================
# 13 — SINGLE-R
# ============================================================

hpca_ref <- HumanPrimaryCellAtlasData()

pred3 <- SingleR(
  test = GetAssayData(
    seurat3,
    assay = "RNA",
    layer = "data"
  ),
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

seurat3$SingleR_label <- pred3$labels

# ============================================================
# 14 — CLUSTER VS SINGLER
# ============================================================

singleR_cluster_table3 <- table(
  seurat3$seurat_clusters,
  seurat3$SingleR_label
)

write.csv(
  as.data.frame(singleR_cluster_table3),
  file = here::here("results", "cluster_vs_SingleR_BRI1458.csv"),
  row.names = FALSE
)

# ============================================================
# 15 — SAVE
# ============================================================

saveRDS(
  seurat3,
  file = here::here("results", "seurat3_primary_processed.rds")
)