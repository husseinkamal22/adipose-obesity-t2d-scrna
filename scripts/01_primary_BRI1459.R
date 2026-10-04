# ============================================================
# 01 — PRIMARY PROCESSING
# BRI-1459
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

sample4_path <- here::here("data", "GSM8548238_BRI-1459")

counts4 <- Read10X(
  data.dir = sample4_path
)

seurat4 <- CreateSeuratObject(
  counts = counts4,
  project = "GSM8548238_BRI-1459"
)

# ============================================================
# 02 — QC
# ============================================================

percent.mt <- PercentageFeatureSet(
  seurat4,
  pattern = "^MT-"
)

seurat4 <- subset(
  seurat4,
  subset = nFeature_RNA > 100 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

# ============================================================
# 03 — DOUBLETS
# ============================================================

sce_d4 <- as.SingleCellExperiment(seurat4)

set.seed(100)

sce_d4 <- scDblFinder(sce_d4)

seurat4$doublet_score <- colData(sce_d4)$scDblFinder.score

seurat4$doublet_class <- colData(sce_d4)$scDblFinder.class

seurat4 <- subset(
  seurat4,
  subset = doublet_class == "singlet"
)

# ============================================================
# 04 — NORMALIZATION
# ============================================================

seurat4 <- NormalizeData(seurat4)

# ============================================================
# 05 — VARIABLE FEATURES
# ============================================================

seurat4 <- FindVariableFeatures(
  seurat4,
  nfeatures = 2000
)

# ============================================================
# 06 — SCALING
# ============================================================

seurat4 <- ScaleData(
  seurat4,
  features = VariableFeatures(seurat4)
)

# ============================================================
# 07 — PCA
# ============================================================

seurat4 <- RunPCA(
  seurat4,
  features = VariableFeatures(seurat4)
)

# ============================================================
# 08 — NEIGHBORS
# ============================================================

seurat4 <- FindNeighbors(
  seurat4,
  dims = 1:10
)

# ============================================================
# 09 — CLUSTERING
# ============================================================

seurat4 <- FindClusters(
  seurat4,
  resolution = 0.5
)

# ============================================================
# 10 — UMAP
# ============================================================

seurat4 <- RunUMAP(
  seurat4,
  dims = 1:10
)

# ============================================================
# 11 — MARKERS
# ============================================================

markers4 <- FindAllMarkers(
  seurat4,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers4,
  file = here::here("results", "cluster_markers_BRI1459.csv"),
  row.names = FALSE
)

# ============================================================
# 12 — TOP 10 MARKERS
# ============================================================

top_markers4 <- markers4 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

write.csv(
  top_markers4,
  file = here::here("results", "top_markers_per_cluster_BRI1459.csv"),
  row.names = FALSE
)

# ============================================================
# 13 — SINGLE-R
# ============================================================

hpca_ref <- HumanPrimaryCellAtlasData()

pred4 <- SingleR(
  test = GetAssayData(
    seurat4,
    assay = "RNA",
    layer = "data"
  ),
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

seurat4$SingleR_label <- pred4$labels

# ============================================================
# 14 — CLUSTER VS SINGLER
# ============================================================

singleR_cluster_table4 <- table(
  seurat4$seurat_clusters,
  seurat4$SingleR_label
)

write.csv(
  as.data.frame(singleR_cluster_table4),
  file = here::here("results", "cluster_vs_SingleR_BRI1459.csv"),
  row.names = FALSE
)

# ============================================================
# 15 — SAVE
# ============================================================

saveRDS(
  seurat4,
  file = here::here("results", "seurat4_primary_processed.rds")
)