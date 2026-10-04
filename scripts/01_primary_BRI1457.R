# ============================================================
# 01 — PRIMARY PROCESSING
# BRI-1457
# ============================================================

# ============================================================
# 01 — LOAD PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(scDblFinder)

# ============================================================
# 02 — READ 10X DATA
# ============================================================

sample2_path <- here::here("data", "GSM8548236_BRI-1457")

counts2 <- Read10X(
  data.dir = sample2_path
)

seurat2 <- CreateSeuratObject(
  counts = counts2,
  project = "GSM8548236_BRI-1457"
)

# ============================================================
# 03 — QC
# ============================================================

percent.mt <- PercentageFeatureSet(
  seurat2,
  pattern = "^MT-"
)

seurat2 <- subset(
  seurat2,
  subset = nFeature_RNA > 100 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

# ============================================================
# 04 — DOUBLETS
# ============================================================

sce_d2 <- as.SingleCellExperiment(seurat2)

set.seed(100)

sce_d2 <- scDblFinder(sce_d2)

seurat2$doublet_score <- colData(sce_d2)$scDblFinder.score

seurat2$doublet_class <- colData(sce_d2)$scDblFinder.class

seurat2 <- subset(
  seurat2,
  subset = doublet_class == "singlet"
)

# ============================================================
# 05 — NORMALIZATION
# ============================================================

seurat2 <- NormalizeData(seurat2)

# ============================================================
# 06 — VARIABLE FEATURES
# ============================================================

seurat2 <- FindVariableFeatures(
  seurat2,
  nfeatures = 2000
)

# ============================================================
# 07 — SCALING
# ============================================================

seurat2 <- ScaleData(
  seurat2,
  features = VariableFeatures(seurat2)
)

# ============================================================
# 08 — PCA
# ============================================================

seurat2 <- RunPCA(
  seurat2,
  features = VariableFeatures(seurat2)
)

# ============================================================
# 09 — NEIGHBORS
# ============================================================

seurat2 <- FindNeighbors(
  seurat2,
  dims = 1:10
)

# ============================================================
# 10 — CLUSTERING
# ============================================================

seurat2 <- FindClusters(
  seurat2,
  resolution = 0.5
)

# ============================================================
# 11 — UMAP
# ============================================================

seurat2 <- RunUMAP(
  seurat2,
  dims = 1:10
)

# ============================================================
# 12 — CLUSTER MARKERS
# ============================================================

markers2 <- FindAllMarkers(
  seurat2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers2,
  file = here::here("results", "cluster_markers_BRI1457.csv"),
  row.names = FALSE
)

# ============================================================
# 13 — TOP 10 MARKERS PER CLUSTER
# ============================================================

top_markers2 <- markers2 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

write.csv(
  top_markers2,
  file = here::here("results", "top_markers_per_cluster_BRI1457.csv"),
  row.names = FALSE
)

# ============================================================
# 14 — SINGLE-R ANNOTATION
# ============================================================

library(SingleR)
library(celldex)

hpca_ref <- HumanPrimaryCellAtlasData()

pred2 <- SingleR(
  test = GetAssayData(
    seurat2,
    assay = "RNA",
    layer = "data"
  ),
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

seurat2$SingleR_label <- pred2$labels

# ============================================================
# 15 — CLUSTER VS SINGLER
# ============================================================

singleR_cluster_table2 <- table(
  seurat2$seurat_clusters,
  seurat2$SingleR_label
)

write.csv(
  as.data.frame(singleR_cluster_table2),
  file = here::here("results", "cluster_vs_SingleR_BRI1457.csv"),
  row.names = FALSE
)

# ============================================================
# 16 — SAVE PROCESSED OBJECT
# ============================================================

saveRDS(
  seurat2,
  file = here::here("results", "seurat2_primary_processed.rds")
)