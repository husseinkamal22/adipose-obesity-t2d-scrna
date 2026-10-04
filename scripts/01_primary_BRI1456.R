# ============================================================
# 01 — LOAD PACKAGES
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)


# ============================================================
# 02 — SET DATA PATH
# ============================================================

data_path <- here::here("data")

list.files(data_path)


# ============================================================
# 03 — SELECT SAMPLE 1
# ============================================================

sample1_path <- file.path(
  data_path,
  "GSM8548235_BRI-1456"
)

list.files(sample1_path)


# ============================================================
# 04 — READ 10X COUNT MATRIX
# ============================================================

counts1 <- Read10X(
  data.dir = sample1_path
)

dim(counts1)
class(counts1)


# ============================================================
# 05 — CREATE SEURAT OBJECT
# ============================================================

seurat1 <- CreateSeuratObject(
  counts = counts1,
  project = "GSM8548235_BRI-1456"
)

seurat1

# ============================================================
# 06 — INSPECT CELL METADATA
# ============================================================

metadata_file <- file.path(
  data_path,
  "GSE278526_cell_barcodes_metadata.tsv.gz"
)

metadata <- read.delim(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

dim(metadata)
head(metadata)
colnames(metadata)

# ============================================================
# 07 — INSPECT METADATA SAMPLE DISTRIBUTION
# ============================================================

table(metadata$sample_name)

table(metadata$disease)

table(metadata$tissue)

table(metadata$status)

# ============================================================
# 08 — CHECK BARCODE MATCHING
# ============================================================

head(colnames(seurat1))

head(metadata$barcode)

sum(colnames(seurat1) %in% metadata$barcode)

# ============================================================
# 09 — CHECK BARCODE FORMAT
# ============================================================

seurat1_barcodes <- gsub(
  "-1$",
  "",
  colnames(seurat1)
)

head(seurat1_barcodes)

sum(seurat1_barcodes %in% metadata$barcode)

# ============================================================
# 10 — CHECK SAMPLE-SPECIFIC BARCODE MATCHING
# ============================================================

metadata_bri1456 <- metadata[
  metadata$sample_name == "BRI-1456",
]

sum(
  seurat1_barcodes %in% metadata_bri1456$barcode
)

dim(metadata_bri1456)

head(metadata_bri1456$barcode)

# ============================================================
# 11 — QC: MITOCHONDRIAL PERCENTAGE
# ============================================================

seurat1[["percent.mt"]] <- PercentageFeatureSet(
  seurat1,
  pattern = "^MT-"
)

head(seurat1@meta.data)

# ============================================================
# 12 — QC: VISUALIZE QUALITY METRICS
# ============================================================

VlnPlot(
  seurat1,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)
FeatureScatter(
  seurat1,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

FeatureScatter(
  seurat1,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

# ============================================================
# 13 — QC: SUMMARY OF QUALITY METRICS
# ============================================================

summary(seurat1$nFeature_RNA)

summary(seurat1$nCount_RNA)

summary(seurat1$percent.mt)

sum(is.na(seurat1$nFeature_RNA))

sum(is.na(seurat1$nCount_RNA))

sum(is.na(seurat1$percent.mt))
# ============================================================
# 14 — QC: INSPECT CELLS WITH DETECTED RNA
# ============================================================

sum(seurat1$nFeature_RNA > 0)

sum(seurat1$nFeature_RNA > 50)

sum(seurat1$nFeature_RNA > 100)

sum(seurat1$nFeature_RNA > 200)

sum(seurat1$nFeature_RNA > 500)

sum(seurat1$nFeature_RNA > 1000)
# ============================================================
# 15 — QC: APPLY FILTERING
# ============================================================

seurat1 <- subset(
  seurat1,
  subset =
    nFeature_RNA > 100 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

seurat1

# ============================================================
# 16 — DOUBLETS
# ============================================================

library(SingleCellExperiment)
library(scDblFinder)

sce_d1 <- as.SingleCellExperiment(seurat1)

set.seed(100)

sce_d1 <- scDblFinder(sce_d1)

seurat1$doublet_score <- colData(sce_d1)$scDblFinder.score

seurat1$doublet_class <- colData(sce_d1)$scDblFinder.class

seurat1 <- subset(
  seurat1,
  subset = doublet_class == "singlet"
)

seurat1

# ============================================================
# 17 — NORMALIZATION
# ============================================================

seurat1 <- NormalizeData(
  seurat1
)

# ============================================================
# 16 — FIND VARIABLE FEATURES
# ============================================================

seurat1 <- FindVariableFeatures(
  seurat1,
  nfeatures = 2000
)
# ============================================================
# 17 — SCALE DATA
# ============================================================

seurat1 <- ScaleData(
  seurat1,
  features = VariableFeatures(seurat1)
)
# ============================================================
# 18 — PCA
# ============================================================

seurat1 <- RunPCA(
  seurat1,
  features = VariableFeatures(seurat1)
)
# ============================================================
# 19 — FIND NEIGHBORS
# ============================================================

seurat1 <- FindNeighbors(
  seurat1,
  dims = 1:10
)
# ============================================================
# 20 — CLUSTERING
# ============================================================

seurat1 <- FindClusters(
  seurat1,
  resolution = 0.5
)

table(seurat1$seurat_clusters)

# ============================================================
# 21 — UMAP
# ============================================================

seurat1 <- RunUMAP(
  seurat1,
  dims = 1:10
)

DimPlot(
  seurat1,
  reduction = "umap",
  label = TRUE
)
# ============================================================
# 22 — FIND CLUSTER MARKERS
# ============================================================

markers <- FindAllMarkers(
  seurat1,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(markers)
# ============================================================
# 23 — SAVE MARKER RESULTS
# ============================================================

write.csv(
  markers,
  file = here::here("results", "cluster_markers.csv"),
  row.names = FALSE
)
# ============================================================
# 24 — TOP MARKERS PER CLUSTER
# ============================================================

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

top_markers
# ============================================================
# 25 — SAVE TOP MARKERS
# ============================================================

write.csv(
  top_markers,
  file = here::here("results", "top_markers_per_cluster.csv"),
  row.names = FALSE
)
# ============================================================
# 26 — SINGLE-R CELL ANNOTATION
# ============================================================

library(SingleR)
library(celldex)

hpca_ref <- HumanPrimaryCellAtlasData()

pred <- SingleR(
  test = GetAssayData(
    seurat1,
    assay = "RNA",
    layer = "data"
  ),
  ref = hpca_ref,
  labels = hpca_ref$label.main
)

table(pred$labels)
# ============================================================
# 27 — ADD SINGLE-R LABELS TO SEURAT
# ============================================================

seurat1$SingleR_label <- pred$labels

table(seurat1$SingleR_label)
# ============================================================
# 28 — COMPARE CLUSTERS WITH SINGLER ANNOTATION
# ============================================================

table(
  seurat1$seurat_clusters,
  seurat1$SingleR_label
)
# ============================================================
# 29 — FINAL CELL TYPE ANNOTATION
# ============================================================

cluster_annotations <- c(
  `0` = "Stromal",
  `1` = "Stromal_mesenchymal",
  `2` = "T_cells",
  `3` = "NK_cells",
  `4` = "Low_complexity_unresolved",
  `5` = "Myeloid_monocyte_macrophage_like",
  `6` = "Fibroblast_stromal",
  `7` = "Vascular_smooth_muscle",
  `8` = "Endothelial",
  `9` = "cDC1_like",
  `10` = "pDC_like"
)

cell_types <- unname(
  cluster_annotations[
    as.character(seurat1$seurat_clusters)
  ]
)

names(cell_types) <- colnames(seurat1)

seurat1$cell_type <- cell_types

table(seurat1$cell_type)
# ============================================================
# 30 — VISUALIZE CELL TYPE ANNOTATION
# ============================================================

DimPlot(
  seurat1,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)
# ============================================================
# 31 — MARKER VALIDATION
# ============================================================

DotPlot(
  seurat1,
  features = c(
    "CD3D", "CD3E",
    "NKG7", "GNLY",
    "FCN1", "C1QB",
    "COL1A1", "COL3A1",
    "MYH11",
    "PTPRB", "ERG",
    "CLEC9A",
    "LILRA4", "CLEC4C"
  ),
  group.by = "cell_type"
) +
  RotatedAxis()
# ============================================================
# 32 — CELL TYPE SUMMARY
# ============================================================

cell_type_summary <- as.data.frame(
  table(seurat1$cell_type)
)

colnames(cell_type_summary) <- c(
  "cell_type",
  "n_cells"
)

cell_type_summary$percentage <- (
  cell_type_summary$n_cells /
    sum(cell_type_summary$n_cells)
) * 100
# ============================================================
# 33 — SAVE PROCESSED SEURAT OBJECT
# ============================================================

saveRDS(
  seurat1,
  file = here::here("results", "seurat1_primary_processed.rds")
)
cell_type_summary

write.csv(
  cell_type_summary,
  file = here::here("results", "cell_type_summary.csv"),
  row.names = FALSE
)