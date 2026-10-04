# Single-cell analysis of human adipose tissue (GSE278526)

Reanalysis of publicly available single-cell RNA-seq data from human
subcutaneous (SQ) and omental adipose tissue, comparing samples collected from
metabolically healthy and metabolically unhealthy donors.

This repository contains the **processing, annotation, descriptive scoring and
ligand–receptor analysis** used to produce the figures and result tables in
`results/`. It is a portfolio / methods-reproducibility project.

> **Scope statement.** No biological analysis was rerun to prepare this
> repository, and no numeric result was altered. The changes made were file
> organisation, path portability, documentation, and — per §1 — correction of
> the tissue and composite sample labels to match the official GEO records.

---

## 1. Data

| Item | Value |
| ---- | ----- |
| GEO accession | [GSE278526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278526) |
| Platform | 10x Genomics scRNA-seq |
| Samples | 4 (one per condition × tissue) |
| Raw files | `GSE278526_RAW.tar` (not tracked in Git — download from GEO) |
| Per-cell metadata | `GSE278526_cell_barcodes_metadata.tsv.gz` (not tracked in Git) |
| Barcodes before QC | 43,887 |
| Cells after QC | 32,851 |

### Sample sheet (as used in this analysis)

| Sample ID | GEO sample | Condition | Tissue |
| --------- | ---------- | --------- | ------ |
| BRI-1456 | GSM8548235 | Healthy | Omentum |
| BRI-1457 | GSM8548236 | Unhealthy | Omentum |
| BRI-1458 | GSM8548237 | Healthy | SQ |
| BRI-1459 | GSM8548238 | Unhealthy | SQ |

> **Tissue labelling.** Tissue labels above are taken from the official
> NCBI GEO sample records for GSE278526 (`!Sample_source_name_ch1`: BRI-1456
> and BRI-1457 are omental, BRI-1458 and BRI-1459 are subcutaneous abdominal).
> The deposited per-cell metadata file
> `GSE278526_cell_barcodes_metadata.tsv.gz` carries the **opposite**
> assignment; that file is excluded by `.gitignore`, is not modified by this
> repository, and was not used to label the output tables. Every result table
> in `results/` and every figure label follows the GEO records. Condition
> (Healthy / Unhealthy) is assigned identically in both sources and is
> unaffected. Relabelling touched tissue and composite sample labels only —
> **no numeric value in any result table was changed**, and no analysis was
> re-executed.

### Design limitation

There is **one sample per condition × tissue combination**. Tissue and
condition are therefore not separable from donor-level variation, there is no
biological replication within any group, and no inferential statistics are
reported. All comparisons below are **descriptive**.

---

## 2. Pipeline

Scripts in `scripts/` are numbered in execution order. Nothing in this section
was re-executed while preparing the repository.

```
01_primary_BRI145*.R   load 10x matrices, QC, doublet removal, clustering
02_annotation_and_composition.R   cell-type annotation + composition tables
03_inflammation_scoring.R         module scores (inflammatory / chemokine / ECM)
04_cellchat_analysis.R            build + run CellChat objects
05_cellchat_results.R             extract CellChat result tables
06_final_figures.R                Figures 1–4
```

`scripts/legacy/` holds earlier, superseded versions of the workflow. They are
**not** part of the pipeline and should not be run — see
[`scripts/legacy/README.md`](scripts/legacy/README.md).

### Processing parameters

| Step | Setting |
| ---- | ------- |
| QC filter | `nFeature_RNA > 100`, `nFeature_RNA < 6000`, `percent.mt < 10` |
| Doublet removal | `scDblFinder` |
| Normalisation | `NormalizeData` |
| Feature selection | `FindVariableFeatures(nfeatures = 2000)` |
| Scaling | `ScaleData` |
| Neighbours | `FindNeighbors(dims = 1:10)` |
| Embedding | `RunUMAP(dims = 1:10)` |
| Clustering | `FindClusters(resolution = 0.5)` |
| Marker detection | `FindAllMarkers(min.pct = 0.25, logfc.threshold = 0.25)` |

### Annotation

Cluster-level labels were assigned from marker expression using an explicit
cluster→label map in `02_annotation_and_composition.R` (functions
`cluster_annotations_*` and `define_broad_cell_type()`), with `SingleR`
predictions used as a cross-check (see `results/cluster_vs_SingleR_*.csv`).
Clusters that did not receive a confident label are retained as
**`Unresolved`** (8,115 cells, 24.7% of the dataset) rather than being
dropped or force-assigned.

### Scores

`03_inflammation_scoring.R` uses `Seurat::AddModuleScore` with three curated
gene sets to produce `InflammatoryScore`, `ChemokineScore` and `ECMScore`
per cell. Gene lists are defined at the top of that script.

### CellChat

`CellChatDB.human` → `subsetData` → `identifyOverExpressedGenes` →
`identifyOverExpressedInteractions` → `computeCommunProb(type = "triMean")`
→ `filterCommunication(min.cells = 10)` → `computeCommunProbPathway` →
`aggregateNet`.

Cells labelled `Unresolved` are **excluded** from all CellChat analyses.

---

## 3. Results

All numbers below are read directly from the files in `results/`. They are
descriptive summaries of these four samples only — no statistical testing was
performed.

### 3.1 Cell composition

32,851 cells after QC: BRI-1456 = 2,528 · BRI-1457 = 4,848 ·
BRI-1458 = 12,839 · BRI-1459 = 12,636.

| Broad cell type | Cells | Share |
| --------------- | ----: | ----: |
| Stromal_fibroblast | 10,266 | 31.3% |
| **Unresolved** | **8,115** | **24.7%** |
| T_cells | 5,385 | 16.4% |
| Myeloid_macrophage | 2,454 | 7.5% |
| Epithelial_mesothelial | 2,029 | 6.2% |
| NK_cells | 1,918 | 5.8% |
| Endothelial | 1,223 | 3.7% |
| Vascular_smooth_muscle | 570 | 1.7% |
| Dendritic_cells | 461 | 1.4% |
| B_cells | 430 | 1.3% |

Descriptively: `Epithelial_mesothelial` and `B_cells` are detected only in
BRI-1458/1459; `Dendritic_cells` are present in BRI-1456/1457/1458 but not
BRI-1459. `Stromal_fibroblast` is the largest annotated compartment in all
four samples. Full counts: `results/broad_cell_type_counts_by_sample.csv`.

### 3.2 Inflammatory / chemokine / ECM scores

Mean `InflammatoryScore` per broad cell type (from
`results/inflammation_score_summary_focused.csv`):

| Cell type | BRI-1456 | BRI-1457 | BRI-1458 | BRI-1459 |
| --------- | -------: | -------: | -------: | -------: |
| Myeloid_macrophage | 0.595 | 0.747 | 0.533 | 0.588 |
| NK_cells | 0.474 | 0.561 | 0.300 | 0.467 |
| Endothelial | 0.476 | 0.363 | 0.103 | 0.100 |
| T_cells | 0.155 | 0.212 | 0.194 | 0.230 |
| Stromal_fibroblast | 0.097 | 0.089 | −0.083 | −0.114 |
| Dendritic_cells | −0.032 | 0.352 | −0.321 | n/a |

Descriptively, `Myeloid_macrophage` and `NK_cells` show the highest mean
inflammatory scores in each sample, while `Stromal_fibroblast` shows the
highest mean `ECMScore` in every sample (0.72–0.90) and negative mean
inflammatory scores in BRI-1458/1459. `Endothelial` mean inflammatory
scores are lower in BRI-1458/1459 than in BRI-1456/1457.
These are four single samples; no group-level inference is supported.

### 3.3 Inferred cell–cell communication (CellChat)

**CellChat is a ligand–receptor co-expression framework.** The values below
are model-derived `communication_probability` scores. They describe which
ligand–receptor pairs *could* be supported by the observed expression; they
are **not** measurements of pathway activity, flux or causality.

Mean communication probability per focused pathway:

| Pathway | Healthy Omentum | Unhealthy Omentum | Healthy SQ | Unhealthy SQ |
| ------- | ---------: | -----------: | --------------: | -----------------: |
| CCL | 1.051 | 1.321 | 0.424 | 0.854 |
| CXCL | 1.723 | 2.289 | 1.029 | 1.178 |
| MHC-II | 1.923 | 1.186 | 2.400 | 0.935 |
| MHC-I | 0.969 | 0.547 | 0.787 | 1.772 |
| COMPLEMENT | 0.115 | 0.183 | 0.356 | 0.313 |
| Prostaglandin | 0.290 | 0.427 | 0.346 | 0.845 |
| IL6 | 0.046 | 0.016 | 0.020 | 0.005 |
| IL1 | 0.045 | 0.032 | 0.013 | 0.016 |
| TNF | 0.031 | 0.127 | not detected | 0.070 |
| TGFb | not detected | not detected | 0.032 | not detected |

Summed across the four samples, the highest-probability pathways are MIF,
COLLAGEN, MHC-II, CXCL, LAMININ and APP
(`results/cellchat_pathway_communication_summary.csv`). MIF and COLLAGEN rank
first and second in every sample.

Descriptively: CCL and CXCL scores are higher in the omental samples than
in the SQ samples in both conditions; MHC-II scores are higher in the healthy
samples than in the unhealthy samples in both tissues; MHC-I and
Prostaglandin are highest in BRI-1459 (Unhealthy SQ). With n = 1 per
group these are observed patterns in these four samples, not group effects.

**Top inferred CXCL interactions** (`results/cellchat_CXCL_top10_sender_receiver.csv`):

| Sample | Top sender → receiver | Score |
| ------ | --------------------- | ----: |
| Healthy Omentum | Myeloid_macrophage → Endothelial | 0.578 |
| Unhealthy Omentum | Myeloid_macrophage → Endothelial | 0.630 |
| Unhealthy SQ | Myeloid_macrophage → Endothelial | 0.430 |
| Healthy SQ | Stromal_fibroblast → Dendritic_cells | 0.209 |

`Myeloid_macrophage → Endothelial` is the highest-scoring inferred CXCL
interaction in three of the four samples; in BRI-1458 the highest-scoring
interaction instead involves `Stromal_fibroblast` as sender.

**Top inferred MHC-II interactions**
(`results/cellchat_MHCII_top10_sender_receiver.csv`): highest-scoring pairs
involve `Dendritic_cells` and `Myeloid_macrophage` as senders and receivers
in all four samples (e.g. Dendritic_cells → Dendritic_cells, 0.527 in
BRI-1458; Myeloid_macrophage → Myeloid_macrophage, 0.512 in BRI-1456).

---

## 4. Repository layout

```
adipose-obesity-t2d-scrna/
├── README.md
├── .gitignore
├── adipose-obesity-t2d-scrna.Rproj
├── data/                  # raw GEO input — NOT tracked, download from GEO
│   ├── GSE278526_RAW.tar
│   ├── GSE278526_cell_barcodes_metadata.tsv.gz
│   └── GSM854823{5,6,7,8}_BRI-145{6,7,8,9}/   # barcodes / features / matrix
├── scripts/
│   ├── 01_primary_BRI1456.R … 01_primary_BRI1459.R
│   ├── 02_annotation_and_composition.R
│   ├── 03_inflammation_scoring.R
│   ├── 04_cellchat_analysis.R
│   ├── 05_cellchat_results.R
│   ├── 06_final_figures.R
│   └── legacy/            # superseded scripts — do not run
├── results/
│   ├── *.rds              # Seurat / CellChat objects — NOT tracked (large)
│   ├── *.csv              # result tables
│   └── figures/           # Figure1–Figure4 PNGs
└── output/
```

Tracked content is ~13 MB. The `data/` directory and `results/*.rds` objects
are excluded by `.gitignore` because they total ~2 GB and are either
re-downloadable from GEO or re-derivable from the scripts.

---

## 5. Reproducing the results

1. Download `GSE278526_RAW.tar` and `GSE278526_cell_barcodes_metadata.tsv.gz`
   from GEO into `data/`.
2. Extract the four sample folders into `data/` as named in §4.
3. Install R packages: `Seurat`, `SeuratObject`, `scDblFinder`, `SingleR`,
   `CellChat`, `here`, `dplyr`, `ggplot2`.
4. Run `scripts/01_*.R` → `06_final_figures.R` in order from the project root.
   All paths are resolved with `here::here()`, so scripts work from any
   working directory inside the project.

**Caveat:** `scripts/04_cellchat_analysis.R` documents the full four-sample
CellChat workflow, including the `saveRDS()` calls that write the four
`results/cellchat_BRI-*.rds` objects. **It was not executed while preparing
this repository.** Those four objects already exist from the original run and
were deliberately not regenerated, because doing so would mean re-running part
of the analysis. The script is therefore an accurate record of how they were
produced, not evidence that they were reproduced here.

---

## 6. Results interpretation limits

- **Descriptive only.** One sample per condition × tissue; no biological
  replication; no p-values, effect sizes or confidence intervals are reported.
- **`Unresolved` is a real, retained population** (8,115 cells, 24.7%),
  excluded from CellChat but kept in all composition tables.
- **CellChat output is inferred communication potential** from ligand–receptor
  co-expression. It does not establish pathway activity, direction of causal
  signalling, or protein-level interaction.
- **No claim of disease mechanism or treatment effect** is made anywhere in
  this repository. Observed patterns are reported as patterns in these four
  samples.
- Tissue labelling follows the official GEO records for GSE278526 — see §1.
  The deposited per-cell metadata file carries the opposite assignment and is
  not used by this repository.

---

## 7. License

This repository is released under the MIT License — see [`LICENSE`](LICENSE).

---

## 8. Data availability

All raw data are public at GEO under accession **GSE278526**. No participant
identifiers, sequencing barcodes beyond those published by the original
submitters, or credentials are included in this repository.
