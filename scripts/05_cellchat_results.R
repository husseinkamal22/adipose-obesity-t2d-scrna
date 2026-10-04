# ============================================================
# 05_cellchat_results.R
# CellChat results extraction
# Project: Adipose obesity / metabolic disease scRNA-seq
# ============================================================

library(CellChat)
library(dplyr)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

results_dir <- here::here("results")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load saved CellChat objects
# ------------------------------------------------------------

cellchat1 <- readRDS(
  file.path(results_dir, "cellchat_BRI-1456.rds")
)

cellchat2 <- readRDS(
  file.path(results_dir, "cellchat_BRI-1457.rds")
)

cellchat3 <- readRDS(
  file.path(results_dir, "cellchat_BRI-1458.rds")
)

cellchat4 <- readRDS(
  file.path(results_dir, "cellchat_BRI-1459.rds")
)

cat("\nAll four CellChat objects loaded successfully.\n")

# ------------------------------------------------------------
# 3. Define sample labels
# ------------------------------------------------------------

sample_labels <- c(
  "Healthy SQ",
  "Unhealthy SQ",
  "Healthy Omentum",
  "Unhealthy Omentum"
)

cellchat_list <- list(
  cellchat1,
  cellchat2,
  cellchat3,
  cellchat4
)

# ------------------------------------------------------------
# 4. Extract pathway communication strength
# ------------------------------------------------------------

extract_pathway_strength <- function(cellchat, sample_label) {
  
  prob <- cellchat@netP$prob
  
  pathway_names <- dimnames(prob)[[3]]
  
  pathway_strength <- apply(
    prob,
    3,
    sum,
    na.rm = TRUE
  )
  
  data.frame(
    sample = sample_label,
    pathway = pathway_names,
    communication_probability = pathway_strength,
    stringsAsFactors = FALSE
  )
}

cellchat_summary <- bind_rows(
  extract_pathway_strength(cellchat1, sample_labels[1]),
  extract_pathway_strength(cellchat2, sample_labels[2]),
  extract_pathway_strength(cellchat3, sample_labels[3]),
  extract_pathway_strength(cellchat4, sample_labels[4])
)

# ------------------------------------------------------------
# 5. Add disease and tissue information
# ------------------------------------------------------------

cellchat_summary <- cellchat_summary %>%
  mutate(
    disease = ifelse(
      grepl("Unhealthy", sample),
      "Unhealthy",
      "Healthy"
    ),
    tissue = ifelse(
      grepl("SQ", sample),
      "SQ",
      "Omentum"
    )
  ) %>%
  select(
    sample,
    disease,
    tissue,
    pathway,
    communication_probability
  )

write.csv(
  cellchat_summary,
  file.path(
    results_dir,
    "cellchat_pathway_communication_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Focused pathways
# ------------------------------------------------------------

focused_pathways <- c(
  "MHC-II",
  "CXCL",
  "CCL",
  "IL6",
  "IL1",
  "TNF",
  "TGFb",
  "MHC-I",
  "COMPLEMENT",
  "Prostaglandin"
)

cellchat_focus <- cellchat_summary %>%
  filter(pathway %in% focused_pathways) %>%
  arrange(pathway, sample)

write.csv(
  cellchat_focus,
  file.path(
    results_dir,
    "cellchat_focused_pathways.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Extract detailed pathway communication
# ------------------------------------------------------------

extract_pathway_table <- function(
    cellchat,
    signaling_name,
    sample_label
) {
  
  result <- subsetCommunication(
    cellchat,
    signaling = signaling_name
  )
  
  if (is.null(result) || nrow(result) == 0) {
    return(data.frame())
  }
  
  result %>%
    mutate(sample = sample_label) %>%
    select(
      sample,
      source,
      target,
      ligand,
      receptor,
      prob
    )
}

# ------------------------------------------------------------
# 8. MHC-II detailed table
# ------------------------------------------------------------

mhcii_table <- bind_rows(
  extract_pathway_table(
    cellchat1,
    "MHC-II",
    sample_labels[1]
  ),
  extract_pathway_table(
    cellchat2,
    "MHC-II",
    sample_labels[2]
  ),
  extract_pathway_table(
    cellchat3,
    "MHC-II",
    sample_labels[3]
  ),
  extract_pathway_table(
    cellchat4,
    "MHC-II",
    sample_labels[4]
  )
)

write.csv(
  mhcii_table,
  file.path(
    results_dir,
    "cellchat_MHCII_detailed.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. MHC-II top sender-receiver pairs
# ------------------------------------------------------------

mhcii_top10 <- mhcii_table %>%
  group_by(
    sample,
    source,
    target
  ) %>%
  summarise(
    total_prob = sum(
      prob,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  slice_max(
    total_prob,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    sample,
    desc(total_prob)
  )

write.csv(
  mhcii_top10,
  file.path(
    results_dir,
    "cellchat_MHCII_top10_sender_receiver.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. CXCL detailed table
# ------------------------------------------------------------

cxcl_table <- bind_rows(
  extract_pathway_table(
    cellchat1,
    "CXCL",
    sample_labels[1]
  ),
  extract_pathway_table(
    cellchat2,
    "CXCL",
    sample_labels[2]
  ),
  extract_pathway_table(
    cellchat3,
    "CXCL",
    sample_labels[3]
  ),
  extract_pathway_table(
    cellchat4,
    "CXCL",
    sample_labels[4]
  )
)

write.csv(
  cxcl_table,
  file.path(
    results_dir,
    "cellchat_CXCL_detailed.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. CXCL top sender-receiver pairs
# ------------------------------------------------------------

cxcl_top10 <- cxcl_table %>%
  group_by(
    sample,
    source,
    target
  ) %>%
  summarise(
    total_prob = sum(
      prob,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  slice_max(
    total_prob,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    sample,
    desc(total_prob)
  )

write.csv(
  cxcl_top10,
  file.path(
    results_dir,
    "cellchat_CXCL_top10_sender_receiver.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("CellChat results extraction completed.\n")
cat("============================================\n\n")

cat("Files created:\n\n")

cat("1. cellchat_pathway_communication_summary.csv\n")
cat("2. cellchat_focused_pathways.csv\n")
cat("3. cellchat_MHCII_detailed.csv\n")
cat("4. cellchat_MHCII_top10_sender_receiver.csv\n")
cat("5. cellchat_CXCL_detailed.csv\n")
cat("6. cellchat_CXCL_top10_sender_receiver.csv\n\n")

cat("No CellChat analysis was recomputed.\n")