# ------------------------------------------------------------------------------
# 04_integrate_annotate.R
# Integrates and annotates the merged object from script 03:
#   1. load the merged object
#   2. remove batch and capture effects from the PCA with Harmony
#   3. UMAP and clustering on the Harmony embedding
#   4. label cell types with Azimuth (PBMC reference, run locally)
#   5. diagnostic UMAPs, cell-count tables and plots
#
# Run:  Rscript new_scripts/04_integrate_annotate.R
#
# Input:  pipeline/merged/seurat_merged_preharmony.rds   (from script 03)
# Output: pipeline/merged/seurat_merged_harmony_azimuth.rds
#         pipeline/qc_plots/*.png, *.tsv
# ------------------------------------------------------------------------------

# Load a local GLPK build before Seurat (needed on Gadi). Edit or remove.
if (file.exists("/path/to/libglpk.so.40")) dyn.load("/path/to/libglpk.so.40")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(harmony)
  library(Azimuth)
  library(SeuratData)
  library(clustree)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(future)
  library(cowplot)
  library(scales)
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# Azimuth passes the whole object to future workers. Raise the size limit.
options(future.globals.maxSize = 32 * 1024^3)

# --- 1. Load merged object ----------------------------------------------------
in_path <- file.path(merged_dir, "seurat_merged_preharmony.rds")
if (!file.exists(in_path)) {
  stop("Missing merged object: ", in_path,
       " - did script 03 finish?")
}
log_msg("== Loading ", in_path)
seu <- readRDS(in_path)

# --- 2. Harmony on batch and capture ------------------------------------------
log_msg("== Running Harmony on PCA (group.by.vars = batch + capture)")
# Harmony has no seed argument, so set one here.
set.seed(42)
# Harmony adjusts the first 30 PCs so cells cluster by biology rather than by
# batch or capture. The counts themselves are not changed.
seu <- RunHarmony(
  seu,
  group.by.vars = c("batch", "capture"),
  reduction     = "pca",
  dims.use      = 1:30,
  reduction.save = "harmony",
  verbose       = FALSE
)

# --- 3. UMAP + neighbours + multi-resolution clustering -----------------------
log_msg("== UMAP + clustering on Harmony embedding")
# UMAP and nearest-neighbour graph from the Harmony embedding.
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
               reduction.name = "umap", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30,
                     verbose = FALSE)

# Clusters at several resolutions for clustree.
res_grid <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in res_grid) {
  seu <- FindClusters(seu, resolution = r, verbose = FALSE)
}
# Default identity for downstream use.
Idents(seu) <- "RNA_snn_res.0.6"
seu$seurat_clusters <- Idents(seu)

# Clustree diagnostic.
log_msg("  saving clustree")
ct_plot <- clustree(seu, prefix = "RNA_snn_res.")
ggsave(file.path(qc_dir, "clustree.png"),
       ct_plot, width = 8, height = 10, dpi = 150)

# --- 4. Local Azimuth annotation ----------------------------------------------
# Uses a local copy of the PBMC reference (no internet on compute nodes).
if (!exists("azimuth_pbmcref_dir") || is.null(azimuth_pbmcref_dir)) {
  stop("azimuth_pbmcref_dir is not set in config.R - required because Gadi ",
       "compute nodes cannot download the Azimuth reference.")
}
# Accept either the folder holding ref.Rds or the package folder above it.
azimuth_ref_path <- if (file.exists(file.path(azimuth_pbmcref_dir, "ref.Rds"))) {
  azimuth_pbmcref_dir
} else if (file.exists(file.path(azimuth_pbmcref_dir, "azimuth", "ref.Rds"))) {
  file.path(azimuth_pbmcref_dir, "azimuth")
} else {
  hit <- list.files(azimuth_pbmcref_dir, pattern = "^ref\\.Rds$",
                    recursive = TRUE, full.names = TRUE)
  if (length(hit) == 0L)
    stop("Could not find ref.Rds under ", azimuth_pbmcref_dir,
         " - check the local Azimuth reference install.")
  dirname(hit[1])
}
log_msg("== Running Azimuth (pbmcref) locally from ", azimuth_ref_path)
# Point Azimuth's homolog lookup at a local copy instead of downloading it.
if (!exists("azimuth_homologs_path") || is.null(azimuth_homologs_path)) {
  stop("azimuth_homologs_path is not set in config.R - required because Gadi ",
       "compute nodes cannot fetch homologs.rds from seurat.nygenome.org.")
}
# patch_azimuth_homologs() (utils.R) makes Azimuth read its gene-name table
# from the local file set in config.R.
patch_azimuth_homologs(azimuth_homologs_path)
# Map every cell onto the Azimuth PBMC reference to predict its cell type.
seu <- Azimuth::RunAzimuth(seu, reference = azimuth_ref_path, verbose = FALSE)
# RunAzimuth adds predicted.celltype.l1/l2, their scores and mapping.score.

# --- 5. Diagnostic plots ------------------------------------------------------
log_msg("== Saving post-Harmony UMAP diagnostics")
# Fixed colours per cell type and capture, so they match across plots
# (make_*_palette() in utils.R).
celltype_l1_pal <- make_celltype_palette(seu$predicted.celltype.l1)
celltype_l2_pal <- make_celltype_palette(seu$predicted.celltype.l2)
capture_pal     <- make_capture_palette(seu$capture)

# Post-Harmony UMAPs coloured by batch, capture, cohort, cluster and cell type.
p_batch <- DimPlot(seu, group.by = "batch", reduction = "umap",
                   pt.size = 1, label.size = 5) +
  ggtitle("Post-Harmony: batch")
p_capture <- DimPlot(seu, group.by = "capture", reduction = "umap",
                     pt.size = 1, cols = capture_pal) +
  ggtitle("Post-Harmony: capture") + theme(legend.position = "none")
p_cohort <- DimPlot(seu, group.by = "cohort", reduction = "umap",
                    pt.size = 1) +
  ggtitle("Post-Harmony: cohort")
p_cluster <- DimPlot(seu, group.by = "seurat_clusters", reduction = "umap",
                     label = TRUE, repel = TRUE,
                     pt.size = 1, label.size = 5) +
  ggtitle("Post-Harmony: clusters (res=0.6)")
p_l1 <- DimPlot(seu, group.by = "predicted.celltype.l1", reduction = "umap",
                label = TRUE, repel = TRUE,
                pt.size = 1, label.size = 5,
                cols = celltype_l1_pal) +
  ggtitle("Azimuth L1")
p_l2 <- DimPlot(seu, group.by = "predicted.celltype.l2", reduction = "umap",
                label = TRUE, repel = TRUE,
                pt.size = 1, label.size = 3,
                cols = celltype_l2_pal) +
  ggtitle("Azimuth L2") + theme(legend.position = "none")

ggsave(file.path(qc_dir, "postharmony_umap_batch.png"),
       p_batch, width = 6, height = 5, dpi = 150)
ggsave(file.path(qc_dir, "postharmony_umap_capture.png"),
       p_capture, width = 6, height = 5, dpi = 150)
ggsave(file.path(qc_dir, "postharmony_umap_cohort.png"),
       p_cohort, width = 6, height = 5, dpi = 150)
ggsave(file.path(qc_dir, "postharmony_umap_cluster.png"),
       p_cluster, width = 7, height = 6, dpi = 150)
ggsave(file.path(qc_dir, "azimuth_umap_l1.png"),
       p_l1, width = 7, height = 6, dpi = 150)
ggsave(file.path(qc_dir, "azimuth_umap_l2.png"),
       p_l2, width = 8, height = 7, dpi = 150)

# --- 5b. Per-cohort UMAPs -----------------------------------------------------
# GEM108 is split into GEM108_pre and GEM108_post. Everyone else keeps their
# cohort label.
log_msg("== Saving per-cohort UMAPs")

# cohort_or_patient_tx: cohort label, except GEM108 cells become GEM108_pre or
# GEM108_post. GEM108 cells with no treatment label keep their cohort label.
# Later scripts use this column.
tx_short <- sub("_treatment$", "", seu$treatment)
seu$cohort_or_patient_tx <- ifelse(
  seu$patient == "GEM108" & !is.na(tx_short) & nzchar(tx_short),
  paste0("GEM108_", tx_short),
  as.character(seu$cohort)
)

cohort_dir <- file.path(qc_dir, "umap_per_cohort")
dir.create(cohort_dir, showWarnings = FALSE, recursive = TRUE)

# One UMAP per group, with that group's cells highlighted.
cohort_groups <- sort(unique(stats::na.omit(seu$cohort_or_patient_tx)))
log_msg("   groups: ", paste(cohort_groups, collapse = ", "))

for (g in cohort_groups) {
  cells_in <- colnames(seu)[which(seu$cohort_or_patient_tx == g)]
  n_in     <- length(cells_in)
  p_g <- DimPlot(
    seu,
    reduction       = "umap",
    cells.highlight = list(group = cells_in),
    cols.highlight  = "firebrick",
    cols            = "grey85",
    sizes.highlight = 0.6,
    pt.size         = 0.4,
    order           = TRUE
  ) +
    ggtitle(paste0(g, "  (n = ", n_in, " cells)")) +
    theme(legend.position = "none",
          plot.title      = element_text(hjust = 0.5))
  out_file <- file.path(
    cohort_dir,
    paste0("umap_cohort_", gsub("[^A-Za-z0-9]+", "_", g), ".png")
  )
  ggsave(out_file, p_g, width = 6, height = 5, dpi = 150)
}

# Faceted UMAPs, one panel per cohort.
n_panels <- length(cohort_groups)
facet_w  <- max(8, 4 * n_panels)
p_split_l2 <- DimPlot(
  seu,
  reduction = "umap",
  group.by  = "predicted.celltype.l2",
  split.by  = "cohort_or_patient_tx",
  label     = TRUE, repel = TRUE, label.size = 3,
  pt.size   = 0.4,
  cols      = celltype_l2_pal
) + ggtitle("Azimuth L2, split by cohort/treatment") +
  theme(legend.position = "none")
ggsave(file.path(qc_dir, "umap_split_cohort_azimuth_l2.png"),
       p_split_l2, width = facet_w, height = 5, dpi = 150)

# Same, coloured by cluster.
p_split_cluster <- DimPlot(
  seu,
  reduction = "umap",
  group.by  = "seurat_clusters",
  split.by  = "cohort_or_patient_tx",
  label     = TRUE, repel = TRUE, label.size = 3,
  pt.size   = 0.4
) + ggtitle("Clusters (res=0.6), split by cohort/treatment")
ggsave(file.path(qc_dir, "umap_split_cohort_clusters.png"),
       p_split_cluster, width = facet_w, height = 5, dpi = 150)

# Cell-count tables for the record.
ct_l1 <- seu@meta.data %>%
  dplyr::count(cohort, predicted.celltype.l1) %>%
  tidyr::pivot_wider(names_from = predicted.celltype.l1,
                     values_from = n, values_fill = 0)
readr::write_tsv(ct_l1, file.path(qc_dir, "celltype_counts_l1.tsv"))

ct_l2 <- seu@meta.data %>%
  dplyr::count(cohort, predicted.celltype.l2) %>%
  tidyr::pivot_wider(names_from = predicted.celltype.l2,
                     values_from = n, values_fill = 0)
readr::write_tsv(ct_l2, file.path(qc_dir, "celltype_counts_l2.tsv"))

# Same tables split by cohort_or_patient_tx.
ct_l1_tx <- seu@meta.data %>%
  dplyr::count(cohort_or_patient_tx, predicted.celltype.l1) %>%
  tidyr::pivot_wider(names_from = predicted.celltype.l1,
                     values_from = n, values_fill = 0)
readr::write_tsv(ct_l1_tx,
                 file.path(qc_dir, "celltype_counts_l1_by_tx.tsv"))

ct_l2_tx <- seu@meta.data %>%
  dplyr::count(cohort_or_patient_tx, predicted.celltype.l2) %>%
  tidyr::pivot_wider(names_from = predicted.celltype.l2,
                     values_from = n, values_fill = 0)
readr::write_tsv(ct_l2_tx,
                 file.path(qc_dir, "celltype_counts_l2_by_tx.tsv"))

# --- 5c. Azimuth QC violins ---------------------------------------------------
log_msg("== Azimuth QC violins")
# Azimuth prediction score, UMIs and genes per cell type.
md <- seu@meta.data
md$predicted.celltype.l2 <- factor(
  md$predicted.celltype.l2,
  levels = sort(unique(stats::na.omit(as.character(md$predicted.celltype.l2))))
)
v_theme <- theme(text = element_text(size = 8),
                 axis.text.x = element_text(angle = 90,
                                            vjust = 0.5, hjust = 1))
v_score <- ggplot(md,
    aes(y = predicted.celltype.l2.score,
        x = predicted.celltype.l2,
        fill = predicted.celltype.l2)) +
  geom_violin(scale = "width") +
  scale_fill_manual(values = celltype_l2_pal) +
  v_theme + NoLegend()
v_count <- ggplot(md,
    aes(y = nCount_RNA,
        x = predicted.celltype.l2,
        fill = predicted.celltype.l2)) +
  geom_violin(scale = "area") +
  scale_y_log10() +
  scale_fill_manual(values = celltype_l2_pal) +
  v_theme + NoLegend()
v_feat  <- ggplot(md,
    aes(y = nFeature_RNA,
        x = predicted.celltype.l2,
        fill = predicted.celltype.l2)) +
  geom_violin(scale = "area") +
  scale_y_log10() +
  scale_fill_manual(values = celltype_l2_pal) +
  v_theme + NoLegend()
v_stack <- v_score / v_count / v_feat
ggsave(file.path(qc_dir, "azimuth_qc_violins.png"),
       v_stack, width = 9, height = 9, dpi = 150)

# --- 5d. Cell-count bar plots -------------------------------------------------
log_msg("== Cell-count stacked bars")
# Cells per cohort/treatment group and per patient, coloured by cell type.
p_cohort_bar <- ggplot(md,
    aes(x = cohort_or_patient_tx,
        fill = predicted.celltype.l2)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            hjust = 1.5, size = 2) +
  scale_fill_manual(values = celltype_l2_pal) +
  coord_flip() +
  labs(title = "Cells per cohort/treatment, fill = Azimuth L2",
       x = "cohort_or_patient_tx", y = "cells") +
  cowplot::theme_cowplot(font_size = 8) +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"))
ggsave(file.path(qc_dir, "cell_count_by_cohort.png"),
       p_cohort_bar, width = 8, height = 5, dpi = 150)

p_patient_bar <- ggplot(md,
    aes(x = patient,
        fill = predicted.celltype.l2)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            hjust = 1.5, size = 2) +
  scale_fill_manual(values = celltype_l2_pal) +
  coord_flip() +
  labs(title = "Cells per patient, fill = Azimuth L2",
       x = "patient", y = "cells") +
  cowplot::theme_cowplot(font_size = 8) +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"))
ggsave(file.path(qc_dir, "cell_count_by_patient.png"),
       p_patient_bar, width = 8, height = 6, dpi = 150)

# --- 6. Save ------------------------------------------------------------------
# Save for scripts 05 to 11.
out_path <- file.path(merged_dir, "seurat_merged_harmony_azimuth.rds")
saveRDS(seu, out_path)
log_msg("== saved ", out_path)
log_msg("== 04_integrate_annotate.R done.")
