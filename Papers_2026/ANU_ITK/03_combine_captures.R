# ------------------------------------------------------------------------------
# 03_combine_captures.R
# Combines the 20 QC'd captures into one Seurat object:
#   1. load the output of script 02 for every capture
#   2. add patient, cohort and treatment from the sample sheet
#   3. keep clean singlets only
#   4. draw per-capture QC and hashtag plots, convert each capture to Seurat
#   5. merge the captures
#   6. normalise, find variable genes and scale
#   7. PCA and a pre-Harmony UMAP, to show the batch effect before script 04
#
# Run:  Rscript new_scripts/03_combine_captures.R
#
# Input:  pipeline/SCEs/dataset_{1..20}_qcd.rds   (from script 02)
#         pipeline/sample_sheet.csv
# Output: pipeline/merged/seurat_merged_preharmony.rds
#         pipeline/qc_plots/*.png
# ------------------------------------------------------------------------------

# Load a local GLPK build before Seurat (needed on Gadi). Edit or remove.
if (file.exists("/path/to/libglpk.so.40")) dyn.load("/path/to/libglpk.so.40")

suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(readr)
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# --- 1. Load per-capture QCd SCEs ---------------------------------------------
log_msg("== Loading per-capture SCEs")
sce_list <- list()
for (i in capture_table$i) {
  # Dataset ID and capture name for capture i (get_capture_info() in utils.R).
  info <- get_capture_info(i)
  p <- file.path(sce_dir, paste0(info$dataset, "_qcd.rds"))
  if (!file.exists(p)) {
    stop("Missing QCd SCE: ", p, " - did script 02 finish for capture ", i, "?")
  }
  sce_list[[info$dataset]] <- readRDS(p)
  log_msg("  ", info$dataset, ": ", ncol(sce_list[[info$dataset]]), " cells")
}

# --- 2. Sample-sheet attach ---------------------------------------------------
log_msg("== Attaching sample-sheet metadata")
# Read pipeline/sample_sheet.csv. load_sample_sheet() (utils.R) checks the
# required columns are present and no (capture, hto) pair appears twice.
ss <- load_sample_sheet()

# Add patient, cohort and treatment to every cell, matched on capture and the
# cell's best hashtag (HTO_best from script 02).
attach_sample_sheet <- function(sce, info) {
  # Sample sheet rows for this capture.
  ss_cap <- ss[ss$capture == info$capture, , drop = FALSE]
  if (nrow(ss_cap) == 0) {
    stop("Sample sheet has no rows for capture ", info$capture)
  }
  # Match on the hashtag name, dropping any suffix Cell Ranger adds.
  hto_key <- sub("[-_ ].*$", "", sce$HTO_best)
  # Row of the sample sheet for each cell.
  idx <- match(paste0(sce$capture, "|", hto_key),
               paste0(ss_cap$capture, "|", ss_cap$hto))
  sce$patient   <- ss_cap$patient[idx]
  sce$cohort    <- ss_cap$cohort[idx]
  sce$treatment <- ss_cap$treatment[idx]
  # Cells whose hashtag is not in the sheet get NA and are removed in step 3.

  # --- Single-HTO captures: rebuild is_singlet --------------------------------
  # With only one HTO loaded, hashedDrops cannot call confidence, so is_singlet
  # from script 02 would drop every cell. For these captures a singlet is a cell
  # whose best HTO is the expected one and that scDblFinder calls a singlet.
  if (nrow(ss_cap) == 1) {
    expected_hto <- ss_cap$hto[1]
    sce$is_singlet <- !is.na(sce$HTO_best) &
      sce$HTO_best == expected_hto &
      !is.na(sce$scDblFinder.class) &
      sce$scDblFinder.class == "singlet"
    log_msg("    [", info$capture, "] single-HTO capture: is_singlet ",
            "rebuilt against expected '", expected_hto, "' (",
            sum(sce$is_singlet, na.rm = TRUE), " singlets)")
  }
  sce
}

for (nm in names(sce_list)) {
  # Attach sample-sheet information to every capture.
  i <- capture_table$i[capture_table$dataset == nm]
  info_obj <- get_capture_info(i)
  sce_list[[nm]] <- attach_sample_sheet(sce_list[[nm]], info_obj)
}

# --- 3. Filter to clean singlets ----------------------------------------------
# Keep singlets that are not mito outliers, not DropletQC empty droplets and
# have a patient and cohort in the sample sheet. NA counts as fail.
filter_clean <- function(sce) {
  # Singlet and mito flags come from script 02.
  is_singlet_safe   <- !is.na(sce$is_singlet)   & sce$is_singlet
  mito_outlier_safe <-  is.na(sce$mito_outlier) | sce$mito_outlier   # NA = exclude
  keep <- is_singlet_safe &
    !mito_outlier_safe &
    !is.na(sce$patient) & !is.na(sce$cohort)
  # Drop DropletQC empty droplets. Cells without a DropletQC call are kept.
  if ("dropletqc_cell_status" %in% colnames(colData(sce))) {
    dqc_ok <- is.na(sce$dropletqc_cell_status) |
      sce$dropletqc_cell_status != "empty_droplet"
    keep <- keep & dqc_ok
  }
  # Drop samples listed in config.R samples_to_exclude.
  if (exists("samples_to_exclude") && !is.null(samples_to_exclude) &&
      nrow(samples_to_exclude) > 0L) {
    hto_key_excl <- sub("[-_ ].*$", "", sce$HTO_best)
    excl_key <- paste0(samples_to_exclude$capture, "|", samples_to_exclude$hto)
    is_excluded <- paste0(sce$capture, "|", hto_key_excl) %in% excl_key
    if (any(is_excluded)) {
      log_msg("    [", unique(sce$capture), "] dropping ",
              sum(is_excluded), " cell(s) via samples_to_exclude")
    }
    keep <- keep & !is_excluded
  }
  sce[, which(keep)]
}

log_msg("== Filtering to clean singlets")
# Apply the filter to every capture and log how many cells remain.
clean_counts <- integer(length(sce_list))
names(clean_counts) <- names(sce_list)
for (nm in names(sce_list)) {
  before <- ncol(sce_list[[nm]])
  sce_list[[nm]] <- filter_clean(sce_list[[nm]])
  after <- ncol(sce_list[[nm]])
  clean_counts[nm] <- after
  log_msg("  ", nm, ": ", before, " -> ", after, " singlets")
}

# Drop captures that ended up empty (shouldn't happen, but guard).
sce_list <- sce_list[clean_counts > 0]
stopifnot(length(sce_list) > 0)

# --- 4a. Per-capture QC violins -----------------------------------------------
log_msg("== Plotting per-capture QC violins")
for (nm in names(sce_list)) {
  sce <- sce_list[[nm]]
  # UMIs, genes detected and mitochondrial % per cell (from script 02).
  qc_df <- data.frame(
    nCount_RNA           = sce$sum,
    nFeature_RNA         = sce$detected,
    subsets_Mito_percent = sce$subsets_Mito_percent
  )
  qc_long <- tidyr::pivot_longer(qc_df, dplyr::everything(),
                                 names_to = "metric", values_to = "value")
  p <- ggplot(qc_long, aes(x = metric, y = value)) +
    geom_violin(fill = "grey80", scale = "width", trim = FALSE) +
    facet_wrap(~ metric, scales = "free", ncol = 3) +
    theme_bw() +
    labs(title = nm, x = NULL, y = NULL) +
    theme(strip.text = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  ggsave(file.path(qc_dir, paste0(nm, "__qc_violin.png")),
         p, width = 9, height = 3.5, dpi = 150)
}

# --- 4a2. HTO and capture demultiplexing bar plots ----------------------------
log_msg("== HTO/capture demultiplexing stacked bars")
# One row per cell across all captures: capture, best hashtag, confident call.
demux_df <- do.call(rbind, lapply(names(sce_list), function(nm) {
  cd <- as.data.frame(colData(sce_list[[nm]]))
  data.frame(
    capture       = as.character(cd$capture),
    HTO_best      = as.character(cd$HTO_best),
    HTO_confident = !is.na(cd$HTO_confident) & cd$HTO_confident,
    stringsAsFactors = FALSE
  )
}))
# Fixed colours for hashtags and captures (utils.R).
hto_pal     <- make_hto_palette(demux_df$HTO_best)
capture_pal <- make_capture_palette(demux_df$capture)

# Cells per hashtag (counts and proportions), coloured by confident call.
p_hto_count <- ggplot(demux_df,
    aes(x = HTO_best, fill = HTO_confident)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  ylab("Number of droplets") +
  cowplot::theme_cowplot(font_size = 7)
p_hto_prop <- ggplot(demux_df,
    aes(x = HTO_best, fill = HTO_confident)) +
  geom_bar(position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Proportion of droplets") +
  cowplot::theme_cowplot(font_size = 7)
ggsave(file.path(qc_dir, "demux_hto_confidence.png"),
       p_hto_count / p_hto_prop, width = 6, height = 6, dpi = 150)

# Hashtag mix in each capture.
p_cap_hto <- ggplot(demux_df,
    aes(x = capture, fill = HTO_best)) +
  geom_bar(position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  scale_fill_manual(values = hto_pal) +
  cowplot::theme_cowplot(font_size = 8)

# Cells per capture.
p_cap_count <- ggplot(demux_df,
    aes(x = capture, fill = capture)) +
  geom_bar() +
  coord_flip() +
  ylab("Number of droplets") +
  scale_fill_manual(values = capture_pal) +
  cowplot::theme_cowplot(font_size = 8) +
  guides(fill = "none")

ggsave(file.path(qc_dir, "demux_capture_hto.png"),
       p_cap_hto, width = 6, height = 5, dpi = 150)
ggsave(file.path(qc_dir, "demux_capture_count.png"),
       p_cap_count, width = 6, height = 5, dpi = 150)

# --- 4b. Convert to Seurat ----------------------------------------------------
# Replace '_' with '-' in gene names, as CreateSeuratObject would do anyway,
# and log how many genes are affected.
log_msg("== Converting each capture to Seurat")
seu_list <- lapply(names(sce_list), function(nm) {
  sce <- sce_list[[nm]]
  cnt <- counts(sce)
  n_underscore <- sum(grepl("_", rownames(cnt)))
  if (n_underscore > 0) {
    rownames(cnt) <- gsub("_", "-", rownames(cnt), fixed = TRUE)
  }
  # Seurat object with all genes and cells, keeping the cell metadata.
  s <- CreateSeuratObject(
    counts    = cnt,
    project   = nm,
    min.cells = 0, min.features = 0,
    meta.data = as.data.frame(colData(sce))
  )
  s$dataset <- nm
  attr(s, "n_renamed_underscore") <- n_underscore
  s
})
names(seu_list) <- names(sce_list)
# Log the rename count.
total_renamed <- vapply(seu_list, function(s) {
  v <- attr(s, "n_renamed_underscore")
  if (is.null(v)) 0L else as.integer(v)
}, integer(1))
log_msg("  '_' → '-' rename per capture: ",
        paste(names(total_renamed), total_renamed, sep = "=",
              collapse = ", "))
rm(sce_list); gc(verbose = FALSE)

# --- 5. Merge all captures ----------------------------------------------------
log_msg("== Merging all captures into a single Seurat")
first <- seu_list[[1]]
rest  <- seu_list[-1]
# Merge all captures. Cell names get the dataset ID as a prefix so barcodes
# stay unique.
seu <- merge(
  first,
  y = rest,
  add.cell.ids = names(seu_list),
  project = "scRNA_ITK_merged",
  merge.data = FALSE
)
# Join the per-capture layers into single layers.
if (inherits(seu[["RNA"]], "Assay5")) {
  seu <- JoinLayers(seu)
}
rm(seu_list, first, rest); gc(verbose = FALSE)
log_msg("  merged object: ", ncol(seu), " cells, ",
        nrow(seu), " genes")

# --- 6. Normalise + variable features + scale ---------------------------------
log_msg("== NormalizeData / FindVariableFeatures / ScaleData")
# Log-normalise counts and pick the 3,000 most variable genes.
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
# All genes stay in the assay. Only the variable genes are scaled.
seu <- ScaleData(seu, features = VariableFeatures(seu),
                 vars.to.regress = NULL, verbose = FALSE)

# --- 7. Pre-integration PCA / UMAP / clustering -------------------------------
log_msg("== Pre-integration PCA + UMAP (Harmony will re-do this in script 04)")
# 50 PCs on the variable genes. The first 30 are used for the UMAP.
seu <- RunPCA(seu, features = VariableFeatures(seu),
              npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
               reduction.name = "umap_preharmony", verbose = FALSE)

# Pre-Harmony UMAPs by batch and capture.
p_batch <- DimPlot(seu, group.by = "batch",
                   reduction = "umap_preharmony") +
  ggtitle("Pre-Harmony: coloured by batch")
p_capture <- DimPlot(seu, group.by = "capture",
                     reduction = "umap_preharmony") +
  ggtitle("Pre-Harmony: coloured by capture") +
  theme(legend.position = "none")
ggsave(file.path(qc_dir, "preharmony_umap_batch.png"),
       p_batch, width = 6, height = 5, dpi = 150)
ggsave(file.path(qc_dir, "preharmony_umap_capture.png"),
       p_capture, width = 6, height = 5, dpi = 150)

# --- 7b. ElbowPlot and DimHeatmap ---------------------------------------------
log_msg("== ElbowPlot + DimHeatmap")
# Plots used to choose the number of PCs.
p_elbow <- ElbowPlot(seu, ndims = 30, reduction = "pca") +
  ggtitle("Pre-Harmony PCA elbow")
ggsave(file.path(qc_dir, "preharmony_elbow.png"),
       p_elbow, width = 6, height = 4, dpi = 150)

# DimHeatmap draws to the open device.
png(file.path(qc_dir, "preharmony_dimheatmap.png"),
    width = 1500, height = 2000, res = 150)
DimHeatmap(seu, dims = 1:30, cells = 500, balanced = TRUE,
           reduction = "pca", fast = TRUE)
dev.off()

# --- 8. Save ------------------------------------------------------------------
# Save for script 04.
out_path <- file.path(merged_dir, "seurat_merged_preharmony.rds")
saveRDS(seu, out_path)
log_msg("== saved ", out_path)
log_msg("== 03_combine_captures.R done.")
