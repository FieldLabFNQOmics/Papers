# ------------------------------------------------------------------------------
# 02_qc_and_demux.R
# Per-capture QC. For each capture:
#   1. keep droplets that emptyDrops called as cells
#   2. score doublets with scDblFinder
#   3. assign each cell to a donor from its hashtag (HTO) with hashedDrops
#   4. drop genes seen in very few cells
#   5. flag cells with high mitochondrial %
#   6. score nuclear fraction with DropletQC (needs the Cell Ranger BAM)
#   7. save QC plots
# Cells are flagged here, not removed. Script 03 does the filtering.
#
# Run:  Rscript new_scripts/02_qc_and_demux.R [i1 i2 ...] [--bam-template TMPL]
#   i1 i2 ...        capture numbers 1-20 (default: all 20)
#   --bam-template   path pattern for each capture's BAM file. May contain
#                    {batch_root}, {capture} and {sample}.
#
# Input:  pipeline/SCEs/dataset_{i}.rds          (from script 01)
# Output: pipeline/SCEs/dataset_{i}_qcd.rds
#         pipeline/qc_plots/per_capture/dataset_{i}/*.png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(scran)
  library(scuttle)
  library(scDblFinder)
  library(DropletQC)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# --- Argument parsing ---------------------------------------------------------
# Command-line arguments after the script name, e.g. c("3", "7").
raw_args <- commandArgs(trailingOnly = TRUE)
# Default BAM path. Cell Ranger was run from a folder named after the
# capture, so the capture name appears twice in the path.
bam_template_default <- file.path(
  "{batch_root}", "{capture}", "{capture}",
  "outs", "per_sample_outs",
  "{sample}", "count", "sample_alignments.bam"
)
# Read capture numbers and an optional --bam-template from the arguments.
bam_template <- bam_template_default
indices <- integer(0)
i <- 1
while (i <= length(raw_args)) {
  a <- raw_args[i]
  if (a == "--bam-template") {
    bam_template <- raw_args[i + 1]
    i <- i + 2
  } else {
    indices <- c(indices, as.integer(a))
    i <- i + 1
  }
}
# No capture numbers given: process all 20.
if (length(indices) == 0) indices <- capture_table$i
# Stop before doing any work if a number is not a valid capture (1-20).
bad <- indices[!(indices %in% capture_table$i) | is.na(indices)]
if (length(bad) > 0L) {
  stop("Bad capture index/indices: ",
       paste(bad, collapse = ", "),
       ". Valid range is 1..", max(capture_table$i),
       " (one per 10x capture). Submit script 02 with `-J 1-10` ",
       "then `-J 11-20`; the 440-element array is for script 05 only.",
       call. = FALSE)
}

# Build the BAM path for one capture by filling the {batch_root}, {capture}
# and {sample} placeholders in bam_template.
resolve_bam_path <- function(info, sample = info$capture) {
  p <- bam_template
  p <- gsub("\\{batch_root\\}", info$batch_root, p, fixed = FALSE)
  p <- gsub("\\{capture\\}",    info$capture,    p, fixed = FALSE)
  p <- gsub("\\{sample\\}",     sample,          p, fixed = FALSE)
  p
}

# --- Process one capture ------------------------------------------------------
# Runs steps 1-7 for capture number i and saves the QC'd object.
process_capture <- function(i) {
  # Dataset ID, capture name, batch folder and file paths for this capture
  # (get_capture_info() in utils.R reads capture_table).
  info <- get_capture_info(i)
  log_msg("== Capture ", info$dataset, " (", info$capture, ")")

  # Load the SingleCellExperiment saved by script 01.
  in_path <- file.path(sce_dir, paste0(info$dataset, ".rds"))
  sce <- readRDS(in_path)

  # --- 1. Filter to called cells ----------------------------------------------
  # Keep droplets whose emptyDrops FDR is below the cutoff in config.R.
  is_cell <- !is.na(sce$FDR) & sce$FDR < emptydrops_fdr
  log_msg("  emptyDrops cells: ", sum(is_cell), " / ", ncol(sce))
  sce <- sce[, is_cell]

  # --- 2. Doublet scoring via scDblFinder -------------------------------------
  log_msg("  running scDblFinder")
  # scDblFinder scores each cell against simulated doublets. Seeded so the
  # result is reproducible.
  set.seed(42)
  sce <- scDblFinder::scDblFinder(sce, clusters = FALSE,
                                  BPPARAM = BiocParallel::SerialParam())
  # scDblFinder adds colData columns: scDblFinder.score, scDblFinder.class
  log_msg("  scDblFinder doublets: ",
          sum(sce$scDblFinder.class == "doublet"), " / ", ncol(sce))

  # --- 3. HTO demux via hashedDrops -------------------------------------------
  # HTO counts were stored as an alternative experiment in script 01.
  hto_altexp <- altExp(sce, "HTO")
  if (is.null(hto_altexp)) {
    stop("No HTO altExp on ", info$dataset,
         " - check script 01 feature-type split.")
  }
  log_msg("  running hashedDrops")
  # hashedDrops picks the most abundant HTO in each cell and says whether the
  # call is confident and whether the cell looks like a doublet.
  hashed <- DropletUtils::hashedDrops(
    counts(hto_altexp),
    confident.min = hashed_confident_min
  )
  # hashed rows align with sce columns (same barcodes). Sanity-check.
  stopifnot(identical(rownames(hashed), colnames(sce)))

  # Store the hashtag calls on each cell.
  sce$HTO_best       <- rownames(hto_altexp)[hashed$Best]
  sce$HTO_confident  <- hashed$Confident
  sce$HTO_doublet    <- hashed$Doublet
  sce$HTO_second     <- rownames(hto_altexp)[hashed$Second]
  sce$HTO_logFC      <- hashed$LogFC

  # --- 3b. Final singlet call -------------------------------------------------
  # Singlet = confident HTO, not an HTO doublet and scDblFinder singlet.
  # hashedDrops returns NA for Confident when only one HTO was loaded in a
  # capture, so NA is treated as FALSE here. Script 03 re-derives singlets for
  # single-HTO captures.
  hto_conf  <- !is.na(sce$HTO_confident) & sce$HTO_confident
  hto_doub  <- !is.na(sce$HTO_doublet)   & sce$HTO_doublet
  scdbl_sing <- !is.na(sce$scDblFinder.class) &
                sce$scDblFinder.class == "singlet"
  sce$is_singlet <- hto_conf & !hto_doub & scdbl_sing
  log_msg("  singlets after consensus demux: ",
          sum(sce$is_singlet), " / ", ncol(sce))
  # Non-singlets are kept here and removed in script 03.

  # --- 4. Gene annotation + min-cells-per-gene filter -------------------------
  # Use gene symbols as row names. Duplicated symbols get the Ensembl ID added.
  rownames(sce) <- scuttle::uniquifyFeatureNames(
    rowData(sce)$ID, rowData(sce)$Symbol
  )
  # Keep genes detected in at least min_cells_per_gene cells (config.R).
  gene_kept <- Matrix::rowSums(counts(sce) > 0) >= min_cells_per_gene
  log_msg("  genes kept (>=", min_cells_per_gene, " cells): ",
          sum(gene_kept), " / ", length(gene_kept))
  sce <- sce[gene_kept, ]

  # --- 5. Flag mito % outliers (removed in script 03) -------------------------
  # Mitochondrial genes start with MT-.
  mito_genes <- grep("^MT-", rowData(sce)$Symbol, value = FALSE)
  if (length(mito_genes) == 0) {
    # Fall back to rownames (after uniquify) for species with different prefix.
    mito_genes <- grep("^MT-", rownames(sce), value = FALSE)
  }
  # Per-cell totals: UMIs (sum), genes detected and mitochondrial %.
  qc <- scuttle::perCellQCMetrics(
    sce, subsets = list(Mito = mito_genes)
  )
  colData(sce) <- cbind(colData(sce), qc)
  # Outlier = mitochondrial % more than mito_nmads MADs above the median.
  mito_outlier <- scuttle::isOutlier(
    qc$subsets_Mito_percent, nmads = mito_nmads, type = "higher"
  )
  sce$mito_outlier <- as.logical(mito_outlier)
  log_msg("  mito outliers (>", mito_nmads, " MAD): ",
          sum(sce$mito_outlier), " / ", ncol(sce))

  # --- 6. DropletQC nuclear fraction ------------------------------------------
  # DropletQC uses the BAM to measure the fraction of reads from unspliced
  # (intronic) RNA in each cell. Empty droplets have a low nuclear fraction and
  # few UMIs. Skipped if the BAM is missing.
  bam_path <- resolve_bam_path(info)
  if (file.exists(bam_path)) {
    log_msg("  DropletQC: reading ", bam_path)
    nf <- tryCatch(
      # Nuclear fraction per cell barcode.
      DropletQC::nuclear_fraction_tags(
        bam = bam_path,
        barcodes = colnames(sce),
        tiles = 1, cores = 1, verbose = FALSE
      ),
      error = function(e) {
        log_msg("  DropletQC failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(nf)) {
      sce$nuclear_fraction <- nf$nuclear_fraction
      # Label each barcode as cell or empty droplet, using the thresholds in config.R.
      ec <- DropletQC::identify_empty_drops(
        nf_umi = data.frame(nf = sce$nuclear_fraction,
                            umi = sce$total),
        nf_rescue = dropletqc_nf_rescue,
        umi_rescue = dropletqc_umi_rescue
      )
      sce$dropletqc_cell_status <- ec$cell_status
    }
  } else {
    log_msg("  DropletQC: BAM not found at ", bam_path, " - skipping")
    sce$nuclear_fraction      <- NA_real_
    sce$dropletqc_cell_status <- NA_character_
  }

  # --- 7. Diagnostic plots ----------------------------------------------------
  # Plots go to pipeline/qc_plots/per_capture/dataset_{i}/.
  capture_qc_dir <- file.path(qc_dir, "per_capture", info$dataset)
  dir.create(capture_qc_dir, showWarnings = FALSE, recursive = TRUE)

  # One row per cell with its hashtag call.
  hto_df <- as.data.frame(colData(sce)) %>%
    dplyr::transmute(
      barcode   = colnames(sce),
      capture   = info$capture,
      HTO_best  = factor(HTO_best),
      HTO_confident = ifelse(is.na(HTO_confident), FALSE, HTO_confident),
      HTO_doublet   = ifelse(is.na(HTO_doublet),   FALSE, HTO_doublet)
    )
  if (nrow(hto_df) > 0 && nlevels(hto_df$HTO_best) >= 1) {
    # Fixed colours per hashtag (make_hto_palette() in utils.R).
    hto_pal <- make_hto_palette(levels(hto_df$HTO_best))
    # Cells per hashtag, coloured by whether the call was confident.
    p_hto_stack <- ggplot(hto_df,
        aes(x = HTO_best, fill = HTO_confident)) +
      geom_bar(position = position_stack(reverse = TRUE)) +
      coord_flip() +
      scale_fill_manual(values = c(`FALSE` = "grey80",
                                   `TRUE`  = "steelblue")) +
      labs(title = paste0(info$dataset, " — HTO call (count)"),
           x = "HTO_best", y = "cells") +
      cowplot::theme_cowplot(font_size = 7)
    ggsave(file.path(capture_qc_dir, "hto_stacked_count.png"),
           p_hto_stack, width = 5, height = 4, dpi = 150)

    # Same, as proportions.
    p_hto_fill <- ggplot(hto_df,
        aes(x = HTO_best, fill = HTO_confident)) +
      geom_bar(position = position_fill(reverse = TRUE)) +
      coord_flip() +
      scale_fill_manual(values = c(`FALSE` = "grey80",
                                   `TRUE`  = "steelblue")) +
      labs(title = paste0(info$dataset, " — HTO call (proportion)"),
           x = "HTO_best", y = "fraction") +
      cowplot::theme_cowplot(font_size = 7)
    ggsave(file.path(capture_qc_dir, "hto_stacked_fill.png"),
           p_hto_fill, width = 5, height = 4, dpi = 150)

    # Hashtag mix within the capture.
    p_cap_hto <- ggplot(hto_df,
        aes(x = capture, fill = HTO_best)) +
      geom_bar() +
      scale_fill_manual(values = hto_pal) +
      labs(title = paste0(info$dataset, " — capture × HTO"),
           x = "capture", y = "cells") +
      cowplot::theme_cowplot(font_size = 7)
    ggsave(file.path(capture_qc_dir, "capture_by_hto.png"),
           p_cap_hto, width = 4, height = 4, dpi = 150)
  }

  # DropletQC scatter. Skipped if the BAM was missing.
  if (any(!is.na(sce$nuclear_fraction))) {
    dq_df <- data.frame(
      nuclear_fraction = sce$nuclear_fraction,
      total            = sce$total,
      cell_status      = if (!is.null(sce$dropletqc_cell_status))
                           sce$dropletqc_cell_status else NA_character_
    )
    p_dq <- ggplot(dq_df,
        aes(x = total, y = nuclear_fraction, colour = cell_status)) +
      geom_point(alpha = 0.4, size = 0.6) +
      scale_x_log10() +
      labs(title = paste0(info$dataset, " — DropletQC"),
           x = "total UMI (log10)", y = "nuclear fraction") +
      cowplot::theme_cowplot(font_size = 8) +
      theme(legend.position = "right")
    ggsave(file.path(capture_qc_dir, "dropletqc_scatter.png"),
           p_dq, width = 6, height = 4, dpi = 150)

    # UMI counts by DropletQC status.
    if (!all(is.na(dq_df$cell_status))) {
      p_dq_v <- ggplot(dq_df,
          aes(x = factor(cell_status), y = total, fill = cell_status)) +
        geom_violin(scale = "width") +
        scale_y_log10() +
        labs(title = paste0(info$dataset, " — UMI by cell_status"),
             x = "cell_status", y = "total UMI (log10)") +
        cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none")
      ggsave(file.path(capture_qc_dir, "dropletqc_status_violin.png"),
             p_dq_v, width = 5, height = 4, dpi = 150)
    }
  }

  # --- Save -------------------------------------------------------------------
  # Save to pipeline/SCEs/dataset_{i}_qcd.rds for script 03.
  out_path <- file.path(sce_dir, paste0(info$dataset, "_qcd.rds"))
  saveRDS(sce, out_path)
  log_msg("  saved ", out_path)
  invisible(out_path)
}

# Run every requested capture in turn.
for (i in indices) process_capture(i)

log_msg("== 02_qc_and_demux.R done.")
