# ------------------------------------------------------------------------------
# 01_create_sce.R
# Reads the Cell Ranger raw matrix for each capture, splits gene expression
# from hashtag (HTO) counts, and runs emptyDrops to score which droplets
# contain cells. Droplets are filtered on this score in script 02.
#
# Run:  Rscript new_scripts/01_create_sce.R [i1 i2 ...]
#   i1 i2 ...   capture numbers 1-20 (default: all 20)
#
# Input:  Cell Ranger raw_feature_bc_matrix for each capture (paths in config.R)
# Output: pipeline/SCEs/dataset_{i}.rds
#         pipeline/emptyDrops/diagnostics/dataset_{i}__barcode_ranks.rds
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(BiocParallel)
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# Captures to process, given as index numbers 1-20 (the i column of
# capture_table in config.R: 1-4 = Aug batch, 5-20 = May batch).
# No arguments = all 20 captures.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  indices <- capture_table$i
} else {
  indices <- suppressWarnings(as.integer(args))
}
# Stop before doing any work if an index is not a valid capture number.
bad <- indices[!(indices %in% capture_table$i) | is.na(indices)]
if (length(bad) > 0L) {
  stop("Bad capture index/indices: ",
       paste(bad, collapse = ", "),
       ". Valid range is 1..", max(capture_table$i),
       " (one per 10x capture). Submit script 01 with `-J 1-10` ",
       "then `-J 11-20`; the 440-element array is for script 05 only.",
       call. = FALSE)
}

# Process one capture: read the Cell Ranger raw matrix, split gene expression
# from HTO counts, run emptyDrops and save the result as an SCE.
process_capture <- function(i) {
  # Dataset ID, capture name, batch and paths for this capture.
  info <- get_capture_info(i)
  log_msg("== Capture ", info$dataset, " (", info$capture, ", batch ",
          info$batch, ")")

  # Read the raw (unfiltered) matrix, which includes empty droplets.
  raw_dir <- info$raw_matrix_dir
  if (!dir.exists(raw_dir)) {
    stop("Raw matrix dir not found: ", raw_dir)
  }
  # Create single cell experiment object
  sce <- DropletUtils::read10xCounts(raw_dir, col.names = TRUE)
  # Tag with canonical dataset ID in colData.
  sce$dataset <- info$dataset
  sce$capture <- info$capture
  sce$batch   <- info$batch

  # Split features: GEX vs HTO. Cellranger 7.1.0 stores feature type in rowData.
  feat_type <- rowData(sce)$Type
  if (is.null(feat_type)) {
    stop("rowData(sce)$Type is NULL for ", info$capture,
         " - cannot split GEX vs HTO.")
  }
  gex_idx <- which(feat_type == "Gene Expression")
  hto_idx <- which(feat_type == "Antibody Capture")
  if (length(gex_idx) == 0) stop("No Gene Expression features in ", info$capture)
  if (length(hto_idx) == 0) stop("No Antibody Capture (HTO) features in ", info$capture)

  # Keep gene expression as the main matrix and store HTO counts as an altExp.
  gex_sce <- sce[gex_idx, ]
  hto_sce <- sce[hto_idx, ]
  altExp(gex_sce, "HTO") <- hto_sce

  # Barcode-rank data for QC plots.
  br <- DropletUtils::barcodeRanks(counts(gex_sce))
  diag_path <- file.path(emptydrops_diag,
                         paste0(info$dataset, "__barcode_ranks.rds"))
  saveRDS(br, diag_path)

  # emptyDrops, seeded for reproducibility.
  log_msg("  running emptyDrops (lower = ", emptydrops_lower, ")")
  ed <- DropletUtils::emptyDrops(
    counts(gex_sce),
    lower  = emptydrops_lower,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 100L)
  )
  # Keep the full emptyDrops result. FDR filtering is done in script 02.
  colData(gex_sce) <- cbind(colData(gex_sce), ed)

  # Save to pipeline/SCEs/dataset_{i}.rds.
  out_path <- file.path(sce_dir, paste0(info$dataset, ".rds"))
  saveRDS(gex_sce, out_path)
  log_msg("  saved ", out_path, " - ",
          ncol(gex_sce), " barcodes, ",
          sum(!is.na(ed$FDR) & ed$FDR < emptydrops_fdr),
          " pass FDR<", emptydrops_fdr)
  invisible(out_path)
}

for (i in indices) {
  process_capture(i)
}

log_msg("== 01_create_sce.R done.")
