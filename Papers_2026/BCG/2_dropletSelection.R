# ==============================================================================
# 2. Droplet selection, doublet annotation, and sample merging
# ==============================================================================
# Loads the per-sample SCEs, inspects barcode-rank distributions, annotates
# within-sample doublets, and merges the samples into a single SCE.
#
# Adapted from:
#   https://oshlacklab.com/paed-cf-cite-seq/03_C133_Neeland.emptyDrops.html
#
# Inputs:
#   data/<sample>/SCEs/<sample>.CellRanger.SCE.rds   (from script 1)
# Outputs:
#   data/merged/SCEs/merged.doublets.rds
#   data/merged/SCEs/merged.emptyDrops.SCE.rds
#
# ==============================================================================

# 1. Load libraries
library(BiocStyle)
library(tidyverse)
library(here)
library(glue)
library(DropletUtils)
library(scran)
library(scater)
library(scuttle)
library(scds)
library(scDblFinder)
library(Matrix)

set.seed(42)
options(scipen = 999)
# Maximum size of globals exported to workers by the future framework.
options(future.globals.maxSize = 6500 * 1024^2)

sample_names <- c("SD_BCG", "SD_PE25SS", "SD_RD1")

# 2. Load the per-sample SCEs
for (i in sample_names) {
  file_path <- here("data", i, "SCEs", glue("{i}.CellRanger.SCE.rds"))
  sce <- readRDS(file_path)
  sce$capture <- factor(sce$Sample)
  assign(glue("sce_{i}"), sce)
}

# Prefix cell IDs with a per-sample letter so barcodes remain unique after
# merging (the same barcode whitelist is reused across captures).
colnames(sce_SD_BCG)    <- paste0("A-", colnames(sce_SD_BCG))
colnames(sce_SD_PE25SS) <- paste0("B-", colnames(sce_SD_PE25SS))
colnames(sce_SD_RD1)    <- paste0("C-", colnames(sce_SD_RD1))

# Merge the per-sample SCEs.
sce_list <- lapply(sample_names, function(i) get(glue("sce_{i}")))
sce <- do.call(cbind, sce_list)

# 3. Examine barcode-rank distributions per sample
par(mfrow = c(2, 2))
lapply(levels(sce$capture), function(s) {
  sce <- sce[, sce$capture == s]
  bcrank <- barcodeRanks(counts(sce))

  # Plot only unique points for speed.
  uniq <- !duplicated(bcrank$rank)
  plot(
    x = bcrank$rank[uniq],
    y = bcrank$total[uniq],
    log = "xy",
    xlab = "Rank",
    ylab = "Total UMI count",
    main = s,
    cex.lab = 1.2,
    xlim = c(1, 500000),
    ylim = c(1, 200000))
  abline(h = metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
  abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
})

# 4. Call within-sample doublets
dir.create(here("data", "merged", "SCEs"), recursive = TRUE, showWarnings = FALSE)
out <- here("data", "merged", "SCEs", "merged.doublets.rds")

if (!file.exists(out)) {

  sceLst <- sapply(levels(sce$capture), function(cap) {
    # Annotate doublets using the scds three-step process (as in Demuxafy).
    sce1 <- bcds(sce[, sce$capture == cap], retRes = TRUE, estNdbl = TRUE)
    sce1 <- cxds(sce1, retRes = TRUE, estNdbl = TRUE)
    sce1 <- cxds_bcds_hybrid(sce1, estNdbl = TRUE)
    # Annotate doublets using scDblFinder. The expected doublet rate follows
    # the ~0.8% per 1,000 cells rule of thumb for 10x data.
    sce1 <- scDblFinder(sce1, dbr = ncol(sce1) / 1000 * 0.008)
    sce1
  })

  lapply(sceLst, function(s) {
    colData(s) %>%
      data.frame() %>%
      rownames_to_column(var = "cell")
  }) %>%
    bind_rows() %>%
    saveRDS(file = out)
}

# 5. Save the merged SCE
out <- here("data", "merged", "SCEs", "merged.emptyDrops.SCE.rds")
if (!file.exists(out)) saveRDS(sce, file = out)

# sessionInfo()
