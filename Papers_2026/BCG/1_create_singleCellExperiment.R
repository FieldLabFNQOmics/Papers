# ==============================================================================
# 1. Create SingleCellExperiment objects from 10x Cell Ranger output
# ==============================================================================
# Reads the per-sample filtered feature-barcode matrices produced by
# `cellranger multi` and saves one SingleCellExperiment (SCE) per sample.
#
# Adapted from a pipeline by Peter Hickey:
#   https://github.com/Oshlack/paed-cf-cite-seq/blob/submission2/code/C133_Neeland-dropletutils.R
#
# Inputs:
#   cellranger/<sample>/outs/per_sample_outs/<sample>/count/
#       sample_filtered_feature_bc_matrix/   (one per sample)
# Outputs:
#   data/<sample>/SCEs/<sample>.CellRanger.SCE.rds
#
# ==============================================================================

library(DropletUtils)
library(here)

# Sample identifiers. These correspond to the `--id` values used for
# `cellranger multi` and to the sub-directories under cellranger/ and data/.
sample_names <- c("SD_BCG", "SD_PE25SS", "SD_RD1")

for (i in sample_names) {

  dir.create(here("data", i, "SCEs"), recursive = TRUE, showWarnings = FALSE)

  # Use the pre-filtered Cell Ranger output (rather than the raw matrix).
  captures <- setNames(
    here(
      "cellranger",
      i,
      "outs",
      "per_sample_outs",
      i,
      "count",
      "sample_filtered_feature_bc_matrix"),
    i)

  sce <- read10xCounts(samples = captures, col.names = TRUE)
  stopifnot(!anyDuplicated(colnames(sce)))

  # Save the SCE built on the filtered data.
  output_sce_file <- paste0(i, ".CellRanger.SCE.rds")
  saveRDS(
    object = sce,
    file = here("data", i, "SCEs", output_sce_file),
    compress = "xz")
}

# sessionInfo()
