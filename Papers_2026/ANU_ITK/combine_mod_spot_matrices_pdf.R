# ------------------------------------------------------------------------------
# combine_mod_spot_matrices_pdf.R
# Combines the single-grid Hallmark spot matrices made by script 06 into one
# PDF, one page per contrast. Script 06 already does this at the end of a run
# unless it is run with --no-pdf, so this is only needed to rebuild the PDF.
#
# Run:  Rscript new_scripts/combine_mod_spot_matrices_pdf.R [options]
#   --output FILE     write the PDF here instead of the default
#   --pattern REGEX   only include contrasts whose name matches REGEX
#
# Input:  pipeline/GSEA/spot_matrices/*_mod_pos_neg.png   (from script 06)
# Output: pipeline/GSEA/spot_matrices/all_mod_spot_matrices.pdf
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
})

# Load settings (config.R) and shared functions (utils.R).
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))

# --- CLI ----------------------------------------------------------------------
# Options typed after the script name.
args <- commandArgs(trailingOnly = TRUE)
# Value that follows an option, or NULL if absent.
get_flag_val <- function(flag) {
  k <- which(args == flag)
  if (length(k) == 1L && length(args) > k) args[k + 1L] else NULL
}
opt_output  <- get_flag_val("--output")
opt_pattern <- get_flag_val("--pattern")

sm_dir <- file.path(gsea_dir, "spot_matrices")

# bundle_mod_spot_matrices_pdf() (utils.R) writes one PNG per page and returns
# the PDF path, or NULL if there were no PNGs.
out <- bundle_mod_spot_matrices_pdf(sm_dir, out_pdf = opt_output,
                                    pattern = opt_pattern)

# Stop if there is nothing to combine.
if (is.null(out)) {
  stop("No *_mod_pos_neg.png files found in ", sm_dir,
       " - has 06_immune_pathway_dotplots.R been run?")
}
