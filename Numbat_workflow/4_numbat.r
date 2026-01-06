library(numbat)
library(here)

# Load custom reference from 2_custom_ref.r
custom_ref <- readRDS("ref.rds")

# Load exp from 1_qc.r
file_path <- here("countsData.RData")
load(file_path)

# Load allele counts created with 3_pileup_and_phase.sh
allele <- read.table(gzfile(here("allele_counts.tsv.gz")), header = TRUE, sep = "\t")

# Run Numbat iteratively with loosening parameters
attempts <- list(
  list(min_LLR = 25, max_entropy = 0.7, name = "First"),
  list(min_LLR = 10, max_entropy = 0.7, name = "Second"),
  list(min_LLR = 5, max_entropy = 0.7, name = "Third"),
  list(min_LLR = 5, max_entropy = 0.5, name = "Fourth")
)

out <- NULL
for (i in seq_along(attempts)) {
  params <- attempts[[i]]
  tryCatch({
    print(paste0(params$name, " attempt (min_LLR=", params$min_LLR, 
                 ", max_entropy=", params$max_entropy, ")..."))
    
    out <- run_numbat(exp, custom_ref, allele, genome = "hg38", t = 1e-05, 
                      ncores = 12, plot = TRUE, out_dir = here("."),
                      max_entropy = params$max_entropy, init_k = 4,
                      min_cells = 10, tau = 0.2, min_LLR = params$min_LLR, 
                      max_iter = 3)
    
    if (!is.null(out)) {  # Check if result is valid
      print(paste0(params$name, " run_numbat succeeded!"))
      break
    }
  }, error = function(e) {
    print(paste0(params$name, " run_numbat failed: ", e$message))
  })
}

if (is.null(out)) {
  stop("All numbat attempts failed")
}