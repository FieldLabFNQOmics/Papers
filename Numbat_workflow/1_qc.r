library(Seurat)
library(here)
library(Azimuth)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(biomaRt)
library(scuttle)
library(glue)
library(scds)
library(scDblFinder)
library(magrittr)
library(dplyr)
library(UpSetR)

# Sample names and paths
paths <- "/path/to/sample_filtered_feature_bc_matrix"

##### Remove doublets and damaged cells #####
#### Doublets ####
### Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(seu)

### Run scds pipeline
sce <- bcds(sce, retRes = TRUE, estNdbl = TRUE)

# Try to run cxds
cat("Running cxds...\n")
cxds_success <- FALSE
try({
  sce <- cxds(sce, retRes = TRUE, estNdbl = TRUE)
  cxds_success <- TRUE
}, silent = TRUE)

# Run hybrid if cxds worked
if (cxds_success && "cxds_score" %in% colnames(colData(sce))) {
  cat("Running scds hybrid method...\n")
  sce <- cxds_bcds_hybrid(sce, estNdbl = TRUE)
} else {
  cat("cxds failed, using bcds scores only...\n")
  sce$hybrid_score <- sce$bcds_score
}

### Run scDblFinder using 10x's standard doublet rate estimate of approximately 0.8% per 1,000 cells loaded
cat("Running scDblFinder...\n")
sce <- scDblFinder(sce, dbr = ncol(sce)/1000*0.008)

# Step 2: Create doublet calls for each method
# Define threshold for scds (top 5% as doublets, following typical practice)
scds_threshold <- quantile(sce$hybrid_score, 0.95)
sce$scds_singlet <- sce$hybrid_score <= scds_threshold
sce$scds_doublet <- sce$hybrid_score > scds_threshold
sce$scDblFinder_singlet <- sce$scDblFinder.class == "singlet"
sce$scDblFinder_doublet <- sce$scDblFinder.class == "doublet"

### Prepare data for upset plot
# Extract metadata
metadata <- colData(sce) %>% data.frame()

# Create binary matrix for UpSetR (1 = singlet, 0 = doublet)
upset_data_singlet <- metadata %>%
  select(scds_singlet, scDblFinder_singlet) %>%
  mutate(
    scds = as.numeric(scds_singlet),
    scDblFinder = as.numeric(scDblFinder_singlet)
  ) %>%
  select(scds, scDblFinder)

upset_data_doublet <- metadata %>%
  select(scds_doublet, scDblFinder_doublet) %>%
  mutate(
    scds = as.numeric(scds_doublet),
    scDblFinder = as.numeric(scDblFinder_doublet)
  ) %>%
  select(scds, scDblFinder)

### Create upset plot
cat("Creating upset plot...\n")

# Create the upset plot
upset_plot_singlet <- upset(
  upset_data_singlet,
  sets = c("scds", "scDblFinder"),
  sets.bar.color = c("lightblue", "lightgreen"),
  main.bar.color = "steelblue",
  matrix.color = "darkblue",
  text.scale = 1.5,
  point.size = 3,
  line.size = 1,
  mainbar.y.label = "Number of Cells",
  sets.x.label = "Total Singlets per Method",
  order.by = "freq",
  keep.order = TRUE
)

upset_plot_doublet <- upset(
  upset_data_doublet,
  sets = c("scds", "scDblFinder"),
  sets.bar.color = c("lightblue", "lightgreen"),
  main.bar.color = "steelblue",
  matrix.color = "darkblue",
  text.scale = 1.5,
  point.size = 3,
  line.size = 1,
  mainbar.y.label = "Number of Cells",
  sets.x.label = "Total Doublets per Method",
  order.by = "freq",
  keep.order = TRUE
)


out <- here("singlets.png")
png(out, width = 3000, height = 3000, res = 300)
print(upset_plot_singlet)
dev.off()

out <- here("doublets.png")
png(out, width = 3000, height = 3000, res = 300)
print(upset_plot_doublet)
dev.off()

### Remove doublets from sce - union approach
# Print summary before filtering
cat("Before filtering:\n")
cat("Total cells:", ncol(sce), "\n")
cat("scds doublets:", sum(sce$scds_doublet), "\n")
cat("scDblFinder doublets:", sum(sce$scDblFinder_doublet), "\n")
cat("Either method calls doublet:", sum(sce$scds_doublet | sce$scDblFinder_doublet), "\n")
cat("Both methods call doublet:", sum(sce$scds_doublet & sce$scDblFinder_doublet), "\n")

# Create filtering criteria - remove if EITHER method calls doublet
doublets_to_remove <- sce$scds_doublet | sce$scDblFinder_doublet

# Filter the SCE object - keep cells that are NOT flagged as doublets by either method
sce_filtered <- sce[, !doublets_to_remove]

# Print summary after filtering
cat("\nAfter filtering:\n")
cat("Remaining cells:", ncol(sce_filtered), "\n")
cat("Removed cells:", ncol(sce) - ncol(sce_filtered), "\n")
cat("Percentage removed:", round((ncol(sce) - ncol(sce_filtered))/ncol(sce) * 100, 2), "%\n")

# Verify the filtering worked correctly
cat("\nVerification:\n")
cat("Remaining cells should all be singlets by both methods...\n")
cat("scds singlets in filtered data:", sum(sce_filtered$scds_singlet), "\n")
cat("scDblFinder singlets in filtered data:", sum(sce_filtered$scDblFinder_singlet), "\n")
cat("Both methods agree singlet:", sum(sce_filtered$scds_singlet & sce_filtered$scDblFinder_singlet), "\n")

##### Remove damaged cells ######
# Get chromosome info for your genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
chr_info <- getBM(attributes = c("hgnc_symbol", "chromosome_name"),
                  filters = "hgnc_symbol", 
                  values = rownames(sce_filtered),
                  mart = mart)

# Match and add to rowData
rowData(sce_filtered)$CHR <- chr_info$chromosome_name[match(rownames(sce_filtered), chr_info$hgnc_symbol)]

# Identify mito set and percentage
mito_set <- rownames(sce_filtered)[which(rowData(sce_filtered)$CHR == "MT")]
is_mito <- rownames(sce_filtered) %in% mito_set
summary(is_mito)
sce_filtered <- addPerCellQC(
  sce_filtered, 
  subsets = list(Mito = which(is_mito)))

# Identify mito outliers
mito_drop <- isOutlier(
  metric = sce_filtered$subsets_Mito_percent, 
  nmads = 3, 
  type = "higher")

cat("Before mito filtering:\n")
cat("Total cells:", ncol(sce_filtered), "\n")
cat("Cells flagged for high mito:", sum(mito_drop), "\n")

# Filter out the high mito cells (keep cells where mito_drop is FALSE)
sce_filtered_filtered <- sce_filtered[, !mito_drop]

cat("\nAfter mito filtering:\n")
cat("Remaining cells:", ncol(sce_filtered_filtered), "\n")
cat("Removed cells:", ncol(sce_filtered) - ncol(sce_filtered_filtered), "\n")
cat("Percentage removed:", round((ncol(sce_filtered) - ncol(sce_filtered_filtered))/ncol(sce_filtered) * 100, 2), "%\n")

# Verify filtering worked
cat("\nVerification - remaining cells should have lower mito percentages:\n")
cat("Max mito % before filtering:", max(sce_filtered$subsets_Mito_percent), "\n")
cat("Max mito % after filtering:", max(sce_filtered_filtered$subsets_Mito_percent), "\n")
cat("Mito threshold was:", attr(mito_drop, "thresholds")["higher"], "\n")

### Checking for removal of biologically relevant subpopulations
lost <- calculateAverage(counts(sce_filtered)[, !mito_drop])
kept <- calculateAverage(counts(sce_filtered)[, mito_drop])
library(edgeR)
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)

is_mito <- rownames(sce_filtered) %in% mito_set

out <- here("mt_content.pdf")
pdf(out, width = 10, height = 10)
par(mfrow = c(1, 1))
plot(
  abundance,
  logFC,
  xlab = "Average count",
  ylab = "Log-FC (lost/kept)",
  pch = 16)
points(
  abundance[is_mito],
  logFC[is_mito],
  col = "dodgerblue",
  pch = 16,
  cex = 1)

legend(
  "topleft", # Position of the legend
  legend = c("High mitochondrial content"),
  col = c("dodgerblue"),
  pch = 16,
  title = "Legend",
  x.intersp = 0.1, # Reduce spacing between legend symbols and text
  text.width = 1 # Reduce the width allocated for legend text
)

##### Get expression matrix for Numbat #####
seu <- as.Seurat(sce_filtered_filtered, 
                 counts = "counts",
                 data = NULL)

exp <- as.matrix(GetAssayData(seu, layer = "counts"))
file_path <- here("countsData.RData")
save(exp, file = file_path)

out <- here("seu.RData")
save(seu, file = out)