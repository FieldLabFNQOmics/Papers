# ==============================================================================
# 3. Gene annotation, quality control, and outlier removal
# ==============================================================================
# Annotates genes (Ensembl + NCBI), computes per-cell QC metrics, removes
# outlier cells and doublets, and writes a QC report.
#
# Adapted from:
#   https://oshlacklab.com/paed-cf-cite-seq/05_C133_Neeland.preprocess.html
#
# Inputs:
#   data/merged/SCEs/merged.emptyDrops.SCE.rds   (from script 2)
#   data/merged/SCEs/merged.doublets.rds         (from script 2)
# Outputs:
#   data/merged/results/QC_metrics_and_plots.pdf
#   data/merged/SCEs/merged.preprocessed.SCE.rds
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
library(janitor)
library(cowplot)
library(patchwork)
library(scales)
library(DropletQC)
library(AnnotationHub)
library(ensembldb)
library(gridExtra)
library(kableExtra)
library(grid)

sce <- readRDS(here("data", "merged", "SCEs", "merged.emptyDrops.SCE.rds"))

# 2. Gene annotation
# Use a project-local AnnotationHub cache so runs are self-contained.
ah_cache <- here(".cache", "R", "AnnotationHub")
dir.create(ah_cache, recursive = TRUE, showWarnings = FALSE)
setAnnotationHubOption("CACHE", ah_cache)
ah <- AnnotationHub(ask = FALSE)

# EnsDb for mouse on GRCm38 (mm10), matching the Cell Ranger reference.
# Record AH89211 was selected as the latest GRCm38 EnsDb at time of analysis;
# re-run the query below if this record is unavailable in your AnnotationHub
# snapshot, and update the accession accordingly.
#   query_result   <- query(ah, "Mus musculus")
#   ensdb_records  <- query_result[grep("EnsDb", query_result$rdataclass)]
#   filtered       <- subset(ensdb_records, genome == "GRCm38")
ensdb_id <- "AH89211"
EnsDb.Mmusculus <- ah[[ensdb_id]]

# Gene-based annotation on the SCE.
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

# Chromosome location (used to flag mitochondrial genes). Gene version numbers
# are stripped from the IDs prior to lookup.
location <- mapIds(
  x = EnsDb.Mmusculus,
  keys = rowData(sce)$ID,
  keytype = "GENEID",
  column = "SEQNAME")
rowData(sce)$CHR <- location

# Additional Ensembl gene metadata.
ensdb_columns <- c(
  "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
stopifnot(all(ensdb_columns %in% columns(EnsDb.Mmusculus)))
ensdb_df <- DataFrame(
  lapply(ensdb_columns, function(column) {
    mapIds(
      x = EnsDb.Mmusculus,
      keys = rowData(sce)$ID,
      keytype = "GENEID",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = rowData(sce)$ID)
# GENEID cannot be looked up with a GENEID key, so add it manually.
ensdb_df$ENSEMBL.GENEID <- rowData(sce)$ID

# Additional NCBI gene metadata.
library(Mus.musculus)
ncbi_columns <- c("ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
stopifnot(all(ncbi_columns %in% columns(Mus.musculus)))
ncbi_df <- DataFrame(
  lapply(ncbi_columns, function(column) {
    mapIds(
      x = Mus.musculus,
      keys = rowData(sce)$ID,
      keytype = "ENSEMBL",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = rowData(sce)$ID)
rowData(sce) <- cbind(rowData(sce), ensdb_df, ncbi_df)

# Useful gene sets.
mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
ribo_set <- grep("^Rp(s|l)", rownames(sce), value = TRUE)
library(msigdbr)
c2_sets <- msigdbr(species = "Mus musculus", category = "C2")
ribo_set <- union(
  ribo_set,
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$gene_symbol)
is_ribo <- rownames(sce) %in% ribo_set
sex_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))]
pseudogene_set <- rownames(sce)[
  any(grepl("pseudogene", rowData(sce)$ENSEMBL.GENEBIOTYPE))]

# 3. Quality control
# Define and visualise per-cell QC metrics.
is_mito <- rownames(sce) %in% mito_set
is_ribo <- rownames(sce) %in% ribo_set
sce <- addPerCellQC(
  sce,
  subsets = list(Mito = which(is_mito), Ribo = which(is_ribo)))

# Flag doublets identified in script 2.
doublets <- readRDS(here("data", "merged", "SCEs", "merged.doublets.rds"))
doublet_indices <- doublets$scDblFinder.class == "doublet"
sce$Sample[doublet_indices] <- "Doublet"
sce$Sample <- factor(
  sce$Sample,
  levels = c(levels(sce$Sample), "SD_BCG", "SD_PE25SS", "SD_RD1", "Doublet"))

p1 <- plotColData(
  sce, "sum", x = "Sample",
  other_fields = c("capture", "Sample"),
  colour_by = "capture", point_size = 0.5) +
  scale_y_log10() +
  labs(x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- plotColData(
  sce, "detected", x = "Sample",
  other_fields = c("capture", "Sample"),
  colour_by = "capture", point_size = 0.5) +
  labs(x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- plotColData(
  sce, "subsets_Mito_percent", x = "Sample",
  other_fields = c("capture", "Sample"),
  colour_by = "capture", point_size = 0.5) +
  labs(x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4 <- plotColData(
  sce, "subsets_Ribo_percent", x = "Sample",
  other_fields = c("capture", "Sample"),
  colour_by = "capture", point_size = 0.5) +
  labs(x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
qc <- p1 + p2 + p3 + p4 + plot_layout(guides = "collect", ncol = 2)

# Identify outliers by each metric (per capture).
mito_drop <- isOutlier(
  metric = sce$subsets_Mito_percent,
  nmads = 3, type = "higher",
  batch = sce$capture,
  subset = !grepl("Unknown", sce$capture))
mito_drop_df <- data.frame(
  sample = factor(
    colnames(attributes(mito_drop)$thresholds),
    levels(sce$capture)),
  lower = attributes(mito_drop)$thresholds["higher", ])
ribo_drop <- isOutlier(
  metric = sce$subsets_Ribo_percent,
  nmads = 3, type = "higher",
  batch = sce$capture,
  subset = !grepl("Unknown", sce$capture))
ribo_drop_df <- data.frame(
  sample = factor(
    colnames(attributes(ribo_drop)$thresholds),
    levels(sce$capture)),
  lower = attributes(ribo_drop)$thresholds["higher", ])

qc_cutoffs_df <- dplyr::inner_join(mito_drop_df, ribo_drop_df, by = "sample")
colnames(qc_cutoffs_df) <- c("capture", "%mito", "%ribo")
T1 <- inner_join(
  qc_cutoffs_df,
  distinct(as.data.frame(colData(sce)[, c("capture"), drop = FALSE])),
  by = "capture") %>%
  dplyr::select(capture, everything()) %>%
  arrange(capture)

sce_pre_QC_outlier_removal <- sce

# Number of cells flagged by each metric.
T2 <- as.data.frame(table(ribo_drop))
T3 <- as.data.frame(table(mito_drop))

# Exclude cells flagged by mito, ribo, or doublet status.
keep <- !mito_drop & !ribo_drop & !doublet_indices
sce_pre_QC_outlier_removal$keep <- keep
sce <- sce[, keep]
T4 <- data.frame(
  ByMito = tapply(mito_drop, sce_pre_QC_outlier_removal$capture, sum, na.rm = TRUE),
  ByRibo = tapply(ribo_drop, sce_pre_QC_outlier_removal$capture, sum, na.rm = TRUE),
  ByDoublet = tapply(doublet_indices, sce_pre_QC_outlier_removal$capture, sum, na.rm = TRUE),
  Remaining = as.vector(unname(table(sce$capture)))) %>%
  tibble::rownames_to_column("capture") %>%
  dplyr::arrange(dplyr::desc(Remaining))

# Check that QC has not preferentially removed a biological subpopulation.
lost <- calculateAverage(counts(sce_pre_QC_outlier_removal)[, !keep])
kept <- calculateAverage(counts(sce_pre_QC_outlier_removal)[, keep])
library(edgeR)
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)

is_mito <- rownames(sce) %in% mito_set
is_ribo <- rownames(sce) %in% ribo_set
sub_pop <- ggplot() +
  geom_point(aes(x = abundance, y = logFC), pch = 16) +
  labs(x = "Average count", y = "Log-FC (lost/kept)") +
  geom_point(
    data = data.frame(abundance = abundance[is_mito], logFC = logFC[is_mito]),
    aes(x = abundance, y = logFC, color = "Mito_Drop"), pch = 16) +
  geom_point(
    data = data.frame(abundance = abundance[is_ribo], logFC = logFC[is_ribo]),
    aes(x = abundance, y = logFC, color = "Ribo_Drop"), pch = 16) +
  geom_hline(yintercept = c(-1, 1), col = "red", linetype = 2) +
  scale_color_manual(
    values = c("Mito_Drop" = "dodgerblue", "Ribo_Drop" = "orange"),
    name = "Groups") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Check whether removed cells derive preferentially from particular captures.
# Some colnames are non-unique (shared barcodes across captures), so make.unique.
colnames(sce_pre_QC_outlier_removal) <- make.unique(colnames(sce_pre_QC_outlier_removal))
dropped_by_capture <- ggcells(sce_pre_QC_outlier_removal) +
  geom_bar(aes(x = capture, fill = keep)) +
  ylab("Number of droplets") +
  theme_cowplot(font_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(capture ~ ., scales = "free_y")

# Compare QC metrics of discarded vs retained droplets.
p1 <- plotColData(
  sce_pre_QC_outlier_removal, "sum", x = "Sample",
  colour_by = "keep", point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"), mid = unit(0.06, "cm"), long = unit(0.09, "cm"))
p2 <- plotColData(
  sce_pre_QC_outlier_removal, "detected", x = "Sample",
  colour_by = "keep", point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- plotColData(
  sce_pre_QC_outlier_removal, "subsets_Mito_percent", x = "Sample",
  colour_by = "keep", point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4 <- plotColData(
  sce_pre_QC_outlier_removal, "subsets_Ribo_percent", x = "Sample",
  colour_by = "keep", point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
qc_discarded <- p1 + p2 + p3 + p4 + plot_layout(guides = "collect")

# NOTE: A DropletQC nuclear-fraction step was evaluated here but is not
# applicable to this dataset: the 10x Flex chemistry is probe-based and yields
# no intronic reads, so nuclear-fraction scores cannot be computed.

# 4. Assemble the QC report
T1_grob <- tableGrob(T1)
T2_grob <- tableGrob(T2)
T3_grob <- tableGrob(T3)
T4_grob <- tableGrob(T4)

heading_T1 <- textGrob("Table 1: Sample-specific QC metric cutoffs")
heading_T2 <- textGrob("Table 2: No. of cells dropped due to high ribosomal content")
heading_T3 <- textGrob("Table 3: No. of cells dropped due to high mitochondrial content")
heading_T4 <- textGrob("Table 4: Summary of cells dropped by each QC step")

table_arrangement <- arrangeGrob(
  heading_T1, T1_grob, heading_T2, T2_grob,
  heading_T3, T3_grob, heading_T4, T4_grob,
  ncol = 1,
  heights = c(0.5, 5, 0.5, 5, 0.5, 5, 0.5, 5))

dir.create(here("data", "merged", "results"), recursive = TRUE, showWarnings = FALSE)
pdf_file <- here("data", "merged", "results", "QC_metrics_and_plots.pdf")
pdf(pdf_file, width = 8.27, height = 11.69)  # A4
grid.draw(table_arrangement)
print(qc)
print(sub_pop)
print(dropped_by_capture)
print(qc_discarded)
dev.off()

# 5. Save the preprocessed SCE
out <- here("data", "merged", "SCEs", "merged.preprocessed.SCE.rds")
if (!file.exists(out)) saveRDS(sce, out)

# sessionInfo()
