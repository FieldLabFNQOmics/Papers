# ==============================================================================
# 4. Clustering, cell-type annotation, differential expression, and GSEA
# ==============================================================================
# Converts the preprocessed SCE to a Seurat object, clusters cells, annotates
# cell types (SingleR + scMCA + marker-based gating), then runs per-cell-type
# differential expression (Sanity + limma) and preranked GSEA, with dot-plot
# and network visualisations.
#
# The culstering and annotation sections were adapted from:
#   https://oshlacklab.com/paed-cf-cite-seq/06_COMBO.clustering_annotation.html
#
# Inputs:
#   data/merged/SCEs/merged.preprocessed.SCE.rds   (from script 3)
# Outputs (under data/merged/):
#   SCEs/merged.clustered.SEU.rds
#   plots/, gsea/, sanity/, visualizations/         (per cell type)
#
# External tools (set the paths below to match your installation):
#   - Sanity            https://github.com/jmbreda/Sanity
#   - GSEA (CLI)        https://www.gsea-msigdb.org/gsea/
#   - MSigDB mouse .gmt collections (m2, m3, m5, mh)
#
# ==============================================================================

# ---- External tool / resource paths (edit these) -----------------------------
SANITY_BIN   <- "/path/to/Sanity"
GSEA_CLI     <- "/path/to/gsea-cli.sh"
GSEA_DB_PATH <- "/path/to/msigdb/mouse"   # directory containing the .gmt files
# ------------------------------------------------------------------------------
# 1. Load libraries
library(here)
library(scater)
library(Seurat)
library(clustree)
library(dplyr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(sctransform)
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(org.Mm.eg.db)
library(SingleR)
library(scRNAseq)
library(scuttle)
library(celldex)
library(scMCA)
library(grid)
library(UpSetR)
library(limma)
library(statmod)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(stringr)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(msigdbr)

set.seed(42)

# 1. Load data
sce <- readRDS(here("data", "merged", "SCEs", "merged.preprocessed.SCE.rds"))

# 2. Identify uninformative gene sets
mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
ribo_set <- grep("^Rp(s|l)", rownames(sce), value = TRUE)
c2_sets <- msigdbr(species = "Mus musculus", category = "C2")
ribo_set <- union(
  ribo_set,
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$gene_symbol)
sex_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))]
pseudogene_set <- rownames(sce)[
  any(grepl("pseudogene", rowData(sce)$ENSEMBL.GENEBIOTYPE))]

# 3. Per-cell QC metrics
is_mito <- rownames(sce) %in% mito_set
is_ribo <- rownames(sce) %in% ribo_set
sce <- addPerCellQC(
  sce,
  subsets = list(Mito = which(is_mito), Ribo = which(is_ribo)))

# Percentage of zero counts per cell.
sce$zero_percent <- colSums(counts(sce) == 0) / nrow(sce)

# Visualise QC metrics.
p1 <- plotColData(sce, "sum", x = "capture", colour_by = "capture", point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(size = 6)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"), mid = unit(0.06, "cm"), long = unit(0.09, "cm"))
p2 <- plotColData(sce, "detected", x = "capture", colour_by = "capture", point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p3 <- plotColData(sce, "subsets_Mito_percent", x = "capture", colour_by = "capture", point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p4 <- plotColData(sce, "subsets_Ribo_percent", x = "capture", colour_by = "capture", point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p5 <- plotColData(sce, x = "capture", y = "zero_percent", colour_by = "capture", point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))

qc <- ((p1 | p2) / (p3 | p4) / p5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Flag uninformative genes (mito, ribo, sex-linked, pseudogenes).
uninformative <- is_mito | is_ribo |
  rownames(sce) %in% sex_set | rownames(sce) %in% pseudogene_set

# 4. Convert to Seurat
counts <- counts(sce)
rownames(counts) <- rowData(sce)$Symbol

# A small number of gene symbols are duplicated; drop the duplicates.
dup_rows <- duplicated(rownames(counts))
counts_clean <- counts[!dup_rows, ]

seu <- CreateSeuratObject(
  counts = counts_clean,
  meta.data = data.frame(colData(sce)))

# Normalise, reduce dimensionality, and embed.
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE, dims = 1:30) %>%
  RunUMAP(verbose = FALSE, dims = 1:30)

# Inspect PCA and dimensionality.
DimHeatmap(seu, dims = 1:50, cells = 1000, balanced = TRUE)
ElbowPlot(seu, ndims = 50)

# 5. Cluster
out <- here("data", "merged", "SCEs", "merged.clustered.SEU.rds")
if (!file.exists(out)) {
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
  seu <- FindClusters(seu, algorithm = 3, resolution = seq(0.1, 1, by = 0.1))
  seu <- RunUMAP(seu, dims = 1:30)
  saveRDS(seu, file = out)
} else {
  seu <- readRDS(out)
}

# Flag cells with expression in <5% of genes (very high zero fraction).
seu@meta.data$group <- ifelse(seu@meta.data$zero_percent > 0.95, "High", "Low")
zero_percent_umap <- DimPlot(
  seu, reduction = "umap", repel = TRUE, pt.size = 1,
  group.by = "group", split.by = "capture", cols = c("blue", "red"))

# Remove high-zero-fraction cells and re-cluster.
seu <- subset(seu, subset = group == "Low")
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
seu <- FindClusters(seu, algorithm = 3, resolution = seq(0.1, 1, by = 0.1))
seu <- RunUMAP(seu, dims = 1:30)
saveRDS(seu, file = out)

# Visualise clustering at the default resolution.
clustered_umap <- DimPlot(
  seu, reduction = "umap", label = TRUE, repel = TRUE,
  label.size = 5, pt.size = 1) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

# Gene-by-cell count matrix (available for external annotation tools).
counts_matrix <- GetAssayData(seu)

# ==============================================================================
# 6. Cell-type annotation
# ==============================================================================
# Reference-based labels are generated with SingleR (celldex mouse reference)
# and scMCA, then combined with marker-based gating to isolate populations of
# interest.

# ---- Reference-based labels --------------------------------------------------
ref <- celldex::MouseRNAseqData()
results <- SingleR(
  test = as.SingleCellExperiment(seu),
  ref = ref,
  labels = ref$label.fine)
seu$singler_labels <- results$labels

singler <- DimPlot(
  seu, reduction = "umap", label = FALSE, repel = TRUE,
  label.size = 5, pt.size = 1, group.by = "singler_labels") +
  NoLegend()

# scMCA labels.
mca_result <- scMCA(scdata = counts_matrix, numbers_plot = 1)
mca_labels <- data.frame(
  Original.index = names(mca_result$scMCA),
  Predicted.cell.types = as.character(mca_result$scMCA))
seu$mca_labels <- mca_labels$Predicted.cell.types

meta <- seu@meta.data

# ==============================================================================
# 7. Marker-based gating of populations of interest
# ==============================================================================
# Each block isolates a Seurat subset for a population, combining cluster
# membership, reference labels, and marker expression. Cluster numbers are
# specific to this dataset's clustering above.

### Neutrophils: cluster 1 granulocytes.
neutrophil_cells <- intersect(
  WhichCells(seu, idents = 1),
  WhichCells(seu, expression = singler_labels == "Granulocytes"))
seu_neutrophils <- subset(seu, cells = neutrophil_cells)

### Alveolar macrophages: cluster 14.
avm_cells <- WhichCells(seu, idents = 14)
seu_lung_macrophage <- subset(seu, cells = avm_cells)

### Macrophages (alveolar + non-alveolar): clusters 14 and 17.
seu_macrophage <- subset(seu, cells = WhichCells(seu, idents = c(14, 17)))

### Non-alveolar macrophages: cluster 17.
navm_cells <- WhichCells(seu, idents = 17)
seu_non_alveolar_macrophages <- subset(seu, cells = navm_cells)

### Dendritic cells: cluster 15, excluding one mislabelled cell type.
dc_cells <- rownames(meta)[
  meta$seurat_clusters %in% c(15) &
    meta$mca_labels != "Dendritic cell_Cst3 high(Mammary-Gland-Virgin)"]
seu_DCs <- subset(seu, cells = dc_cells)

### B cells: selected clusters, SingleR B cells, excluding contaminating types.
b_exclude <- paste(c("Macrophage", "Dividing", "Dendritic", "Stem"), collapse = "|")
bcell_cells <- rownames(meta)[
  meta$seurat_clusters %in% c(0, 3, 4, 11, 12) &
    meta$singler_labels %in% c("B cells") &
    !grepl(b_exclude, meta$mca_labels)]
seu_Bcells <- subset(seu, cells = bcell_cells)

### Kit+ central cluster.
central_clusters <- c(6, 8, 9, 10, 13, 16, 20, 22, 25)
kit_expr <- FetchData(seu, vars = "Kit", slot = "data")$Kit > 0
kit_exclude <- grepl("Dendritic|dendritic|gdT|T-cells|NK", seu$mca_labels)
kit_cells <- rownames(meta)[
  meta$seurat_clusters %in% central_clusters & kit_expr & !kit_exclude]
seu_kit_positive <- subset(seu, cells = kit_cells)

### Kit+ Sca1+ central cluster.
kit_sca1_cells <- rownames(meta)[
  meta$seurat_clusters %in% central_clusters &
    FetchData(seu, vars = "Kit", slot = "data")$Kit > 0 &
    FetchData(seu, vars = "Ly6a", slot = "data")$Ly6a > 0 &
    !grepl("Dendritic|dendritic|gdT|T-cells|NK", seu$mca_labels)]
seu_kit_sca1_positive <- subset(seu, cells = kit_sca1_cells)

### Sca1+ (Kit-independent) central cluster.
sca1_expr <- FetchData(seu, vars = "Ly6a", slot = "data")$Ly6a > 0
sca1_exclude <- grepl(
  "Dendritic|dendritic|gdT|T-cells|T cell|thymocyte|NK|B cell|AT1|AT2|Ciliated",
  seu$mca_labels)
sca1_cells <- rownames(meta)[
  meta$seurat_clusters %in% central_clusters & sca1_expr & !sca1_exclude]
seu_sca1_positive <- subset(seu, cells = sca1_cells)

# ==============================================================================
# 8. Combined custom labels for plotting
# ==============================================================================
seu$custom_labels <- "Other"
seu$custom_labels[rownames(meta) %in% neutrophil_cells]  <- "Neutrophils"
seu$custom_labels[rownames(meta) %in% avm_cells]         <- "Alveolar Macrophages"
seu$custom_labels[rownames(meta) %in% navm_cells]        <- "Non-alveolar Macrophages"
seu$custom_labels[rownames(meta) %in% dc_cells]          <- "DCs"
seu$custom_labels[rownames(meta) %in% kit_sca1_cells]    <- "Kit+, Sca1+"

# Colour palette with "Other" greyed out.
unique_labels <- unique(seu$custom_labels)
colors <- rainbow(length(unique_labels) - 1)
names(colors) <- setdiff(unique_labels, "Other")
colors["Other"] <- "lightgrey"

umap_dir <- here("data", "merged", "results")
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(umap_dir, "umap_custom_labels.pdf"), width = 12, height = 12)
print(DimPlot(seu, group.by = "custom_labels", label = TRUE, cols = colors))
dev.off()

pdf(file.path(umap_dir, "umap_custom_noLabels.pdf"), width = 12, height = 12)
print(DimPlot(seu, group.by = "custom_labels", label = FALSE, cols = colors))
dev.off()

# ==============================================================================
# 9. Differential expression (Sanity + limma) and preranked GSEA
# ==============================================================================
# For a given cell-type subset: estimate denoised expression with Sanity, fit a
# linear model across samples with limma, then run preranked GSEA against the
# mouse MSigDB collections and save volcano plots and rank files.

run_differential_expression_and_gsea <- function(seu_subset,
                                                  cell_type_name,
                                                  gsea_db_path = GSEA_DB_PATH,
                                                  output_base_dir = NULL) {
  if (is.null(output_base_dir)) output_base_dir <- here("data", "merged")

  sanity_dir <- file.path(output_base_dir, "sanity", cell_type_name)
  plot_dir   <- file.path(output_base_dir, "plots", cell_type_name)
  gsea_dir   <- file.path(output_base_dir, "gsea", cell_type_name)
  for (d in c(sanity_dir, plot_dir, gsea_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # Restrict to the 3,000 most variable genes.
  seu_subset <- FindVariableFeatures(seu_subset, selection.method = "vst", nfeatures = 3000)
  features <- SelectIntegrationFeatures(object.list = list(seu_subset), nfeatures = 3000)

  raw_counts <- GetAssayData(seu_subset, layer = "counts")
  filtered_counts <- raw_counts[intersect(rownames(raw_counts), features), ]
  countmatrix <- as.matrix(filtered_counts)

  # Write the count matrix in Sanity's expected format.
  forsanity <- rbind(colnames(countmatrix), countmatrix)
  forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)
  count_out <- file.path(sanity_dir, "countmatrix.txt")
  if (!file.exists(count_out)) {
    write.table(forsanity, file = count_out, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  # Run Sanity if not already done.
  likelihood_out <- file.path(sanity_dir, "likelihood.txt")
  if (!file.exists(likelihood_out)) {
    system(paste(SANITY_BIN, "-n 28 -e 1 -f", count_out, "-d", sanity_dir))
  }

  # Read Sanity output.
  expr_file <- file.path(sanity_dir, "log_transcription_quotients.txt")
  expr <- read.table(expr_file, row.names = 1, header = TRUE, stringsAsFactors = FALSE)
  expr <- data.matrix(expr)
  expr <- expr - min(expr)
  colnames(expr) <- gsub("\\.", "-", colnames(expr))

  # Define contrasts from the Sample metadata.
  seu_subset$contrasts <- ifelse(
    seu_subset$Sample %in% c("SD_BCG"), "SD_BCG",
    ifelse(seu_subset$Sample %in% c("SD_PE25SS"), "SD_PE25SS",
           ifelse(seu_subset$Sample %in% c("SD_RD1"), "SD_RD1", NA)))

  contrasts <- as.character(seu_subset$contrasts)
  expr <- expr[, match(colnames(expr), colnames(seu_subset))]
  keep <- !is.na(contrasts)
  expr <- expr[, keep]
  contrasts <- contrasts[keep]

  design <- model.matrix(~0 + contrasts)
  colnames(design) <- gsub("contrasts", "", colnames(design))

  fit <- lmFit(expr, design)
  contrast_matrix <- makeContrasts(
    SD_PE25SS_vs_SD_BCG  = SD_PE25SS - SD_BCG,
    SD_RD1_vs_SD_BCG     = SD_RD1 - SD_BCG,
    SD_RD1_vs_SD_PE25SS  = SD_RD1 - SD_PE25SS,
    levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  contrasts_list <- colnames(contrast_matrix)
  results_list <- lapply(contrasts_list, function(contrast) {
    topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH")
  })
  names(results_list) <- contrasts_list

  make_volcano_plot <- function(results, contrast_name) {
    results$Significance <- "Not Significant"
    results$Significance[results$adj.P.Val < 0.05 & abs(results$logFC) > 1] <- "Significant"
    label_subset <- subset(results, abs(logFC) > 1.5)

    p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
      geom_point(size = 1.5, alpha = 0.6) +
      scale_color_manual(values = c("grey", "red")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
      geom_text_repel(
        data = label_subset, aes(label = rownames(label_subset)),
        max.overlaps = 10, size = 3) +
      labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value",
           title = paste("Volcano Plot:", contrast_name)) +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.background = element_rect(fill = "white", color = NA))

    ggsave(file.path(plot_dir, paste0("volcano_", contrast_name, ".png")),
           p, width = 8, height = 6, dpi = 300)
    p
  }

  # MSigDB mouse collections (edit the filenames if you use a different release).
  gsea_dbs <- c(
    m2 = "m2.all.v2024.1.Mm.symbols.gmt",
    m3 = "m3.all.v2024.1.Mm.symbols.gmt",
    m5 = "m5.all.v2024.1.Mm.symbols.gmt",
    mh = "mh.all.v2024.1.Mm.symbols.gmt")

  for (contrast in contrasts_list) {
    p <- make_volcano_plot(results_list[[contrast]], contrast)
    print(p)

    result_table <- results_list[[contrast]]
    write.csv(result_table, file = file.path(gsea_dir, paste0(contrast, "_topTable.csv")))
    message(paste("Significant genes in", contrast, ":",
                  sum(result_table$adj.P.Val < 0.05)))

    # Extract gene symbols (handle "ID_Symbol" style rownames).
    if (grepl("_", rownames(result_table)[1])) {
      result_table$gene <- sub(".*_", "", rownames(result_table))
    } else {
      result_table$gene <- rownames(result_table)
    }

    # Build the ranked list (by moderated t-statistic).
    gseatt <- data.frame(gene = result_table$gene, t = result_table$t)
    gseatt <- gseatt[!duplicated(gseatt$gene), ]
    gseatt <- gseatt[order(gseatt$t, decreasing = TRUE), ]
    rank_file <- file.path(gsea_dir, paste0(contrast, "_gsea.rnk"))
    write.table(gseatt, file = rank_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Run preranked GSEA against each collection.
    for (db_name in names(gsea_dbs)) {
      db_file <- file.path(gsea_db_path, gsea_dbs[db_name])
      if (!file.exists(db_file)) {
        warning(paste("Database file not found:", db_file))
        next
      }
      label <- paste0(db_name, ".", contrast)
      gsea_cmd <- paste0(
        GSEA_CLI, " GSEAPreranked ",
        "-gmx ", db_file, " ",
        "-norm meandiv -nperm 10000 ",
        "-rnk ", rank_file, " ",
        "-scoring_scheme weighted ",
        "-rpt_label ", label, " ",
        "-create_svgs false -make_sets true -plot_top_x 40 ",
        "-rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false ",
        "-out ", gsea_dir)
      message(paste("Running GSEA for", contrast, "with database", db_name))
      system(gsea_cmd)
    }
  }

  list(results = results_list, expr = expr, design = design,
       gsea_dir = gsea_dir, plot_dir = plot_dir)
}

# Batch driver: run DE + GSEA over a named list of cell-type subsets.
batch_run_de_and_gsea <- function(subset_list,
                                  gsea_db_path = GSEA_DB_PATH,
                                  output_base_dir = NULL) {
  all_results <- list()
  for (cell_type in names(subset_list)) {
    tryCatch({
      message(paste("===== Processing", cell_type, "====="))
      seu_subset <- subset_list[[cell_type]]
      if (ncol(seu_subset) < 10) {
        warning(paste("Skipping", cell_type, "- too few cells:", ncol(seu_subset)))
        next
      }
      all_results[[cell_type]] <- run_differential_expression_and_gsea(
        seu_subset = seu_subset,
        cell_type_name = cell_type,
        gsea_db_path = gsea_db_path,
        output_base_dir = output_base_dir)
      message(paste("Completed processing for", cell_type))
    }, error = function(e) {
      message(paste("ERROR in", cell_type, ":", e$message))
    })
  }
  all_results
}

subsets <- list(
  neutrophils              = seu_neutrophils,
  lung_macrophage          = seu_lung_macrophage,
  DCs                      = seu_DCs,
  Bcells                   = seu_Bcells,
  kit_positive             = seu_kit_positive,
  kit_sca1_positive        = seu_kit_sca1_positive,
  macrophage               = seu_macrophage,
  sca1_positive            = seu_sca1_positive,
  non_alveolar_macrophages = seu_non_alveolar_macrophages)

all_results <- batch_run_de_and_gsea(subsets)

# ==============================================================================
# 10. GSEA result visualisation (dot plots and network plots)
# ==============================================================================

visualize_gsea_results <- function(gsea_dir, cell_type, output_dir = NULL) {
  if (is.null(output_dir)) output_dir <- gsea_dir
  viz_dir <- file.path(output_dir, cell_type)
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  message("Saving visualizations to: ", viz_dir)

  read_gsea_results <- function(gsea_dir) {
    gsea_result_dirs <- list.files(
      gsea_dir, pattern = "GseaPreranked", full.names = TRUE, include.dirs = TRUE)
    all_results <- list()

    for (result_dir in gsea_result_dirs) {
      report_files <- list.files(
        result_dir, pattern = "gsea_report.*\\.tsv$", full.names = TRUE, recursive = TRUE)

      for (report_file in report_files) {
        if (!file.exists(report_file) || file.size(report_file) == 0) next
        is_pos <- grepl("pos", basename(report_file))
        db_contrast <- str_extract(basename(result_dir), "(m[2-5h]\\.[^.]+)")
        if (!is.na(db_contrast)) {
          parts <- strsplit(db_contrast, "\\.")[[1]]
          db <- parts[1]
          contrast <- paste(parts[-1], collapse = ".")
          result <- tryCatch(
            read.delim(report_file, check.names = FALSE),
            error = function(e) {
              message("Error reading file ", report_file, ": ", e$message)
              NULL
            })
          if (!is.null(result) && nrow(result) > 0) {
            result$database <- db
            result$contrast <- contrast
            result$direction <- ifelse(is_pos, "positive", "negative")
            key <- paste(db, contrast, ifelse(is_pos, "pos", "neg"), sep = "_")
            all_results[[key]] <- result
          }
        }
      }
    }

    if (length(all_results) == 0) {
      message("No GSEA results found in ", gsea_dir)
      return(NULL)
    }

    # Coerce reported columns to numeric, dropping any non-numeric NES rows.
    numeric_columns <- c("NES", "NOM p-val", "FDR q-val", "FWER p-val",
                         "ES", "SIZE", "RANK AT MAX")
    for (i in seq_along(all_results)) {
      if ("NES" %in% names(all_results[[i]]) && is.character(all_results[[i]]$NES)) {
        problematic_rows <- which(is.na(suppressWarnings(as.numeric(all_results[[i]]$NES))))
        if (length(problematic_rows) > 0) {
          all_results[[i]] <- all_results[[i]][-problematic_rows, ]
        }
      }
      for (col in numeric_columns) {
        if (col %in% names(all_results[[i]])) {
          all_results[[i]][[col]] <- as.numeric(as.character(all_results[[i]][[col]]))
        }
      }
    }
    dplyr::bind_rows(all_results)
  }

  gsea_results <- read_gsea_results(gsea_dir)
  if (is.null(gsea_results) || nrow(gsea_results) == 0) {
    message("No GSEA results to visualize.")
    return(NULL)
  }

  gsea_results_clean <- gsea_results %>%
    filter(abs(NES) > 1, `FDR q-val` < 0.25) %>%
    arrange(`FDR q-val`) %>%
    mutate(
      NAME = ifelse(nchar(NAME) > 40, paste0(substr(NAME, 1, 37), "..."), NAME),
      contrast_db = paste(contrast, database, sep = "_"))

  if (nrow(gsea_results_clean) == 0) {
    message("No significant GSEA results to visualize.")
    return(NULL)
  }

  dotplot_data <- gsea_results_clean %>%
    dplyr::select(NAME, NES, `FDR q-val`, contrast, database, direction) %>%
    dplyr::mutate(
      neg_log10_fdr = -log10(`FDR q-val`),
      neg_log10_fdr = pmin(neg_log10_fdr, 10)) %>%
    filter(neg_log10_fdr >= 1)

  create_dotplot <- function(data, contrast_value, db_value) {
    plot_data <- data %>% filter(contrast == contrast_value, database == db_value)
    if (nrow(plot_data) == 0) return(NULL)

    plot_data <- plot_data %>%
      mutate(`FDR q-val` = ifelse(`FDR q-val` == 0, 1e-5, `FDR q-val`),
             NAME_display = stringr::str_trunc(NAME, 50, "right")) %>%
      arrange(NES) %>%
      group_by(NAME_display) %>%
      mutate(NAME_display = if (n() > 1)
        paste0(NAME_display, " (", row_number(), ")") else NAME_display) %>%
      ungroup()
    plot_data$NAME_display <- factor(plot_data$NAME_display, levels = unique(plot_data$NAME_display))

    p <- ggplot(plot_data, aes(x = NES, y = NAME_display,
                               size = neg_log10_fdr, color = `FDR q-val`)) +
      geom_point() +
      scale_color_gradientn(
        colors = colorRampPalette(c("darkred", "red", "orange", "yellow"))(100),
        trans = "log10", name = "FDR q-val",
        guide = guide_colorbar(reverse = TRUE)) +
      scale_size_continuous(range = c(2, 8)) +
      labs(title = paste("GSEA Dotplot:", contrast_value, "-", db_value),
           subtitle = paste("Cell type:", cell_type),
           x = "Normalized Enrichment Score (NES)", y = "Pathway") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = "right")

    filename <- file.path(
      viz_dir, paste0("gsea_dotplot_", contrast_value, "_", db_value, "_", cell_type, ".png"))
    ggsave(filename, p, width = 10, height = 8, dpi = 300)
    message("Saved dotplot: ", filename)
  }

  unique_combinations <- unique(dotplot_data %>% dplyr::select(contrast, database))
  for (i in seq_len(nrow(unique_combinations))) {
    create_dotplot(dotplot_data,
                   unique_combinations$contrast[i],
                   unique_combinations$database[i])
  }

  gsea_results_clean
}

create_network_plot <- function(data, contrast_value, db_filter = NULL,
                                 cell_type, viz_dir, gsea_dir) {
  if (is.null(data)) {
    message("No data provided for network plot")
    return(NULL)
  }

  filtered_data <- data %>% dplyr::filter(contrast == contrast_value)
  if (!is.null(db_filter)) {
    filtered_data <- filtered_data %>% dplyr::filter(database == db_filter)
  }
  if (nrow(filtered_data) == 0) {
    message("No data for network plot of contrast: ", contrast_value)
    return(NULL)
  }

  # Top pathways by |NES| (capped for readability/performance).
  top_pathways <- filtered_data %>% dplyr::slice_max(abs(NES), n = 30)
  pathway_names <- top_pathways$NAME

  # Read gene sets for these pathways from the GMT files.
  gmt_files <- list.files(path = gsea_dir, pattern = "\\.gmt$",
                          recursive = TRUE, full.names = TRUE)
  gene_sets <- list()
  for (gmt_file in gmt_files) {
    for (line in readLines(gmt_file)) {
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) < 3) next
      if (parts[1] %in% pathway_names) gene_sets[[parts[1]]] <- parts[-c(1, 2)]
    }
  }
  if (length(gene_sets) == 0) {
    message("No matching gene sets found. Skipping network plot.")
    return(NULL)
  }

  # Jaccard similarity between pathways (edges above a threshold).
  edges <- expand.grid(names(gene_sets), names(gene_sets)) %>%
    dplyr::filter(Var1 != Var2) %>%
    rowwise() %>%
    mutate(weight = length(intersect(gene_sets[[Var1]], gene_sets[[Var2]])) /
             length(union(gene_sets[[Var1]], gene_sets[[Var2]]))) %>%
    dplyr::filter(weight > 0.1) %>%
    rename(from = Var1, to = Var2)
  if (nrow(edges) == 0) {
    message("No significant pathway connections found. Skipping network plot.")
    return(NULL)
  }

  nodes <- top_pathways %>%
    dplyr::mutate(
      id = NAME,
      neg_log10_fdr = -log10(ifelse(`FDR q-val` == 0, 1e-10, `FDR q-val`)),
      neg_log10_fdr = pmin(neg_log10_fdr, 10)) %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::select(id, NAME, NES, `FDR q-val`, database)

  edges <- edges %>% dplyr::filter(from %in% nodes$id & to %in% nodes$id)
  if (nrow(edges) == 0) {
    message("No matching pathway connections found after validation. Skipping network plot.")
    return(NULL)
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  V(g)$color <- ifelse(V(g)$NES > 0, "red", "blue")
  V(g)$size <- 5 + 2 * abs(V(g)$NES)
  if (ecount(g) > 0) E(g)$width <- 1 + 5 * E(g)$weight

  set.seed(42)
  if (vcount(g) == 0) {
    message("No vertices in the graph. Skipping network plot.")
    return(NULL)
  }

  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(linewidth = weight), alpha = 0.3, colour = "grey") +
    geom_node_point(aes(size = size, color = NES)) +
    geom_node_text(aes(label = NAME), repel = TRUE, size = 3) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(3, 10), name = "-log10(FDR)") +
    scale_edge_width(range = c(0.5, 3)) +
    labs(title = paste("GSEA Network Plot:", contrast_value,
                       if (!is.null(db_filter)) paste("-", db_filter) else ""),
         subtitle = paste("Cell type:", cell_type)) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10))

  tryCatch({
    filename <- if (is.null(db_filter)) {
      paste0("gsea_network_", contrast_value, "_", cell_type, ".png")
    } else {
      paste0("gsea_network_", contrast_value, "_", db_filter, "_", cell_type, ".png")
    }
    ggsave(file.path(viz_dir, filename), p, width = 12, height = 10, dpi = 300)
    message("Saved network plot: ", filename)
  }, error = function(e) {
    message("Error saving network plot: ", e$message)
  })

  p
}

batch_visualize_gsea <- function(base_dir, cell_types) {
  all_visualizations <- list()
  main_viz_dir <- file.path(base_dir, "visualizations")
  dir.create(main_viz_dir, recursive = TRUE, showWarnings = FALSE)

  for (cell_type in cell_types) {
    gsea_dir <- file.path(base_dir, "gsea", cell_type)
    viz_dir <- file.path(main_viz_dir, cell_type)
    dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

    if (!dir.exists(gsea_dir)) {
      message(paste("GSEA directory not found for", cell_type))
      next
    }
    message(paste("Processing GSEA results for", cell_type))

    tryCatch({
      gsea_results_data <- visualize_gsea_results(gsea_dir, cell_type, output_dir = main_viz_dir)
      if (is.null(gsea_results_data)) {
        message(paste("No significant GSEA results found for", cell_type))
        next
      }
      all_visualizations[[cell_type]] <- gsea_results_data

      gsea_dirs <- list.dirs(gsea_dir, full.names = TRUE, recursive = FALSE)
      gsea_dirs <- gsea_dirs[grepl("GseaPreranked", gsea_dirs)]
      if (length(gsea_dirs) == 0) gsea_dirs <- gsea_dir

      unique_combinations <- unique(gsea_results_data %>% dplyr::select(contrast, database))
      for (i in seq_len(nrow(unique_combinations))) {
        contrast_value <- unique_combinations$contrast[i]
        db_value <- unique_combinations$database[i]

        specific_gsea_dir <- gsea_dir
        for (dir in gsea_dirs) {
          if (grepl(paste0(db_value, ".", contrast_value), basename(dir), fixed = TRUE)) {
            specific_gsea_dir <- dir
            break
          }
        }

        tryCatch({
          create_network_plot(
            data = gsea_results_data, contrast_value = contrast_value,
            db_filter = db_value, cell_type = cell_type,
            viz_dir = viz_dir, gsea_dir = specific_gsea_dir)
        }, error = function(e) {
          message(paste("Error creating network plot for",
                        contrast_value, "-", db_value, ":", e$message))
        })
      }
      message(paste("Completed all visualizations for", cell_type))
    }, error = function(e) {
      message(paste("Error processing", cell_type, ":", e$message))
    })
  }
  all_visualizations
}

cell_types <- c("lung_macrophage", "DCs", "Bcells", "kit_positive",
                "kit_sca1_positive", "macrophage", "sca1_positive",
                "non_alveolar_macrophages", "neutrophils")

base_dir <- here("data", "merged")
visualizations <- batch_visualize_gsea(base_dir, cell_types)

# ==============================================================================
# 11. TPM values by sample (example: alveolar macrophages)
# ==============================================================================
# Transcript lengths are retrieved from Ensembl via biomaRt; for genes with
# multiple transcripts the longest is used.

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_names <- rownames(seu_lung_macrophage)

transcript_info <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id",
                 "ensembl_transcript_id", "transcript_length"),
  filters = "external_gene_name",
  values = gene_names,
  mart = ensembl)

gene_transcript_lengths <- transcript_info %>%
  group_by(external_gene_name) %>%
  summarise(max_transcript_length = max(transcript_length, na.rm = TRUE))

gene_lengths <- gene_transcript_lengths$max_transcript_length[
  match(gene_names, gene_transcript_lengths$external_gene_name)]
names(gene_lengths) <- gene_names

valid_genes <- !is.na(gene_lengths)
gene_lengths_clean <- gene_lengths[valid_genes]
gene_names_clean <- names(gene_lengths_clean)
message(paste("Found transcript lengths for", length(gene_lengths_clean),
              "out of", length(gene_names), "genes"))

calculate_tpm <- function(counts, gene_lengths) {
  rate <- counts / (gene_lengths / 1000)       # counts per kb
  t(t(rate) / colSums(rate)) * 1e6             # per million
}

raw_counts <- LayerData(seu_lung_macrophage, assay = "RNA", layer = "counts")[gene_names_clean, ]
tpm_matrix <- calculate_tpm(raw_counts, gene_lengths_clean)

captures_char <- as.character(unique(seu_lung_macrophage$capture))
tpm_by_capture <- matrix(0, nrow = length(gene_names_clean), ncol = length(captures_char))
rownames(tpm_by_capture) <- gene_names_clean
colnames(tpm_by_capture) <- captures_char

for (capture in captures_char) {
  cells_in_capture <- colnames(seu_lung_macrophage)[seu_lung_macrophage$capture == capture]
  cells_in_capture <- intersect(cells_in_capture, colnames(tpm_matrix))
  if (length(cells_in_capture) > 0) {
    tpm_by_capture[, capture] <- rowMeans(tpm_matrix[, cells_in_capture, drop = FALSE])
  }
}

tpm_table <- as.data.frame(tpm_by_capture)
tpm_table$gene <- rownames(tpm_table)
tpm_table <- tpm_table[, c("gene", captures_char)]

write.csv(tpm_table,
          file = here("data", "merged", "results", "alveolar_macrophage_tpm_by_capture.csv"),
          row.names = FALSE)

# sessionInfo()
