# ------------------------------------------------------------------------------
# umaps_by_cohort.R
# UMAPs of each cohort, drawn separately for T cells and non-T cells and
# coloured by Azimuth L2 cell type. Each panel is saved as its own PDF (with
# and without a legend) and all panels are collated into one PDF.
# Also writes Supplementary Tables 6 and 7: cell counts and % per Azimuth L2
# cell type for the four groups shown in the paper.
#
# Run:  Rscript new_scripts/umaps_by_cohort.R
#
# Input:  pipeline/merged/seurat_merged_harmony_azimuth.rds   (from script 04)
# Output: pipeline/qc_plots/umaps_by_cohort/
# ------------------------------------------------------------------------------

# Load a local GLPK build before Seurat (needed on Gadi). Edit or remove.
if (file.exists("/path/to/libglpk.so.40")) dyn.load("/path/to/libglpk.so.40")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(grid)
})

# Load settings (config.R).
source(here::here("new_scripts", "config.R"))

# --- Settings -----------------------------------------------------------------
# Harmony UMAP from script 04.
reduction_name <- "umap"

# Output directory for the PDFs.
out_dir <- file.path(qc_dir, "umaps_by_cohort")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save each panel with and without a legend.
save_legend_versions <- TRUE

# 8 x 8 in plots.
pdf_w <- 8; pdf_h <- 8; pdf_dpi <- 300

# --- Cell-type palette --------------------------------------------------------
# Fixed colour for each Azimuth L2 cell type.
my_palette <- c(
  "CD4 TCM" = "#0000FF", "NK" = "#FF0000", "B naive" = "#00FF00", "B intermediate" = "#000033",
  "CD14 Mono" = "#201A01", "CD16 Mono" = "#DD00FF", "CD8 TEM" = "#005300", "CD8 TCM" = "#009FFF",
  "MAIT" = "#9A4D42", "NK_CD56bright" = "#00FFBE", "CD4 TEM" = "#783FC1", "B memory" = "#1F9698",
  "NK Proliferating" = "#FFACFD", "CD8 Naive" = "#B1CC71", "gdT" = "#F1085C", "Platelet" = "#FE8F42",
  "CD4 Naive" = "#FF00B6", "pDC" = "#720055", "cDC2" = "#766C95", "CD4 CTL" = "#02AD24",
  "Plasmablast" = "#C8FF00", "ILC" = "#886C00", "ASDC" = "#FFB79F", "dnT" = "#858567", "HSPC" = "#A10300",
  "Eryth" = "#14F9FF", "cDC1" = "#00479E", "Treg" = "#FFD300",
  # Proliferating T subtypes added to the palette.
  "CD4 Proliferating" = "#B15928", "CD8 Proliferating" = "#6A3D9A"
)

# Azimuth L2 labels on the T-cell panel (config.R umap_tcell_labels). MAIT, dnT
# and gdT are on the T-cell panel here, unlike the heatmaps and violins.
t_cell_types <- umap_tcell_labels

# --- Load annotated object ----------------------------------------------------
# Path to the annotated object from script 04 (annotated_rds_path() in config.R).
in_path <- annotated_rds_path()
if (!file.exists(in_path)) {
  stop("Annotated object not found: ", in_path,
       "\n  Did 04_integrate_annotate.R finish? (use_cellsweep = ",
       use_cellsweep, ")")
}
message("== Loading ", in_path)
seu <- readRDS(in_path)

if (!reduction_name %in% Reductions(seu)) {
  stop("Reduction '", reduction_name, "' not in object. Available: ",
       paste(Reductions(seu), collapse = ", "))
}

# Drop unused factor levels so labels/colours only reflect cells present.
seu$predicted.celltype.l2 <- forcats::fct_drop(
  factor(seu$predicted.celltype.l2))

# T-cell vs non-T-cell tag.
seu$cell_category <- ifelse(
  as.character(seu$predicted.celltype.l2) %in% t_cell_types,
  "T_cells", "Non_T_cells")

# --- Define per-cohort subsets ------------------------------------------------
# GEM108 is shown separately from E42K_affected.
md <- seu@meta.data
is_gem108 <- md$patient == "GEM108"
tx        <- md$treatment

# Cell names in each cohort panel.
subset_cells <- list(
  HBD              = rownames(md)[md$cohort == "HBD"],
  E42K_affected    = rownames(md)[md$cohort == "E42K_affected" & !is_gem108],
  # E42K_carriers = E42K_affected + E42K_unaffected, GEM108 pre-treatment only.
  E42K_carriers    = rownames(md)[
    (md$cohort %in% c("E42K_affected", "E42K_unaffected") & !is_gem108) |
    (is_gem108 & !is.na(tx) & tx == "pre_treatment")
  ],
  E42K_unaffected  = rownames(md)[md$cohort == "E42K_unaffected"],
  T504S            = rownames(md)[md$cohort == "T504S"],
  D135Y            = rownames(md)[md$cohort == "D135Y"],
  # GEM108 pre- and post-treatment cells. The May-batch GEM108 samples have no
  # treatment label and are left out.
  GEM108           = rownames(md)[is_gem108 & !is.na(tx) &
                                    tx %in% c("pre_treatment", "post_treatment")],
  GEM108_pre       = rownames(md)[is_gem108 & !is.na(tx) & tx == "pre_treatment"],
  GEM108_post      = rownames(md)[is_gem108 & !is.na(tx) & tx == "post_treatment"],
  GEM108_prePost   = rownames(md)[is_gem108 & !is.na(tx) &
                                    tx %in% c("pre_treatment", "post_treatment")]
)

# --- Supplementary Tables 6 and 7 ---------------------------------------------
# Cell counts (Table 6) and % of each group's cells (Table 7) per Azimuth L2
# label, from the same cells as the UMAP panels.
#   A.I.1                           = GEM108 pre + post
#   A.II.1, A.II.3, A.II.4, A.III.1 = PMAI0017, PMAI0018, PMAI0023, PMAI0024
# Stop if any of the four E42K family members has no cells.
e42k_family <- c("PMAI0017", "PMAI0018", "PMAI0023", "PMAI0024")
missing_fam <- setdiff(e42k_family, unique(md$patient))
if (length(missing_fam)) {
  stop("Supplementary tables: no cells for ", paste(missing_fam, collapse = ", "))
}
# The four table columns.
supp_groups <- list(
  "HBD"                              = subset_cells$HBD,
  "A.I.1"                            = subset_cells$GEM108,
  "A.II.1, A.II.3, A.II.4, A.III.1"  = rownames(md)[md$patient %in% e42k_family],
  "ITK T504S"                        = subset_cells$T504S
)
# Count cells per cell type in each group, then convert to % of the group.
all_types <- sort(levels(md$predicted.celltype.l2))
supp_counts <- sapply(supp_groups, function(cells)
  as.integer(table(factor(md[cells, "predicted.celltype.l2"],
                          levels = all_types))))
rownames(supp_counts) <- all_types
supp_pct <- round(100 * sweep(supp_counts, 2, colSums(supp_counts), "/"), 2)

# Write a table as TSV, optionally with a Total row.
write_supp <- function(m, file, total_row) {
  df <- data.frame("Cell type" = rownames(m), m, check.names = FALSE,
                   stringsAsFactors = FALSE)
  if (total_row) {
    df <- rbind(df, data.frame("Cell type" = "Total", t(colSums(m)),
                               check.names = FALSE))
  }
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  saved ", file)
}
write_supp(supp_counts,
           file.path(out_dir, "supp_table6_celltype_counts.tsv"), TRUE)
write_supp(supp_pct,
           file.path(out_dir, "supp_table7_celltype_percent.tsv"), FALSE)

# --- Plot helper --------------------------------------------------------------
# One UMAP panel for a cohort and cell category.
make_umap <- function(cells, title, with_legend = FALSE) {
  if (length(cells) == 0) return(NULL)
  obj <- subset(seu, cells = cells)
  obj$predicted.celltype.l2 <- forcats::fct_drop(obj$predicted.celltype.l2)
  p <- DimPlot(
    obj,
    reduction  = reduction_name,
    order      = c("Treg"),  # Treg drawn last (on top)
    shuffle    = TRUE,
    seed       = 42,
    label      = TRUE,
    repel      = TRUE,
    label.size = 3,
    pt.size    = 1,
    group.by   = "predicted.celltype.l2",
    cols       = my_palette
  ) + ggtitle(title)
  if (!with_legend) p <- p + NoLegend()
  p
}

# --- Generate and save --------------------------------------------------------
# Every panel is also kept for the collated PDF.
collated <- list()

# For each cohort: a T-cell panel and a non-T-cell panel.
for (cohort_name in names(subset_cells)) {
  all_cells <- subset_cells[[cohort_name]]
  if (length(all_cells) == 0) {
    message("  [skip] ", cohort_name, " - no cells")
    next
  }

  cats <- list(
    Tcells    = intersect(all_cells, rownames(md)[md$cell_category == "T_cells"]),
    NonTcells = intersect(all_cells, rownames(md)[md$cell_category == "Non_T_cells"])
  )

  for (cat_name in names(cats)) {
    cells <- cats[[cat_name]]
    if (length(cells) == 0) {
      message("  [skip] ", cohort_name, " ", cat_name, " - no cells")
      next
    }
    label <- paste0(cohort_name, " - ",
                    ifelse(cat_name == "Tcells", "T cells", "Non-T cells"),
                    " (n = ", length(cells), ")")

    p <- make_umap(cells, label, with_legend = FALSE)
    f <- file.path(out_dir, paste0(cohort_name, "_UMAP_", cat_name, ".pdf"))
    ggsave(f, p, width = pdf_w, height = pdf_h, dpi = pdf_dpi, device = cairo_pdf)
    message("  saved ", f)
    collated[[length(collated) + 1]] <- list(plot = p, label = label)

    if (save_legend_versions) {
      pl <- make_umap(cells, label, with_legend = TRUE)
      fl <- file.path(out_dir,
                      paste0(cohort_name, "_UMAP_", cat_name, "_legend.pdf"))
      ggsave(fl, pl, width = pdf_w, height = pdf_h, dpi = pdf_dpi, device = cairo_pdf)
      message("  saved ", fl)
      collated[[length(collated) + 1]] <- list(plot = pl,
                                               label = paste0(label, " [legend]"))
    }
  }
}

# --- Collated multi-page PDF --------------------------------------------------
# One UMAP per page, with page numbers.
if (length(collated) > 0) {
  collated_pdf <- file.path(out_dir, "umaps_by_cohort_all.pdf")
  n_pages <- length(collated)
  cairo_pdf(collated_pdf, width = pdf_w, height = pdf_h, onefile = TRUE)
  for (i in seq_along(collated)) {
    print(collated[[i]]$plot)
    grid.text(
      sprintf("Page %d of %d", i, n_pages),
      x = unit(0.5, "npc"), y = unit(4, "mm"),
      gp = gpar(fontsize = 8, col = "grey40")
    )
  }
  invisible(dev.off())
  message("  saved ", collated_pdf, " (", n_pages, " pages)")
}

message("== umaps_by_cohort.R done. Output in ", out_dir)
