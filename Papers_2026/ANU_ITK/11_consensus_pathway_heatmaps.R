# ------------------------------------------------------------------------------
# 11_consensus_pathway_heatmaps.R
# Heatmaps of the consensus pathways (rows) across cell types (columns), one
# per contrast, for Hallmark (script 05 fgsea) and Disease Ontology (script 08
# DO GSEA cache). Tile colour = NES, a dot = not significant (BH padj > 0.05),
# grey = not tested in that cell type.
#
# Two versions of each heatmap:
#   full     every consensus row and every consensus cell type (PNG + PDF)
#   figure   the trimmed rows used in the paper, split into a T-cell panel and
#            an "other" panel (config.R heatmap_figure_panels), as PDFs
# NES is flipped for some contrasts so the patient group reads as red
# (heatmap_orientation() in utils.R). Run script 08 first.
#
# Run:  Rscript new_scripts/11_consensus_pathway_heatmaps.R [options]
#   --contrast NAME                  one contrast only
#   --panel NAME                     one figure panel (tcell or other)
#   --collection cc                  Hallmark fgsea collection (default h)
#   --skip-hallmark / --skip-do      skip one heatmap type
#   --skip-full / --skip-figure      skip the full heatmaps / the figure panels
#   --no-methods                     no methods page in the multi-page PDFs
#   --no-pdf / --no-png              skip the multi-page PDFs / the PNGs
#   --force                          overwrite existing outputs
#   --list                           print the available contrasts and exit
#
# Input:  pipeline/GSEA/*__gsea_h.tsv                  (from script 05)
#         pipeline/dea_plots/cache/*__do_gsea_full.rds (from script 08)
# Output: pipeline/consensus_heatmaps/
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(DOSE)  # needed to read script 08's cached gseDO results
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()
# Consensus pathway lists and figure panels are set in config.R.

# --- 1. CLI -------------------------------------------------------------------
# Options typed after the script name.
args <- commandArgs(trailingOnly = TRUE)
opt_list_only     <- "--list"          %in% args
opt_force         <- "--force"         %in% args
opt_no_pdf        <- "--no-pdf"        %in% args
opt_no_png        <- "--no-png"        %in% args
opt_no_methods    <- "--no-methods"    %in% args
opt_skip_hallmark <- "--skip-hallmark" %in% args
opt_skip_do       <- "--skip-do"       %in% args
opt_skip_full     <- "--skip-full"     %in% args   # skip the full-row heatmaps
opt_skip_figure   <- "--skip-figure"   %in% args   # skip the trimmed figure PDFs

# Value that follows an option, e.g. --panel tcell, or the default.
get_opt <- function(flag, default = NULL) {
  k <- which(args == flag)
  if (length(k) == 1L && length(args) > k) args[k + 1L] else default
}
opt_contrast   <- get_opt("--contrast")
opt_collection <- get_opt("--collection", "h")
opt_panel      <- get_opt("--panel")

# Which column panels to draw for the figure variant.
panels_to_run <- heatmap_figure_panels
if (!is.null(opt_panel)) {
  if (!opt_panel %in% names(heatmap_figure_panels)) {
    stop("--panel ", opt_panel, " not in heatmap_figure_panels (",
         paste(names(heatmap_figure_panels), collapse = ", "), ")")
  }
  panels_to_run <- heatmap_figure_panels[opt_panel]
}

# --- 2. Output paths ----------------------------------------------------------
out_root <- file.path(pipeline_root,
                      if (isTRUE(use_cellsweep)) "consensus_heatmaps_cellsweep"
                                                else "consensus_heatmaps")
hd_dir   <- file.path(out_root, "hallmark")           # full-row PNGs
do_dir   <- file.path(out_root, "do")
# Figure panels are written as PDFs under figure/.
hd_fig_dir <- file.path(out_root, "figure", "hallmark")
do_fig_dir <- file.path(out_root, "figure", "do")
for (d in c(hd_dir, do_dir, hd_fig_dir, do_fig_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

hallmark_pdf <- file.path(out_root, "consensus_hallmark_heatmaps.pdf")
do_pdf       <- file.path(out_root, "consensus_do_heatmaps.pdf")
# One paged figure PDF per kind x panel.
fig_pdf_path <- function(kind, panel_key)
  file.path(out_root, sprintf("consensus_%s_heatmaps_figure_%s.pdf",
                              kind, panel_key))

# NES colour scale runs from -4 to 4.
NES_LIMIT <- 4

# --- 2b. Canvas geometry ------------------------------------------------------
# Canvas size is derived from the grid so the panel fills the page.
CELL_IN_FIGURE <- 0.60  # tile width (in), figure panels
CELL_IN_FULL   <- 0.32  # tile width (in), full heatmaps
# Tile height = width * Y_SQUISH (full heatmaps).
Y_SQUISH <- 0.52
# Tile height = width * Y_SQUISH_FIGURE (figure panels).
Y_SQUISH_FIGURE <- 0.80
# Pale border around every tile.
TILE_BORDER    <- "grey85"
TILE_BORDER_LW <- 0.75
# Space for titles and legends, in inches.
TITLE_IN  <- 1.00           # title + 3-line subtitle above the panel
LEGEND_IN <- 0.85           # NES colourbar + its margin, to the right
AXIS_PAD  <- 0.15
FS_Y      <- 9              # axis.text.y size (pt), must match the theme below
FS_X      <- 8              # axis.text.x size (pt), must match the theme below
# Minimum width so the subtitle is not clipped.
MIN_W_SUBTITLE <- 5.6
# Minimum page size for the multi-page PDFs (room for the methods page).
MIN_W_PAGED <- 8.0
MIN_H_PAGED <- 6.0

# Approximate width in inches of the longest label at a given font size.
text_in <- function(labels, fontsize_pt) {
  if (!length(labels)) return(0)
  AXIS_PAD + max(nchar(labels)) * fontsize_pt * 0.5 / 72
}

# Canvas size from the number and length of row and column labels.
heatmap_canvas <- function(row_labels, col_labels, cell_in, y_squish = Y_SQUISH) {
  w <- text_in(row_labels, FS_Y) + (length(col_labels) + 0.2) * cell_in + LEGEND_IN
  h <- TITLE_IN + (length(row_labels) + 0.2) * cell_in * y_squish +
       text_in(col_labels, FS_X)
  list(width = max(MIN_W_SUBTITLE, w), height = h)
}

# --- 3. Heatmap builder -------------------------------------------------------
# variant "full" = all consensus rows and cell types; "figure" = trimmed rows and
# the columns of one panel. NES sign is set by heatmap_orientation() (utils.R).
# Returns the plot and its width and height.
make_consensus_heatmap <- function(contrast_full, kind, variant = "full",
                                   panel = NULL) {
  stopifnot(variant %in% c("full", "figure"))
  if (variant == "figure" && is.null(panel))
    stop("make_consensus_heatmap(variant='figure') requires a `panel`")

  levels_vec <-
    if (kind == "hallmark")
      (if (variant == "figure") hallmark_figure_levels else hallmark_levels)
    else
      (if (variant == "figure") do_figure_levels       else do_levels)

  # Column keys are de_cell_subsets names; labels are what is displayed.
  if (variant == "figure") {
    col_keys   <- unname(panel$columns)
    col_labels <- names(panel$columns)
  } else {
    col_keys   <- consensus_cell_types
    col_labels <- celltype_label(consensus_cell_types)
  }

  # NES and padj for every consensus row and column (consensus_long() in
  # utils.R). Missing results are NA and show as grey tiles.
  long <- consensus_long(contrast_full, kind, opt_collection,
                         cell_types = col_keys)   # label, cell_type, NES, padj
  if (is.null(long) || !nrow(long)) return(NULL)
  long <- long[long$label %in% levels_vec, , drop = FALSE]      # trim to variant rows
  if (!nrow(long)) return(NULL)

  # Orient NES so the intended cohort reads as positive/red.
  ori <- heatmap_orientation(contrast_full)
  long$NES <- ori$sign * long$NES

  long$row <- factor(long$label, levels = levels_vec)
  long$col <- factor(col_labels[match(long$cell_type, col_keys)],
                     levels = col_labels)

  ttl_kind <- if (kind == "hallmark") "Hallmark consensus"
              else                    "Disease Ontology consensus"
  # contrast_pretty() marks batch-2-only contrasts.
  title <- sprintf("%s — %s", ttl_kind, contrast_pretty(contrast_full))
  if (variant == "figure") title <- sprintf("%s\n%s", title, panel$title)

  # Log any tiles whose |NES| exceeds the colour limit (clamped, not hidden).
  off <- long[!is.na(long$NES) & abs(long$NES) > NES_LIMIT, , drop = FALSE]
  if (nrow(off)) {
    for (i in seq_len(nrow(off))) {
      log_msg("   [colour-clamp] ", title, " :: ", off$label[i], " / ",
              off$cell_type[i], " NES=", round(off$NES[i], 3),
              " (|NES|>", NES_LIMIT, ", clamped to scale end)")
    }
  }

  # Non-significant tiles (have a result but BH padj > 0.05) get a dot.
  ns <- long[!is.na(long$NES) & !is.na(long$padj) & long$padj > 0.05, ,
             drop = FALSE]

  p <- ggplot(long, aes(x = col, y = row, fill = NES)) +
    # Thin pale border on every tile.
    geom_tile(colour = TILE_BORDER, linewidth = TILE_BORDER_LW) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         limits = c(-NES_LIMIT, NES_LIMIT),
                         oob = scales::squish, na.value = "grey92",
                         name = "NES") +
    scale_x_discrete(limits = col_labels, drop = FALSE) +  # labels at the bottom
    scale_y_discrete(limits = rev(levels_vec), drop = FALSE) +
    # Fixed tile shape, flatter in y.
    coord_fixed(ratio = if (variant == "figure") Y_SQUISH_FIGURE else Y_SQUISH,
                clip = "off") +
    labs(title = title,
         subtitle = paste0(
           sprintf(paste0("Positive NES = higher expression in %s    |    ",
                          "Negative NES = higher expression in %s"),
                   ori$pos, ori$neg),
           "\nDot = not significant (BH padj > 0.05); ",
           "grey = not testable. Colour = GSEA NES (not an IPA z-score)."),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title     = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle  = element_text(hjust = 0.5, size = 7, colour = "grey25"),
          # Cell-type labels at the bottom, printed vertically.
          axis.text.x    = element_text(size = FS_X, angle = 90,
                                        hjust = 1, vjust = 0.5),
          axis.text.y    = element_text(size = FS_Y),
          panel.grid     = element_blank(),
          legend.title   = element_text(size = 9),
          legend.text    = element_text(size = 8))

  if (nrow(ns)) {
    p <- p + geom_point(data = ns, aes(x = col, y = row),
                        inherit.aes = FALSE, size = 1.4, colour = "black")
  }

  # Figure panels: fix the absolute tile size (0.60 x 0.48 in).
  if (variant == "figure") {
    p <- patchwork::wrap_plots(
      p,
      widths  = grid::unit((length(col_labels) + 0.2) * CELL_IN_FIGURE, "in"),
      heights = grid::unit((length(levels_vec) + 0.2) * CELL_IN_FIGURE *
                             Y_SQUISH_FIGURE, "in"))
  }

  dims <- if (variant == "figure") {
    heatmap_canvas(levels_vec, col_labels, CELL_IN_FIGURE, Y_SQUISH_FIGURE)
  } else {
    heatmap_canvas(levels_vec, col_labels, CELL_IN_FULL, Y_SQUISH)
  }
  list(plot = p, width = dims$width, height = dims$height)
}

# --- 3b. Methods front page ---------------------------------------------------
# Two-column text page drawn on the open PDF device.
# Left column of the methods page.
methods_lines_left <- function(kind, variant, panel) {
  gs <- if (kind == "hallmark")
    c("GENE SETS — MSigDB HALLMARK",
      "Rows are MSigDB Hallmark gene sets (Liberzon et al. 2015",
      "Cell Systems), matched to the fgsea output by exact gene-set",
      "ID. The trimmed figure panels show the 7 immune-themed sets;",
      "the full diagnostic heatmaps show all 13 consensus sets.")
  else
    c("GENE SETS — DISEASE ONTOLOGY",
      "Rows are Disease Ontology terms (Schriml et al. 2019 Nucleic",
      "Acids Res) scored by DOSE::gseDO (Yu et al. 2015",
      "Bioinformatics). Terms are matched to the consensus list by",
      "exact, case-insensitive description match — never substring.",
      "Values are reused from the shared script-08 cache (computed",
      "at pvalueCutoff = 1) so this figure, the script-08 plots and",
      "the consensus dotplots report identical numbers.",
      "",
      paste0("Gene sets outside ", do_gsea_min_size, "-", do_gsea_max_size,
             " genes are not tested and render"),
      "grey. Size is the term's overlap with the genes tested in",
      "that cell type, so it varies by column. The upper bound was",
      "raised from the DOSE default of 500 because Disease Ontology",
      "is a hierarchy: a broad parent term inherits every",
      "descendant's annotations, and the two broadest consensus",
      "rows were being dropped from the largest cell subsets purely",
      "because those subsets have more genes to overlap with.")

  cols <- if (variant == "figure") {
    c("CELL-TYPE COLUMNS",
      # Wrap the panel's column list.
      strwrap(paste0("This panel: ",
                     paste(names(panel$columns), collapse = ", "), "."),
              width = 62),
      "",
      "Each column is an INDEPENDENT analysis of the cells carrying",
      "the corresponding Azimuth level-2 labels. Activated CD4 =",
      "CD4 TCM + CD4 TEM + CD4 CTL; Activated CD8 = CD8 TCM +",
      "CD8 TEM. Cells are pooled BEFORE pseudobulk aggregation, so",
      "each column is a single limma fit and a single GSEA run — not",
      "an average of the individual subtype results.",
      "",
      "Total CD4 / Total CD8 are the whole-lineage analyses and",
      "include the naive cells, so a Total column is NOT the sum or",
      "mean of the Naive and Activated columns beside it; the three",
      "are separate fits on overlapping cell populations.")
  } else {
    c("CELL-TYPE COLUMNS",
      "Every consensus cell-type subset, including the individual",
      "memory/effector subtypes that the trimmed figure panels fold",
      "into the Activated CD4 / Activated CD8 aggregates. Each",
      "column is an independent limma fit and GSEA run on the cells",
      "carrying the corresponding Azimuth level-2 labels.")
  }

  paste(c(
    "OVERVIEW",
    "Pathway enrichment across immune cell types for each cohort",
    "comparison in the ITK gain-of-function single-cell RNA-seq",
    "cohort (E42K, T504S, D135Y, healthy blood donors, and the",
    "GEM108 pre-/post-tacrolimus patient). One page per comparison.",
    "",
    "DIFFERENTIAL EXPRESSION",
    "Counts are aggregated to pseudobulk per donor x cell type and",
    "tested with limma-voom (Law et al. 2014 Genome Biol; Ritchie",
    "et al. 2015 Nucleic Acids Res) under ~0 + cohort + batch, the",
    "batch term absorbing the Aug-2022 / May-2023 capture rounds.",
    "Pseudobulk rather than per-cell testing follows Squair et al.",
    "2021 Nat Commun and Zimmerman et al. 2021 Nat Commun, which",
    "show per-cell tests treat cells as independent replicates and",
    "badly inflate false positives. p-values are adjusted by",
    "Benjamini-Hochberg (1995 JRSS-B).",
    "",
    "GENE SET ENRICHMENT",
    "Preranked GSEA via fgsea (Korotkevich et al. 2021 bioRxiv;",
    "Subramanian et al. 2005 PNAS) on genes ranked by the limma",
    "moderated t-statistic. NES is the enrichment score normalised",
    "for gene-set size; padj is BH across all sets in the",
    "collection.",
    "",
    gs,
    "",
    cols
  ), collapse = "\n")
}

# Right column of the methods page.
methods_lines_right <- function(kind, variant, panel) {
  paste(c(
    "READING THE HEATMAP",
    "",
    "FILL = normalised enrichment score (NES). Red = higher in the",
    "cohort named in the per-page subtitle; blue = higher in the",
    "other cohort; white = no net enrichment. The scale is fixed at",
    paste0("+/-", NES_LIMIT, " across every page so pages are directly",
           " comparable;"),
    "the largest NES observed anywhere in this dataset is ~3.75, so",
    "in practice nothing is clamped.",
    "",
    "DOT = the gene set was tested in that cell type but is NOT",
    "significant (BH padj > 0.05). This is the de Cevins et al.",
    "2023 Cell Rep Med convention. An undotted tile is significant.",
    "",
    "GREY, NO DOT = no result: the gene set was not testable in",
    "that cell type, usually because too few of its genes were",
    "detected, or the subset had too few cells in enough donors.",
    "",
    "TILE SHAPE AND BORDERS carry no meaning. Tiles are fixed",
    "geometry (coord_fixed) and slightly compressed in the row",
    "direction, with a thin pale border, so that neighbouring",
    "cells of similar colour stay separable. Layout only.",
    "",
    "COLOUR DIRECTION",
    "fgsea returns NES positive for the second cohort in the",
    "contrast name. For the figures the sign is oriented so the",
    "disease / experimental sample reads as positive (red):",
    "GEM108 vs HBD orients to GEM108; GEM108 pre vs post orients to",
    "pre; all other comparisons keep the default orientation. The",
    "subtitle on each page states the resulting direction",
    "explicitly, so it is always authoritative.",
    "",
    "CAVEATS",
    "",
    "1. NES is not an IPA activation z-score. It is a rank-based",
    "   enrichment statistic; it reports coordinated shift of a",
    "   gene set, not inferred pathway activity.",
    "",
    "2. Descriptive contrasts (GEM108 vs HBD) are single-donor on",
    "   the GEM108 side. Those pages rank genes on log fold change",
    "   rather than a moderated t-statistic, because no per-gene",
    "   variance is estimable from n = 1. Their NES is therefore",
    "   not on the same footing as the pseudobulk pages despite",
    "   sharing this colour scale.",
    "",
    "3. Columns are separate analyses, not a partition. Total,",
    "   Naive and Activated columns of the same lineage overlap in",
    "   cells and their significance calls are not independent.",
    "",
    "4. Gene-set size cuts both ways. A set of many hundreds of",
    "   genes approaches the overall distribution shift, so a small",
    "   padj on a broad Disease Ontology parent term carries less",
    "   information per unit of significance than the same padj on",
    "   a narrow term. At the other end, a term sitting near the",
    paste0("   ", do_gsea_min_size, "-gene floor has very little null distribution"),
    "   behind its NES. Compare NES and set size, not padj alone.",
    "",
    "5. Disease Ontology rows are not independent of one another.",
    "   The ontology is a hierarchy and a parent term inherits",
    "   every descendant's gene annotations, so wherever a parent",
    "   and one of its descendants both appear as rows, the two",
    "   tiles are driven by overlapping genes and are not two",
    "   separate pieces of evidence.",
    "",
    "6. Some rows are clinically motivated rather than data-",
    "   driven: they were chosen because patients in these",
    "   families carry the diagnosis, not because the analysis",
    "   surfaced them. They are shown on every panel so the",
    "   panels stay comparable.",
    "",
    "7. Pseudobulk samples are keyed on donor x cell type. Where a",
    "   column pools several Azimuth labels, unequal per-donor cell",
    "   -type composition can contribute to a fold change; treat",
    "   marginal calls in aggregate columns accordingly."
  ), collapse = "\n")
}

METHODS_TOP  <- 0.93   # npc y of the first text line
METHODS_FS   <- 5.4    # preferred font size (pt); shrunk if the text is long
METHODS_LH   <- 1.15   # lineheight multiplier, must match the gpar below

# Draw the methods page.
draw_methods_page <- function(kind, variant, panel = NULL) {
  grid::grid.newpage()
  left  <- methods_lines_left(kind, variant, panel)
  right <- methods_lines_right(kind, variant, panel)

  # Shrink the font if the text would not fit on the page.
  h_in     <- grid::convertHeight(grid::unit(1, "npc"), "in", valueOnly = TRUE)
  avail_in <- METHODS_TOP * h_in - 0.12          # leave room for the page stamp
  measure  <- function(txt) grid::convertHeight(
    grid::grobHeight(grid::textGrob(
      txt, gp = grid::gpar(fontsize = METHODS_FS, fontfamily = "mono",
                           lineheight = METHODS_LH))),
    "in", valueOnly = TRUE)
  h_at_max <- max(vapply(c(left, right), measure, numeric(1)))
  fs       <- if (h_at_max <= 0) METHODS_FS else
    max(3.6, min(METHODS_FS, METHODS_FS * avail_in / h_at_max * 0.98))

  ttl_kind <- if (kind == "hallmark") "Hallmark" else "Disease Ontology"
  ttl_var  <- if (variant == "figure") "figure panels" else "full diagnostic set"
  grid::grid.text(sprintf("Consensus Pathway Heatmaps (%s, %s) — Methods",
                          ttl_kind, ttl_var),
                  x = 0.5, y = 0.975, just = c("centre", "top"),
                  gp = grid::gpar(fontsize = 12, fontface = "bold"))
  for (i in seq_along(c(left, right))) {
    grid::grid.text(c(left, right)[i],
                    x = c(0.025, 0.515)[i], y = METHODS_TOP,
                    just = c("left", "top"),
                    gp = grid::gpar(fontsize = fs, fontfamily = "mono",
                                    lineheight = METHODS_LH))
  }
  invisible(NULL)
}

# --- 4. Discovery -------------------------------------------------------------
# Contrasts with a Hallmark fgsea table (discover_consensus_contrasts() in utils.R).
avail <- discover_consensus_contrasts(opt_collection)

if (opt_list_only) {
  log_msg("Discovered contrasts (Hallmark collection = ", opt_collection, "):")
  cons <- sort(unique(avail$contrast_full))
  if (!length(cons)) log_msg("  (none)") else for (cf in cons) log_msg("  ", cf)
  log_msg("Figure column panels:")
  for (nm in names(heatmap_figure_panels)) {
    p <- heatmap_figure_panels[[nm]]
    log_msg("  ", nm, ": ", paste(names(p$columns), collapse = ", "))
  }
  quit(save = "no", status = 0L)
}

if (!nrow(avail)) {
  stop("No Hallmark GSEA TSVs found for cell types {",
       paste(consensus_cell_types, collapse = ", "), "} under ", gsea_dir,
       ". Run script 05 for these subsets first.")
}

# Warn if any figure-panel subset has no results from script 05.
missing_subsets <- setdiff(unlist(lapply(heatmap_figure_panels,
                                         function(p) unname(p$columns))),
                           unique(avail$subset))
if (length(missing_subsets)) {
  log_msg("!! No GSEA results for figure-panel subset(s): ",
          paste(missing_subsets, collapse = ", "),
          " - those columns will render all-grey. Run script 05 (--list to",
          " regenerate tasks.tsv, then the new task ids) and script 08 first.")
}

# Contrasts to plot: all found, or the one given with --contrast.
contrasts_to_run <- sort(unique(avail$contrast_full))
if (!is.null(opt_contrast)) {
  if (!opt_contrast %in% contrasts_to_run) {
    stop("--contrast ", opt_contrast, " not found among discovered contrasts (",
         paste(contrasts_to_run, collapse = ", "), ")")
  }
  contrasts_to_run <- opt_contrast
}

# --- 5. Drive -----------------------------------------------------------------
log_msg("== 11_consensus_pathway_heatmaps.R: ", length(contrasts_to_run),
        " contrast(s); full columns = ",
        paste(celltype_label(consensus_cell_types), collapse = ", "))

plot_records <- list()
# Keep a heatmap and its size for the multi-page PDFs.
push_record <- function(kind, variant, panel_key, contrast, res) {
  plot_records[[length(plot_records) + 1L]] <<- list(
    kind = kind, variant = variant, panel_key = panel_key,
    contrast = contrast, plot = res$plot,
    width = res$width, height = res$height)
}

# Save a PNG unless --no-png, or unless it exists and --force is not set.
png_save <- function(file, res) {
  if (opt_no_png) return(invisible(NULL))
  if (!opt_force && file.exists(file)) return(invisible(NULL))
  ggplot2::ggsave(file, res$plot, width = res$width, height = res$height,
                  dpi = 200)
  invisible(file)
}
# Single-page PDF per contrast and panel (no methods page).
pdf_save <- function(file, res) {
  if (!opt_force && file.exists(file)) return(invisible(NULL))
  ggplot2::ggsave(file, res$plot, width = res$width, height = res$height,
                  device = grDevices::cairo_pdf)
  invisible(file)
}

# The two heatmap types and where their files go.
kinds <- list(
  list(kind = "hallmark", label = "Hallmark",
       skip = opt_skip_hallmark, full_dir = hd_dir, fig_dir = hd_fig_dir),
  list(kind = "do",       label = "Disease Ontology",
       skip = opt_skip_do,       full_dir = do_dir, fig_dir = do_fig_dir))

# For each contrast: full heatmap and figure panels, for Hallmark and DO.
for (cf in contrasts_to_run) {
  log_msg("== Contrast: ", cf, " (", contrast_pretty(cf), ")")
  for (k in kinds) {
    if (isTRUE(k$skip)) next

    # Full-row diagnostic heatmap (PNG + paged PDF): all consensus columns.
    if (!opt_skip_full) {
      res <- tryCatch(make_consensus_heatmap(cf, k$kind, "full"),
                      error = function(e) {
                        log_msg("   [", k$label, "] ERROR: ",
                                conditionMessage(e)); NULL })
      if (!is.null(res)) {
        out <- file.path(k$full_dir, paste0(cf, ".png"))
        if (!is.null(png_save(out, res)))
          log_msg("   [", k$label, "] wrote ", basename(out))
        push_record(k$kind, "full", NA_character_, cf, res)
      }
    }

    # Figure panels, one per column panel.
    if (!opt_skip_figure) {
      for (pn in names(panels_to_run)) {
        panel <- panels_to_run[[pn]]
        res <- tryCatch(make_consensus_heatmap(cf, k$kind, "figure", panel),
                        error = function(e) {
                          log_msg("   [", k$label, " figure/", pn, "] ERROR: ",
                                  conditionMessage(e)); NULL })
        if (is.null(res)) next
        out <- file.path(k$fig_dir, paste0(cf, "__", pn, ".pdf"))
        if (!is.null(pdf_save(out, res)))
          log_msg("   [", k$label, " figure] wrote ",
                  file.path("figure", k$kind, basename(out)))
        push_record(k$kind, "figure", pn, cf, res)
      }
    }
  }
}

# --- 6. Paged PDFs ------------------------------------------------------------
# All pages in one PDF share a size. The methods page is page 1.
write_pdf <- function(records, path, label, kind, variant, panel = NULL) {
  if (!length(records)) {
    log_msg("== no ", label, " heatmaps - skipping ", basename(path))
    return(invisible(NULL))
  }
  if (!opt_force && file.exists(path)) {
    log_msg("== exists, skipping (use --force): ", basename(path))
    return(invisible(NULL))
  }
  w <- max(MIN_W_PAGED, vapply(records, function(r) r$width,  numeric(1)))
  h <- max(MIN_H_PAGED, vapply(records, function(r) r$height, numeric(1)))
  log_msg("== Assembling ", label, " PDF: ", path,
          sprintf(" (%.2f x %.2f in)", w, h))
  grDevices::cairo_pdf(path, width = w, height = h, onefile = TRUE)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  pg_tot <- length(records) + (if (opt_no_methods) 0L else 1L)
  pg_n   <- 0L
  stamp  <- function() {
    pg_n <<- pg_n + 1L
    grid::grid.text(sprintf("Page %d of %d", pg_n, pg_tot),
                    x = 0.01, y = 0.01, just = c("left", "bottom"),
                    gp = grid::gpar(fontsize = 8, col = "grey50"))
  }

  if (!opt_no_methods) {
    draw_methods_page(kind, variant, panel)
    stamp()
  }
  for (rec in records) {
    print(rec$plot)
    grid::upViewport(0)
    stamp()
  }
  grDevices::dev.off()
  log_msg("== PDF done: ", path)
}

if (!opt_no_pdf) {
  # Heatmaps of one type, version and panel.
  sel <- function(kind, variant, panel_key = NA_character_)
    Filter(function(r) r$kind == kind && r$variant == variant &&
                       identical(r$panel_key, panel_key), plot_records)

  write_pdf(sel("hallmark", "full"), hallmark_pdf, "Hallmark",
            "hallmark", "full")
  write_pdf(sel("do",       "full"), do_pdf,       "Disease Ontology",
            "do",       "full")

  for (pn in names(panels_to_run)) {
    panel <- panels_to_run[[pn]]
    write_pdf(sel("hallmark", "figure", pn), fig_pdf_path("hallmark", pn),
              paste0("Hallmark (figure / ", pn, ")"), "hallmark", "figure", panel)
    write_pdf(sel("do", "figure", pn), fig_pdf_path("do", pn),
              paste0("Disease Ontology (figure / ", pn, ")"), "do", "figure", panel)
  }
}

log_msg("== 11_consensus_pathway_heatmaps.R complete. (",
        length(plot_records), " heatmap(s))")
