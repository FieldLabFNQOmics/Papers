# ------------------------------------------------------------------------------
# 10_consensus_pathway_dotplots.R
# Dotplots of a fixed set of pathways for every cell type and contrast, so the
# same rows appear in every plot whether or not they are significant:
#   Hallmark   consensus Hallmark pathways, from the script 05 fgsea results
#   DO         consensus Disease Ontology terms, from the script 08 DO GSEA cache
# The pathway lists are set in config.R. x = NES, dot size = gene-set size,
# colour = BH-adjusted p. All plots of one kind share the same colour and size
# scales. Run script 08 first.
#
# Run:  Rscript new_scripts/10_consensus_pathway_dotplots.R [options]
#   --contrast NAME / --celltype NAME   restrict to one contrast / cell type
#   --collection cc                     Hallmark fgsea collection (default h)
#   --skip-hallmark / --skip-do         skip one plot type
#   --no-pdf / --no-png                 skip the PDFs / the PNGs
#   --force                             overwrite existing PNGs
#   --list                              print what is available and exit
#
# Input:  pipeline/GSEA/*__gsea_h.tsv                (from script 05)
#         pipeline/dea_plots/cache/*__do_gsea_full.rds (from script 08)
# Output: pipeline/consensus_dotplots/
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
# Consensus pathway lists (hallmark_consensus, do_consensus) are in config.R.

# --- 1. CLI -------------------------------------------------------------------
# Options typed after the script name.
args <- commandArgs(trailingOnly = TRUE)
opt_list_only     <- "--list"          %in% args
opt_force         <- "--force"         %in% args
opt_no_pdf        <- "--no-pdf"        %in% args
opt_no_png        <- "--no-png"        %in% args
opt_skip_hallmark <- "--skip-hallmark" %in% args
opt_skip_do       <- "--skip-do"       %in% args

# Value that follows an option, e.g. --contrast NAME, or the default.
get_opt <- function(flag, default = NULL) {
  k <- which(args == flag)
  if (length(k) == 1L && length(args) > k) args[k + 1L] else default
}
opt_contrast   <- get_opt("--contrast")
opt_celltype   <- get_opt("--celltype")
opt_collection <- get_opt("--collection", "h")

# Cell types from config.R.
dea_cell_types <- consensus_cell_types
if (!is.null(opt_celltype)) {
  if (!opt_celltype %in% dea_cell_types) {
    stop("--celltype ", opt_celltype, " not in supported list: ",
         paste(dea_cell_types, collapse = ", "))
  }
  dea_cell_types <- opt_celltype
}

# --- 2. Output paths ----------------------------------------------------------
out_root  <- file.path(pipeline_root,
                       if (isTRUE(use_cellsweep)) "consensus_dotplots_cellsweep"
                                                 else "consensus_dotplots")
hd_dir    <- file.path(out_root, "hallmark")
do_dir    <- file.path(out_root, "do")
for (d in c(hd_dir, do_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# One PDF for Hallmark, one for DO.
hallmark_pdf <- file.path(out_root, "consensus_hallmark_dotplots.pdf")
do_pdf       <- file.path(out_root, "consensus_do_dotplots.pdf")

# 8 x 8 in plots.
PLOT_W <- 8; PLOT_H <- 8

# --- 3. Dotplot builders ------------------------------------------------------
# d needs columns label, NES, padj, setSize. Rows with no data stay on the axis
# without a dot.
consensus_dotplot <- function(d, y_levels, title, subtitle,
                              padj_limits, size_limits) {
  d <- d[!is.na(d$NES), , drop = FALSE]
  # x-axis fixed at NES -4..4. Any dot outside that range is logged.
  xmax <- if (nrow(d) && max(abs(d$NES)) > 3) 4 else 3
  off  <- d[abs(d$NES) > xmax, , drop = FALSE]
  if (nrow(off)) {
    for (i in seq_len(nrow(off))) {
      log_msg("   [off-scale] ", title, " :: ", off$label[i],
              " NES=", round(off$NES[i], 3), " (|NES|>", xmax, ", clipped)")
    }
  }
  # padj can be exactly 0. Floor it so the dot is still drawn.
  if (nrow(d)) d$padj <- pmax(d$padj, padj_limits[1])
  ggplot(d, aes(x = NES, y = label)) +
    geom_point(aes(size = setSize, colour = padj)) +
    scale_y_discrete(limits = rev(y_levels), drop = FALSE) +
    scale_colour_gradient(low = "red", high = "blue", name = "p.adjust",
                          trans = "log10", limits = padj_limits,
                          guide = guide_colourbar(reverse = TRUE)) +
    scale_size_continuous(range = c(3, 10), name = "Set size",
                          limits = size_limits) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    coord_cartesian(xlim = c(-xmax, xmax)) +
    labs(title = title, subtitle = subtitle, x = "NES", y = NULL) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey25"),
          axis.text.y   = element_text(size = 9),
          axis.text.x   = element_text(size = 9),
          axis.title.x  = element_text(size = 10),
          legend.title  = element_text(size = 9),
          legend.text   = element_text(size = 8))
}

# Hallmark plot data for one cell type and contrast: NES, padj and set size for
# each consensus pathway, plus title and caption. Only the data is built here so
# all plots can later share one colour and size scale. contrast_pretty(),
# celltype_label() and direction_caption() are in utils.R.
hallmark_dotplot_data <- function(gsea_df, contrast_full, cell_type) {
  # Take the gene-set IDs in y-axis label order.
  ids <- unname(hallmark_consensus[hallmark_levels])
  idx <- match(ids, gsea_df$pathway)
  d <- data.frame(
    label   = hallmark_levels,
    NES     = gsea_df$NES[idx],
    padj    = gsea_df$padj[idx],
    setSize = if ("size" %in% colnames(gsea_df)) gsea_df$size[idx] else NA_real_,
    stringsAsFactors = FALSE)
  d$label <- factor(d$label, levels = hallmark_levels)
  list(d = d, y_levels = hallmark_levels,
       title = sprintf("Hallmark consensus — %s   %s",
                       contrast_pretty(contrast_full), celltype_label(cell_type)),
       subtitle = direction_caption(contrast_full))
}

# Same for the DO terms. match_do_consensus() (utils.R) finds each consensus
# term in the gseDO result by exact name. Terms not tested stay empty.
do_dotplot_data <- function(do_df, contrast_full, cell_type) {
  matched <- match_do_consensus(do_df)
  d <- data.frame(label = do_levels, NES = NA_real_, padj = NA_real_,
                  setSize = NA_real_, stringsAsFactors = FALSE)
  if (nrow(matched)) {
    mi <- match(matched$label, d$label)
    d$NES[mi]     <- matched$NES
    d$padj[mi]    <- matched$padj
    d$setSize[mi] <- matched$setSize
  }
  d$label <- factor(d$label, levels = do_levels)
  list(d = d, y_levels = do_levels,
       title = sprintf("Disease Ontology consensus — %s   %s",
                       contrast_pretty(contrast_full), celltype_label(cell_type)),
       subtitle = direction_caption(contrast_full))
}

# Shared padj and set-size limits across all plots of one kind.
kind_scale_limits <- function(data_items) {
  padj <- unlist(lapply(data_items, function(it) it$d$padj))
  size <- unlist(lapply(data_items, function(it) it$d$setSize))
  padj <- padj[!is.na(padj)]; size <- size[!is.na(size)]
  pos  <- padj[padj > 0]
  list(padj = c(if (length(pos)) min(pos) else 1e-300,
                if (length(padj)) max(padj) else 1),
       size = c(if (length(size)) min(size) else 0,
                if (length(size)) max(size) else 1))
}

# --- 4. Discovery -------------------------------------------------------------
# Subset/contrast pairs that have a Hallmark fgsea table
# (discover_consensus_contrasts() in utils.R).
avail <- discover_consensus_contrasts(opt_collection, dea_cell_types)

if (opt_list_only) {
  log_msg("Discovered cell_type x contrast pairs (Hallmark collection = ",
          opt_collection, "):")
  if (!nrow(avail)) log_msg("  (none)") else print(avail, row.names = FALSE)
  quit(save = "no", status = 0L)
}

if (!nrow(avail)) {
  stop("No Hallmark GSEA TSVs found for cell types {",
       paste(dea_cell_types, collapse = ", "), "} under ", gsea_dir,
       ". Run script 05 for these subsets first.")
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
log_msg("== 10_consensus_pathway_dotplots.R: ",
        length(contrasts_to_run), " contrast(s) x ",
        length(dea_cell_types), " cell type(s)")
log_msg("   Hallmark collection: ", opt_collection,
        "   output: ", out_root)

plot_records <- list()
# Keep a plot for the PDFs.
push_record <- function(kind, contrast, cell_type, plot) {
  plot_records[[length(plot_records) + 1L]] <<- list(
    kind = kind, contrast = contrast, cell_type = cell_type, plot = plot)
}

# Save a PNG unless --no-png, or unless it exists and --force is not set.
png_save <- function(file, plot) {
  if (opt_no_png) return(invisible(NULL))
  if (!opt_force && file.exists(file)) return(invisible(NULL))
  ggplot2::ggsave(file, plot, width = PLOT_W, height = PLOT_H, dpi = 200)
  invisible(file)
}

# --- Pass 1: plot data for every contrast and cell type -----------------------
hallmark_items <- list()   # each: list(cf, ct, data)
do_items       <- list()
for (cf in contrasts_to_run) {
  log_msg("== Contrast: ", cf, " (", contrast_pretty(cf), ")")
  for (ct in dea_cell_types) {
    if (!opt_skip_hallmark) {
      # read_gsea() (utils.R) reads the fgsea results table.
      gsea_df <- read_gsea(ct, cf, opt_collection)
      if (is.null(gsea_df) || !nrow(gsea_df)) {
        log_msg("   [", ct, "] [Hallmark] no GSEA TSV - skipping")
      } else {
        dat <- tryCatch(hallmark_dotplot_data(gsea_df, cf, ct),
                        error = function(e) {
                          log_msg("   [", ct, "] [Hallmark] ERROR: ",
                                  conditionMessage(e)); NULL })
        if (!is.null(dat))
          hallmark_items[[length(hallmark_items) + 1L]] <-
            list(cf = cf, ct = ct, data = dat)
      }
    }
    if (!opt_skip_do) {
      # read_do_gsea_full() (utils.R) reads script 08's cached DO GSEA result.
      do_df <- tryCatch(read_do_gsea_full(ct, cf),
                        error = function(e) {
                          log_msg("   [", ct, "] [DO] ERROR: ",
                                  conditionMessage(e)); NULL })
      dat <- tryCatch(do_dotplot_data(do_df, cf, ct),
                      error = function(e) {
                        log_msg("   [", ct, "] [DO] ERROR: ",
                                conditionMessage(e)); NULL })
      if (!is.null(dat))
        do_items[[length(do_items) + 1L]] <- list(cf = cf, ct = ct, data = dat)
    }
  }
}

# --- Compute the shared padj + set-size scales per kind -----------------------
hall_lim <- kind_scale_limits(lapply(hallmark_items, function(x) x$data))
do_lim   <- kind_scale_limits(lapply(do_items,       function(x) x$data))
log_msg("== Shared Hallmark scale: padj [", signif(hall_lim$padj[1], 3), ", ",
        signif(hall_lim$padj[2], 3), "]  set size [", hall_lim$size[1], ", ",
        hall_lim$size[2], "]")
log_msg("== Shared DO scale:       padj [", signif(do_lim$padj[1], 3), ", ",
        signif(do_lim$padj[2], 3), "]  set size [", do_lim$size[1], ", ",
        do_lim$size[2], "]")

# --- Pass 2: render each plot with the shared scales, save PNG ----------------
render_items <- function(items, kind, lim, out_dir) {
  for (it in items) {
    p <- tryCatch(consensus_dotplot(it$data$d, it$data$y_levels, it$data$title,
                                    it$data$subtitle, lim$padj, lim$size),
                  error = function(e) {
                    log_msg("   [", it$ct, "] [", kind, "] render ERROR: ",
                            conditionMessage(e)); NULL })
    if (is.null(p)) next
    out <- file.path(out_dir, paste0(it$cf, "__", it$ct, ".png"))
    saved <- png_save(out, p)
    if (!is.null(saved)) log_msg("   [", it$ct, "] [", kind, "] wrote ",
                                 basename(out))
    push_record(kind, it$cf, it$ct, p)
  }
}
render_items(hallmark_items, "hallmark", hall_lim, hd_dir)
render_items(do_items,       "do",       do_lim,   do_dir)

# --- 6. Paged PDFs ------------------------------------------------------------
# One PDF per plot kind, with page numbers.
write_pdf <- function(records, path, label) {
  if (!length(records)) {
    log_msg("== no ", label, " plots - skipping ", basename(path))
    return(invisible(NULL))
  }
  log_msg("== Assembling ", label, " PDF: ", path)
  grDevices::cairo_pdf(path, width = PLOT_W, height = PLOT_H, onefile = TRUE)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  pg_n  <- 0L
  pg_tot <- length(records)
  for (rec in records) {
    print(rec$plot)
    pg_n <- pg_n + 1L
    grid::upViewport(0)   # reset to device root so the stamp hits the page
    grid::grid.text(sprintf("Page %d of %d", pg_n, pg_tot),
                    x = 0.01, y = 0.01, just = c("left", "bottom"),
                    gp = grid::gpar(fontsize = 8, col = "grey50"))
  }
  grDevices::dev.off()
  log_msg("== PDF done: ", path)
}

if (!opt_no_pdf) {
  is_kind <- function(k) Filter(function(r) r$kind == k, plot_records)
  write_pdf(is_kind("hallmark"), hallmark_pdf, "Hallmark")
  write_pdf(is_kind("do"),       do_pdf,       "Disease Ontology")
}

log_msg("== 10_consensus_pathway_dotplots.R complete. (",
        length(plot_records), " plot record(s))")
