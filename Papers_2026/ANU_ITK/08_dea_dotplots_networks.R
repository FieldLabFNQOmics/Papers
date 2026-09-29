# ------------------------------------------------------------------------------
# 08_dea_dotplots_networks.R
# Enrichment plots for each cell type and contrast:
#   1. Hallmark GSEA dotplot        top pathways from the script 05 fgsea results
#   2. Disease Ontology GSEA        DOSE::gseDO on genes ranked by t (or logFC)
#   3. Disease Ontology ORA         DOSE::enrichDO on genes with P < 0.01
#   4. DO ORA network plot          genes linked to the top DO terms (cnetplot)
# All plots also go into one PDF report. The full DO GSEA result (all terms,
# not only significant ones) is cached and reused by scripts 10 and 11.
#
# Run:  Rscript new_scripts/08_dea_dotplots_networks.R [options]
#   --contrast NAME / --celltype NAME   restrict to one contrast / cell type
#   --collection cc                     fgsea collection to plot (default h)
#   --top N                             pathways per dotplot (default 20)
#   --top-network K                     pathways per network plot (default 5)
#   --max-genes-per-pathway M           genes per pathway in networks (default 25)
#   --ora-p-cutoff p                    P.Value cutoff for ORA genes (default 0.01)
#   --do-min-size / --do-max-size       DO gene-set size bounds (default config.R)
#   --skip-hallmark / --skip-do-gsea / --skip-do-ora
#   --no-pdf / --no-png                 skip the PDF report / the PNGs
#   --force                             overwrite outputs and ignore the caches
#   --list                              print what is available and exit
#
# Input:  pipeline/DE/*.tsv, pipeline/GSEA/*.tsv   (from script 05)
# Output: pipeline/dea_plots/  (PNGs, DEA_plots_report.pdf, cache/)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(ggrepel)
  library(enrichplot)
  library(DOSE)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(rlang)  # provides %||%
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# --- Use a non-forking BiocParallel backend -----------------------------------
# fgsea otherwise forks one worker per core, which can exhaust memory.
suppressPackageStartupMessages(library(BiocParallel))
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

# --- 1. CLI -------------------------------------------------------------------
# Options typed after the script name.
args <- commandArgs(trailingOnly = TRUE)
opt_list_only        <- "--list"          %in% args
opt_force            <- "--force"         %in% args
opt_no_pdf           <- "--no-pdf"        %in% args
opt_no_png           <- "--no-png"        %in% args
opt_skip_hallmark    <- "--skip-hallmark" %in% args
opt_skip_do_gsea     <- "--skip-do-gsea"  %in% args
opt_skip_do_ora      <- "--skip-do-ora"   %in% args

# Value that follows an option, e.g. --top 20, or the default if absent.
get_opt <- function(flag, default = NULL) {
  k <- which(args == flag)
  if (length(k) == 1L && length(args) > k) args[k + 1L] else default
}
opt_contrast    <- get_opt("--contrast")
opt_celltype    <- get_opt("--celltype")
opt_collection  <- get_opt("--collection", "h")
opt_top         <- as.integer(get_opt("--top",                   "20"))
opt_top_network <- as.integer(get_opt("--top-network",            "5"))
opt_max_genes   <- as.integer(get_opt("--max-genes-per-pathway", "25"))
opt_ora_p_cut   <- as.numeric(get_opt("--ora-p-cutoff",        "0.01"))
# DO gene-set size bounds, shared with scripts 10 and 11 via config.R.
opt_do_min_size <- as.numeric(get_opt("--do-min-size", do_gsea_min_size))
opt_do_max_size <- as.numeric(get_opt("--do-max-size", do_gsea_max_size))

# Cell types: the same set the consensus figures (scripts 10 and 11) use.
dea_cell_types <- consensus_cell_types
if (!is.null(opt_celltype)) {
  if (!opt_celltype %in% dea_cell_types) {
    stop("--celltype ", opt_celltype, " not in supported list: ",
         paste(dea_cell_types, collapse = ", "))
  }
  dea_cell_types <- opt_celltype
}

# --- 2. Output paths ----------------------------------------------------------
dea_dir   <- file.path(pipeline_root,
                       if (isTRUE(use_cellsweep)) "dea_plots_cellsweep"
                                                 else "dea_plots")
hd_dir    <- file.path(dea_dir, "hallmark_dotplots")
dg_dir    <- file.path(dea_dir, "do_gsea_dotplots")
do_dir    <- file.path(dea_dir, "do_ora_dotplots")
on_dir    <- file.path(dea_dir, "do_ora_networks")
cache_dir <- file.path(dea_dir, "cache")
for (d in c(hd_dir, dg_dir, do_dir, on_dir, cache_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
pdf_path <- file.path(dea_dir, "DEA_plots_report.pdf")

# --- 3. Helpers ---------------------------------------------------------------

# Functions from utils.R used below: celltype_label() turns a subset key into a
# plot label, display_cohorts() splits a contrast name into its two cohorts,
# contrast_pretty() gives a readable contrast name and read_gsea() reads an
# fgsea results table.

# Caption stating which cohort a positive NES favours.
direction_caption <- function(contrast_full, with_nes = TRUE) {
  cohorts <- display_cohorts(contrast_full)
  if (with_nes) {
    sprintf(
      "Positive NES = higher expression in %s    |    Negative NES = higher expression in %s",
      cohorts[2], cohorts[1]
    )
  } else {
    sprintf(
      "DE-gene selection: P.Value < %g (either direction) | logFC sign in network plot: red = higher in %s, blue = higher in %s",
      opt_ora_p_cut, cohorts[2], cohorts[1]
    )
  }
}

# HALLMARK_INTERFERON_GAMMA_RESPONSE -> INTERFERON GAMMA RESPONSE
prettify_pathway <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  gsub("_", " ", x)
}

# leadingEdge as a pipe-separated string or a list.
parse_leading_edge <- function(x) {
  if (is.list(x)) return(lapply(x, as.character))
  if (is.character(x)) {
    return(lapply(x, function(s) {
      if (is.na(s) || s == "") return(character(0))
      strsplit(s, "|", fixed = TRUE)[[1]]
    }))
  }
  rep(list(character(0)), length(x))
}

# HGNC symbols to Entrez IDs. Unmapped symbols are dropped.
hgnc_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) {
    return(data.frame(SYMBOL = character(0),
                      ENTREZID = character(0), stringsAsFactors = FALSE))
  }
  res <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = symbols, keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    )),
    error = function(e) {
      log_msg("   hgnc_to_entrez() failed: ", conditionMessage(e))
      data.frame(SYMBOL = character(0), ENTREZID = character(0),
                 stringsAsFactors = FALSE)
    }
  )
  res <- res[!is.na(res$ENTREZID), , drop = FALSE]
  res$ENTREZID <- as.character(res$ENTREZID)
  res
}


# Read a DE table, or NULL if missing. Drops any ENSG prefix from gene names.
read_de <- function(cell_type, contrast_full) {
  f <- file.path(de_dir, paste0(cell_type, "__", contrast_full, ".tsv"))
  if (!file.exists(f)) return(NULL)
  df <- as.data.frame(data.table::fread(f))
  if (!"gene" %in% colnames(df) || !"logFC" %in% colnames(df)) return(NULL)
  df$gene <- sub(".*_", "", df$gene)
  df
}

# TRUE if the DE table has no per-gene p-values (descriptive contrast).
is_descriptive_de <- function(de_df) {
  !("P.Value" %in% colnames(de_df))
}

# --- 4. Hallmark GSEA dotplot -------------------------------------------------
# Top N pathways by padj. x = NES, size = set size, colour = padj.
make_hallmark_dotplot <- function(gsea_df, contrast_full, cell_type,
                                  top_n = 20) {
  d <- gsea_df
  d <- d[!is.na(d$padj), , drop = FALSE]
  d <- d[order(d$padj), , drop = FALSE]
  d <- head(d, top_n)
  if (!nrow(d)) return(NULL)

  d$pathway_pretty <- prettify_pathway(d$pathway)
  d$pathway_pretty <- factor(d$pathway_pretty,
                             levels = d$pathway_pretty[order(d$NES)])

  cohorts <- display_cohorts(contrast_full)
  title <- sprintf("Hallmark GSEA — %s vs %s %s",
                   cohorts[1], cohorts[2], celltype_label(cell_type))

  # Symmetric x-axis, at least +/-3, wider if any |NES| is larger.
  x_lim <- max(3, max(abs(d$NES), na.rm = TRUE) + 0.5)

  ggplot(d, aes(x = NES, y = pathway_pretty)) +
    geom_point(aes(size = size, colour = padj)) +
    scale_colour_gradient(low = "red", high = "blue",
                          name = "p.adjust", trans = "log10",
                          guide = guide_colourbar(reverse = TRUE)) +
    scale_size_continuous(range = c(3, 10), name = "Set size") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    labs(title    = title,
         subtitle = direction_caption(contrast_full, with_nes = TRUE),
         x        = "NES",
         y        = NULL) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10,
                                       colour = "grey25"),
          axis.text.y   = element_text(size = 10),
          axis.text.x   = element_text(size = 10),
          axis.title.x  = element_text(size = 11),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 9))
}

# --- 5. Ranked Entrez statistic for DO GSEA -----------------------------------
# Genes ranked by t-statistic (logFC for descriptive contrasts), named by
# Entrez ID as DOSE needs.
ranked_entrez_stats <- function(de_df) {
  stat_col <- if ("t" %in% colnames(de_df)) "t" else "logFC"
  ok <- !is.na(de_df[[stat_col]]) & !is.na(de_df$gene) & nzchar(de_df$gene)
  d <- de_df[ok, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  map <- hgnc_to_entrez(d$gene)
  if (!nrow(map)) return(NULL)
  d <- merge(d[, c("gene", stat_col), drop = FALSE],
             map, by.x = "gene", by.y = "SYMBOL")
  if (!nrow(d)) return(NULL)

  # Keep the most extreme statistic per Entrez ID.
  d <- d[order(-abs(d[[stat_col]])), ]
  d <- d[!duplicated(d$ENTREZID), , drop = FALSE]

  stats <- setNames(d[[stat_col]], d$ENTREZID)
  stats <- stats[!is.na(stats)]
  if (!length(stats)) return(NULL)
  sort(stats, decreasing = TRUE)
}

# --- 6. DO GSEA ---------------------------------------------------------------
# Full, seeded gseDO result (pvalueCutoff = 1), cached so scripts 10 and 11 use
# the same numbers. Plots in this script show only the significant terms.
do_gsea_full_cached <- function(cell_type, contrast_full, force = FALSE,
                                min_size = opt_do_min_size,
                                max_size = opt_do_max_size) {
  cache_file <- file.path(cache_dir,
    paste0(cell_type, "__", contrast_full, "__do_gsea_full.rds"))
  de_file <- file.path(de_dir, paste0(cell_type, "__", contrast_full, ".tsv"))
  if (!file.exists(de_file)) return(NULL)
  if (!force && file.exists(cache_file) &&
      file.info(cache_file)$mtime > file.info(de_file)$mtime) {
    cached <- readRDS(cache_file)
    # The cache is also invalid if the DO size bounds changed.
    bounds <- attr(cached, "do_gsea_size_bounds")
    if (!is.null(cached) && identical(bounds, c(min_size, max_size)))
      return(cached)
    if (!is.null(cached))
      log_msg("   [DO GSEA] cache size bounds ",
              if (is.null(bounds)) "(unstamped)" else paste(bounds, collapse = "-"),
              " != requested ", min_size, "-", max_size, " - recomputing ",
              cell_type, " :: ", contrast_full)
  }
  de_df <- read_de(cell_type, contrast_full)
  if (is.null(de_df)) return(NULL)
  stats <- ranked_entrez_stats(de_df)
  if (is.null(stats) || length(stats) < min_size) return(NULL)
  set.seed(42)   # make the (otherwise stochastic) fgsea result reproducible
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      # Every DO term is kept (pvalueCutoff = 1); BH-adjusted p-values.
      DOSE::gseDO(geneList = stats,
                  organism = "hsa",
                  minGSSize = min_size, maxGSSize = max_size,
                  pvalueCutoff = 1, pAdjustMethod = "BH",
                  eps = 0, verbose = FALSE, by = "fgsea")
    )),
    error = function(e) {
      log_msg("   [DO GSEA] failed for ", cell_type, " :: ",
              contrast_full, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0L) {
    if (file.exists(cache_file)) unlink(cache_file)   # never persist NULL
    return(NULL)
  }
  res_readable <- tryCatch(
    DOSE::setReadable(res, "org.Hs.eg.db", "ENTREZID"),
    error = function(e) res
  )
  # Record the size bounds on the cached object.
  attr(res_readable, "do_gsea_size_bounds") <- c(min_size, max_size)
  saveRDS(res_readable, cache_file)
  res_readable
}

# DO terms with adjusted p <= 0.2 from the full result, for this script's plots.
do_gsea_for_subset_contrast <- function(cell_type, contrast_full,
                                        force = FALSE,
                                        min_size = opt_do_min_size,
                                        max_size = opt_do_max_size,
                                        pcut = 0.2) {
  full <- do_gsea_full_cached(cell_type, contrast_full, force = force,
                              min_size = min_size, max_size = max_size)
  if (is.null(full)) return(NULL)
  keep <- !is.na(full@result$p.adjust) & full@result$p.adjust <= pcut
  full@result <- full@result[keep, , drop = FALSE]
  if (!nrow(full@result)) return(NULL)
  full
}

# DO GSEA dotplot.
make_do_gsea_dotplot <- function(gse_res, contrast_full, cell_type,
                                 top_n = 20) {
  if (is.null(gse_res) || nrow(as.data.frame(gse_res)) == 0L) return(NULL)
  cohorts <- display_cohorts(contrast_full)
  title <- sprintf("Disease Ontology GSEA — %s vs %s %s",
                   cohorts[1], cohorts[2], celltype_label(cell_type))
  # Symmetric x-axis, as for the Hallmark dotplot.
  x_lim <- max(3, max(abs(as.data.frame(gse_res)$NES), na.rm = TRUE) + 0.5)
  p <- enrichplot::dotplot(gse_res, showCategory = top_n,
                           color = "p.adjust", x = "NES",
                           orderBy = "NES",
                           label_format = 40)
  p +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    # Most significant = red.
    scale_colour_gradient(low = "red", high = "blue", name = "p.adjust",
                          guide = guide_colourbar(reverse = TRUE)) +
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    labs(title    = title,
         subtitle = direction_caption(contrast_full, with_nes = TRUE)) +
    theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10,
                                       colour = "grey25"))
}

# --- 7. DO ORA ----------------------------------------------------------------
# enrichDO on genes with P.Value below the cutoff.
do_ora_for_subset_contrast <- function(cell_type, contrast_full,
                                       p_cutoff = opt_ora_p_cut,
                                       force = FALSE,
                                       min_size = opt_do_min_size,
                                       max_size = opt_do_max_size,
                                       qcut = 0.2) {
  cache_file <- file.path(cache_dir,
    paste0(cell_type, "__", contrast_full, "__do_ora.rds"))
  de_file <- file.path(de_dir, paste0(cell_type, "__", contrast_full, ".tsv"))
  if (!file.exists(de_file)) return(NULL)
  if (!force && file.exists(cache_file) &&
      file.info(cache_file)$mtime > file.info(de_file)$mtime) {
    return(readRDS(cache_file))
  }
  de_df <- read_de(cell_type, contrast_full)
  if (is.null(de_df)) return(NULL)
  if (is_descriptive_de(de_df)) {
    # Descriptive contrasts have no p-values, so no ORA.
    saveRDS(NULL, cache_file)
    return(NULL)
  }
  # Genes with raw P below the cutoff, in either direction.
  sig <- de_df[de_df$P.Value < p_cutoff & !is.na(de_df$P.Value), , drop = FALSE]
  if (!nrow(sig)) {
    saveRDS(NULL, cache_file); return(NULL)
  }
  map <- hgnc_to_entrez(sig$gene)
  if (!nrow(map)) {
    saveRDS(NULL, cache_file); return(NULL)
  }
  de_entrez <- unique(map$ENTREZID)
  # Universe = all genes tested in this DE table.
  universe_map <- hgnc_to_entrez(de_df$gene)
  universe <- unique(universe_map$ENTREZID)

  res <- tryCatch(
    suppressMessages(suppressWarnings(
      DOSE::enrichDO(gene = de_entrez,
                     ont = "DO",
                     pvalueCutoff = 1,         # keep everything; filter later
                     pAdjustMethod = "BH",
                     universe = universe,
                     minGSSize = min_size, maxGSSize = max_size,
                     qvalueCutoff = qcut,
                     readable = FALSE)
    )),
    error = function(e) {
      log_msg("   [DO ORA] failed for ", cell_type, " :: ",
              contrast_full, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0L) {
    saveRDS(NULL, cache_file); return(NULL)
  }
  res_readable <- tryCatch(
    DOSE::setReadable(res, "org.Hs.eg.db", "ENTREZID"),
    error = function(e) res
  )

  # Signed logFC per Entrez ID for the network plot node colours.
  lfc_map <- merge(
    de_df[, c("gene", "logFC")], map, by.x = "gene", by.y = "SYMBOL")
  lfc_map <- lfc_map[order(-abs(lfc_map$logFC)), ]
  lfc_map <- lfc_map[!duplicated(lfc_map$ENTREZID), , drop = FALSE]
  lfc_entrez <- setNames(lfc_map$logFC, lfc_map$ENTREZID)

  out <- list(enrich = res_readable, foldChange = lfc_entrez)
  saveRDS(out, cache_file)
  out
}

# DO ORA dotplot.
make_do_ora_dotplot <- function(ora_res, contrast_full, cell_type,
                                top_n = 20) {
  if (is.null(ora_res) ||
      is.null(ora_res$enrich) ||
      nrow(as.data.frame(ora_res$enrich)) == 0L) return(NULL)
  cohorts <- display_cohorts(contrast_full)
  title <- sprintf("Disease Ontology ORA — %s vs %s %s",
                   cohorts[1], cohorts[2], celltype_label(cell_type))
  p <- enrichplot::dotplot(ora_res$enrich, showCategory = top_n,
                           color = "p.adjust",
                           label_format = 40)
  p +
    # Most significant = red.
    scale_colour_gradient(low = "red", high = "blue", name = "p.adjust",
                          guide = guide_colourbar(reverse = TRUE)) +
    labs(title    = title,
         subtitle = direction_caption(contrast_full, with_nes = FALSE)) +
    theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10,
                                       colour = "grey25"))
}

# DO ORA network: linear and circular cnetplot side by side.
make_do_ora_network <- function(ora_res, contrast_full, cell_type,
                                top_k = 5) {
  if (is.null(ora_res) ||
      is.null(ora_res$enrich) ||
      nrow(as.data.frame(ora_res$enrich)) == 0L) return(NULL)
  edox <- ora_res$enrich
  fc   <- ora_res$foldChange

  pA <- tryCatch(
    enrichplot::cnetplot(edox, showCategory = top_k, foldChange = fc),
    error = function(e) {
      log_msg("   [DO ORA cnetplot linear] failed: ", conditionMessage(e))
      NULL
    }
  )
  pB <- tryCatch(
    enrichplot::cnetplot(edox, showCategory = top_k, foldChange = fc,
                         circular = TRUE, colorEdge = TRUE),
    error = function(e) {
      log_msg("   [DO ORA cnetplot circular] failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(pA) && is.null(pB)) return(NULL)

  cohorts    <- display_cohorts(contrast_full)
  page_title <- sprintf("Disease Ontology ORA — Gene-pathway network — %s vs %s %s",
                        cohorts[1], cohorts[2], celltype_label(cell_type))

  panels <- if (!is.null(pA) && !is.null(pB)) {
    cowplot::plot_grid(pA, pB, ncol = 2,
                       labels = c("A", "B"), label_size = 16,
                       rel_widths = c(0.8, 1.2))
  } else if (!is.null(pA)) {
    cowplot::plot_grid(pA, ncol = 1, labels = "A", label_size = 16)
  } else {
    cowplot::plot_grid(pB, ncol = 1, labels = "B", label_size = 16)
  }

  header <- cowplot::ggdraw() +
    cowplot::draw_label(page_title, fontface = "bold", size = 14,
                        hjust = 0.5, y = 0.7) +
    cowplot::draw_label(direction_caption(contrast_full, with_nes = FALSE),
                        size = 10, colour = "grey25",
                        hjust = 0.5, y = 0.25)

  cowplot::plot_grid(header, panels, ncol = 1,
                     rel_heights = c(0.08, 0.92))
}

# --- 8. Discovery -------------------------------------------------------------
# Subset/contrast pairs with an fgsea results table for the chosen collection.
discover_contrasts <- function() {
  pat <- paste0("__gsea_", opt_collection, "\\.tsv$")
  files <- list.files(gsea_dir, pattern = pat, full.names = FALSE)
  if (!length(files)) return(data.frame())
  stems <- sub(pat, "", files)
  m <- regexpr("__", stems, fixed = TRUE)
  ok <- m > 0L
  data.frame(
    subset        = substring(stems[ok], 1L, m[ok] - 1L),
    contrast_full = substring(stems[ok], m[ok] + 2L),
    stringsAsFactors = FALSE
  )
}

avail <- discover_contrasts()
avail <- avail[avail$subset %in% dea_cell_types, , drop = FALSE]

if (opt_list_only) {
  log_msg("Discovered cell_type x contrast pairs (Hallmark collection = ",
          opt_collection, "):")
  if (!nrow(avail)) {
    log_msg("  (none)")
  } else {
    cov <- avail %>%
      dplyr::mutate(present = TRUE) %>%
      tidyr::pivot_wider(id_cols = contrast_full, names_from = subset,
                         values_from = present, values_fill = FALSE) %>%
      as.data.frame()
    print(cov, row.names = FALSE)
  }
  quit(save = "no", status = 0L)
}

if (!nrow(avail)) {
  stop("No Hallmark GSEA TSVs found for cell types {",
       paste(dea_cell_types, collapse = ", "),
       "} under ", gsea_dir,
       ". Run script 05 for these subsets first.")
}

# Contrasts to plot: all found, or the one given with --contrast.
contrasts_to_run <- sort(unique(avail$contrast_full))
if (!is.null(opt_contrast)) {
  if (!opt_contrast %in% contrasts_to_run) {
    stop("--contrast ", opt_contrast, " not found among discovered ",
         "contrasts (", paste(contrasts_to_run, collapse = ", "), ")")
  }
  contrasts_to_run <- opt_contrast
}

# --- 9. Drive -----------------------------------------------------------------
log_msg("== 08_dea_dotplots_networks.R: ",
        length(contrasts_to_run), " contrast(s) x ",
        length(dea_cell_types), " cell type(s)")
log_msg("   Hallmark collection: ", opt_collection,
        "   top (dotplot): ", opt_top,
        "   top (network):", opt_top_network,
        "   ORA P<",       opt_ora_p_cut)

# Plots for the PDF, in contrast / cell type / plot order.
plot_records <- list()

# Keep a plot for the PDF report.
push_record <- function(kind, contrast, cell_type, plot) {
  plot_records[[length(plot_records) + 1L]] <<- list(
    kind = kind, contrast = contrast, cell_type = cell_type, plot = plot)
}

# Save a PNG unless --no-png, or unless it exists and --force is not set.
png_save <- function(file, plot, w, h) {
  if (opt_no_png) return(invisible(NULL))
  if (!opt_force && file.exists(file)) return(invisible(NULL))
  ggplot2::ggsave(file, plot, width = w, height = h, dpi = 200)
  invisible(file)
}

# For every contrast and cell type, make the four plots.
for (cf in contrasts_to_run) {
  log_msg("== Contrast: ", cf, " (", contrast_pretty(cf), ")")

  for (ct in dea_cell_types) {

    # --- Hallmark GSEA dotplot ------------------------------------------------
    if (!opt_skip_hallmark) {
      gsea_df <- read_gsea(ct, cf, opt_collection)
      if (is.null(gsea_df) || !nrow(gsea_df)) {
        log_msg("   [", ct, "] [Hallmark] no GSEA TSV - skipping")
      } else {
        p <- tryCatch(make_hallmark_dotplot(gsea_df, cf, ct, opt_top),
                      error = function(e) {
                        log_msg("   [", ct, "] [Hallmark] ERROR: ",
                                conditionMessage(e)); NULL
                      })
        if (!is.null(p)) {
          out <- file.path(hd_dir, paste0(cf, "__", ct, "_hallmark_dotplot.png"))
          saved <- png_save(out, p, 11, 8)
          if (!is.null(saved)) log_msg("   [", ct, "] [Hallmark] wrote ",
                                       basename(out))
          push_record("hallmark_dotplot", cf, ct, p)
        }
      }
    }

    # --- DO GSEA dotplot ------------------------------------------------------
    if (!opt_skip_do_gsea) {
      gse_do <- tryCatch(
        do_gsea_for_subset_contrast(ct, cf, force = opt_force),
        error = function(e) {
          log_msg("   [", ct, "] [DO GSEA] ERROR: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(gse_do)) {
        log_msg("   [", ct, "] [DO GSEA] no result - skipping")
      } else {
        p <- tryCatch(make_do_gsea_dotplot(gse_do, cf, ct, opt_top),
                      error = function(e) {
                        log_msg("   [", ct, "] [DO GSEA dotplot] ERROR: ",
                                conditionMessage(e)); NULL
                      })
        if (!is.null(p)) {
          out <- file.path(dg_dir, paste0(cf, "__", ct, "_do_gsea_dotplot.png"))
          saved <- png_save(out, p, 11, 8)
          if (!is.null(saved)) log_msg("   [", ct, "] [DO GSEA] wrote ",
                                       basename(out))
          push_record("do_gsea_dotplot", cf, ct, p)
        }
      }
    }

    # --- DO ORA: dotplot + network --------------------------------------------
    if (!opt_skip_do_ora) {
      ora_res <- tryCatch(
        do_ora_for_subset_contrast(ct, cf, p_cutoff = opt_ora_p_cut,
                                   force = opt_force),
        error = function(e) {
          log_msg("   [", ct, "] [DO ORA] ERROR: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(ora_res)) {
        log_msg("   [", ct, "] [DO ORA] no result (descriptive contrast?) ",
                "- skipping ORA dotplot + network")
      } else {
        # ORA dotplot
        p_d <- tryCatch(make_do_ora_dotplot(ora_res, cf, ct, opt_top),
                        error = function(e) {
                          log_msg("   [", ct, "] [DO ORA dotplot] ERROR: ",
                                  conditionMessage(e)); NULL
                        })
        if (!is.null(p_d)) {
          out <- file.path(do_dir, paste0(cf, "__", ct, "_do_ora_dotplot.png"))
          saved <- png_save(out, p_d, 11, 8)
          if (!is.null(saved)) log_msg("   [", ct, "] [DO ORA dot] wrote ",
                                       basename(out))
          push_record("do_ora_dotplot", cf, ct, p_d)
        }

        # ORA network (top 5)
        p_n <- tryCatch(make_do_ora_network(ora_res, cf, ct,
                                            top_k = opt_top_network),
                        error = function(e) {
                          log_msg("   [", ct, "] [DO ORA network] ERROR: ",
                                  conditionMessage(e)); NULL
                        })
        if (!is.null(p_n)) {
          out <- file.path(on_dir, paste0(cf, "__", ct, "_do_ora_network.png"))
          saved <- png_save(out, p_n, 22, 10)
          if (!is.null(saved)) log_msg("   [", ct, "] [DO ORA net] wrote ",
                                       basename(out))
          push_record("do_ora_network", cf, ct, p_n)
        }
      }
    }
  }
}

# --- 10. Multi-page PDF report ------------------------------------------------
if (!opt_no_pdf && length(plot_records)) {
  log_msg("== Assembling multi-page PDF: ", pdf_path)
  pdf(pdf_path, width = 12, height = 12, onefile = TRUE)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  # Page number in the bottom-right corner.
  .pg_n <- 0L
  page_number <- function() {
    .pg_n <<- .pg_n + 1L
    grid::upViewport(0)
    grid::grid.text(sprintf("Page %d", .pg_n), x = 0.99, y = 0.01,
                    just = c("right", "bottom"),
                    gp = grid::gpar(fontsize = 8, col = "grey50"))
  }

  # Cover page.
  cover <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Differential Expression — Enrichment Plots",
      fontface = "bold", size = 22, y = 0.78) +
    cowplot::draw_label(
      paste0(
        "Cell types: ", paste(dea_cell_types, collapse = ", "), "\n",
        "Contrasts: ", length(contrasts_to_run),
        "    Hallmark collection: ", opt_collection,
        "    ORA P-cutoff: ", opt_ora_p_cut, "\n",
        "Dotplot top-N: ", opt_top,
        "    Network top-K: ", opt_top_network, "\n\n",
        "Per (cell type x contrast):\n",
        "  1. Hallmark GSEA dotplot   (fgsea, from script 05)\n",
        "  2. Disease Ontology GSEA dotplot  (DOSE::gseDO)\n",
        "  3. Disease Ontology ORA dotplot   (DOSE::enrichDO)\n",
        "  4. Disease Ontology ORA network   (enrichplot::cnetplot,\n",
        "     linear + circular, top-", opt_top_network, " pathways)"
      ),
      size = 11, y = 0.40, lineheight = 1.25) +
    cowplot::draw_label(
      paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      size = 10, y = 0.06, colour = "grey40")
  print(cover)
  page_number()

  kind_label <- c(
    hallmark_dotplot = "Hallmark GSEA dotplot",
    do_gsea_dotplot  = "Disease Ontology GSEA dotplot",
    do_ora_dotplot   = "Disease Ontology ORA dotplot",
    do_ora_network   = "Disease Ontology ORA network plot"
  )
  # One page per plot, with a header naming contrast, cell type and plot type.
  for (rec in plot_records) {
    page_title <- sprintf("%s — %s — %s",
                          contrast_pretty(rec$contrast),
                          celltype_label(rec$cell_type),
                          kind_label[[rec$kind]] %||% rec$kind)
    header <- cowplot::ggdraw() +
      cowplot::draw_label(page_title, fontface = "bold",
                          size = 14, hjust = 0.5)
    page <- cowplot::plot_grid(header, rec$plot, ncol = 1,
                               rel_heights = c(0.04, 0.96))
    print(page)
    page_number()
  }

  grDevices::dev.off()
  log_msg("== PDF report done: ", pdf_path)
} else if (!opt_no_pdf) {
  log_msg("== no plots produced - skipping PDF assembly")
}

log_msg("== 08_dea_dotplots_networks.R complete. (",
        length(plot_records), " plot record(s))")
