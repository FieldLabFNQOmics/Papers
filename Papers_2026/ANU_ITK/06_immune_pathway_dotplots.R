# ------------------------------------------------------------------------------
# 06_immune_pathway_dotplots.R
# Spot matrices (dot plots) of pathway and disease enrichment across the CD4
# and CD8 T-cell subsets, one set per contrast:
#   DisGeNET      disease-term enrichment (enrichDGN) of DE genes, split into
#                 genes up and genes down
#   Hallmark      seven immune Hallmark pathways from the script 05 GSEA, drawn
#                 two ways: split by direction, and as a single grid
# The single-grid Hallmark plots are also combined into one PDF.
#
# Run:  Rscript new_scripts/06_immune_pathway_dotplots.R [options]
#   --contrast NAME    one contrast only
#   --hallmark-only    skip the DisGeNET plots
#   --dgn-only         skip the Hallmark plots
#   --force            ignore the enrichDGN cache
#   --list             print the available DE and GSEA results and exit
#   --no-pdf           skip the combined PDF
#
# Input:  pipeline/DE/*.tsv, pipeline/GSEA/*__gsea_h.tsv   (from script 05)
# Output: pipeline/GSEA/spot_matrices/
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(reshape2)
  library(stringr)
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# --- 1. CLI -------------------------------------------------------------------
# Options typed after the script name.
args <- commandArgs(trailingOnly = TRUE)
opt_list_only     <- "--list"          %in% args
opt_skip_dgn      <- "--hallmark-only" %in% args
opt_skip_hallmark <- "--dgn-only"      %in% args
opt_force         <- "--force"         %in% args
opt_no_pdf        <- "--no-pdf"        %in% args
opt_contrast <- {
  k <- which(args == "--contrast")
  if (length(k) == 1L && length(args) > k) args[k + 1L] else NULL
}

# --- 2. Output paths ----------------------------------------------------------
# Outputs, and a cache of enrichDGN results so re-runs are fast.
sm_dir    <- file.path(gsea_dir, "spot_matrices")
cache_dir <- file.path(sm_dir, "cache")
dir.create(sm_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Helpers ---------------------------------------------------------------

# display_cohorts() (utils.R) splits a contrast name such as "HBD_vs_T504S"
# into its two cohorts, formatted for plot labels.

# List result files named {subset}__{contrast}{suffix} in a folder and return
# the subset, contrast and file path for each.
discover_stems <- function(dir_path, suffix) {
  pat <- paste0(suffix, "$")
  files <- list.files(dir_path, pattern = pat, full.names = FALSE)
  if (length(files) == 0L) return(data.frame())
  stems <- sub(pat, "", files)
  m <- regexpr("__", stems, fixed = TRUE)
  ok <- m > 0L
  data.frame(
    subset        = substring(stems[ok], 1L, m[ok] - 1L),
    contrast_full = substring(stems[ok], m[ok] + 2L),
    file          = file.path(dir_path, files[ok]),
    stringsAsFactors = FALSE
  )
}

# Map HGNC symbols to Entrez IDs via org.Hs.eg.db (fast, no network).
hgnc_to_entrez <- function(hgnc_symbols) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db",  quietly = TRUE)) {
    stop("HGNC->Entrez mapping needs AnnotationDbi + org.Hs.eg.db. ",
         "Install them or run with --hallmark-only.")
  }
  hgnc_symbols <- unique(hgnc_symbols[!is.na(hgnc_symbols) & nzchar(hgnc_symbols)])
  if (length(hgnc_symbols) == 0L) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0),
                      stringsAsFactors = FALSE))
  }
  res <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = hgnc_symbols, keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    )),
    error = function(e) {
      log_msg("   [DGN] HGNC->Entrez lookup failed: ", conditionMessage(e))
      data.frame(SYMBOL = character(0), ENTREZID = character(0),
                 stringsAsFactors = FALSE)
    }
  )
  res <- res[!is.na(res$ENTREZID), , drop = FALSE]
  res$ENTREZID <- as.character(res$ENTREZID)
  res
}

# DisGeNET enrichment for one subset and contrast. Genes with P < 0.01 and
# |B| > 2 in the DE table are split into up (logFC > 0) and down, and each set
# is tested with DOSE::enrichDGN. Returns the top terms for each direction.
# The result is cached and reused unless the DE table is newer or --force.
dgn_for_subset_contrast <- function(subset_name, contrast_full,
                                    p_cutoff = 0.01, b_cutoff = 2,
                                    top_n = 20, force = FALSE) {
  cache_file <- file.path(cache_dir,
    paste0(subset_name, "__", contrast_full, "__dgn.rds"))
  de_file <- file.path(de_dir,
    paste0(subset_name, "__", contrast_full, ".tsv"))
  if (!file.exists(de_file)) return(NULL)
  if (!force && file.exists(cache_file) &&
      file.info(cache_file)$mtime > file.info(de_file)$mtime) {
    return(readRDS(cache_file))
  }
  if (!requireNamespace("DOSE", quietly = TRUE)) {
    stop("DGN spot matrix requires DOSE. Install it or run --hallmark-only.")
  }

  tt <- as.data.frame(data.table::fread(de_file))
  required <- c("gene", "logFC", "P.Value", "B")
  if (length(setdiff(required, colnames(tt)))) {
    # Descriptive contrast: no p-values.
    return(NULL)
  }
  # Drop any ENSG prefix from gene names.
  tt$gene <- sub(".*_", "", tt$gene)

  # Genes up and down in the second cohort.
  pos_idx <- tt$P.Value < p_cutoff & tt$logFC > 0
  neg_idx <- tt$P.Value < p_cutoff & tt$logFC < 0

  # enrichDGN on one gene set. Adds -log10(adjusted p) and gene ratio (%).
  run_one_side <- function(genes, b_vals) {
    if (length(genes) == 0L) return(NULL)
    map <- hgnc_to_entrez(genes)
    if (nrow(map) == 0L) return(NULL)
    df <- data.frame(SYMBOL = genes, B = b_vals, stringsAsFactors = FALSE)
    df <- merge(df, map, by = "SYMBOL")
    df <- df[abs(df$B) > b_cutoff, , drop = FALSE]
    if (nrow(df) == 0L) return(NULL)
    de_entrez <- unique(df$ENTREZID)
    edo <- tryCatch(
      suppressMessages(DOSE::enrichDGN(de_entrez)),
      error = function(e) {
        log_msg("   [DGN] enrichDGN failed for ", subset_name, " :: ",
                contrast_full, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(edo) || nrow(edo@result) == 0L) return(NULL)
    res <- edo@result[, c("Description", "p.adjust", "GeneRatio")]
    res <- res[seq_len(min(top_n, nrow(res))), , drop = FALSE]
    res$cell_type <- subset_name
    res[["-log10(adjusted_p_value)"]] <- -log10(res$p.adjust)
    parts <- strsplit(res$GeneRatio, "/", fixed = TRUE)
    res[["Gene Ratio"]] <- vapply(parts, function(p) {
      p <- suppressWarnings(as.numeric(p))
      if (length(p) >= 2L && !is.na(p[2]) && p[2] > 0) {
        100 * p[1] / p[2]
      } else 0
    }, numeric(1))
    res
  }

  out <- list(
    pos = run_one_side(tt$gene[pos_idx], tt$B[pos_idx]),
    neg = run_one_side(tt$gene[neg_idx], tt$B[neg_idx])
  )
  saveRDS(out, cache_file)
  out
}

# Long table for the DGN spot matrix. NULL if nothing is left.
build_dgn_long_df <- function(contrast_full, cell_types, force = FALSE) {
  per_ct <- setNames(vector("list", length(cell_types)), cell_types)
  for (ct in cell_types) {
    per_ct[[ct]] <- dgn_for_subset_contrast(ct, contrast_full, force = force)
  }

  # Disease terms to show: the curated list for this contrast in config.R, or
  # terms found in at least 4 cell types if there is no list.
  curated <- disgenet_curated_terms[[contrast_full]]
  if (is.null(curated)) {
    base <- sub("_batch2only$", "", contrast_full)
    curated <- disgenet_curated_terms[[base]]
  }
  if (is.null(curated)) {
    # Fallback: terms appearing in >=4 cell types (any direction).
    ct_per_desc <- list()
    for (ct in names(per_ct)) {
      x <- per_ct[[ct]]
      if (is.null(x)) next
      seen <- unique(c(if (!is.null(x$pos)) x$pos$Description,
                       if (!is.null(x$neg)) x$neg$Description))
      for (d in seen) ct_per_desc[[d]] <- unique(c(ct_per_desc[[d]], ct))
    }
    n_ct <- vapply(ct_per_desc, length, integer(1))
    curated <- names(n_ct)[n_ct >= 4L]
    if (length(curated) == 0L) {
      log_msg("   [DGN] no curated list and fallback (>=4 cell types) ",
              "yielded none for ", contrast_full)
      return(NULL)
    }
    log_msg("   [DGN] using auto-derived terms for ", contrast_full,
            ": ", paste(curated, collapse = ", "))
  }

  # Pad with zeros so every cell type and direction has a row.
  rows <- list()
  for (ct in cell_types) {
    x <- per_ct[[ct]]
    if (is.null(x)) next
    for (sign in c("pos", "neg")) {
      df <- x[[sign]]
      if (is.null(df)) next
      f <- df[df$Description %in% curated, , drop = FALSE]
      if (nrow(f) == 0L) next
      f$variable <- sign
      rows[[length(rows) + 1L]] <- f
    }
  }
  if (length(rows) == 0L) return(NULL)
  long <- do.call(rbind, rows)

  full_grid <- expand.grid(
    Description = curated,
    cell_type   = cell_types,
    variable    = c("pos", "neg"),
    stringsAsFactors = FALSE
  )
  long_key <- paste(long$Description, long$cell_type, long$variable, sep = "||")
  full_key <- paste(full_grid$Description, full_grid$cell_type,
                    full_grid$variable, sep = "||")
  miss <- setdiff(full_key, long_key)
  if (length(miss)) {
    add <- full_grid[match(miss, full_key), , drop = FALSE]
    add$p.adjust <- NA_real_
    add$GeneRatio <- NA_character_
    add[["-log10(adjusted_p_value)"]] <- 0
    add[["Gene Ratio"]] <- 0
    long <- rbind(long, add[, colnames(long), drop = FALSE])
  }

  long$Description <- factor(long$Description, levels = curated)
  long$cell_type   <- factor(long$cell_type,   levels = spot_matrix_y_order)
  long
}

# Hallmark fgsea results for one subset and contrast, 7 immune pathways only.
hallmark_for_subset_contrast <- function(subset_name, contrast_full) {
  f <- file.path(gsea_dir,
                 paste0(subset_name, "__", contrast_full, "__gsea_h.tsv"))
  if (!file.exists(f)) return(NULL)
  df <- as.data.frame(data.table::fread(f))
  cols <- colnames(df)
  for (need in c("pathway", "pval", "padj", "NES")) {
    j <- match(tolower(need), tolower(cols))
    if (!is.na(j)) cols[j] <- need
  }
  colnames(df) <- cols
  df <- df[df$pathway %in% hallmark_immune_pathways,
           c("pathway", "padj", "NES"), drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  data.frame(
    description = df$pathway,
    NES         = df$NES,
    FDR.q.val   = df$padj,
    cell_type   = subset_name,
    variable    = ifelse(df$NES >= 0, "pos", "neg"),
    stringsAsFactors = FALSE
  )
}

# Long table for the Hallmark spot matrices. Pathways with FDR > 0.25 are blanked.
build_hallmark_long_df <- function(contrast_full, cell_types) {
  # Read the Hallmark results for each cell type and pad missing combinations.
  rows <- lapply(cell_types, function(ct) {
    hallmark_for_subset_contrast(ct, contrast_full)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(NULL)
  long <- do.call(rbind, rows)

  full_grid <- expand.grid(
    description = hallmark_immune_pathways,
    cell_type   = cell_types,
    variable    = c("pos", "neg"),
    stringsAsFactors = FALSE
  )
  full_grid$NES <- 0
  full_grid$FDR.q.val <- NA_real_
  long_key <- paste(long$description, long$cell_type, long$variable, sep = "||")
  full_key <- paste(full_grid$description, full_grid$cell_type,
                    full_grid$variable, sep = "||")
  miss <- setdiff(full_key, long_key)
  if (length(miss)) {
    add <- full_grid[match(miss, full_key), , drop = FALSE]
    long <- rbind(long, add[, colnames(long), drop = FALSE])
  }

  # FDR > 0.25: set NES to 0 and FDR to NA.
  long$FDR.q.val <- ifelse(is.na(long$FDR.q.val) | long$FDR.q.val > 0.25,
                           NA_real_, long$FDR.q.val)
  long$NES <- ifelse(is.na(long$FDR.q.val), 0, long$NES)

  # Axis labels without the HALLMARK_ prefix.
  display_levels <- sub("^HALLMARK_", "", hallmark_immune_pathways)
  long$description <- sub("^HALLMARK_", "", long$description)
  long$description <- factor(long$description, levels = display_levels)
  long$cell_type   <- factor(long$cell_type,   levels = spot_matrix_y_order)
  long
}

# --- 4. Plot helpers ----------------------------------------------------------

# DGN spot matrix, faceted by direction. Positive logFC = higher in the
# second cohort of the contrast.
plot_dgn_spot_matrix <- function(long_df, contrast_full) {
  cohorts <- display_cohorts(contrast_full)
  pos_label <- paste("Higher in", cohorts[2])
  neg_label <- paste("Higher in", cohorts[1])

  # Dot size = -log10(adjusted p), fill = gene ratio (%).
  ggplot(long_df, aes(Description, cell_type)) +
    geom_point(aes(size = 0.81 * .data[["-log10(adjusted_p_value)"]],
                   fill = .data[["Gene Ratio"]]), shape = 21) +
    scale_fill_gradient(low = "blue", high = "red") +
    facet_grid(. ~ factor(variable, levels = c("pos", "neg")),
               labeller = as_labeller(c("pos" = pos_label,
                                        "neg" = neg_label)),
               scales = "free_y") +
    labs(size = "-log10_adj_p_value", x = NULL, y = NULL) +
    ggtitle(contrast_full) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5),
          axis.title    = element_text(size = 16),
          axis.text     = element_text(size = 14),
          axis.ticks.x  = element_blank(),
          axis.text.x   = element_text(angle = 45, hjust = 0),
          axis.ticks.y  = element_blank(),
          strip.text    = element_text(size = 16)) +
    scale_x_discrete(position = "top") +
    scale_size_continuous(range = c(1, 20))
}

# Hallmark spot matrix, faceted by direction. Size = |NES|, fill = BH q.
# Blanked (non-significant or untested) results are not drawn.
plot_hallmark_facet_spot_matrix <- function(long_df, contrast_full) {
  cohorts <- display_cohorts(contrast_full)
  pos_label <- paste("Higher in", cohorts[2])
  neg_label <- paste("Higher in", cohorts[1])

  draw <- long_df[!is.na(long_df$FDR.q.val), , drop = FALSE]
  n_drop <- nrow(long_df) - nrow(draw)
  if (n_drop) {
    log_msg("   [hallmark facet] ", contrast_full, ": dropped ", n_drop,
            " padded/blanked cell(s) (FDR NA) that would have rendered as ",
            "size-1 grey dots")
  }
  if (nrow(draw) == 0L) {
    log_msg("   [hallmark facet] no significant points (FDR<=0.25) for ",
            contrast_full)
    return(NULL)
  }

  ggplot(draw, aes(description, cell_type)) +
    geom_point(aes(size = abs(NES), fill = FDR.q.val), shape = 21) +
    scale_fill_gradient(low = "red", high = "blue",
                        name = "BH q",
                        guide = guide_colourbar(reverse = TRUE)) +
    # Keep both facets even if all results point one way.
    facet_grid(. ~ factor(variable, levels = c("pos", "neg")),
               labeller = as_labeller(c("pos" = pos_label,
                                        "neg" = neg_label)),
               scales = "free_y", drop = FALSE) +
    labs(size = "NES", x = NULL, y = NULL,
         caption = paste0("Fill = Benjamini-Hochberg adjusted p (fgsea), ",
                          "over the full Hallmark collection.",
                          "  Pathways with q > 0.25 are not plotted.")) +
    ggtitle(paste0("hallmark_", contrast_full)) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.caption  = element_text(hjust = 0.5, size = 9,
                                       colour = "grey30"),
          axis.title    = element_text(size = 16),
          axis.text     = element_text(size = 14),
          axis.ticks.x  = element_blank(),
          axis.text.x   = element_text(angle = 45, hjust = 0),
          axis.ticks.y  = element_blank(),
          strip.text    = element_text(size = 16),
          legend.title  = element_text(hjust = 0.5)) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_size_continuous(range = c(1, 20), breaks = c(0, 1, 2, 3))
}

# Hallmark spot matrix in a single grid. Size = log10(1/q), fill = NES.
# Blanked results are not drawn.
plot_hallmark_mod_spot_matrix <- function(long_df, contrast_full) {
  cohorts <- display_cohorts(contrast_full)
  draw <- long_df[!is.na(long_df$FDR.q.val), , drop = FALSE]
  if (nrow(draw) == 0L) {
    log_msg("   [hallmark mod] no significant points (FDR<=0.25) for ",
            contrast_full)
    return(NULL)
  }
  # padj can be exactly 0. Floor it to the smallest positive value in the plot.
  zero_fdr <- which(draw$FDR.q.val <= 0)
  if (length(zero_fdr)) {
    pos_fdr   <- draw$FDR.q.val[draw$FDR.q.val > 0]
    fdr_floor <- if (length(pos_fdr)) min(pos_fdr) else 1e-10
    for (i in zero_fdr) {
      log_msg("   [padj==0] ", contrast_full, " :: ",
              as.character(draw$cell_type[i]), " / ",
              as.character(draw$description[i]),
              " - FDR floored to ", format(fdr_floor, digits = 3),
              " for the size scale (log10(1/0) = Inf would be dropped)")
    }
    draw$FDR.q.val[zero_fdr] <- fdr_floor
  }

  # Colour scale symmetric around NES = 0.
  rng <- max(abs(draw$NES), na.rm = TRUE)
  if (!is.finite(rng) || rng == 0) rng <- 1

  ggplot(draw, aes(description, cell_type)) +
    geom_point(aes(size = log10(1 / FDR.q.val), fill = NES), shape = 21) +
    scale_fill_gradient(low = "blue", high = "red", limits = c(-rng, rng)) +
    # FDR.q.val is fgsea's BH-adjusted p-value.
    labs(size = "log10(1/BH q)",
         fill = "NES", x = NULL, y = NULL,
         caption = paste0("Blue = higher in ", cohorts[1],
                          "    Red = higher in ", cohorts[2],
                          "\nq = Benjamini-Hochberg adjusted p (fgsea), ",
                          "over the full Hallmark collection.",
                          "\nPathways with q > 0.25 are not plotted.")) +
    ggtitle(paste0("hallmark_", contrast_full)) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.caption  = element_text(hjust = 0.5, size = 11,
                                       colour = "grey30"),
          axis.title    = element_text(size = 16),
          axis.text     = element_text(size = 14),
          axis.ticks.x  = element_blank(),
          axis.text.x   = element_text(angle = 45, hjust = 0),
          axis.ticks.y  = element_blank(),
          strip.text    = element_text(size = 16),
          legend.title  = element_text(hjust = 0.5)) +
    scale_x_discrete(position = "top") +
    scale_size_continuous(range = c(1, 19),
                          breaks = c(1, 2, 4, 6, 8, 10))
}

# --- 5. Discover what's available ---------------------------------------------
# Find which DE and Hallmark GSEA results exist for the spot-matrix cell types.
avail_de   <- discover_stems(de_dir,   "\\.tsv")
avail_gsea <- discover_stems(gsea_dir, "__gsea_h\\.tsv")

avail_de_sm   <- avail_de[avail_de$subset     %in% spot_matrix_cell_types, ,
                          drop = FALSE]
avail_gsea_sm <- avail_gsea[avail_gsea$subset %in% spot_matrix_cell_types, ,
                            drop = FALSE]

# Print a contrast x cell type table of available results (--list).
coverage_table <- function(df, label) {
  if (nrow(df) == 0L) {
    log_msg("== ", label, ": NONE")
    return(invisible(NULL))
  }
  cov <- df %>%
    dplyr::mutate(present = TRUE) %>%
    tidyr::pivot_wider(id_cols = contrast_full, names_from = subset,
                       values_from = present, values_fill = FALSE)
  log_msg("== ", label, ":")
  print(as.data.frame(cov), row.names = FALSE)
}

if (opt_list_only) {
  coverage_table(avail_de_sm,   "DE coverage (rows=contrast, cols=cell_type)")
  coverage_table(avail_gsea_sm, "Hallmark GSEA coverage")
  quit(save = "no", status = 0L)
}

# Contrasts to plot: all found, or the one given with --contrast.
all_contrasts <- unique(c(avail_de_sm$contrast_full, avail_gsea_sm$contrast_full))
if (length(all_contrasts) == 0L) {
  stop("No DE or Hallmark GSEA TSVs found for any of the spot-matrix ",
       "cell types. Check that script 05 has been run for the new ",
       "CD4_CTL/CD4_*/CD8_* sub-subset tasks.")
}
if (!is.null(opt_contrast)) {
  if (!opt_contrast %in% all_contrasts) {
    stop("--contrast ", opt_contrast, " not found among discovered contrasts.")
  }
  contrasts_to_plot <- opt_contrast
} else {
  contrasts_to_plot <- all_contrasts
}

# --- 6. Drive -----------------------------------------------------------------
log_msg("== 06_immune_pathway_dotplots.R: ",
        length(contrasts_to_plot), " contrast(s) to process")

# For each contrast: build the data table, save it as TSV, draw and save the plots.
for (cf in contrasts_to_plot) {
  log_msg("== Contrast: ", cf)

  # --- DGN spot matrix --------------------------------------------------------
  if (!opt_skip_dgn) {
    dgn_long <- tryCatch(
      build_dgn_long_df(cf, spot_matrix_cell_types, force = opt_force),
      error = function(e) {
        log_msg("   [DGN] ERROR for ", cf, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(dgn_long) && nrow(dgn_long) > 0L) {
      tsv_file <- file.path(sm_dir, paste0(cf, "_dgn_long.tsv"))
      readr::write_tsv(dgn_long, tsv_file)

      p <- plot_dgn_spot_matrix(dgn_long, cf)
      width <- spot_matrix_dgn_widths[[cf]]
      if (is.null(width)) width <- 15
      out_file <- file.path(sm_dir, paste0(cf, "_spotMatrix_pos_neg.png"))
      ggsave(out_file, p, width = width, height = 8, dpi = 200)
      log_msg("   [DGN] wrote ", out_file)
    } else {
      log_msg("   [DGN] no data for ", cf, " - skipping")
    }
  }

  # --- Hallmark spot matrices -------------------------------------------------
  if (!opt_skip_hallmark) {
    hm_long <- tryCatch(
      build_hallmark_long_df(cf, spot_matrix_cell_types),
      error = function(e) {
        log_msg("   [hallmark] ERROR for ", cf, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(hm_long) && nrow(hm_long) > 0L) {
      tsv_file <- file.path(sm_dir, paste0(cf, "_hallmark_long.tsv"))
      readr::write_tsv(hm_long, tsv_file)

      # The facet plot is NULL when nothing passes FDR <= 0.25.
      p_facet <- plot_hallmark_facet_spot_matrix(hm_long, cf)
      if (!is.null(p_facet)) {
        out_file <- file.path(sm_dir,
          paste0(cf, "_hallmark_FDR0.25_reduced_spotMatrix_pos_neg.png"))
        ggsave(out_file, p_facet, width = 13, height = 8, dpi = 200)
        log_msg("   [hallmark facet] wrote ", out_file)
      }

      p_mod <- plot_hallmark_mod_spot_matrix(hm_long, cf)
      if (!is.null(p_mod)) {
        out_file <- file.path(sm_dir,
          paste0(cf, "_hallmark_FDR0.25_reduced_spotMatrix_mod_pos_neg.png"))
        ggsave(out_file, p_mod, width = 8, height = 8, dpi = 200)
        log_msg("   [hallmark mod] wrote ", out_file)
      }
    } else {
      log_msg("   [hallmark] no GSEA TSVs for ", cf,
              " - skipping hallmark plots")
    }
  }
}

# Combine the single-grid Hallmark plots into one PDF (skip with --no-pdf).
if (!opt_no_pdf) {
  # bundle_mod_spot_matrices_pdf() is in utils.R.
  bundle_mod_spot_matrices_pdf(sm_dir)
}

log_msg("== 06_immune_pathway_dotplots.R done.")
