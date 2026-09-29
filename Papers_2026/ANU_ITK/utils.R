# ------------------------------------------------------------------------------
# utils.R
# Shared helper functions. Sourced after config.R by every script.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
})

# --- Logging ------------------------------------------------------------------
# Print a message with a timestamp.
log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
}

# --- Directory setup ----------------------------------------------------------
# Creates any missing pipeline folders.
ensure_pipeline_dirs <- function() {
  dirs <- c(pipeline_root, sce_dir, emptydrops_dir, emptydrops_diag,
            doublets_dir, qc_dir, merged_dir, de_dir, gsea_dir, sanity_dir)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

# --- Capture lookup -----------------------------------------------------------
# Dataset ID, capture name and paths for capture index i (1..20).
get_capture_info <- function(i) {
  if (!isTRUE(i %in% capture_table$i)) {
    stop(
      "Capture index i = ", i, " is out of range. ",
      "Valid range: 1..", max(capture_table$i),
      " (one per 10x capture). ",
      "Did the PBS array driver pass an index intended for script 05? ",
      "Scripts 01 / 01b / 02 must be submitted with `-J 1-10` and ",
      "`-J 11-20` only - see README.md (run order).",
      call. = FALSE
    )
  }
  row <- capture_table[capture_table$i == i, , drop = FALSE]
  # Cell Ranger was run from a folder named after the capture, so outputs are at
  # <batch_root>/<capture>/<capture>/outs/.
  cellranger_run_dir <- file.path(row$batch_root, row$capture, row$capture)
  raw_matrix_dir <- file.path(cellranger_run_dir,
                              "outs", "multi", "count",
                              "raw_feature_bc_matrix")
  # BAM folder. Script 02 fills in the sample name from its BAM template.
  bam_dir <- file.path(cellranger_run_dir, "outs", "per_sample_outs")
  list(
    i                  = row$i,
    dataset            = row$dataset,
    capture            = row$capture,
    batch              = row$batch,
    batch_root         = row$batch_root,
    cellranger_run_dir = cellranger_run_dir,
    raw_matrix_dir     = raw_matrix_dir,
    bam_dir            = bam_dir
  )
}

# --- Sample sheet -------------------------------------------------------------
# Read and check the sample sheet. Each (capture, hto) pair must appear once.
load_sample_sheet <- function(path = sample_sheet_path) {
  ss <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  req <- c("capture", "hto", "patient", "cohort")
  missing_cols <- setdiff(req, colnames(ss))
  if (length(missing_cols)) {
    stop("Sample sheet missing columns: ", paste(missing_cols, collapse = ", "))
  }
  # Duplicate (capture, hto) would break the join.
  dup <- ss[duplicated(ss[, c("capture", "hto")]), ]
  if (nrow(dup)) {
    stop("Sample sheet has duplicate (capture, hto) rows:\n",
         paste(capture, dup$capture, dup$hto, sep = " / ", collapse = "\n"))
  }
  ss
}

# --- Barcode-keyed metadata merge ---------------------------------------------
# Join on barcode and check every barcode matches before writing columns.
merge_by_barcode <- function(dest_df, src_df, key = "barcode", cols) {
  stopifnot(key %in% colnames(dest_df), key %in% colnames(src_df))
  stopifnot(!anyDuplicated(src_df[[key]]))
  stopifnot(all(dest_df[[key]] %in% src_df[[key]]))
  src_idx <- match(dest_df[[key]], src_df[[key]])
  for (col in cols) {
    stopifnot(col %in% colnames(src_df))
    dest_df[[col]] <- src_df[[col]][src_idx]
  }
  dest_df
}

# --- Sanity input writer ------------------------------------------------------
# Writes counts as a text matrix for Sanity, in blocks of genes. Not used by the
# included scripts.
write_sanity_input <- function(sce, out_path, assay_name = "counts",
                               chunk_cells = 500L) {
  m <- SummarizedExperiment::assay(sce, assay_name)
  stopifnot(inherits(m, "sparseMatrix"))
  gene_ids <- rownames(m)
  cell_ids <- colnames(m)
  con <- file(out_path, open = "wt")
  on.exit(close(con), add = TRUE)
  # Header: "GeneID\tcell1\tcell2\t..."
  writeLines(paste(c("GeneID", cell_ids), collapse = "\t"), con)
  n <- ncol(m)
  for (start in seq(1L, n, by = chunk_cells)) {
    end <- min(start + chunk_cells - 1L, n)
    block <- as.matrix(m[, start:end, drop = FALSE])
    # Collect rows in blocks and write each block with write.table.
    if (start == 1L) {
      buf <- block
    } else {
      buf <- cbind(buf, block)
    }
  }
  utils::write.table(
    data.frame(GeneID = gene_ids, buf, check.names = FALSE),
    file = out_path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE
  )
  invisible(out_path)
}

# --- DE: pseudobulk limma-voom ------------------------------------------------
# Sums counts per patient x cell type and fits limma-voom with ~0 + cohort
# (+ batch if include_batch). contrast = c(level1, level2); positive logFC =
# higher in level2. Samples with fewer than min_cells cells are dropped.
# Returns a BH-adjusted topTable.
pseudobulk_de <- function(counts, col_meta, contrast,
                          include_batch = TRUE, min_cells = 10L) {
  stopifnot(ncol(counts) == nrow(col_meta))
  required <- c("patient", "cohort", "batch", "cell_type")
  missing_cols <- setdiff(required, colnames(col_meta))
  if (length(missing_cols)) {
    stop("col_meta missing columns: ", paste(missing_cols, collapse = ", "))
  }
  keep_cells <- col_meta$cohort %in% contrast
  if (!any(keep_cells)) {
    stop("No cells matching contrast levels: ", paste(contrast, collapse = " vs "))
  }
  counts <- counts[, keep_cells, drop = FALSE]
  col_meta <- col_meta[keep_cells, , drop = FALSE]

  # One pseudobulk sample per patient x cell type.
  sample_key <- paste(col_meta$patient, col_meta$cell_type, sep = "__")
  sample_tbl <- col_meta %>%
    dplyr::mutate(.sample_key = sample_key) %>%
    dplyr::group_by(.sample_key, patient, cohort, batch, cell_type) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_cells >= min_cells)

  if (nrow(sample_tbl) < 2) {
    stop("Fewer than 2 pseudobulk samples survive min_cells filter.")
  }
  # Sum counts with a sparse cells x samples indicator matrix.
  cell_to_sample <- match(sample_key, sample_tbl$.sample_key)
  keep_cell <- !is.na(cell_to_sample)
  counts <- counts[, keep_cell, drop = FALSE]
  cell_to_sample <- cell_to_sample[keep_cell]
  ind <- Matrix::sparseMatrix(
    i = seq_along(cell_to_sample),
    j = cell_to_sample,
    x = 1,
    dims = c(length(cell_to_sample), nrow(sample_tbl))
  )
  pb <- counts %*% ind  # genes x samples
  pb <- as.matrix(pb)
  colnames(pb) <- sample_tbl$.sample_key

  # Filter genes expressed in enough samples.
  keep_gene <- rowSums(pb >= 5) >= 2
  pb <- pb[keep_gene, , drop = FALSE]

  # Design matrix.
  sample_tbl$cohort <- factor(sample_tbl$cohort, levels = contrast)
  if (include_batch && length(unique(sample_tbl$batch)) > 1L) {
    design <- model.matrix(~0 + cohort + batch, data = sample_tbl)
  } else {
    design <- model.matrix(~0 + cohort, data = sample_tbl)
  }
  colnames(design) <- make.names(colnames(design))

  dge <- edgeR::DGEList(counts = pb)
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  ctr_name <- paste0("cohort", contrast[2], "-cohort", contrast[1])
  # Re-derive names after make.names.
  colA <- make.names(paste0("cohort", contrast[1]))
  colB <- make.names(paste0("cohort", contrast[2]))
  cm <- limma::makeContrasts(contrasts = paste0(colB, "-", colA),
                             levels = design)
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2, robust = TRUE)
  tt <- limma::topTable(fit2, number = Inf, adjust.method = "BH",
                        sort.by = "P")
  tt <- tibble::rownames_to_column(tt, var = "gene")
  attr(tt, "n_samples") <- nrow(sample_tbl)
  attr(tt, "sample_tbl") <- sample_tbl
  tt
}

# --- DE: single-cell limma with duplicateCorrelation --------------------------
# Used for GEM108 pre vs post. Blocks on capture.
single_cell_de_paired <- function(logcounts, col_meta, group_col, levels,
                                  block_col = "capture") {
  stopifnot(ncol(logcounts) == nrow(col_meta))
  stopifnot(group_col %in% colnames(col_meta))
  stopifnot(block_col %in% colnames(col_meta))

  keep <- col_meta[[group_col]] %in% levels
  logcounts <- logcounts[, keep, drop = FALSE]
  col_meta <- col_meta[keep, , drop = FALSE]

  grp <- factor(col_meta[[group_col]], levels = levels)
  block <- factor(col_meta[[block_col]])
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels

  dupcor <- limma::duplicateCorrelation(logcounts, design, block = block)
  fit <- limma::lmFit(logcounts, design,
                      block = block, correlation = dupcor$consensus.correlation)
  cm <- limma::makeContrasts(contrasts = paste0(levels[2], "-", levels[1]),
                             levels = design)
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2, robust = TRUE)
  tt <- limma::topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
  tt <- tibble::rownames_to_column(tt, var = "gene")
  attr(tt, "rho") <- dupcor$consensus.correlation
  tt
}

# --- DE: descriptive (no p-values) --------------------------------------------
# For GEM108 vs HBD, where one side is a single patient. Pseudobulk per patient;
# returns logFC and mean expression only.
descriptive_de <- function(counts, col_meta, patient_col = "patient",
                           group_col = "cohort_or_patient", levels) {
  stopifnot(group_col %in% colnames(col_meta))
  stopifnot(patient_col %in% colnames(col_meta))
  keep <- col_meta[[group_col]] %in% levels
  counts <- counts[, keep, drop = FALSE]
  col_meta <- col_meta[keep, , drop = FALSE]

  # Pseudobulk per patient.
  sample_key <- col_meta[[patient_col]]
  sample_tbl <- data.frame(
    patient = unique(sample_key),
    stringsAsFactors = FALSE
  )
  sample_tbl$group <- col_meta[[group_col]][match(sample_tbl$patient, sample_key)]
  cell_to_sample <- match(sample_key, sample_tbl$patient)
  ind <- Matrix::sparseMatrix(
    i = seq_along(cell_to_sample),
    j = cell_to_sample,
    x = 1,
    dims = c(length(cell_to_sample), nrow(sample_tbl))
  )
  pb <- as.matrix(counts %*% ind)
  colnames(pb) <- sample_tbl$patient

  dge <- edgeR::DGEList(counts = pb)
  dge <- edgeR::calcNormFactors(dge)
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

  grp <- sample_tbl$group
  mean_A <- rowMeans(logcpm[, grp == levels[1], drop = FALSE])
  mean_B <- rowMeans(logcpm[, grp == levels[2], drop = FALSE])
  tt <- data.frame(
    gene         = rownames(pb),
    logFC        = mean_B - mean_A,
    mean_logCPM  = (mean_A + mean_B) / 2,
    mean_A       = mean_A,
    mean_B       = mean_B,
    stringsAsFactors = FALSE
  )
  # Sort by |logFC|. No p-values.
  tt <- tt[order(-abs(tt$logFC)), ]
  rownames(tt) <- NULL
  attr(tt, "levels") <- levels
  attr(tt, "note")   <- "Descriptive only; no p-values (n=1 on one side)."
  tt
}

# --- GSEA with fgsea ----------------------------------------------------------
# Runs every collection in gsea_collections on a ranked vector (t or logFC).
# Writes per collection a results TSV, a summary table PNG and per-pathway
# enrichment PNGs, plus one HTML index per contrast.
run_gsea <- function(ranked_stats, out_prefix,
                     collections = gsea_collections,
                     gmt_dir = gsea_gmt_dir, version = gsea_version,
                     species = "Hs", min_size = 15, max_size = 500,
                     nperm = 10000,
                     top_n_plots = 20,
                     plot_width = 6, plot_height = 4, dpi = 150,
                     seed = gsea_seed) {
  suppressPackageStartupMessages({
    library(fgsea)
    library(ggplot2)
  })
  stopifnot(!is.null(names(ranked_stats)))
  ranked_stats <- sort(ranked_stats, decreasing = TRUE)
  results <- list()
  html_blocks <- list()

  # Make pathway names safe for filenames.
  safe_name <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    substr(x, 1, 120)
  }

  for (cc in collections) {
    gmt_file <- file.path(gmt_dir, sprintf("%s.all.%s.gmt", cc, version))
    if (!file.exists(gmt_file)) {
      log_msg("  GSEA skipping ", cc, " - no gmt file at ", gmt_file)
      next
    }
    pathways <- fgsea::gmtPathways(gmt_file)
    # eps = 0 allows p-values below 1e-50. Tie warnings are muffled.
    # Seeded per collection so any subset of collections gives the same result.
    set.seed(seed + sum(utf8ToInt(cc)))
    res <- withCallingHandlers(
      fgsea::fgsea(pathways = pathways, stats = ranked_stats,
                   minSize = min_size, maxSize = max_size,
                   nPermSimple = nperm, eps = 0),
      warning = function(w) {
        if (grepl("ties in the preranked stats",
                  conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    res <- res[order(res$padj), ]
    out_file <- paste0(out_prefix, "__gsea_", cc, ".tsv")
    data.table::fwrite(res, file = out_file, sep = "\t")
    results[[cc]] <- res

    # --- pick top N up + top N down by padj for plotting ----------------------
    res_up <- res[!is.na(res$NES) & res$NES > 0, ]
    res_dn <- res[!is.na(res$NES) & res$NES < 0, ]
    top_up <- head(res_up[order(res_up$padj), ], top_n_plots)
    top_dn <- head(res_dn[order(res_dn$padj), ], top_n_plots)
    top_combined <- rbind(top_up, top_dn)

    # --- summary gseaTable PNG ------------------------------------------------
    table_png <- paste0(out_prefix, "__gsea_", cc, "_table.png")
    table_ok <- FALSE
    if (nrow(top_combined) > 0) {
      tryCatch({
        png(table_png, width = 1400,
            height = max(300, 32 * nrow(top_combined)), res = 120)
        on.exit(try(dev.off(), silent = TRUE), add = TRUE)
        gt <- fgsea::plotGseaTable(pathways[top_combined$pathway],
                                   ranked_stats, top_combined,
                                   gseaParam = 0.5)
        # Newer fgsea versions return the table instead of drawing it.
        if (!is.null(gt)) {
          if (inherits(gt, c("gtable", "grob"))) {
            grid::grid.draw(gt)
          } else if (inherits(gt, c("gg", "ggplot"))) {
            print(gt)
          }
        }
        dev.off()
        table_ok <- TRUE
      }, error = function(e) {
        log_msg("    plotGseaTable failed for ", cc, ": ",
                conditionMessage(e))
        try(dev.off(), silent = TRUE)
      })
    }

    # --- per-pathway running-enrichment PNGs ----------------------------------
    plots_dir <- paste0(out_prefix, "__gsea_", cc, "_plots")
    if (nrow(top_combined) > 0) {
      dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
    }
    plot_rel_paths <- character(0)
    for (i in seq_len(nrow(top_combined))) {
      pw <- top_combined$pathway[i]
      direction <- if (top_combined$NES[i] > 0) "up" else "dn"
      png_name <- file.path(
        plots_dir,
        sprintf("%02d_%s_%s.png", i, direction, safe_name(pw))
      )
      tryCatch({
        p <- fgsea::plotEnrichment(pathways[[pw]], ranked_stats) +
          ggplot2::labs(
            title = pw,
            subtitle = sprintf("NES=%.2f  padj=%.2g  size=%d",
                               top_combined$NES[i],
                               top_combined$padj[i],
                               top_combined$size[i])
          )
        ggplot2::ggsave(png_name, p, width = plot_width,
                        height = plot_height, dpi = dpi)
        plot_rel_paths <- c(plot_rel_paths,
                            file.path(basename(plots_dir), basename(png_name)))
      }, error = function(e) {
        log_msg("    plotEnrichment failed for ", pw, ": ",
                conditionMessage(e))
      })
    }

    html_blocks[[cc]] <- list(
      collection   = cc,
      tsv          = basename(out_file),
      table_png    = if (table_ok) basename(table_png) else NA_character_,
      plot_paths   = plot_rel_paths,
      top_combined = top_combined
    )
    log_msg("    ", cc, ": ", nrow(res), " pathways tested, ",
            nrow(top_combined), " plotted (",
            nrow(top_up), " up / ", nrow(top_dn), " down)")
  }

  # --- HTML report ------------------------------------------------------------
  if (length(html_blocks)) {
    html_path <- paste0(out_prefix, "__gsea_report.html")
    write_gsea_html(html_path, basename(out_prefix), html_blocks)
    log_msg("    wrote ", html_path)
  }

  invisible(results)
}

# Write the per-contrast GSEA HTML index. Image links are relative.
write_gsea_html <- function(path, title, blocks) {
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(c(
    "<!doctype html>",
    "<html><head><meta charset='utf-8'>",
    paste0("<title>GSEA — ", esc(title), "</title>"),
    "<style>",
    "body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;",
    "  margin:2rem;color:#222;max-width:1400px;}",
    "h1{font-size:1.4rem;}",
    "h2{margin-top:2rem;border-bottom:1px solid #ddd;padding-bottom:0.25rem;}",
    "nav a{margin-right:1rem;}",
    ".plots{display:grid;grid-template-columns:repeat(auto-fill,minmax(420px,1fr));",
    "  gap:0.75rem;}",
    ".plots figure{margin:0;border:1px solid #eee;padding:0.25rem;}",
    ".plots img{width:100%;height:auto;display:block;}",
    ".plots figcaption{font-size:0.78rem;color:#444;padding:0.25rem;",
    "  word-break:break-word;}",
    "a{color:#06c;text-decoration:none;} a:hover{text-decoration:underline;}",
    "</style></head><body>",
    paste0("<h1>GSEA report — ", esc(title), "</h1>"),
    "<p>Preranked GSEA via <code>fgsea</code>. Top pathways shown per ",
    "collection (by BH-adjusted p-value, split into up- and down-regulated).</p>",
    "<nav><strong>Collections:</strong> ",
    paste(vapply(names(blocks), function(cc)
                 sprintf("<a href='#%s'>%s</a>", esc(cc), esc(cc)),
                 character(1)), collapse = " "),
    "</nav>"
  ), con)

  for (cc in names(blocks)) {
    b <- blocks[[cc]]
    writeLines(c(
      paste0("<h2 id='", esc(cc), "'>", esc(cc), "</h2>"),
      paste0("<p>Full results table: <a href='", esc(b$tsv), "'>",
             esc(b$tsv), "</a></p>")
    ), con)
    if (!is.na(b$table_png)) {
      writeLines(c(
        "<p><em>Top pathways summary (running-enrichment minigraphs):</em></p>",
        paste0("<p><img src='", esc(b$table_png),
               "' alt='gseaTable' style='max-width:100%;'></p>")
      ), con)
    }
    if (length(b$plot_paths)) {
      writeLines("<div class='plots'>", con)
      tc <- b$top_combined
      for (i in seq_along(b$plot_paths)) {
        cap <- sprintf("%s — NES=%.2f, padj=%.2g, size=%d",
                       tc$pathway[i], tc$NES[i], tc$padj[i], tc$size[i])
        writeLines(c(
          "<figure>",
          paste0("<img src='", esc(b$plot_paths[i]),
                 "' alt='", esc(tc$pathway[i]), "'>"),
          paste0("<figcaption>", esc(cap), "</figcaption>"),
          "</figure>"
        ), con)
      }
      writeLines("</div>", con)
    } else {
      writeLines("<p><em>No plottable pathways for this collection.</em></p>",
                 con)
    }
  }
  writeLines("</body></html>", con)
}

# --- Cohort display label -----------------------------------------------------
# Shows "E42K_affected" as "E42K carriers affected" in plots. File and contrast
# names are unchanged.
relabel_cohort <- function(x) {
  ifelse(x == "E42K_affected", "E42K carriers affected", x)
}

# --- Consensus figure helpers (scripts 10 and 11) -----------------------------
# The pathway lists and cell-type columns are defined in config.R.

# Display labels for cell-type keys, e.g. "CD4_Naive" -> "CD4 Naive".
# Used by the full heatmaps, script 10 and script 08.
.CELLTYPE_LABELS <- c(
  Bcell         = "B cells",
  NK_ILC        = "NK + ILC",
  DC            = "DCs",
  Hematopoietic = "Hematopoietic"
)
# Label for a subset key: the map above, otherwise "_" replaced by a space.
celltype_label <- function(ct) {
  hit <- unname(.CELLTYPE_LABELS[ct])
  ifelse(is.na(hit), gsub("_", " ", ct), hit)
}

# --- Violin cell-type groups (scripts 07 and 07a) -----------------------------

# Group labels in panel order.
celltype_group_levels <- function(panels = violin_celltype_groups) {
  unlist(lapply(panels, function(p) names(p$groups)), use.names = FALSE)
}

# One row per (group, L2 label) pair.
celltype_group_membership <- function(panels = violin_celltype_groups) {
  do.call(rbind, lapply(names(panels), function(nm) {
    p <- panels[[nm]]
    do.call(rbind, lapply(names(p$groups), function(g)
      data.frame(celltype_group = g,
                 l2             = as.character(p$groups[[g]]),
                 stringsAsFactors = FALSE)))
  }))
}

# Repeat each cell once per group it belongs to (T-cell groups overlap) and
# add celltype_group. Cells in no group are dropped and logged.
expand_celltype_groups <- function(md, celltype_col = "celltype",
                                   groups  = NULL,
                                   panels  = violin_celltype_groups,
                                   verbose = TRUE) {
  if (!celltype_col %in% names(md))
    stop("expand_celltype_groups(): no column '", celltype_col, "' in md")
  memb <- celltype_group_membership(panels)
  l2   <- as.character(md[[celltype_col]])

  if (isTRUE(verbose)) {
    miss <- !is.na(l2) & !(l2 %in% memb$l2)
    if (any(miss)) {
      tab <- sort(table(l2[miss]), decreasing = TRUE)
      log_msg("   [celltype_group] ", sum(miss), " cell(s) in ", length(tab),
              " label(s) covered by no group - excluded from the violins: ",
              paste(sprintf("%s (n=%d)", names(tab), as.integer(tab)),
                    collapse = ", "))
    }
  }

  if (!is.null(groups))
    memb <- memb[memb$celltype_group %in% groups, , drop = FALSE]

  pieces <- lapply(seq_len(nrow(memb)), function(i) {
    sel <- which(l2 == memb$l2[i])
    if (!length(sel)) return(NULL)
    d <- md[sel, , drop = FALSE]
    d$celltype_group <- memb$celltype_group[i]
    d
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) {
    out <- md[0, , drop = FALSE]
    out$celltype_group <- character(0)
    return(out)
  }
  do.call(rbind, pieces)
}

# "HBD_vs_E42K_affected" -> c("HBD", "E42K_affected"). Positive NES = higher in
# the second cohort.
parse_contrast_direction <- function(contrast_full) {
  cf <- sub("_batch2only$", "", contrast_full)
  parts <- strsplit(cf, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) == 1L) return(c(parts[1], "(unknown)"))
  if (length(parts) > 2L) parts <- c(parts[1], paste(parts[-1], collapse = "_vs_"))
  if (parts[1] == "GEM108_pre"  && parts[2] == "post") parts[2] <- "GEM108_post"
  if (parts[1] == "GEM108_post" && parts[2] == "pre")  parts[2] <- "GEM108_pre"
  parts
}

# Cohorts as shown to a reader: parsed names with the E42K_affected rename.
display_cohorts <- function(contrast_full) {
  relabel_cohort(parse_contrast_direction(contrast_full))
}

# Readable contrast name, e.g. "HBD vs T504S (batch 2 only)".
contrast_pretty <- function(contrast_full) {
  ch <- display_cohorts(contrast_full)
  base <- paste0(ch[1], " vs ", ch[2])
  if (grepl("_batch2only$", contrast_full)) paste0(base, " (batch 2 only)") else base
}

# NES direction caption for the consensus figures.
direction_caption <- function(contrast_full) {
  ch <- display_cohorts(contrast_full)
  sprintf(
    "Positive NES = higher expression in %s    |    Negative NES = higher expression in %s",
    ch[2], ch[1])
}

# NES orientation for the script 11 heatmaps. Positive = higher in GEM108 for
# GEM108 vs HBD, higher in GEM108_pre for pre vs post, and higher in the second
# cohort otherwise. Returns list(sign, pos, neg).
heatmap_orientation <- function(contrast_full) {
  raw  <- parse_contrast_direction(contrast_full)   # c(cohort1, cohort2), raw
  disp <- display_cohorts(contrast_full)            # same order, reader labels
  pos_idx <- 2L                                      # default: positive = cohort 2
  gem <- grepl("^GEM108", raw)
  if (any(gem) && "HBD" %in% raw) {
    pos_idx <- which(gem)[1]                         # the GEM108 side is positive
  } else if (setequal(raw, c("GEM108_pre", "GEM108_post"))) {
    pos_idx <- which(raw == "GEM108_pre")            # pre is positive
  }
  sign <- if (pos_idx == 2L) 1 else -1               # read NES is positive-for-2
  neg_idx <- if (pos_idx == 2L) 1L else 2L
  list(sign = sign, pos = disp[pos_idx], neg = disp[neg_idx])
}

# Subset/contrast pairs with a Hallmark fgsea table.
discover_consensus_contrasts <- function(collection = "h",
                                         cell_types = consensus_cell_types) {
  pat <- paste0("__gsea_", collection, "\\.tsv$")
  files <- list.files(gsea_dir, pattern = pat, full.names = FALSE)
  if (!length(files)) return(data.frame())
  stems <- sub(pat, "", files)
  m <- regexpr("__", stems, fixed = TRUE)
  ok <- m > 0L
  out <- data.frame(subset        = substring(stems[ok], 1L, m[ok] - 1L),
                    contrast_full = substring(stems[ok], m[ok] + 2L),
                    stringsAsFactors = FALSE)
  out[out$subset %in% cell_types, , drop = FALSE]
}

# Read a Hallmark fgsea table with standard column names, or NULL if missing.
read_gsea <- function(cell_type, contrast_full, collection = "h") {
  f <- file.path(gsea_dir,
                 paste0(cell_type, "__", contrast_full, "__gsea_", collection, ".tsv"))
  if (!file.exists(f)) return(NULL)
  df <- as.data.frame(data.table::fread(f))
  if (!nrow(df)) return(NULL)
  cols <- colnames(df)
  for (need in c("pathway", "pval", "padj", "NES", "size", "leadingEdge")) {
    j <- match(tolower(need), tolower(cols))
    if (!is.na(j)) cols[j] <- need
  }
  colnames(df) <- cols
  df
}

# Read script 08's full DO GSEA result for one subset and contrast, or NULL.
read_do_gsea_full <- function(cell_type, contrast_full) {
  f <- file.path(dea08_cache_dir,
                 paste0(cell_type, "__", contrast_full, "__do_gsea_full.rds"))
  if (!file.exists(f)) {
    log_msg("   [", cell_type, "] [DO] no shared full cache (",
            basename(f), ") - run script 08 first")
    return(NULL)
  }
  obj <- tryCatch(readRDS(f), error = function(e) {
    log_msg("   [", cell_type, "] [DO] unreadable cache: ",
            conditionMessage(e)); NULL })
  if (is.null(obj)) return(NULL)
  df <- tryCatch(as.data.frame(obj), error = function(e) {
    log_msg("   [", cell_type, "] [DO] could not coerce cache: ",
            conditionMessage(e)); NULL })
  if (is.null(df) || !nrow(df)) return(NULL)
  if (!all(c("Description", "NES", "p.adjust", "setSize") %in% colnames(df))) {
    log_msg("   [", cell_type, "] [DO] cache missing expected columns - skipping")
    return(NULL)
  }
  df
}

# Match consensus DO labels to a gseDO result by exact, case-insensitive term
# name. Unmatched labels are left out.
match_do_consensus <- function(do_df) {
  empty <- data.frame(label = character(0), NES = numeric(0),
                      padj = numeric(0), setSize = numeric(0),
                      stringsAsFactors = FALSE)
  if (is.null(do_df) || !nrow(do_df)) return(empty)
  desc_lc <- tolower(trimws(do_df$Description))
  rows <- lapply(do_consensus, function(entry) {
    hit <- which(desc_lc %in% entry$aliases)   # exact, case-insensitive
    if (!length(hit)) {
      log_msg("   [DO] consensus term not found (left blank): ", entry$label)
      return(NULL)
    }
    hit <- hit[which.min(do_df$p.adjust[hit])]   # in the (rare) tie, most sig.
    data.frame(label = entry$label, NES = do_df$NES[hit],
               padj = do_df$p.adjust[hit], setSize = do_df$setSize[hit],
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}

# Long table (label, cell_type, NES, padj) for one contrast across all
# consensus cell types. kind = "hallmark" or "do". Missing results are NA.
consensus_long <- function(contrast_full, kind, collection = "h",
                           cell_types = consensus_cell_types) {
  levels_vec <- if (kind == "hallmark") hallmark_levels else do_levels
  per_ct <- lapply(cell_types, function(ct) {
    base <- data.frame(label = levels_vec, cell_type = ct,
                       NES = NA_real_, padj = NA_real_,
                       stringsAsFactors = FALSE)
    if (kind == "hallmark") {
      gsea_df <- read_gsea(ct, contrast_full, collection)
      if (!is.null(gsea_df) && nrow(gsea_df)) {
        ids <- unname(hallmark_consensus[levels_vec])   # IDs in row order
        idx <- match(ids, gsea_df$pathway)
        base$NES  <- gsea_df$NES[idx]
        base$padj <- gsea_df$padj[idx]
      }
    } else {
      do_df   <- read_do_gsea_full(ct, contrast_full)
      matched <- match_do_consensus(do_df)
      if (nrow(matched)) {
        mi <- match(matched$label, base$label)
        base$NES[mi]  <- matched$NES
        base$padj[mi] <- matched$padj
      }
    }
    base
  })
  do.call(rbind, per_ct)
}

# --- Palettes -----------------------------------------------------------------
# Cell-type colours from ggplot's default hues, keyed on sorted names.
make_celltype_palette <- function(cell_types) {
  ct <- sort(unique(stats::na.omit(as.character(cell_types))))
  setNames(scales::hue_pal()(length(ct)), ct)
}

# HTO and capture palettes (palette.colors, recycled if needed).
make_hto_palette <- function(htos) {
  lv <- sort(unique(stats::na.omit(as.character(htos))))
  cols <- grDevices::palette.colors(length(lv), hto_palette_name,
                                    recycle = TRUE)
  setNames(cols, lv)
}
# Same for captures.
make_capture_palette <- function(captures) {
  lv <- sort(unique(stats::na.omit(as.character(captures))))
  cols <- grDevices::palette.colors(length(lv), capture_palette_name,
                                    recycle = TRUE)
  setNames(cols, lv)
}

# --- Plots --------------------------------------------------------------------
# Volcano plot. Significant = BH p < 0.05 and |logFC| > 1. Top 50 genes by
# p-value are labelled. contrast_levels = c(neg_label, pos_label) adds direction
# labels.
plot_volcano <- function(tt, title = "", logfc_col = "logFC",
                         p_col = "adj.P.Val", pad_col = "P.Value",
                         logfc_thresh = 1, padj_thresh = 0.05,
                         label_top = 50,
                         contrast_levels = NULL) {
  df <- tt
  df$sig <- df[[p_col]] < padj_thresh & abs(df[[logfc_col]]) > logfc_thresh
  df$neglog10p <- -log10(df[[pad_col]])
  # Top 50 genes by p-value.
  top <- df[order(df[[pad_col]]), ]
  top <- head(top, label_top)

  # Direction labels: first level = negative logFC side.
  if (!is.null(contrast_levels) && length(contrast_levels) == 2L) {
    neg_lab <- contrast_levels[1]
    pos_lab <- contrast_levels[2]
    subtitle <- sprintf(
      "Log2FC > 0: higher in %s   •   Log2FC < 0: higher in %s",
      pos_lab, neg_lab
    )
    xlab <- sprintf("Log2 Fold Change  (positive = higher in %s)", pos_lab)
  } else {
    subtitle <- NULL
    xlab <- "Log2 Fold Change"
  }

  ggplot(df, aes(x = .data[[logfc_col]], y = .data[["neglog10p"]],
                 colour = .data[["sig"]])) +
    geom_point(size = 6) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick")) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh),
               linetype = 2, alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 1, alpha = 0.3) +
    geom_hline(yintercept = -log10(padj_thresh),
               linetype = 2, alpha = 0.4) +
    ggrepel::geom_text_repel(data = top, aes(label = gene),
                             box.padding = 0.5, point.padding = 0.5,
                             size = 3, max.overlaps = Inf,
                             colour = "black", show.legend = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = xlab, y = "-log10(P.Value)") +
    theme_minimal() +
    theme(plot.title       = element_text(hjust = 0.5),
          plot.subtitle    = element_text(hjust = 0.5, size = 10,
                                          colour = "grey25"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank(),
          axis.line        = element_line(color = "black"),
          legend.position  = "none")
}

# MA plot (mean expression vs logFC), used for the descriptive contrasts.
# contrast_levels works as in plot_volcano.
plot_ma <- function(tt, title = "", logfc_col = "logFC",
                    mean_col = "mean_logCPM", p_col = NULL,
                    logfc_thresh = 1, padj_thresh = 0.05,
                    label_top = 20,
                    contrast_levels = NULL) {
  df <- tt
  if (!is.null(p_col) && p_col %in% colnames(df)) {
    df$sig <- df[[p_col]] < padj_thresh & abs(df[[logfc_col]]) > logfc_thresh
  } else {
    df$sig <- abs(df[[logfc_col]]) > logfc_thresh
  }
  top <- df[df$sig, ]
  top <- top[order(-abs(top[[logfc_col]])), ]
  top <- head(top, label_top)

  if (!is.null(contrast_levels) && length(contrast_levels) == 2L) {
    neg_lab <- contrast_levels[1]
    pos_lab <- contrast_levels[2]
    subtitle <- sprintf(
      "Log2FC > 0: higher in %s   •   Log2FC < 0: higher in %s",
      pos_lab, neg_lab
    )
    ylab <- sprintf("log2 fold-change  (positive = higher in %s)", pos_lab)
  } else {
    subtitle <- NULL
    ylab <- "log2 fold-change"
  }

  ggplot(df, aes(x = .data[[mean_col]], y = .data[[logfc_col]],
                 colour = .data[["sig"]])) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
    geom_hline(yintercept = c(-logfc_thresh, logfc_thresh), linetype = 2, alpha = 0.4) +
    geom_hline(yintercept = 0, alpha = 0.3) +
    ggrepel::geom_text_repel(data = top, aes(label = gene),
                             size = 3, max.overlaps = Inf) +
    labs(title = title, subtitle = subtitle,
         x = "mean log2 CPM", y = ylab) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle   = element_text(hjust = 0.5, size = 10,
                                         colour = "grey25"),
          legend.position = "none")
}

# --- Offline Azimuth homolog table --------------------------------------------
# RunAzimuth downloads a homolog table, which fails without internet. This makes
# it read a local copy instead. Pass a plain file path, not a file:// URL. Call
# once before RunAzimuth.
patch_azimuth_homologs <- function(homologs_path) {
  if (!file.exists(homologs_path)) {
    stop("azimuth_homologs_path does not exist: ", homologs_path,
         " - required because Gadi compute nodes cannot fetch homologs.rds.")
  }
  # Check the file is a readable RDS.
  smoke <- tryCatch(readRDS(homologs_path), error = function(e) e)
  if (inherits(smoke, "error")) {
    head_bytes <- tryCatch(
      paste(readBin(homologs_path, "raw", 32), collapse = " "),
      error = function(e) "(could not read)")
    stop("homologs file at ", homologs_path, " is not a valid RDS:\n  ",
         conditionMessage(smoke),
         "\n  file size:      ", file.size(homologs_path), " bytes",
         "\n  first 32 bytes: ", head_bytes,
         "\nDownload the canonical file from https://seurat.nygenome.org/azimuth/references/homologs.rds ",
         "on a machine with internet and copy it to this path.")
  }
  rm(smoke)
  local_path <- normalizePath(homologs_path, mustWork = TRUE)
  ns   <- asNamespace("Azimuth")
  orig <- get("ConvertGeneNames", envir = ns)
  patched <- function(object, reference.names, homolog.table) {
    # Call the original function with the local path.
    orig(object = object, reference.names = reference.names,
         homolog.table = local_path)
  }
  unlockBinding("ConvertGeneNames", ns)
  on.exit(tryCatch(lockBinding("ConvertGeneNames", ns), error = function(e) NULL),
          add = TRUE)
  assign("ConvertGeneNames", patched, envir = ns)
  lockBinding("ConvertGeneNames", ns)
  invisible(TRUE)
}

# --- Spot-matrix PDF ----------------------------------------------------------
# Combine the single-grid Hallmark spot matrices into one PDF, one page per
# contrast. Used by script 06 and combine_mod_spot_matrices_pdf.R.
bundle_mod_spot_matrices_pdf <- function(sm_dir, out_pdf = NULL, pattern = NULL) {
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("bundle_mod_spot_matrices_pdf() needs the 'png' package (png::readPNG).")
  }
  mod_tail <- "_hallmark_FDR0\\.25_reduced_spotMatrix_mod_pos_neg\\.png$"
  files <- list.files(sm_dir, pattern = "_mod_pos_neg\\.png$", full.names = TRUE)
  if (length(files) == 0L) {
    log_msg("   [mod PDF] no *_mod_pos_neg.png in ", sm_dir, " - nothing to bundle")
    return(invisible(NULL))
  }
  stems <- sub(mod_tail, "", basename(files))
  if (!is.null(pattern)) {
    keep <- grepl(pattern, stems, perl = TRUE)
    if (!any(keep)) {
      stop("pattern '", pattern, "' matched no files. Available stems: ",
           paste(stems, collapse = ", "))
    }
    files <- files[keep]; stems <- stems[keep]
  }
  ord <- order(stems); files <- files[ord]; stems <- stems[ord]
  if (is.null(out_pdf)) out_pdf <- file.path(sm_dir, "all_mod_spot_matrices.pdf")
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

  log_msg("== Bundling ", length(files), " mod spot matrices -> ", out_pdf)
  for (s in stems) log_msg("   - ", s)

  # Match the 8 x 8 in PNGs.
  grDevices::pdf(out_pdf, width = 8, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_along(files)) {
    img <- png::readPNG(files[i])
    grid::grid.newpage()
    grid::grid.raster(img, interpolate = FALSE)
    grid::grid.text(sprintf("Page %d", i), x = 0.99, y = 0.01,
                    just = c("right", "bottom"),
                    gp = grid::gpar(fontsize = 8, col = "grey50"))
  }
  log_msg("== done: ", out_pdf)
  invisible(out_pdf)
}
