# ------------------------------------------------------------------------------
# 05_differential_expression.R
# Differential expression (DE) and GSEA. The work is split into tasks, one per
# cell subset and contrast (subsets and contrasts are set in config.R). Each
# task writes a DE table, a volcano or MA plot and GSEA results.
#
# Three kinds of contrast:
#   pseudobulk      cohort vs cohort. Counts are summed per donor and cell type,
#                   then tested with limma-voom (design ~0 + cohort + batch).
#   within_patient  GEM108 pre vs post treatment. Single-cell limma, with
#                   cells from the same capture treated as correlated.
#   descriptive     GEM108 vs HBD. One patient on one side, so fold changes and
#                   GSEA only, no per-gene p-values.
#
# Run:  Rscript new_scripts/05_differential_expression.R [options]
#   (no args)           run every task in one process (slow, high memory)
#   --list              write the task table to pipeline/DE/tasks.tsv and exit
#   --task K [K2 ...]   run the given task ids (used by run_script_05.pbs)
#   --gsea-only         redo GSEA from existing DE tables
#   --pdf-report        build pipeline/DE/DE_plots_report.pdf from DE tables
#   --force             re-run tasks whose output already exists
#
# Input:  pipeline/merged/seurat_merged_harmony_azimuth.rds   (from script 04)
# Output: pipeline/DE/{subset}__{contrast}.tsv and plots
#         pipeline/GSEA/{subset}__{contrast}__gsea_{collection}.tsv and plots
# ------------------------------------------------------------------------------

# Load a local GLPK build before Seurat (needed on Gadi). Edit or remove.
if (file.exists("/path/to/libglpk.so.40")) dyn.load("/path/to/libglpk.so.40")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(ggrepel)
  library(data.table)
  library(cowplot)        # plot_grid + ggdraw for the PDF report header
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

in_path <- annotated_rds_path()  # raw or CellSweep object, see config.R

# --- Argument parsing ---------------------------------------------------------
# Turn the command-line options into a list: mode, task ids and flags.
parse_args <- function(argv) {
  out <- list(
    mode       = "serial",   # "serial" | "list" | "task"
    task_ids   = integer(),
    force      = FALSE,
    gsea_only  = FALSE,
    pdf_report = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--list") {
      out$mode <- "list"
      i <- i + 1L
    } else if (a == "--force") {
      out$force <- TRUE
      i <- i + 1L
    } else if (a == "--gsea-only") {
      out$gsea_only <- TRUE
      i <- i + 1L
    } else if (a == "--pdf-report") {
      out$pdf_report <- TRUE
      i <- i + 1L
    } else if (a == "--task") {
      out$mode <- "task"
      i <- i + 1L
      # Consume every following argv item that parses as a positive integer.
      while (i <= length(argv) && grepl("^[0-9]+$", argv[i])) {
        out$task_ids <- c(out$task_ids, as.integer(argv[i]))
        i <- i + 1L
      }
      if (length(out$task_ids) == 0L) {
        stop("--task requires at least one integer task id")
      }
    } else {
      stop("Unrecognised argument: ", a)
    }
  }
  out
}
# commandArgs(trailingOnly = TRUE) returns the options typed after the script name.
opts <- parse_args(commandArgs(trailingOnly = TRUE))

# Cell-type column (Azimuth L2).
cell_type_col <- "predicted.celltype.l2"

# --- Task table ---------------------------------------------------------------
# One row per subset x contrast. task_id is stable as long as the order of
# de_cell_subsets and the contrast lists in config.R does not change.
build_tasks <- function() {
  rows <- list()
  add_row <- function(...) {
    rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
  }
  # Every subset gets every contrast. HBD_vs_* contrasts also get a repeat using
  # the May batch only (batch2_only = TRUE).
  for (subset_name in names(de_cell_subsets)) {
    # A. Pseudobulk cross-cohort
    for (cname in names(pseudobulk_contrasts)) {
      add_row(subset = subset_name, kind = "pseudobulk",
              contrast = cname, batch2_only = FALSE)
      if (grepl("^HBD_vs_", cname)) {
        add_row(subset = subset_name, kind = "pseudobulk",
                contrast = cname, batch2_only = TRUE)
      }
    }
    # B. Within-patient
    for (cname in names(within_patient_contrasts)) {
      add_row(subset = subset_name, kind = "within_patient",
              contrast = cname, batch2_only = FALSE)
    }
    # C. Descriptive
    for (cname in names(descriptive_contrasts)) {
      add_row(subset = subset_name, kind = "descriptive",
              contrast = cname, batch2_only = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  # contrast_full is the name used in output files.
  out$contrast_full <- ifelse(out$batch2_only,
                              paste0(out$contrast, "_batch2only"),
                              out$contrast)
  out$output_tsv <- file.path(
    de_dir,
    paste0(out$subset, "__", out$contrast_full, ".tsv")
  )
  out$task_id <- seq_len(nrow(out))
  out[, c("task_id", "subset", "kind", "contrast",
          "batch2_only", "contrast_full", "output_tsv")]
}

# --- Output paths -------------------------------------------------------------
# File names for one subset and contrast.
de_paths <- function(subset_name, contrast_name) {
  stub <- file.path(de_dir, paste0(subset_name, "__", contrast_name))
  list(
    tt_file       = paste0(stub, ".tsv"),
    volcano_file  = paste0(stub, "__volcano.png"),
    ma_file       = paste0(stub, "__ma.png"),
    gsea_prefix   = file.path(gsea_dir, paste0(subset_name, "__", contrast_name))
  )
}

# --- Per-task data loader -----------------------------------------------------
# Loads the object, keeps only the cells and assay layer this task needs, then
# frees the rest.
load_task_data <- function(task) {
  log_msg("    loading ", in_path)
  seu <- readRDS(in_path)
  stopifnot(cell_type_col %in% colnames(seu@meta.data))

  # Cell metadata, with helper columns used to pick cells below.
  meta <- seu@meta.data
  meta$cell_type <- meta[[cell_type_col]]
  # cohort_or_patient: cohort, except GEM108 cells are labelled GEM108.
  meta$cohort_or_patient <- ifelse(meta$patient == "GEM108",
                                   "GEM108", meta$cohort)
  # GEM108 split into GEM108_pre / GEM108_post, others keep their cohort.
  tx_short <- sub("_treatment$", "", meta$treatment)   # pre / post / NA
  meta$cohort_or_patient_tx <- ifelse(
    meta$patient == "GEM108" & !is.na(tx_short) & nzchar(tx_short),
    paste0("GEM108_", tx_short),
    meta$cohort
  )

  # Azimuth L2 labels in this subset (NULL = all cells).
  cell_types <- de_cell_subsets[[task$subset]]

  # --- Choose the cells for this task -----------------------------------------
  # Pseudobulk: cells from the two cohorts in the contrast.
  if (task$kind == "pseudobulk") {
    levels <- pseudobulk_contrasts[[task$contrast]]
    # Collapse group levels such as E42K_carriers into their member cohorts. A
    # member cohort that is also a level in this contrast keeps its own label.
    eff_cohort <- meta$cohort
    if (exists("cohort_groups") && length(cohort_groups)) {
      for (grp_name in names(cohort_groups)) {
        if (grp_name %in% levels) {
          members <- setdiff(cohort_groups[[grp_name]], levels)
          if (length(members)) {
            eff_cohort[!is.na(eff_cohort) & eff_cohort %in% members] <- grp_name
          }
        }
      }
    }
    meta$cohort <- eff_cohort
    # GEM108 contributes only its pre-treatment cells, as one donor, to contrasts
    # using E42K_affected or E42K_carriers.
    gem108_pre <- meta$patient == "GEM108" &
            !is.na(tx_short) & tx_short == "pre"
    keep <- (meta$patient != "GEM108" | gem108_pre) &
            meta$cohort %in% levels &
            !is.na(meta$patient) &
            !is.na(meta$cohort) &
            !is.na(meta$cell_type)
    # Batch-2-only repeat: May batch cells only.
    if (isTRUE(task$batch2_only)) keep <- keep & meta$batch == "May"
  } else if (task$kind == "within_patient") {
    spec <- within_patient_contrasts[[task$contrast]]
    # Within-patient: GEM108 cells labelled pre or post treatment.
    keep <- meta$patient == spec$patient &
            meta[[spec$group_col]] %in% spec$levels &
            !is.na(meta[[spec$group_col]])
  } else if (task$kind == "descriptive") {
    spec <- descriptive_contrasts[[task$contrast]]
    # Descriptive: GEM108 (all, pre or post) and HBD cells.
    grp <- meta[[spec$group_col]]
    keep <- grp %in% spec$levels &
            !is.na(grp) &
            !is.na(meta$patient)
  } else {
    stop("Unknown task kind: ", task$kind)
  }
  # Keep only cells of this subset's cell types.
  if (!is.null(cell_types)) keep <- keep & meta$cell_type %in% cell_types

  # Drop the GEM108 May-batch samples (no treatment label) from every task.
  keep <- keep & !(meta$patient %in% "GEM108" & meta$batch %in% "May")

  log_msg("    keep mask: ", sum(keep), " / ", length(keep), " cells")
  # Skip the task if fewer than 50 cells are left.
  if (sum(keep) < 50L) {
    rm(seu); gc(verbose = FALSE)
    return(NULL)
  }

  # Return the metadata and expression matrix for the kept cells only.
  meta_sub <- meta[keep, , drop = FALSE]

  # within_patient uses log-normalised data; the other kinds use raw counts.
  if (task$kind == "within_patient") {
    m_sub <- GetAssayData(seu, assay = "RNA", layer = "data")[, keep, drop = FALSE]
    out <- list(meta = meta_sub, logcounts = m_sub)
  } else {
    m_sub <- GetAssayData(seu, assay = "RNA", layer = "counts")[, keep, drop = FALSE]
    out <- list(meta = meta_sub, counts = m_sub)
  }

  # Drop the full Seurat object before any DE work runs.
  rm(seu); gc(verbose = FALSE)
  out
}

# --- A. Pseudobulk between cohorts --------------------------------------------
run_pseudobulk_contrast <- function(task, td) {
  contrast_name <- task$contrast_full
  levels <- pseudobulk_contrasts[[task$contrast]]
  log_msg("-- [pseudobulk", if (task$batch2_only) " batch2-only" else "",
          "] ", task$subset, " :: ", contrast_name)

  paths <- de_paths(task$subset, contrast_name)
  # pseudobulk_de() (utils.R): sum counts per donor and cell type, then
  # limma-voom with ~0 + cohort (+ batch). Returns the limma topTable with
  # BH-adjusted p-values. The batch2_only repeat has one batch, so no batch term.
  tt <- pseudobulk_de(
    counts        = td$counts,
    col_meta      = td$meta,
    contrast      = levels,
    include_batch = !task$batch2_only
  )
  readr::write_tsv(tt, paths$tt_file)
  log_msg("   wrote ", paths$tt_file, " (", nrow(tt), " genes, ",
          sum(tt$adj.P.Val < 0.05, na.rm = TRUE), " at BH<0.05)")

  v <- plot_volcano(
    tt,
    title = paste0(task$subset, ": ", contrast_name,
                   if (task$batch2_only) " (batch-2-only sensitivity)" else ""),
    p_col = "adj.P.Val", pad_col = "P.Value",
    # Positive logFC = higher in the second level.
    contrast_levels = levels
  )
  ggsave(paths$volcano_file, v, width = 12, height = 12, dpi = 150)

  # GSEA on genes ranked by the limma t-statistic. run_gsea() (utils.R) runs
  # fgsea on each MSigDB collection in config.R and writes tables, plots and
  # an HTML summary.
  if ("t" %in% colnames(tt) && !all(is.na(tt$t))) {
    ranked <- setNames(tt$t, tt$gene)
    ranked <- ranked[!is.na(ranked)]
    run_gsea(ranked, out_prefix = paths$gsea_prefix)
  }
  invisible(NULL)
}

# --- B. Within-patient single-cell: GEM108 pre vs post ------------------------
run_within_patient_contrast <- function(task, td) {
  contrast_name <- task$contrast_full
  spec <- within_patient_contrasts[[task$contrast]]
  log_msg("-- [within-patient] ", task$subset, " :: ", contrast_name)

  # Drop genes expressed in fewer than 10 cells.
  keep_gene <- Matrix::rowSums(td$logcounts > 0) >= 10
  logc_sub <- td$logcounts[keep_gene, , drop = FALSE]
  log_msg("    ", sum(keep_gene), " genes survive expression filter")

  paths <- de_paths(task$subset, contrast_name)
  # single_cell_de_paired() (utils.R): limma on log-normalised expression per
  # cell, with duplicateCorrelation so cells from the same capture are not
  # treated as independent.
  tt <- single_cell_de_paired(
    logcounts = as.matrix(logc_sub),
    col_meta  = td$meta,
    group_col = spec$group_col,
    levels    = spec$levels,
    block_col = "capture"
  )
  rm(logc_sub); gc(verbose = FALSE)

  readr::write_tsv(tt, paths$tt_file)
  log_msg("   wrote ", paths$tt_file, " (", nrow(tt), " genes, rho=",
          signif(attr(tt, "rho"), 3), ")")

  v <- plot_volcano(
    tt,
    title = paste0(task$subset, ": ", contrast_name,
                   "  (single-cell limma, blocked on capture)"),
    p_col = "adj.P.Val", pad_col = "P.Value",
    # Positive logFC = higher in post_treatment.
    contrast_levels = spec$levels
  )
  ggsave(paths$volcano_file, v, width = 12, height = 12, dpi = 150)

  if ("t" %in% colnames(tt)) {
    ranked <- setNames(tt$t, tt$gene)
    ranked <- ranked[!is.na(ranked)]
    run_gsea(ranked, out_prefix = paths$gsea_prefix)
  }
  invisible(NULL)
}

# --- C. Descriptive: GEM108 vs HBD (no per-gene p-values) ---------------------
run_descriptive_contrast <- function(task, td) {
  contrast_name <- task$contrast_full
  spec <- descriptive_contrasts[[task$contrast]]
  log_msg("-- [descriptive]   ", task$subset, " :: ", contrast_name)

  paths <- de_paths(task$subset, contrast_name)
  # descriptive_de() (utils.R): pseudobulk per patient, then the difference in
  # mean log-CPM between the two groups. No p-values.
  tt <- descriptive_de(
    counts      = td$counts,
    col_meta    = td$meta,
    patient_col = "patient",
    group_col   = spec$group_col,
    levels      = spec$levels
  )
  readr::write_tsv(tt, paths$tt_file)
  log_msg("   wrote ", paths$tt_file, " (", nrow(tt),
          " genes - NO p-values, descriptive only)")

  # MA plot: mean expression vs fold change (plot_ma() in utils.R).
  ma <- plot_ma(
    tt,
    title = paste0(task$subset, ": ", contrast_name,
                   "  (n=1 on GEM108 side — descriptive only)"),
    logfc_col = "logFC", mean_col = "mean_logCPM",
    p_col = NULL, logfc_thresh = 1,
    # Positive logFC = higher in the second level (HBD).
    contrast_levels = spec$levels
  )
  ggsave(paths$ma_file, ma, width = 6, height = 5, dpi = 150)

  ranked <- setNames(tt$logFC, tt$gene)
  ranked <- ranked[!is.na(ranked)]
  run_gsea(ranked, out_prefix = paths$gsea_prefix)
  invisible(NULL)
}

# --- GSEA-only re-run from an existing DE table -------------------------------
run_gsea_only_task <- function(task) {
  log_msg("=== Task ", task$task_id, " [GSEA-only]: ", task$subset,
          " :: ", task$kind, " :: ", task$contrast_full)

  if (!file.exists(task$output_tsv)) {
    log_msg("    SKIPPING - DE TSV missing: ", task$output_tsv)
    return(invisible(NULL))
  }

  tt <- tryCatch(readr::read_tsv(task$output_tsv, show_col_types = FALSE),
                 error = function(e) {
                   log_msg("    READ FAILED: ", conditionMessage(e)); NULL
                 })
  if (is.null(tt) || !"gene" %in% colnames(tt)) {
    log_msg("    SKIPPING - DE TSV unreadable or missing 'gene' column")
    return(invisible(NULL))
  }

  # Rank by t-statistic, or by logFC for descriptive contrasts.
  rank_col <- if (task$kind == "descriptive") "logFC" else "t"
  if (!rank_col %in% colnames(tt)) {
    log_msg("    SKIPPING - expected rank column '", rank_col,
            "' not in TSV (have: ", paste(colnames(tt), collapse = ","), ")")
    return(invisible(NULL))
  }
  ranked <- setNames(tt[[rank_col]], tt$gene)
  ranked <- ranked[!is.na(ranked) & !is.na(names(ranked))]
  if (length(ranked) == 0L) {
    log_msg("    SKIPPING - no usable ranked stats")
    return(invisible(NULL))
  }

  paths <- de_paths(task$subset, task$contrast_full)
  run_gsea(ranked, out_prefix = paths$gsea_prefix)
  invisible(NULL)
}

# --- DE plots PDF report ------------------------------------------------------
# One page per subset x contrast, rebuilt from the saved DE tables.
# The two groups compared in a task.
contrast_levels_for_task <- function(task) {
  switch(task$kind,
    pseudobulk     = pseudobulk_contrasts[[task$contrast]],
    within_patient = within_patient_contrasts[[task$contrast]]$levels,
    descriptive    = descriptive_contrasts[[task$contrast]]$levels,
    NULL
  )
}

# Contrast name for the page header.
contrast_pretty_05 <- function(task) {
  lv <- contrast_levels_for_task(task)
  base <- if (!is.null(lv) && length(lv) == 2L) {
    paste0(lv[1], " vs ", lv[2])
  } else {
    task$contrast
  }
  if (isTRUE(task$batch2_only)) paste0(base, " (batch 2 only)") else base
}

# Page header text for each kind of contrast.
kind_label_05 <- c(
  pseudobulk     = "Volcano — pseudobulk limma-voom",
  within_patient = "Volcano — single-cell limma (paired within-patient)",
  descriptive    = "MA plot — descriptive (n=1, no p-values)"
)

# One PDF page: header plus volcano or MA plot.
build_de_plot_page <- function(task) {
  if (!file.exists(task$output_tsv)) {
    log_msg("    [pdf] skip - TSV missing: ", basename(task$output_tsv))
    return(NULL)
  }
  tt <- tryCatch(readr::read_tsv(task$output_tsv, show_col_types = FALSE),
                 error = function(e) {
                   log_msg("    [pdf] read failed: ", conditionMessage(e))
                   NULL
                 })
  if (is.null(tt) || !"gene" %in% colnames(tt) ||
      !"logFC" %in% colnames(tt)) {
    log_msg("    [pdf] skip - TSV unreadable or missing required columns")
    return(NULL)
  }

  lv <- contrast_levels_for_task(task)

  plot_title <- paste0(task$subset, " :: ", task$contrast_full)

  body <- tryCatch({
    if (task$kind == "descriptive") {
      plot_ma(tt,
              title           = plot_title,
              logfc_col       = "logFC",
              mean_col        = "mean_logCPM",
              p_col           = NULL,
              logfc_thresh    = 1,
              contrast_levels = lv)
    } else {
      plot_volcano(tt,
                   title           = plot_title,
                   p_col           = "adj.P.Val",
                   pad_col         = "P.Value",
                   contrast_levels = lv)
    }
  }, error = function(e) {
    log_msg("    [pdf] plot build failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(body)) return(NULL)

  kind_text <- if (task$kind %in% names(kind_label_05)) {
    kind_label_05[[task$kind]]
  } else {
    task$kind
  }
  page_header <- sprintf("%s  —  %s  —  %s",
                         contrast_pretty_05(task),
                         task$subset,
                         kind_text)

  header <- cowplot::ggdraw() +
    cowplot::draw_label(page_header, fontface = "bold",
                        size = 13, hjust = 0.5)
  cowplot::plot_grid(header, body, ncol = 1,
                     rel_heights = c(0.04, 0.96))
}

# Write every selected task's plot to one PDF.
build_de_pdf_report <- function(tasks, selected_ids) {
  out_path <- file.path(de_dir, "DE_plots_report.pdf")
  sel <- tasks[tasks$task_id %in% selected_ids, , drop = FALSE]
  # Order pages by contrast, then subset, then kind.
  subset_levels <- names(de_cell_subsets)
  sel <- sel[order(sel$contrast_full,
                   match(sel$subset, subset_levels),
                   sel$kind), , drop = FALSE]

  log_msg("== Building DE plots PDF: ", out_path,
          " (", nrow(sel), " task(s) selected)")
  pdf(out_path, width = 12, height = 12, onefile = TRUE)
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
      "Differential Expression — Volcano / MA plots",
      fontface = "bold", size = 22, y = 0.82) +
    cowplot::draw_label(
      paste0(
        "Subsets: ", paste(names(de_cell_subsets), collapse = ", "), "\n",
        "Contrasts: ", length(unique(sel$contrast_full)),
        "    Pages: ", nrow(sel), "\n\n",
        "DE kinds (one plot per page):\n",
        "  • Cross-cohort (pseudobulk limma-voom):  volcano plot\n",
        "  • GEM108 pre vs post (within-patient limma + duplicateCorrelation):  volcano plot\n",
        "  • GEM108 vs HBD (descriptive, n=1):  MA plot, no p-values\n\n",
        "Direction convention (all DE kinds):\n",
        "  positive logFC = higher in the SECOND level of the contrast\n",
        "  (e.g. HBD_vs_E42K_affected → positive = higher in E42K_affected;\n",
        "        GEM108_pre_vs_post → positive = higher in post_treatment;\n",
        "        GEM108_pre_vs_HBD → positive = higher in HBD).\n",
        "Each plot's subtitle and axis label spells out the cohort names."
      ),
      size = 11, y = 0.42, lineheight = 1.35) +
    cowplot::draw_label(
      paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      size = 10, y = 0.05, colour = "grey40")
  print(cover)
  page_number()

  built <- 0L
  for (k in seq_len(nrow(sel))) {
    t <- sel[k, , drop = FALSE]
    log_msg("    [pdf] page ", k, "/", nrow(sel), ": ",
            t$subset, " :: ", t$contrast_full, " (", t$kind, ")")
    pg <- build_de_plot_page(t)
    if (!is.null(pg)) {
      print(pg)
      page_number()
      built <- built + 1L
    }
  }

  grDevices::dev.off()
  log_msg("== DE plots PDF report done: ", out_path,
          "  (", built, " of ", nrow(sel), " requested page(s) drawn)")
  invisible(out_path)
}

# --- Run one task -------------------------------------------------------------
# Skip if the output exists (unless --force), load the cells, run the DE and
# GSEA. Errors are logged and the task is skipped.
run_task <- function(task, force = FALSE) {
  log_msg("=== Task ", task$task_id, ": ", task$subset, " :: ",
          task$kind, " :: ", task$contrast_full)

  if (!force && file.exists(task$output_tsv)) {
    log_msg("    skipping (output exists): ", task$output_tsv)
    return(invisible(NULL))
  }

  td <- tryCatch(load_task_data(task),
                 error = function(e) {
                   log_msg("    LOAD FAILED: ", conditionMessage(e))
                   NULL
                 })
  if (is.null(td)) {
    log_msg("    skipping (no cells / load failed)")
    return(invisible(NULL))
  }

  tryCatch({
    switch(task$kind,
      pseudobulk      = run_pseudobulk_contrast(task, td),
      within_patient  = run_within_patient_contrast(task, td),
      descriptive     = run_descriptive_contrast(task, td),
      stop("Unknown task kind: ", task$kind)
    )
  }, error = function(e) {
    log_msg("    FAILED: ", conditionMessage(e))
  })

  rm(td); gc(verbose = FALSE)
  invisible(NULL)
}

# --- Main ---------------------------------------------------------------------
# Build the task table from config.R.
tasks <- build_tasks()

# --list: write the task table and stop.
if (opts$mode == "list") {
  out_path <- file.path(de_dir, "tasks.tsv")
  readr::write_tsv(tasks, out_path)
  log_msg("Wrote task table: ", out_path, " (", nrow(tasks), " tasks)")
  cat("\nTask table preview:\n")
  print(tasks, row.names = FALSE)
  quit(save = "no", status = 0)
}

# --gsea-only and --pdf-report only read DE tables.
if (!opts$gsea_only && !opts$pdf_report && !file.exists(in_path)) {
  stop("Missing annotated object: ", in_path,
       " - did script 04 finish?")
}

# Task ids to run.
selected_ids <- if (opts$mode == "task") opts$task_ids else tasks$task_id

# --gsea-only: redo GSEA for the selected tasks and stop.
if (opts$gsea_only) {
  bad <- setdiff(selected_ids, tasks$task_id)
  if (length(bad)) {
    stop("Unknown task id(s): ", paste(bad, collapse = ", "),
         "  (valid range: 1..", nrow(tasks), ")")
  }
  log_msg("== GSEA-only mode: regenerating GSEA outputs for ",
          length(selected_ids), " task(s) from existing DE TSVs.")
  for (id in selected_ids) {
    run_gsea_only_task(tasks[tasks$task_id == id, , drop = FALSE])
  }
  log_msg("== 05_differential_expression.R --gsea-only done.")
  quit(save = "no", status = 0)
}

# --pdf-report: build the PDF and stop.
if (opts$pdf_report) {
  bad <- setdiff(selected_ids, tasks$task_id)
  if (length(bad)) {
    stop("Unknown task id(s): ", paste(bad, collapse = ", "),
         "  (valid range: 1..", nrow(tasks), ")")
  }
  build_de_pdf_report(tasks, selected_ids)
  log_msg("== 05_differential_expression.R --pdf-report done.")
  quit(save = "no", status = 0)
}

# --task: run the given task ids and stop.
if (opts$mode == "task") {
  bad <- setdiff(opts$task_ids, tasks$task_id)
  if (length(bad)) {
    stop("Unknown task id(s): ", paste(bad, collapse = ", "),
         "  (valid range: 1..", nrow(tasks), ")")
  }
  for (id in opts$task_ids) {
    run_task(tasks[tasks$task_id == id, , drop = FALSE], force = opts$force)
  }
  log_msg("== 05_differential_expression.R --task done.")
  quit(save = "no", status = 0)
}

# Serial mode: every task in one process. Uses more memory than --task.
log_msg("== Serial mode: running all ", nrow(tasks), " tasks in one R process.")
log_msg("   For better memory behaviour use --task <id> via a PBS array.")
for (id in tasks$task_id) {
  run_task(tasks[tasks$task_id == id, , drop = FALSE], force = opts$force)
}

log_msg("== 05_differential_expression.R done.")
log_msg("   DE tables under ", de_dir)
log_msg("   GSEA under ",     gsea_dir)
