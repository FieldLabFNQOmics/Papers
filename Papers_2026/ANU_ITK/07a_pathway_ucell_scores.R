# ------------------------------------------------------------------------------
# 07a_pathway_ucell_scores.R
# The same analysis as script 07 (same six gene sets, comparisons, violins,
# stats and PDF report), but pathway activity is scored with UCell instead of
# AddModuleScore. UCell scores each cell from the ranks of the pathway genes
# within that cell, so scores run from 0 to 1 and do not depend on the other
# cells in the object.
#
# Run:  Rscript new_scripts/07a_pathway_ucell_scores.R [options]
#   Same options as script 07: --list, --comparison NAME, --pathway NAME,
#   --no-pdf, --no-png, --force
#
# Input:  pipeline/merged/seurat_merged_harmony_azimuth.rds   (from script 04)
# Output: pipeline/pathway_reports_ucell/  (violins/, stats/, PDF report)
# ------------------------------------------------------------------------------

# Load a local GLPK build before Seurat (needed on Gadi). Edit or remove.
if (file.exists("/path/to/libglpk.so.40")) {
  dyn.load("/path/to/libglpk.so.40")
}

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(msigdbr)
  library(grid)
  library(gridExtra)
})

# Load settings (config.R) and shared functions (utils.R), then create any
# missing output folders.
source(here::here("new_scripts", "config.R"))
source(here::here("new_scripts", "utils.R"))
ensure_pipeline_dirs()

# --- 1. CLI -------------------------------------------------------------------
# Options typed after the script name.
args <- commandArgs(trailingOnly = TRUE)
opt_list_only <- "--list"   %in% args
opt_no_pdf    <- "--no-pdf" %in% args
opt_no_png    <- "--no-png" %in% args
opt_force     <- "--force"  %in% args
opt_comparison <- {
  k <- which(args == "--comparison")
  if (length(k) == 1L && length(args) > k) args[k + 1L] else NULL
}
opt_pathway <- {
  k <- which(args == "--pathway")
  if (length(k) == 1L && length(args) > k) args[k + 1L] else NULL
}

# --- 2. Output paths ----------------------------------------------------------
pathway_dir   <- file.path(pipeline_root,
                           if (isTRUE(use_cellsweep)) "pathway_reports_ucell_cellsweep"
                                                     else "pathway_reports_ucell")
pathway_cache <- file.path(pathway_dir, "cache")
violin_dir            <- file.path(pathway_dir, "violins")
violin_dir_per_sample <- file.path(pathway_dir, "violins_per_sample")
stats_dir             <- file.path(pathway_dir, "stats")
for (d in c(pathway_dir, pathway_cache, violin_dir,
            violin_dir_per_sample, stats_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# --- 3. Display constants -----------------------------------------------------
# Cell-type order on the x-axis.
celltype_order <- c(
  "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating",
  "Treg",
  "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating",
  "MAIT", "dnT", "gdT",
  "NK", "NK Proliferating", "NK_CD56bright",
  "ILC",
  "B naive", "B intermediate", "B memory",
  "Plasmablast",
  "CD14 Mono", "CD16 Mono",
  "cDC1", "cDC2", "pDC", "ASDC",
  "Platelet", "Eryth", "HSPC",
  "Doublet"
)

# Violin groups from config.R violin_celltype_groups, matching the heatmap
# panels in script 11.
celltype_split <- lapply(violin_celltype_groups, function(p)
  list(suffix = p$suffix, label = p$label, types = names(p$groups)))

# X-axis order of the violin groups (celltype_group_levels() in utils.R).
celltype_group_order <- celltype_group_levels()

# Fixed colour per cell type (make_celltype_palette() in utils.R).
celltype_palette <- make_celltype_palette(celltype_order)

# --- 4. Pathway gene sets -----------------------------------------------------
log_msg("== Building pathway gene sets")
# Same gene sets as script 07.
# NFAT (all C2): every gene in any MSigDB C2 set with NFAT in its name.
c2 <- msigdbr(species = "Homo sapiens", category = "C2")
c2_nfat <- c2 %>% dplyr::filter(grepl("NFAT", gs_name))
nfat_all_genes <- c2_nfat %>% dplyr::pull(gene_symbol) %>% unique()

# Apoptosis (de Cevins et al. 2023), with outdated symbols updated below.
apoptosis_genes <- c(
  "ASM1", "BAD", "BAK1", "BAX", "BCL2", "BCL10", "BCL2L1", "BIK",
  "CARD19", "BIRC8", "CARD8", "CASP1", "CASP2", "CASP3", "CASP4",
  "CASP5", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10", "CASP12",
  "CFLAR", "CRADD", "DIABLO", "AIMP1", "ERO1A", "FADD", "FASL",
  "MEOX2", "GZMA", "GZMB", "HBXIP", "SYVN1", "LCN2", "LTBR",
  "MAPT", "MFN2", "MLKL", "NAIP1", "NAIP5", "NFIL3", "PMAIP1",
  "OPTN", "CDK5R1", "PCNA", "PDCD4", "PDCD8", "PIDD", "PRF1",
  "PTPN6", "PUMA", "RFC4", "SARP2", "SERPINB9", "BIRC5", "TGFB1",
  "TGFB2", "TNFAIP8", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C",
  "TNFRSF10D", "TNFRSF11B", "TRADD", "TNFSF10", "XIAP"
)
apoptosis_symbol_fixes <- c(
  "ASM1" = "SMPD1", "BIRC8" = "NLRP2", "FASL" = "FASLG",
  "HBXIP" = "LAMTOR5", "NAIP1" = "NAIP", "NAIP5" = "NAIP",
  "PDCD8" = "AIFM1", "PIDD" = "PIDD1", "PUMA" = "BBC3",
  "SARP2" = "SFRP1", "XIAP" = "BIRC4"
)
for (old in names(apoptosis_symbol_fixes)) {
  apoptosis_genes[apoptosis_genes == old] <- apoptosis_symbol_fixes[old]
}
apoptosis_genes <- unique(apoptosis_genes)

# Tacrolimus-sensitive NFAT targets.
tacro_nfat_targets <- c(
  "IL2", "IFNG", "TNF", "IL10",
  "IL2RA", "CTLA4", "ICOS", "PDCD1", "FASLG", "LAG3", "HAVCR2", "TIGIT",
  "NFATC1", "TBX21", "GATA3", "RORC", "BCL6", "BATF", "IRF4", "EGR2", "EGR3",
  "GZMB", "PRF1",
  "CXCR5", "SLAMF1"
)

# Hallmark IL2-STAT5, TNFa-NFkB and apoptosis gene sets.
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
il2_stat5_genes <- hallmark %>%
  dplyr::filter(gs_name == "HALLMARK_IL2_STAT5_SIGNALING") %>%
  dplyr::pull(gene_symbol) %>% unique()
tnfa_nfkb_genes <- hallmark %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  dplyr::pull(gene_symbol) %>% unique()
hallmark_apoptosis_genes <- hallmark %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  dplyr::pull(gene_symbol) %>% unique()
rm(c2, c2_nfat, hallmark)

log_msg("   NFAT (all C2):           ", length(nfat_all_genes), " genes")
log_msg("   Tacrolimus NFAT targets: ", length(tacro_nfat_targets), " genes")
log_msg("   Apoptosis (de Cevins):   ", length(apoptosis_genes), " genes")
log_msg("   Apoptosis (Hallmark):    ", length(hallmark_apoptosis_genes), " genes")
log_msg("   IL2-STAT5:               ", length(il2_stat5_genes), " genes")
log_msg("   TNFa-NFkB:               ", length(tnfa_nfkb_genes), " genes")

# --- 5. Pathway registry ------------------------------------------------------
# UCell appends "_UCell" to each name. Pathways with fewer than min_genes genes
# present are skipped.
pathways <- list(
  list(name = "NFAT_allC2",   label = "NFAT (all C2)",             genes = nfat_all_genes,           min_genes = 10),
  list(name = "TacroNFAT",    label = "Tacrolimus-sensitive NFAT", genes = tacro_nfat_targets,       min_genes = 5),
  list(name = "Apoptosis",    label = "Apoptosis (de Cevins)",     genes = apoptosis_genes,          min_genes = 5),
  list(name = "Apoptosis_HM", label = "Apoptosis (Hallmark)",      genes = hallmark_apoptosis_genes, min_genes = 10),
  list(name = "IL2_STAT5",    label = "IL2-STAT5 signaling",       genes = il2_stat5_genes,          min_genes = 10),
  list(name = "TNFa_NFkB",    label = "TNFa-NFkB signaling",       genes = tnfa_nfkb_genes,          min_genes = 10)
)

# --- 6. Comparisons -----------------------------------------------------------
# Same comparisons as script 07. single_donor = TRUE gives effect sizes only.
comparisons <- list(
  list(cohorts = c("GEM108_pre", "GEM108_post"),                  label = "GEM108 pre vs post",                       single_donor = TRUE),
  list(cohorts = c("HBD", "GEM108_pre", "GEM108_post"),           label = "HBD vs GEM108 pre vs post",                single_donor = TRUE),
  list(cohorts = c("HBD", "E42K_affected", "E42K_unaffected"),    label = "HBD vs E42K_affected vs E42K_unaffected",  single_donor = TRUE),
  list(cohorts = c("HBD", "E42K_affected", "T504S"),              label = "HBD vs E42K_affected vs T504S",            single_donor = FALSE),
  list(cohorts = c("HBD", "GEM108_pre"),                          label = "HBD vs GEM108_pre",                        single_donor = TRUE),
  list(cohorts = c("HBD", "GEM108_post"),                         label = "HBD vs GEM108_post",                       single_donor = TRUE),
  list(cohorts = c("HBD", "E42K_carriers"),                       label = "HBD vs E42K_carriers",                     single_donor = FALSE),
  list(cohorts = c("HBD", "E42K_affected"),                       label = "HBD vs E42K_affected",                     single_donor = FALSE),
  list(cohorts = c("HBD", "E42K_unaffected"),                     label = "HBD vs E42K_unaffected",                   single_donor = TRUE),
  list(cohorts = c("HBD", "T504S"),                               label = "HBD vs T504S",                             single_donor = FALSE)
)

if (opt_list_only) {
  log_msg("== Pathways:")
  for (pw in pathways) {
    log_msg("   - ", pw$name, " :: ", pw$label,
            " (", length(pw$genes), " genes)")
  }
  log_msg("== Comparisons:")
  for (cmp in comparisons) {
    log_msg("   - ", cmp$label, " :: ",
            paste(cmp$cohorts, collapse = ", "),
            if (cmp$single_donor) "  [single-donor mode]" else "  [multi-donor]")
  }
  quit(save = "no", status = 0L)
}

# --- 7. Filter by --comparison / --pathway ------------------------------------
if (!is.null(opt_comparison)) {
  hit <- vapply(comparisons, function(c) c$label == opt_comparison, logical(1))
  if (!any(hit)) stop("--comparison ", opt_comparison, " not found.")
  comparisons <- comparisons[hit]
}
if (!is.null(opt_pathway)) {
  hit <- vapply(pathways, function(p) p$name == opt_pathway, logical(1))
  if (!any(hit)) stop("--pathway ", opt_pathway, " not found.")
  pathways <- pathways[hit]
}

# --- 8. Load merged Seurat ----------------------------------------------------
# Path to the annotated object from script 04 (annotated_rds_path() in utils.R).
seu_path <- annotated_rds_path()
if (!file.exists(seu_path)) {
  stop("Merged Seurat not found at ", seu_path,
       " - has script 04 (and 04b, if use_cellsweep=TRUE) been run?")
}
log_msg("== Loading ", seu_path)
seu_combo <- readRDS(seu_path)
DefaultAssay(seu_combo) <- "RNA"

# Check the RNA assay holds the full gene set.
n_genes <- nrow(seu_combo[["RNA"]])
if (n_genes < 10000L) {
  stop("RNA assay has only ", n_genes, " genes - expected the full matrix.")
}
log_msg("   RNA assay: ", n_genes, " genes  x  ", ncol(seu_combo), " cells")

if (!"cohort_or_patient_tx" %in% colnames(seu_combo@meta.data)) {
  stop("`cohort_or_patient_tx` column missing - re-run script 04.")
}

# --- 9. Cohort subsets --------------------------------------------------------
log_msg("== Building cohort subsets")
# One Seurat object per cohort group. Each entry gives the rule for picking
# its cells from the full object (combo = all cells).
group_lookup <- list(
  combo           = NULL,
  # E42K_affected = PMAI0017/0018/0023 + GEM108 pre-treatment cells.
  E42K_affected   = list(filter = function(s)
    (s$cohort == "E42K_affected" & s$patient != "GEM108") |
    (s$patient == "GEM108" & s$cohort_or_patient_tx == "GEM108_pre")),
  E42K_unaffected = list(filter = function(s) s$cohort == "E42K_unaffected"),
  # GEM108 contributes pre-treatment cells only.
  E42K_carriers   = list(filter = function(s)
    (s$cohort %in% c("E42K_affected", "E42K_unaffected") & s$patient != "GEM108") |
    (s$patient == "GEM108" & s$cohort_or_patient_tx == "GEM108_pre")),
  GEM108_pre      = list(filter = function(s) s$cohort_or_patient_tx == "GEM108_pre"),
  GEM108_post     = list(filter = function(s) s$cohort_or_patient_tx == "GEM108_post"),
  HBD             = list(filter = function(s) s$cohort_or_patient_tx == "HBD"),
  T504S           = list(filter = function(s) s$cohort_or_patient_tx == "T504S")
)

seu_list <- list()
for (nm in names(group_lookup)) {
  if (is.null(group_lookup[[nm]])) {
    seu_list[[nm]] <- seu_combo
  } else {
    keep <- group_lookup[[nm]]$filter(seu_combo)
    keep[is.na(keep)] <- FALSE
    if (sum(keep) == 0L) {
      log_msg("   ", nm, ": no cells - skipping")
      next
    }
    seu_list[[nm]] <- seu_combo[, keep]
  }
  log_msg("   ", nm, ": ", ncol(seu_list[[nm]]), " cells")
}

# --- 10. UCell scoring --------------------------------------------------------
# Scores are rank-based and lie in [0, 1]. Run on one core.
# Add a UCell score for each pathway, using the genes present in the data.
# Scores are stored as metadata columns named <pathway>_UCell.
add_ucell <- function(seu, pathways, min_genes_lookup = NULL) {
  prev_assay <- DefaultAssay(seu)
  DefaultAssay(seu) <- "RNA"
  features <- list()
  for (pw in pathways) {
    present <- intersect(unique(pw$genes), rownames(seu))
    if (length(present) < pw$min_genes) {
      warning("Skipping ", pw$name, " - only ", length(present),
              " genes present (min ", pw$min_genes, ")")
      next
    }
    features[[pw$name]] <- present
  }
  if (length(features) == 0L) {
    DefaultAssay(seu) <- prev_assay
    return(seu)
  }
  seu <- UCell::AddModuleScore_UCell(
    seu,
    features = features,
    assay    = "RNA",
    slot     = "data",
    name     = "_UCell",
    ncores   = 1L
  )
  DefaultAssay(seu) <- prev_assay
  seu
}

# Scores for the per-cohort objects are cached. The cache is reused unless the
# input object is newer or --force is given.
cache_file <- file.path(pathway_cache, "seu_ucell_scored.rds")
score_cache_valid <- file.exists(cache_file) && !opt_force &&
                     file.info(cache_file)$mtime > file.info(seu_path)$mtime

if (score_cache_valid && is.null(opt_pathway)) {
  log_msg("== Loading cached UCell-scored Seurat (use --force to redo)")
  seu_list <- readRDS(cache_file)
} else {
  log_msg("== Scoring pathways with UCell across cohort subsets")
  for (nm in names(seu_list)) {
    log_msg("   ", nm)
    seu_list[[nm]] <- add_ucell(seu_list[[nm]], pathways)
  }
  if (is.null(opt_pathway)) {
    saveRDS(seu_list, cache_file)
    log_msg("   cached -> ", cache_file)
  }
}
gc(verbose = FALSE)

# --- 11. Helpers --------------------------------------------------------------
# Merge the per-cohort objects of one comparison and label each cell with its
# group in compare_group.
merge_for_comparison <- function(seu_list, names) {
  subs <- seu_list[names]
  for (nm in names) subs[[nm]]$compare_group <- nm
  m <- merge(subs[[1]], subs[-1])
  m <- tryCatch(
    SeuratObject::JoinLayers(m, assay = "RNA"),
    error = function(e) {
      message("JoinLayers(RNA) skipped: ", e$message)
      m
    }
  )
  m
}

# For each cell type (or violin group when grouped = TRUE) and each pair of
# groups: means, Cohen's d and a Wilcoxon test with BH correction. Same as in
# script 07. mode = "effect_size" gives Cohen's d only.
run_celltype_tests_pairwise <- function(seu, feature, group_col = "compare_group",
                                        celltype_col = "predicted.celltype.l2",
                                        min_cells = 20, min_d = 0.2,
                                        mode = c("pvalue", "effect_size"),
                                        grouped = FALSE) {
  mode <- match.arg(mode)
  md <- seu@meta.data %>%
    tibble::as_tibble() %>%
    dplyr::select(celltype = dplyr::all_of(celltype_col),
                  group    = dplyr::all_of(group_col),
                  score    = dplyr::all_of(feature)) %>%
    dplyr::filter(!is.na(score), !is.na(group))
  # Map L2 labels to violin groups. expand_celltype_groups() (utils.R) repeats
  # a cell once for each group it belongs to.
  if (isTRUE(grouped)) {
    md <- expand_celltype_groups(md, celltype_col = "celltype")
    md$celltype <- md$celltype_group
    md$celltype_group <- NULL
    md <- tibble::as_tibble(md)
  }
  if (!is.factor(md$group)) md$group <- factor(md$group)
  groups <- levels(md$group)
  if (length(groups) < 2L) {
    return(tibble::tibble(celltype = character(), grp1 = character(),
                          grp2 = character(), n_grp1 = integer(),
                          n_grp2 = integer(), mean_grp1 = double(),
                          mean_grp2 = double(), cohens_d = double(),
                          p_value = double(), p_adj = double(),
                          sig = character(), comparison = character()))
  }
  pairs <- combn(groups, 2, simplify = FALSE)
  results <- purrr::map_dfr(pairs, function(pair) {
    md_pair <- md %>% dplyr::filter(group %in% pair) %>%
      dplyr::mutate(group = droplevels(group))
    md_pair %>%
      dplyr::group_by(celltype) %>%
      dplyr::filter(dplyr::n_distinct(group) == 2L,
                    all(table(group) >= min_cells)) %>%
      dplyr::summarise(
        grp1      = pair[1],
        grp2      = pair[2],
        n_grp1    = sum(group == pair[1]),
        n_grp2    = sum(group == pair[2]),
        mean_grp1 = mean(score[group == pair[1]]),
        mean_grp2 = mean(score[group == pair[2]]),
        cohens_d  = (mean_grp1 - mean_grp2) /
          sqrt((var(score[group == pair[1]]) +
                var(score[group == pair[2]])) / 2),
        p_value   = tryCatch(wilcox.test(score ~ group)$p.value,
                             error = function(e) NA_real_),
        .groups   = "drop"
      )
  })
  if (mode == "pvalue") {
    results <- results %>%
      dplyr::mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        sig = dplyr::case_when(
          p_adj < 0.001 & abs(cohens_d) > min_d ~ "***",
          p_adj < 0.01  & abs(cohens_d) > min_d ~ "**",
          p_adj < 0.05  & abs(cohens_d) > min_d ~ "*",
          TRUE                                  ~ "ns"
        ),
        comparison = paste(grp1, "v", grp2)
      )
  } else {
    results <- results %>%
      dplyr::mutate(
        p_adj = NA_real_,
        sig = dplyr::case_when(
          abs(cohens_d) > 0.8 ~ "d>0.8",
          abs(cohens_d) > 0.5 ~ "d>0.5",
          abs(cohens_d) > 0.2 ~ "d>0.2",
          TRUE                ~ "ns"
        ),
        comparison = paste(grp1, "v", grp2)
      )
  }
  results %>% dplyr::arrange(dplyr::desc(abs(cohens_d)))
}

# Violin y-axis range from the data (UCell scores lie in [0, 1]).
violin_ylim_ucell <- function(scores) {
  m <- suppressWarnings(max(scores, na.rm = TRUE))
  if (!is.finite(m)) return(c(0, 0.4))
  if      (m > 0.6) c(0, 1.0)
  else if (m > 0.4) c(0, 0.8)
  else if (m > 0.2) c(0, 0.6)
  else              c(0, 0.4)
}

# Violin plot of one score by cell type, split by group. stats_df is not drawn.
plot_comparison_violin <- function(seu, feature, stats_df = NULL,
                                   ct_order = NULL, restrict = NULL,
                                   celltype_col = "predicted.celltype.l2",
                                   group_col = "compare_group",
                                   grouped = FALSE) {
  md <- seu@meta.data %>%
    tibble::as_tibble() %>%
    dplyr::select(celltype = dplyr::all_of(celltype_col),
                  group    = dplyr::all_of(group_col),
                  score    = dplyr::all_of(feature)) %>%
    dplyr::filter(!is.na(score), !is.na(group))
  # restrict names the groups of one panel.
  if (isTRUE(grouped)) {
    md <- expand_celltype_groups(md, celltype_col = "celltype",
                                 groups = restrict)
    md$celltype <- md$celltype_group
    md$celltype_group <- NULL
    md <- tibble::as_tibble(md)
    restrict <- if (is.null(restrict)) unique(md$celltype) else restrict
  } else if (!is.null(restrict)) {
    md <- md %>% dplyr::filter(celltype %in% restrict)
  }
  if (!is.factor(md$group)) md$group <- factor(md$group)
  if (is.null(ct_order)) {
    ct_order <- md %>% dplyr::group_by(celltype) %>%
      dplyr::summarise(med = median(score)) %>%
      dplyr::arrange(dplyr::desc(med)) %>% dplyr::pull(celltype)
  } else {
    ct_order <- ct_order[ct_order %in% unique(md$celltype)]
    if (is.null(restrict)) {
      ct_order <- c(ct_order, setdiff(unique(md$celltype), ct_order))
    }
  }
  md$celltype <- factor(md$celltype, levels = ct_order)
  dodge_width <- 0.8
  ggplot(md, aes(x = celltype, y = score, fill = group)) +
    geom_violin(scale = "width", position = position_dodge(dodge_width),
                trim = TRUE, linewidth = 0.3) +
    coord_cartesian(ylim = violin_ylim_ucell(md$score)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(x = NULL, y = paste0(feature, "  (UCell)"), fill = NULL)
}

# Run the stats and draw the T-cell and other violin panels for one
# comparison and one score.
compare_cohorts <- function(merged_cache, cohort_names, feature,
                            ct_order = NULL, title = NULL, min_cells = 20,
                            single_donor = FALSE,
                            group_col = "compare_group") {
  key <- paste(cohort_names, collapse = "_")
  seu_merged <- merged_cache[[key]]
  mode <- if (single_donor) "effect_size" else "pvalue"
  # Two stats tables: one per violin group, one per Azimuth L2 label.
  stats <- run_celltype_tests_pairwise(seu_merged, feature,
                                       min_cells = min_cells, mode = mode,
                                       group_col = group_col, grouped = TRUE)
  stats_l2 <- run_celltype_tests_pairwise(seu_merged, feature,
                                          min_cells = min_cells, mode = mode,
                                          group_col = group_col, grouped = FALSE)
  base_title <- title %||% paste(cohort_names, collapse = " vs ")
  plots <- lapply(celltype_split, function(sp) {
    plot_comparison_violin(seu_merged, feature, stats_df = stats,
                           ct_order = sp$types, restrict = sp$types,
                           group_col = group_col, grouped = TRUE) +
      ggtitle(paste0(base_title, " — ", sp$label))
  })
  names(plots) <- vapply(celltype_split, function(sp) sp$suffix, character(1))
  invisible(list(seu = seu_merged, stats = stats, stats_l2 = stats_l2,
                 plots = plots))
}

# Split GEM108_pre / GEM108_post into their individual samples.
expand_per_sample_levels <- function(cohort_levels) {
  out <- character(0)
  for (lvl in cohort_levels) {
    if (identical(lvl, "GEM108_pre")) {
      out <- c(out, "GEM108_pre_s1", "GEM108_pre_s2")
    } else if (identical(lvl, "GEM108_post")) {
      out <- c(out, "GEM108_post_s1", "GEM108_post_s2")
    } else {
      out <- c(out, lvl)
    }
  }
  out
}

# Per-sample group: GEM108 Aug cells get their sample id, others keep cohort.
assign_compare_group_per_sample <- function(seu, cohort_levels) {
  md <- seu@meta.data
  cg <- as.character(md$compare_group)
  hto_short <- if ("HTO_best" %in% colnames(md))
                 sub("[-_ ].*$", "", as.character(md$HTO_best))
               else rep(NA_character_, nrow(md))
  is_gem <- (md$patient == "GEM108") & (md$batch == "Aug")
  sid <- rep(NA_character_, nrow(md))
  sid[is_gem & md$capture == "Maurice2" & hto_short == "Hashtag2"] <- "GEM108_pre_s2"
  sid[is_gem & md$capture == "Maurice3" & hto_short == "Hashtag1"] <- "GEM108_pre_s1"
  sid[is_gem & md$capture == "Maurice2" & hto_short == "Hashtag1"] <- "GEM108_post_s2"
  sid[is_gem & md$capture == "Maurice4" & hto_short == "Hashtag1"] <- "GEM108_post_s1"
  cg_per <- ifelse(!is.na(sid), sid, cg)
  per_levels <- expand_per_sample_levels(cohort_levels)
  per_levels <- per_levels[per_levels %in% unique(cg_per)]
  seu$compare_group_per_sample <- factor(cg_per, levels = per_levels)
  seu
}

# --- 12. Build merged-cache per comparison + UCell-score them -----------------
log_msg("== Building merged objects for each comparison")
# Merge the cohorts of each comparison into one object.
merged_cache <- list()
for (cmp in comparisons) {
  key <- paste(cmp$cohorts, collapse = "_")
  if (!all(cmp$cohorts %in% names(seu_list))) {
    log_msg("   skipping ", cmp$label, " - missing cohort(s): ",
            paste(setdiff(cmp$cohorts, names(seu_list)), collapse = ", "))
    next
  }
  log_msg("   merging: ", key)
  seu_merged <- merge_for_comparison(seu_list, cmp$cohorts)
  seu_merged$compare_group <- factor(seu_merged$compare_group,
                                     levels = cmp$cohorts)
  merged_cache[[key]] <- seu_merged
}
rm(seu_merged); gc(verbose = FALSE)

# Score the merged objects.
log_msg("== Scoring pathways on merged comparison objects (UCell)")
for (key in names(merged_cache)) {
  log_msg("   ", key)
  merged_cache[[key]] <- add_ucell(merged_cache[[key]], pathways)
}

# Add compare_group_per_sample, which splits GEM108 into its samples.
log_msg("== Assigning per-sample compare_group on merged cache")
for (cmp in comparisons) {
  key <- paste(cmp$cohorts, collapse = "_")
  if (!key %in% names(merged_cache)) next
  merged_cache[[key]] <- assign_compare_group_per_sample(
    merged_cache[[key]], cmp$cohorts
  )
}

# --- 13. Pooled comparison violins --------------------------------------------
log_msg("== Generating comparison violins + stats (pooled)")
# One set of violins and stats per comparison and pathway.
all_results <- list()
for (cmp in comparisons) {
  key0 <- paste(cmp$cohorts, collapse = "_")
  if (!key0 %in% names(merged_cache)) next
  for (pw in pathways) {
    feat  <- paste0(pw$name, "_UCell")     # UCell column suffix
    key   <- paste0(cmp$label, " | ", pw$label)
    title <- paste0(pw$label, " — ", cmp$label)
    log_msg("   violin: ", key)
    res <- tryCatch(
      compare_cohorts(merged_cache, cmp$cohorts, feat,
                      ct_order = celltype_group_order, title = title,
                      single_donor = cmp$single_donor),
      error = function(e) {
        log_msg("      ERROR: ", e$message)
        NULL
      }
    )
    all_results[[key]] <- res
    if (!is.null(res)) {
      # Stats per violin group.
      readr::write_tsv(
        res$stats,
        file.path(stats_dir,
                  paste0(gsub("[^A-Za-z0-9]+", "_", cmp$label),
                         "__", pw$name, ".tsv"))
      )
      # Stats per Azimuth L2 label.
      readr::write_tsv(
        res$stats_l2,
        file.path(stats_dir,
                  paste0(gsub("[^A-Za-z0-9]+", "_", cmp$label),
                         "__", pw$name, "__by_l2.tsv"))
      )
      if (!opt_no_png) {
        for (sfx in names(res$plots)) {
          out_file <- file.path(violin_dir,
            paste0(gsub("[^A-Za-z0-9]+", "_", cmp$label),
                   "__", pw$name, "__", sfx, ".png"))
          tryCatch(
            ggsave(out_file, res$plots[[sfx]], width = 12, height = 6,
                   dpi = 200, bg = "white"),
            error = function(e) log_msg("      ggsave failed: ", e$message)
          )
        }
      }
    }
  }
}

# --- 14. Per-sample violins for GEM108 comparisons ----------------------------
log_msg("== Generating per-sample comparison violins + stats")
# Only for comparisons that include GEM108_pre or GEM108_post.
involves_gem108 <- function(cohort_names) {
  any(cohort_names %in% c("GEM108_pre", "GEM108_post"))
}

all_results_per_sample <- list()
for (cmp in comparisons) {
  if (!involves_gem108(cmp$cohorts)) next
  key0 <- paste(cmp$cohorts, collapse = "_")
  if (!key0 %in% names(merged_cache)) next
  for (pw in pathways) {
    feat  <- paste0(pw$name, "_UCell")
    key   <- paste0(cmp$label, " | ", pw$label)
    title <- paste0(pw$label, " — ", cmp$label, " (per sample)")
    log_msg("   per-sample violin: ", key)
    res <- tryCatch(
      compare_cohorts(merged_cache, cmp$cohorts, feat,
                      ct_order = celltype_group_order, title = title,
                      single_donor = cmp$single_donor,
                      group_col = "compare_group_per_sample"),
      error = function(e) {
        log_msg("      ERROR: ", e$message)
        NULL
      }
    )
    all_results_per_sample[[key]] <- res
    if (!is.null(res)) {
      readr::write_tsv(
        res$stats,
        file.path(stats_dir,
                  paste0(gsub("[^A-Za-z0-9]+", "_", cmp$label),
                         "__", pw$name, "__per_sample.tsv"))
      )
      readr::write_tsv(
        res$stats_l2,
        file.path(stats_dir,
                  paste0(gsub("[^A-Za-z0-9]+", "_", cmp$label),
                         "__", pw$name, "__per_sample__by_l2.tsv"))
      )
      if (!opt_no_png) {
        for (sfx in names(res$plots)) {
          out_file <- file.path(violin_dir_per_sample,
            paste0(gsub("[^A-Za-z0-9]+", "_", cmp$label),
                   "__", pw$name, "__", sfx, ".png"))
          tryCatch(
            ggsave(out_file, res$plots[[sfx]], width = 12, height = 6,
                   dpi = 200, bg = "white"),
            error = function(e) log_msg("      ggsave failed: ", e$message)
          )
        }
      }
    }
  }
}

# --- 15. PDF report -----------------------------------------------------------
# One violin panel ("Tcells" or "other"), or a placeholder if missing.
safe_plot <- function(key, suffix) {
  obj <- all_results[[key]]
  if (!is.null(obj) && !is.null(obj$plots) && !is.null(obj$plots[[suffix]]))
    return(obj$plots[[suffix]])
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                      label = paste("Not available:\n", key, "\n", suffix)) +
    theme_void()
}

# Page number in the bottom-right corner.
make_page_numberer <- function() {
  n <- 0L
  function() {
    n <<- n + 1L
    grid::upViewport(0)
    grid::grid.text(sprintf("Page %d", n), x = 0.99, y = 0.01,
                    just = c("right", "bottom"),
                    gp = grid::gpar(fontsize = 8, col = "grey50"))
  }
}
page_number <- make_page_numberer()

# A PDF page of plain text.
text_page <- function(title, body_lines) {
  grid::grid.newpage()
  grid::grid.text(title, x = 0.05, y = 0.95, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 16, fontface = "bold"))
  grid::grid.text(paste(body_lines, collapse = "\n"),
                  x = 0.05, y = 0.88, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 9, fontfamily = "sans",
                                  lineheight = 1.3))
  page_number()
}

if (!opt_no_pdf) {
  out_pdf <- file.path(pathway_dir, "pathway_analysis_report_ucell.pdf")
  log_msg("== Writing PDF report -> ", out_pdf)
  grDevices::pdf(out_pdf, width = 16, height = 10)

  grid::grid.newpage()
  grid::grid.text("Pathway Analysis Report — UCell",
                  x = 0.5, y = 0.97, just = c("centre", "top"),
                  gp = grid::gpar(fontsize = 16, fontface = "bold"))
  grid::grid.text(paste(c(
    "OVERVIEW",
    "Sibling of 07_pathway_module_scores.R that uses UCell instead of",
    "Seurat::AddModuleScore. Same six gene sets, same comparisons, same",
    "T-cell / other panel split. Different scoring method, different",
    "y-axis range, different sensitivity profile.",
    "",
    "SCORING",
    "UCell::AddModuleScore_UCell computes a Mann-Whitney U statistic per",
    "cell over the pathway genes' ranks within that cell's expression",
    "vector. Score range is [0, 1]; 0.5 indicates the pathway genes are",
    "evenly distributed through the cell's rank distribution. UCell is",
    "rank-based, so it is robust to dropouts and depth differences and",
    "comparable across datasets without re-normalisation. Andreatta &",
    "Carmona 2021 Comput Struct Biotechnol J.",
    "",
    "VIOLIN Y-AXIS",
    "Per-panel y-limit chosen from the panel's max score:",
    "  max > 0.6  -> [0, 1.0]",
    "  max > 0.4  -> [0, 0.8]",
    "  max > 0.2  -> [0, 0.6]",
    "  max <= 0.2 -> [0, 0.4]",
    "",
    "STATISTICAL TESTING",
    "Pairwise Wilcoxon rank-sum (BH-adjusted) for multi-donor comparisons,",
    "Cohen's d only for single-donor (GEM108) comparisons. Stats are in",
    "stats/ TSVs; brackets are not drawn on violins."
  ), collapse = "\n"),
  x = 0.05, y = 0.92, just = c("left", "top"),
  gp = grid::gpar(fontsize = 9, fontfamily = "mono", lineheight = 1.3))
  page_number()

  for (cmp in comparisons) {
    text_page(paste("Comparison:", cmp$label), c(
      "",
      paste("Cohorts:", paste(cmp$cohorts, collapse = ", ")),
      "",
      if (cmp$single_donor) {
        "NOTE: single-donor comparison; Cohen's d only in stats TSVs."
      } else {
        "Multi-donor comparison; BH-adj p + Cohen's d in stats TSVs."
      }
    ))
    for (pw in pathways) {
      key <- paste0(cmp$label, " | ", pw$label)
      p_t <- safe_plot(key, "Tcells")
      p_o <- safe_plot(key, "other")
      tryCatch({
        print(p_t / p_o +
                patchwork::plot_annotation(
                  title = paste(pw$label, "—", cmp$label, "(UCell)"),
                  theme = theme(plot.title = element_text(size = 14,
                                                          face = "bold"))))
        page_number()
      }, error = function(e) {
        log_msg("Violin page error: ", pw$label, " ", e$message)
        print(p_t + ggtitle(paste(pw$label, "—", cmp$label, "(UCell)")))
        page_number()
      })
    }
  }
  grDevices::dev.off()
  log_msg("== PDF report complete: ", out_pdf)
}

log_msg("== 07a_pathway_ucell_scores.R done.")
