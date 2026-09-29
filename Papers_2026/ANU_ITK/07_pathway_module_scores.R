# ------------------------------------------------------------------------------
# 07_pathway_module_scores.R
# Pathway activity per cell with Seurat::AddModuleScore, for six gene sets:
# NFAT (all C2), tacrolimus-sensitive NFAT targets, apoptosis (de Cevins),
# Hallmark apoptosis, Hallmark IL2-STAT5 and Hallmark TNFa-NFkB.
#
# For each comparison of cohorts (section 6) the script draws violins of the
# scores per cell-type group and tests each group with a Wilcoxon test and
# Cohen's d. Comparisons with a single donor on one side (GEM108) report
# Cohen's d only. It also draws UMAPs of each score and writes a PDF report.
#
# Run:  Rscript new_scripts/07_pathway_module_scores.R [options]
#   --list               print pathways and comparisons and exit
#   --comparison NAME    one comparison only (label as in section 6)
#   --pathway NAME       one pathway only (e.g. TacroNFAT)
#   --no-pdf             skip the PDF report
#   --no-png             skip the individual PNGs
#   --force              ignore the module-score cache
#
# Input:  pipeline/merged/seurat_merged_harmony_azimuth.rds   (from script 04)
# Output: pipeline/pathway_reports/  (violins/, umaps/, stats/, PDF report)
# ------------------------------------------------------------------------------

# Load a local GLPK build before Seurat (needed on Gadi). Edit or remove.
if (file.exists("/path/to/libglpk.so.40")) {
  dyn.load("/path/to/libglpk.so.40")
}

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(msigdbr)
  library(ggsignif)
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
opt_list_only   <- "--list"   %in% args
opt_no_pdf      <- "--no-pdf" %in% args
opt_no_png      <- "--no-png" %in% args
opt_force       <- "--force"  %in% args
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
                           if (isTRUE(use_cellsweep)) "pathway_reports_cellsweep"
                                                     else "pathway_reports")
pathway_cache <- file.path(pathway_dir, "cache")
violin_dir            <- file.path(pathway_dir, "violins")
violin_dir_per_sample <- file.path(pathway_dir, "violins_per_sample")
umap_out_dir          <- file.path(pathway_dir, "umaps")
stats_dir             <- file.path(pathway_dir, "stats")
for (d in c(pathway_dir, pathway_cache, violin_dir, violin_dir_per_sample,
            umap_out_dir, stats_dir)) {
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

# Each violin is drawn as a T-cell panel and an "other" panel. The groups come
# from config.R violin_celltype_groups and match the heatmap panels in script 11.
# T-cell groups overlap (Total CD4 = Naive CD4 + Activated CD4), so a cell can
# appear in more than one violin.
celltype_split <- lapply(violin_celltype_groups, function(p)
  list(suffix = p$suffix, label = p$label, types = names(p$groups)))

# X-axis order of the violin groups (celltype_group_levels() in utils.R).
celltype_group_order <- celltype_group_levels()

# Cell types for the per-gene tacrolimus table.
t_cell_types <- c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL",
                  "CD8 Naive", "CD8 TCM", "CD8 TEM",
                  "Treg", "MAIT", "gdT", "NK")

# Fixed colour per cell type (make_celltype_palette() in utils.R).
celltype_palette <- make_celltype_palette(celltype_order)

# --- 4. Pathway gene sets -----------------------------------------------------
log_msg("== Building pathway gene sets")

# 4a. NFAT (all C2): every gene in any MSigDB C2 set with NFAT in its name.
c2 <- msigdbr(species = "Homo sapiens", category = "C2")
c2_nfat <- c2 %>% dplyr::filter(grepl("NFAT", gs_name))
nfat_all_genes <- c2_nfat %>% dplyr::pull(gene_symbol) %>% unique()

# 4b. Apoptosis (de Cevins et al. 2023, Cell Rep Med). HGNC fixes applied.
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
# Replace outdated gene symbols with current HGNC symbols.
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

# 4c. Tacrolimus-sensitive NFAT targets (literature-curated).
tacro_nfat_targets <- c(
  "IL2", "IFNG", "TNF", "IL10",
  "IL2RA", "CTLA4", "ICOS", "PDCD1", "FASLG", "LAG3", "HAVCR2", "TIGIT",
  "NFATC1", "TBX21", "GATA3", "RORC", "BCL6", "BATF", "IRF4", "EGR2", "EGR3",
  "GZMB", "PRF1",
  "CXCR5", "SLAMF1"
)

# 4d. Hallmark IL2-STAT5, TNFa-NFkB and apoptosis gene sets.
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
# name = score column name, ctrl = control genes per bin for AddModuleScore,
# min_genes = skip the pathway if fewer of its genes are in the data.
pathways <- list(
  list(name = "NFAT_allC2",   label = "NFAT (all C2)",              genes = nfat_all_genes,           ctrl = 100, min_genes = 10),
  list(name = "TacroNFAT",    label = "Tacrolimus-sensitive NFAT",  genes = tacro_nfat_targets,       ctrl = 100, min_genes = 5),
  list(name = "Apoptosis",    label = "Apoptosis (de Cevins)",      genes = apoptosis_genes,          ctrl = 100, min_genes = 5),
  list(name = "Apoptosis_HM", label = "Apoptosis (Hallmark)",       genes = hallmark_apoptosis_genes, ctrl = 100, min_genes = 10),
  list(name = "IL2_STAT5",    label = "IL2-STAT5 signaling",        genes = il2_stat5_genes,          ctrl = 100, min_genes = 10),
  list(name = "TNFa_NFkB",    label = "TNFa-NFkB signaling",        genes = tnfa_nfkb_genes,          ctrl = 100, min_genes = 10)
)

# Pathway pairs shown together in the PDF.
pathway_groups <- list(
  list(label = "NFAT pathway activity",
       keys  = c("NFAT (all C2)", "Tacrolimus-sensitive NFAT")),
  list(label = "Apoptosis",
       keys  = c("Apoptosis (de Cevins)", "Apoptosis (Hallmark)")),
  list(label = "IL2-STAT5 signaling",
       keys  = c("IL2-STAT5 signaling")),
  list(label = "TNFa-NFkB signaling",
       keys  = c("TNFa-NFkB signaling"))
)

# --- 6. Comparisons -----------------------------------------------------------
# Cohort labels are values of cohort_or_patient_tx (set in script 04).
# single_donor = TRUE gives effect sizes only (no p-values).
comparisons <- list(
  # Three-way and paired comparisons
  list(cohorts = c("GEM108_pre", "GEM108_post"),                  label = "GEM108 pre vs post",                       single_donor = TRUE),
  list(cohorts = c("HBD", "GEM108_pre", "GEM108_post"),           label = "HBD vs GEM108 pre vs post",                single_donor = TRUE),
  list(cohorts = c("HBD", "E42K_affected", "E42K_unaffected"),    label = "HBD vs E42K_affected vs E42K_unaffected",  single_donor = TRUE),
  list(cohorts = c("HBD", "E42K_affected", "T504S"),              label = "HBD vs E42K_affected vs T504S",            single_donor = FALSE),
  # Pairwise HBD vs X comparisons
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
  stop("RNA assay has only ", n_genes, " genes - expected the full ~",
       "20000+ matrix. Check that script 03 used NormalizeData() (not ",
       "SCTransform) and script 04 didn't subset variable features ",
       "into the RNA assay.")
}
log_msg("   RNA assay: ", n_genes, " genes  x  ", ncol(seu_combo), " cells")

# Verify cohort_or_patient_tx is present.
if (!"cohort_or_patient_tx" %in% colnames(seu_combo@meta.data)) {
  stop("`cohort_or_patient_tx` column missing - re-run script 04 ",
       "(it is built in section 5b of that script).")
}

# --- 9. Cohort subsets --------------------------------------------------------
log_msg("== Building cohort subsets")
# One Seurat object per cohort group. Each entry gives the rule for picking
# its cells from the full object (combo = all cells).
group_lookup <- list(
  combo           = NULL,                                            # full
  # E42K_affected = PMAI0017/0018/0023 + GEM108 pre-treatment cells.
  E42K_affected   = list(filter = function(s)
    (s$cohort == "E42K_affected" & s$patient != "GEM108") |
    (s$patient == "GEM108" & s$cohort_or_patient_tx == "GEM108_pre")),
  E42K_unaffected = list(filter = function(s) s$cohort == "E42K_unaffected"),
  # E42K_carriers = E42K_affected + E42K_unaffected, with GEM108 pre-treatment
  # cells only.
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

# --- 10. AddModuleScore -------------------------------------------------------
# Add one module score to a Seurat object. The score is the mean expression
# of the pathway genes minus the mean of randomly chosen control genes with
# similar expression. Stored as a metadata column named <module_name>1.
add_module <- function(seu, genes, module_name,
                       ctrl = 100, min_genes = 10, verbose = TRUE) {
  prev_assay <- DefaultAssay(seu)
  DefaultAssay(seu) <- "RNA"
  genes_present <- intersect(unique(genes), rownames(seu))
  if (verbose) {
    log_msg("      ", module_name, ": ", length(genes_present), "/",
            length(unique(genes)), " genes present")
  }
  if (length(genes_present) < min_genes) {
    warning("Skipping ", module_name, " - only ", length(genes_present),
            " genes (min ", min_genes, ")")
    DefaultAssay(seu) <- prev_assay
    return(seu)
  }
  seu <- AddModuleScore(seu, features = list(genes_present),
                        name = module_name, ctrl = ctrl)
  DefaultAssay(seu) <- prev_assay
  seu
}

# Scores for the per-cohort objects are cached. The cache is reused unless the
# input object is newer or --force is given.
cache_file <- file.path(pathway_cache, "seu_module_scored.rds")
score_cache_valid <- file.exists(cache_file) && !opt_force &&
                     file.info(cache_file)$mtime > file.info(seu_path)$mtime

if (score_cache_valid && is.null(opt_pathway)) {
  log_msg("== Loading cached module-scored Seurat (use --force to redo)")
  seu_list <- readRDS(cache_file)
} else {
  log_msg("== Scoring pathways across cohort subsets")
  for (pw in pathways) {
    # Seed per pathway. AddModuleScore samples control genes at random, so this
    # keeps each pathway's score the same whether it runs alone or with others.
    set.seed(gsea_seed + sum(utf8ToInt(pw$name)))
    log_msg("   pathway: ", pw$name)
    seu_list <- purrr::map(seu_list, add_module,
                           genes = pw$genes, module_name = pw$name,
                           ctrl = pw$ctrl, min_genes = pw$min_genes)
  }
  if (is.null(opt_pathway)) {
    saveRDS(seu_list, cache_file)
    log_msg("   cached -> ", cache_file)
  }
}
gc(verbose = FALSE)

# --- 11. Reusable helpers -----------------------------------------------------
# Merge the per-cohort objects of one comparison and label each cell with its
# group in compare_group.
merge_for_comparison <- function(seu_list, names) {
  subs <- seu_list[names]
  for (nm in names) subs[[nm]]$compare_group <- nm
  m <- merge(subs[[1]], subs[-1])
  # Join per-sample layers so AddModuleScore can read the merged assay.
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
# groups, compare the score: means, Cohen's d and a Wilcoxon test with BH
# correction. Cell types with fewer than min_cells cells in either group are
# skipped. mode = "effect_size" gives Cohen's d only (single-donor comparisons).
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
  # Every pair of groups.
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

  # Significance stars need BH p < cutoff and |d| > min_d.
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

# Violin y-axis range from the maximum score.
violin_ylim <- function(scores) {
  m <- suppressWarnings(max(scores, na.rm = TRUE))
  if (!is.finite(m)) return(c(-0.2, 0.4))
  if (m > 0.6)      c(-0.2, 0.8)
  else if (m > 0.4) c(-0.2, 0.6)
  else              c(-0.2, 0.4)
}

# Violin plot of one score by cell type, split by group. stats_df is not
# drawn. restrict limits the plot to one panel's cell types or groups.
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

  # X-axis order: as given, or by median score.
  if (is.null(ct_order)) {
    ct_order <- md %>% dplyr::group_by(celltype) %>%
      dplyr::summarise(med = median(score)) %>%
      dplyr::arrange(dplyr::desc(med)) %>% dplyr::pull(celltype)
  } else {
    # Keep only requested cell types that are actually present, in order.
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
    coord_cartesian(ylim = violin_ylim(md$score)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(x = NULL, y = feature, fill = NULL)
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
  # Build one violin panel per cell-type split (T cells / other).
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

  # Sample id for GEM108 Aug cells from capture and HTO.
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

# UMAP coloured by cell type next to a UMAP coloured by the score.
plot_umap_pair <- function(seu, feature, reduction = "umap",
                           palette = celltype_palette, title = NULL) {
  p_dim <- DimPlot(seu, reduction = reduction,
                   label = TRUE, repel = TRUE, label.size = 3, pt.size = 1,
                   group.by = "predicted.celltype.l2", cols = palette) +
    NoLegend() +
    ggtitle(title %||% "Cell types")
  p_feat <- FeaturePlot(seu, features = feature, reduction = reduction,
                        cols = c("#FFFFCC", "#CC0000")) +
    ggtitle(feature)
  p_dim + p_feat
}

# --- 12. Build merged-cache for each comparison + score it --------------------
log_msg("== Building merged objects for each comparison")
# Merge the cohorts of each comparison into one object.
merged_cache <- list()
for (cmp in comparisons) {
  key <- paste(cmp$cohorts, collapse = "_")
  if (!all(cmp$cohorts %in% names(seu_list))) {
    log_msg("   skipping ", cmp$label,
            " - missing cohort(s): ",
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

# Score each merged object again. AddModuleScore picks control genes from the
# cells in the object, so the scores on the merged objects are the ones used
# for the violins and stats.
log_msg("== Scoring pathways on merged comparison objects")
for (key in names(merged_cache)) {
  log_msg("   ", key)
  for (pw in pathways) {
    # Seed per comparison and pathway (see above).
    set.seed(gsea_seed + sum(utf8ToInt(paste0(key, "|", pw$name))))
    merged_cache[[key]] <- add_module(
      merged_cache[[key]],
      genes = pw$genes, module_name = pw$name,
      ctrl = pw$ctrl, min_genes = pw$min_genes,
      verbose = FALSE
    )
  }
}

# Add compare_group_per_sample, which splits GEM108 into its individual
# samples, for section 14b.
log_msg("== Assigning per-sample compare_group on merged cache")
for (cmp in comparisons) {
  key <- paste(cmp$cohorts, collapse = "_")
  if (!key %in% names(merged_cache)) next
  merged_cache[[key]] <- assign_compare_group_per_sample(
    merged_cache[[key]], cmp$cohorts
  )
}

# --- 13. UMAPs per cohort and pathway -----------------------------------------
log_msg("== Generating UMAP plots")
# Score UMAPs for each cohort group and pathway.
umap_cohorts <- intersect(c("combo", "GEM108_pre", "GEM108_post",
                            "E42K_affected", "E42K_unaffected",
                            "E42K_carriers", "HBD", "T504S"),
                          names(seu_list))
umap_plots <- list()
for (cohort in umap_cohorts) {
  for (pw in pathways) {
    feat <- paste0(pw$name, "1")  # AddModuleScore appends "1"
    key  <- paste0(cohort, " | ", pw$label)
    log_msg("   UMAP: ", key)
    umap_plots[[key]] <- tryCatch(
      plot_umap_pair(seu_list[[cohort]], feat, title = cohort),
      error = function(e) {
        log_msg("      ERROR: ", e$message)
        NULL
      }
    )
    if (!opt_no_png && !is.null(umap_plots[[key]])) {
      out_file <- file.path(umap_out_dir,
        paste0(cohort, "__", pw$name, ".png"))
      tryCatch(
        ggsave(out_file, umap_plots[[key]], width = 10, height = 5,
               dpi = 200, bg = "white"),
        error = function(e) log_msg("      ggsave failed: ", e$message)
      )
    }
  }
}

# --- 14. Generate comparison violins ------------------------------------------
log_msg("== Generating comparison violins + stats")
# One set of violins and stats per comparison and pathway.
all_results <- list()
for (cmp in comparisons) {
  key0 <- paste(cmp$cohorts, collapse = "_")
  if (!key0 %in% names(merged_cache)) next
  for (pw in pathways) {
    feat  <- paste0(pw$name, "1")
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

# --- 14b. Per-sample violins (GEM108 samples shown separately) ----------------
# Only for comparisons that include GEM108_pre or GEM108_post.
log_msg("== Generating per-sample comparison violins + stats")
involves_gem108 <- function(cohort_names) {
  any(cohort_names %in% c("GEM108_pre", "GEM108_post"))
}

all_results_per_sample <- list()
for (cmp in comparisons) {
  if (!involves_gem108(cmp$cohorts)) next
  key0 <- paste(cmp$cohorts, collapse = "_")
  if (!key0 %in% names(merged_cache)) next
  for (pw in pathways) {
    feat  <- paste0(pw$name, "1")
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

# --- 15. Per-gene tacrolimus descriptive table --------------------------------
# For each tacrolimus-sensitive NFAT gene: mean expression and % of cells
# expressing it in GEM108 pre vs post, per T-cell type, then summarised per gene.
gem_key <- "GEM108_pre_GEM108_post"
if (gem_key %in% names(merged_cache)) {
  log_msg("== Per-gene tacrolimus targets (GEM108 pre vs post)")
  seu_gem <- merged_cache[[gem_key]]
  DefaultAssay(seu_gem) <- "RNA"
  genes_in <- intersect(tacro_nfat_targets, rownames(seu_gem))
  expr <- GetAssayData(seu_gem, layer = "data")[genes_in, , drop = FALSE]
  expr <- t(as.matrix(expr))
  md_gene <- seu_gem@meta.data %>%
    tibble::as_tibble() %>%
    dplyr::select(celltype = predicted.celltype.l2,
                  group    = compare_group) %>%
    dplyr::bind_cols(as.data.frame(expr))

  summary_all_ct_tacro <- md_gene %>%
    dplyr::filter(celltype %in% t_cell_types) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(genes_in),
                        names_to = "gene", values_to = "expr") %>%
    dplyr::group_by(gene, celltype, group) %>%
    dplyr::summarise(
      mean_expr = mean(expr),
      pct_expr  = 100 * mean(expr > 0),
      n_cells   = dplyr::n(),
      .groups   = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = group,
                       values_from = c(mean_expr, pct_expr, n_cells)) %>%
    dplyr::mutate(
      delta_pct_expressing = pct_expr_GEM108_pre - pct_expr_GEM108_post,
      log2fc_mean_expr = log2(pmax(mean_expr_GEM108_pre, 0.001) /
                              pmax(mean_expr_GEM108_post, 0.001))
    )

  gene_summary <- summary_all_ct_tacro %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_celltypes      = dplyr::n(),
      n_lower_in_post  = sum(delta_pct_expressing > 0),
      n_higher_in_post = sum(delta_pct_expressing < 0),
      mean_delta_pct   = mean(delta_pct_expressing),
      mean_log2fc      = mean(log2fc_mean_expr),
      .groups          = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_delta_pct))

  readr::write_tsv(summary_all_ct_tacro,
    file.path(pathway_dir, "tacro_per_celltype_GEM108_pre_vs_post.tsv"))
  readr::write_tsv(gene_summary,
    file.path(pathway_dir, "tacro_per_gene_GEM108_pre_vs_post.tsv"))
} else {
  log_msg("== Skipping per-gene tacrolimus table - GEM108 pre+post merge unavailable")
}

# --- 16. PDF report -----------------------------------------------------------
# One violin panel ("Tcells" or "other") for a comparison x pathway.
safe_plot <- function(key, suffix) {
  obj <- all_results[[key]]
  if (!is.null(obj) && !is.null(obj$plots) && !is.null(obj$plots[[suffix]]))
    return(obj$plots[[suffix]])
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                      label = paste("Not available:\n", key, "\n", suffix)) +
    theme_void()
}
# Stored UMAP pair, or a placeholder.
safe_umap <- function(key) {
  obj <- umap_plots[[key]]
  if (!is.null(obj)) return(obj)
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                      label = paste("Not available:\n", key)) +
    theme_void()
}

# Page number in the bottom-right corner.
make_page_numberer <- function() {
  n <- 0L
  function() {
    n <<- n + 1L
    # Draw in the root viewport so the number lands in the page corner.
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
  body_text <- paste(body_lines, collapse = "\n")
  grid::grid.text(body_text, x = 0.05, y = 0.88, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 9, fontfamily = "sans",
                                  lineheight = 1.3))
  page_number()
}

if (!opt_no_pdf) {
  out_pdf <- file.path(pathway_dir, "pathway_analysis_report.pdf")
  log_msg("== Writing PDF report -> ", out_pdf)
  grDevices::pdf(out_pdf, width = 16, height = 10)

  # --- Methods page -----------------------------------------------------------
  grid::grid.newpage()
  grid::grid.text("Pathway Analysis Report — Methods",
                  x = 0.5, y = 0.97, just = c("centre", "top"),
                  gp = grid::gpar(fontsize = 16, fontface = "bold"))

  left_text <- paste(c(
    "OVERVIEW",
    "Module score analysis of immune signalling pathways in ITK",
    "gain-of-function patients (E42K, T504S) and healthy blood",
    "donors (HBD). GEM108 is an E42K patient sampled before and",
    "after tacrolimus treatment (calcineurin/NFAT inhibitor).",
    "ITK E42K and T504S are novel gain-of-function variants that",
    "enhance TCR signalling with increased calcium flux and",
    "augmented NFAT activity.",
    "",
    "GENE SETS",
    "",
    "1. NFAT (all C2): All NFAT-related gene sets from MSigDB C2",
    "   curated collection, merged. Source: msigdbr, category='C2',",
    "   filtered for 'NFAT' in gene set name.",
    "",
    "2. Tacrolimus-sensitive NFAT targets: 25 literature-curated",
    "   genes representing direct NFAT transcriptional targets",
    "   suppressed by calcineurin inhibition. Includes cytokines",
    "   (IL2, IFNG, TNF, IL10), surface receptors (IL2RA, CTLA4,",
    "   ICOS, PDCD1, FASLG, LAG3, HAVCR2, TIGIT), transcription",
    "   factors (NFATC1, TBX21, GATA3, RORC, BCL6, BATF, IRF4,",
    "   EGR2, EGR3), effector molecules (GZMB, PRF1), and Tfh",
    "   markers (CXCR5, SLAMF1). Sources: Martinez et al. 2018",
    "   Front Immunol; Vafadari et al. 2013 PLOS ONE;",
    "   Klein-Hessling et al. 2017 Nat Commun; Martinez et al.",
    "   2015 Immunity.",
    "",
    "3. Apoptosis (de Cevins): 67-gene custom list from de Cevins",
    "   et al. 2023 Cell Rep Med. 11 gene aliases corrected to",
    "   current HGNC symbols.",
    "",
    "4. Apoptosis (Hallmark): MSigDB Hallmark APOPTOSIS.",
    "   Liberzon et al. 2015 Cell Systems.",
    "",
    "5. IL2-STAT5 signaling: MSigDB Hallmark IL2_STAT5_SIGNALING.",
    "   Genes upregulated by STAT5 in response to IL-2 stimulation.",
    "   Liberzon et al. 2015 Cell Systems.",
    "",
    "6. TNFa-NFkB signaling: MSigDB Hallmark",
    "   TNFA_SIGNALING_VIA_NFKB. Genes regulated by NF-kB in",
    "   response to TNF-alpha. Liberzon et al. 2015 Cell Systems."
  ), collapse = "\n")
  grid::grid.text(left_text, x = 0.03, y = 0.93, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 7.5, fontfamily = "mono",
                                  lineheight = 1.2))

  right_text <- paste(c(
    "MODULE SCORING",
    "Scores computed using Seurat AddModuleScore (Tirosh et al.",
    paste0("2016 Science) on the RNA assay (~", n_genes,
           " genes, log-"),
    "normalised). ctrl=100 control genes per expression bin.",
    "Seed=42 for reproducibility. Unlike the original",
    "natComs_reviews.R, this pipeline integrates with Harmony",
    "(not SCTransform), so the RNA assay retains the full",
    "transcriptome — no separate RNA_full reconstruction is",
    "required.",
    "",
    "Y-AXIS: VIOLIN PLOTS",
    "The y-axis shows the AddModuleScore value: average log-",
    "normalised expression of pathway genes minus the average",
    "expression of expression-matched control genes, per cell.",
    "Positive = pathway genes more highly expressed than expected;",
    "negative = lower than expected. The score is unitless and",
    "relative.",
    "",
    "STATISTICAL TESTING",
    "",
    "Per-cell-type statistics are still computed and written to",
    "the stats/ TSV files (pairwise Wilcoxon rank-sum with BH",
    "adjustment for multi-donor comparisons; Cohen's d effect",
    "size for single-donor/GEM108 comparisons). They are NO",
    "LONGER drawn on the violins: significance brackets and",
    "Cohen's d annotations were removed (2026-05-29) for a",
    "cleaner figure and to avoid clipping under the fixed",
    "y-axis ranges below. Consult the TSVs for significance.",
    "",
    "VIOLIN LAYOUT",
    "Each comparison x pathway is split into two panels:",
    "  - T cells: CD4 (Naive/TCM/TEM/CTL/Proliferating), Treg,",
    "    CD8 (Naive/TCM/TEM/CTL/Proliferating), MAIT, dnT, gdT",
    "  - Other cell types: NK, ILC, B, monocytes/DCs, etc.",
    "",
    "Y-axis range is set per panel from the data:",
    "  max score > 0.6        -> [-0.2, 0.8]",
    "  max score in (0.4,0.6] -> [-0.2, 0.6]",
    "  max score <= 0.4       -> [-0.2, 0.4]",
    "",
    "BATCH NOTE",
    "RNA assay is not batch-corrected. AddModuleScore uses",
    "within-cell control subtraction which is inherently batch-",
    "robust. Single-gene cross-batch comparisons should be",
    "interpreted cautiously."
  ), collapse = "\n")
  grid::grid.text(right_text, x = 0.52, y = 0.93, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 7.5, fontfamily = "mono",
                                  lineheight = 1.2))
  grid::grid.lines(x = c(0.50, 0.50), y = c(0.05, 0.93),
                   gp = grid::gpar(col = "grey70", lwd = 0.5))
  page_number()

  # --- Per-comparison sections ------------------------------------------------
  for (cmp in comparisons) {
    cmp_label <- cmp$label

    text_page(paste("Comparison:", cmp_label), c(
      "",
      paste("Cohorts:", paste(cmp$cohorts, collapse = ", ")),
      "",
      if (cmp$single_donor) {
        paste("NOTE: This comparison involves a single donor (n=1).",
              "Effect sizes (Cohen's d) are in the stats/ TSVs;",
              "results are descriptive and exploratory.")
      } else {
        paste("This comparison has multiple donors per group.",
              "BH-adjusted p-values and Cohen's d are in the stats/ TSVs.")
      },
      "",
      "Violins are split into a T-cell panel and an 'other cell",
      "types' panel (one page per pathway, T cells on top)."
    ))

    # One page per pathway: T-cell panel above the other panel.
    for (pw in pathways) {
      key <- paste0(cmp_label, " | ", pw$label)
      p_t <- safe_plot(key, "Tcells")
      p_o <- safe_plot(key, "other")
      tryCatch({
        print(p_t / p_o +
                patchwork::plot_annotation(
                  title = paste(pw$label, "—", cmp_label),
                  theme = theme(plot.title = element_text(size = 14,
                                                          face = "bold"))))
        page_number()
      }, error = function(e) {
        log_msg("Violin page error: ", pw$label, " ", e$message)
        print(p_t + ggtitle(paste(pw$label, "—", cmp_label)))
        page_number()
      })
    }
  }

  grDevices::dev.off()
  log_msg("== PDF report complete: ", out_pdf)
}

log_msg("== 07_pathway_module_scores.R done.")
