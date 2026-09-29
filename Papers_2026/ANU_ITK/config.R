# ------------------------------------------------------------------------------
# config.R
# Paths, capture table, thresholds, cell-type subsets, contrasts and figure
# settings. Sourced by every other script. Definitions only.
# ------------------------------------------------------------------------------

# --- Project root -------------------------------------------------------------
# All outputs go under this directory.
project_root <- "/path/to/project_root"

# --- Raw data locations (read-only) -------------------------------------------
# Batch 1: Aug 2022, 4 captures.
batch1_root <- "/path/to/cellranger/batch1_Aug2022"
batch1_captures <- paste0("Maurice", 1:4)

# Batch 2: May 2023, 16 captures.
batch2_root <- "/path/to/cellranger/batch2_May2023"
batch2_captures <- paste0("Maurice_GEX_Feature_", 1:16)

# --- Capture to dataset_{i} mapping -------------------------------------------
# Index i = 1..20. Rows 1-4 are batch 1, rows 5-20 are batch 2.
capture_table <- data.frame(
  i         = 1:20,
  dataset   = paste0("dataset_", 1:20),
  capture   = c(batch1_captures, batch2_captures),
  batch     = c(rep("Aug", 4), rep("May", 16)),
  batch_root = c(rep(batch1_root, 4), rep(batch2_root, 16)),
  stringsAsFactors = FALSE
)

# --- CellSweep toggle ---------------------------------------------------------
# Published results use FALSE (raw counts). TRUE reads a CellSweep-denoised
# object and writes to *_cellsweep output folders. The CellSweep scripts are not
# included in this release.
use_cellsweep <- FALSE

# --- Output directory layout --------------------------------------------------
# Everything is written under project_root/pipeline/.
pipeline_root   <- file.path(project_root, "pipeline")
sce_dir         <- file.path(pipeline_root, "SCEs")
emptydrops_dir  <- file.path(pipeline_root, "emptyDrops")
emptydrops_diag <- file.path(emptydrops_dir, "diagnostics")
doublets_dir    <- file.path(pipeline_root, "doublets")
qc_dir          <- file.path(pipeline_root, "qc_plots")
merged_dir      <- file.path(pipeline_root, "merged")
de_dir          <- file.path(pipeline_root,
                             if (isTRUE(use_cellsweep)) "DE_cellsweep"   else "DE")
gsea_dir        <- file.path(pipeline_root,
                             if (isTRUE(use_cellsweep)) "GSEA_cellsweep" else "GSEA")
sanity_dir      <- file.path(pipeline_root, "sanity")

# --- Annotated object path ----------------------------------------------------
# Raw or CellSweep object, depending on use_cellsweep.
annotated_rds_filename <- function() {
  if (isTRUE(use_cellsweep)) {
    "seurat_merged_harmony_azimuth_cellsweep.rds"
  } else {
    "seurat_merged_harmony_azimuth.rds"
  }
}
# Full path to the annotated object written by script 04.
annotated_rds_path <- function() {
  file.path(merged_dir, annotated_rds_filename())
}

# --- External tools and references --------------------------------------------
# Edit these paths for your system. sanity_bin and gsea_cli are not used by the
# included scripts.
sanity_bin   <- "/path/to/Sanity/bin/Sanity"
gsea_cli     <- "/path/to/GSEA_Linux_4.3.2/gsea-cli.sh"
gsea_gmt_dir <- "/path/to/msigdb_gmt_dir"
gsea_version <- "v2023.1.Hs.symbols"

# Local Azimuth PBMC reference: the folder holding ref.Rds and idx.annoy.
# Compute nodes have no internet, so Azimuth must not download it.
azimuth_pbmcref_dir <- "/path/to/pbmcref.SeuratData/inst/azimuth"

# Local copy of the Azimuth homolog table (used by script 04).
azimuth_homologs_path <- "/path/to/azimuth/homologs.rds"

# CellSweep install (not used when use_cellsweep is FALSE).
cellsweep_root   <- "/path/to/cellsweep"
cellsweep_python <- file.path(cellsweep_root, ".venv/bin/python")
cellsweep_cli    <- file.path(cellsweep_root, ".venv/bin/cellsweep")

cellsweep_io_dir <- file.path(pipeline_root, "cellsweep_io")
cellsweep_log_dir <- file.path(pipeline_root, "cellsweep_io", "logs")

# --- Parallelism --------------------------------------------------------------
# One BLAS / OpenMP thread per process.
suppressPackageStartupMessages(library(RhpcBLASctl))
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
num_cores <- max(1, parallel::detectCores() - 1)

# --- QC thresholds / pipeline parameters --------------------------------------
emptydrops_lower     <- 100     # UMI threshold below which barcodes form ambient profile
emptydrops_fdr       <- 0.001   # FDR cutoff for retaining non-empty droplets
hashed_confident_min <- 2       # log2 fold-change between best and 2nd-best HTO
mito_nmads           <- 3       # MAD cutoff for mitochondrial % outlier detection
min_cells_per_gene   <- 20      # genes detected in fewer than this many cells are dropped
dropletqc_umi_rescue <- 10^2.5
dropletqc_nf_rescue  <- 0.02

# --- Sample sheet -------------------------------------------------------------
# Columns: capture, hto, patient, cohort, treatment, batch, notes. Read in script 03.
sample_sheet_path <- file.path(pipeline_root, "sample_sheet.csv")

# --- Samples excluded from the analysis ---------------------------------------
# Cells matching a (capture, HTO) row here are removed in script 03.
# GEM108 pre-treatment sample 2 (Maurice2 / Hashtag2) is excluded.
# PMAI0025 is not listed here: its cohort is blank in the sample sheet, so script
# 03 drops it.
samples_to_exclude <- data.frame(
  capture = c("Maurice2"),
  hto     = c("Hashtag2"),
  reason  = c("GEM108 pre_s2 — confounded NF-kB/inflammation artefact (dropped 2026-05-29)"),
  stringsAsFactors = FALSE
)

# --- Cell-type subsets for DE -------------------------------------------------
# Script 05 runs every contrast on each subset (Azimuth L2 labels).
de_cell_subsets <- list(
  all_cells = NULL,
  CD4       = c("CD4 T", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL",
                "CD4 Proliferating"),
  CD8       = c("CD8 T", "CD8 Naive", "CD8 TCM", "CD8 TEM",
                "CD8 Proliferating"),
  Bcell     = c("B naive", "B intermediate", "B memory", "Plasmablast"),
  NK        = c("NK", "NK_CD56bright", "NK Proliferating"),
  Monocyte  = c("CD14 Mono", "CD16 Mono"),
  Treg      = "Treg",
  # T-cell subtypes used by the script 06 spot matrices.
  CD4_CTL   = "CD4 CTL",
  CD4_Naive = "CD4 Naive",
  CD4_TCM   = "CD4 TCM",
  CD4_TEM   = "CD4 TEM",
  CD8_Naive = "CD8 Naive",
  CD8_TCM   = "CD8 TCM",
  CD8_TEM   = "CD8 TEM",
  # Activated (memory / effector) subsets. Cells are pooled before pseudobulk,
  # so each is its own fit. Task ids follow list order, so new subsets must be
  # added at the end.
  CD4_Activated = c("CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating"),
  CD8_Activated = c("CD8 TCM", "CD8 TEM", "CD8 Proliferating"),
  # Non-conventional T cells and other lineages for the "other" figure panel.
  # NK_ILC = NK + NK Proliferating + NK_CD56bright + ILC. These subsets are rare,
  # so some contrasts have one donor on a side.
  MAIT          = "MAIT",
  dnT           = "dnT",
  gdT           = "gdT",
  NK_ILC        = c("NK", "NK_CD56bright", "NK Proliferating", "ILC"),
  DC            = c("cDC1", "cDC2", "ASDC", "pDC"),
  Hematopoietic = c("Platelet", "Eryth", "HSPC")
)

# Cell-type columns for scripts 08, 10 and 11, in display order. Taken from
# de_cell_subsets. all_cells is not included.
consensus_cell_types <- intersect(
  c("CD4", "CD4_Naive", "CD4_Activated", "CD4_TCM", "CD4_TEM", "CD4_CTL",
    "Treg",
    "CD8", "CD8_Naive", "CD8_Activated", "CD8_TCM", "CD8_TEM",
    "MAIT", "gdT",  # dnT not plotted
    "NK", "NK_ILC", "Bcell", "Monocyte", "DC", "Hematopoietic"),
  setdiff(names(de_cell_subsets), "all_cells"))

# --- Consensus pathway sets (scripts 10 and 11) -------------------------------
# Hallmark: display label -> MSigDB gene-set ID. KRAS signaling has two rows
# (UP and DN).
hallmark_consensus <- c(
  "Allograft Rejection"        = "HALLMARK_ALLOGRAFT_REJECTION",
  "Apoptosis"                  = "HALLMARK_APOPTOSIS",
  "Complement"                 = "HALLMARK_COMPLEMENT",
  "IL2 STAT5 Signaling"        = "HALLMARK_IL2_STAT5_SIGNALING",
  "IL6 JAK STAT3 Signaling"    = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "TNFA Signaling via NFkB"    = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "Inflammatory Response"      = "HALLMARK_INFLAMMATORY_RESPONSE",
  "Interferon Alpha Response"  = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "Interferon Gamma Response"  = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "KRAS Signaling Up"          = "HALLMARK_KRAS_SIGNALING_UP",
  "KRAS Signaling Dn"          = "HALLMARK_KRAS_SIGNALING_DN",
  "PI3K AKT MTOR Signaling"    = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "Oxidative Phosphorylation"  = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

# Disease Ontology: display label -> DO term description(s). Matched exactly
# (case-insensitive) by match_do_consensus() in utils.R.
do_consensus <- list(
  list(label = "Bacterial Infectious Disease",
       aliases = c("bacterial infectious disease",
                   "primary bacterial infectious disease")),
  list(label = "Rheumatic Disease",
       aliases = c("rheumatic disease", "rheumatism", "rheumatic fever")),
  list(label = "Systemic Scleroderma",
       aliases = c("systemic scleroderma", "systemic sclerosis", "scleroderma")),
  list(label = "Pancreas Disease",
       aliases = c("pancreas disease", "pancreatic disease")),
  list(label = "Pneumonia",
       aliases = c("pneumonia")),
  list(label = "COVID-19",
       aliases = c("covid-19", "covid 19", "coronavirus disease 2019")),
  list(label = "Myeloid Leukemia",
       aliases = c("myeloid leukemia", "myeloid leukaemia")),
  list(label = "Rheumatoid Arthritis",
       aliases = c("rheumatoid arthritis", "rheumatic arthritis")),
  list(label = "Atherosclerosis",
       aliases = c("atherosclerosis", "arteriosclerosis")),
  list(label = "Autoimmune Disease",
       aliases = c("autoimmune disease")),
  list(label = "Inflammatory Bowel Disease",
       aliases = c("inflammatory bowel disease")),
  # Rows added for the clinical features of the E42K and T504S families.
  list(label = "Sarcoidosis",
       aliases = c("sarcoidosis")),
  list(label = "Type 1 Diabetes",
       aliases = c("type 1 diabetes mellitus")),
  list(label = "Viral Infections",
       aliases = c("viral infectious disease"))
)

# Rows are alphabetical, with the two KRAS rows last.
.hallmark_kras  <- c("KRAS Signaling Up", "KRAS Signaling Dn")
hallmark_levels <- c(sort(setdiff(names(hallmark_consensus), .hallmark_kras)),
                     .hallmark_kras)
do_levels       <- sort(vapply(do_consensus, function(x) x$label, character(1)))

# --- Figure rows for the script 11 figure heatmaps ----------------------------
# Subsets of the consensus rows, kept in the same order.
.hallmark_figure_want <- c(
  "Allograft Rejection", "Apoptosis", "IL2 STAT5 Signaling",
  "IL6 JAK STAT3 Signaling", "Inflammatory Response",
  "Interferon Gamma Response", "TNFA Signaling via NFkB")
.do_figure_want <- c(
  "Autoimmune Disease", "Inflammatory Bowel Disease", "Rheumatic Disease",
  "Rheumatoid Arthritis", "Sarcoidosis", "Systemic Scleroderma",
  "Type 1 Diabetes", "Viral Infections")
stopifnot(all(.hallmark_figure_want %in% hallmark_levels),
          all(.do_figure_want       %in% do_levels))
hallmark_figure_levels <- intersect(hallmark_levels, .hallmark_figure_want)
do_figure_levels       <- intersect(do_levels,       .do_figure_want)

# --- Column panels for the script 11 figure heatmaps --------------------------
# Names = axis labels, values = de_cell_subsets keys, in left-to-right order.
# The full heatmaps keep all consensus_cell_types columns.
heatmap_figure_panels <- list(
  tcell = list(
    key      = "tcell",
    title    = "T cell subsets",
    # Total CD4 / CD8 include naive cells. Activated excludes them. Each column is
    # its own fit.
    columns  = c("Total CD4"     = "CD4",
                 "Naive CD4"     = "CD4_Naive",
                 "Activated CD4" = "CD4_Activated",
                 "Tregs"         = "Treg",
                 "Total CD8"     = "CD8",
                 "Naive CD8"     = "CD8_Naive",
                 "Activated CD8" = "CD8_Activated")),
  other = list(
    key      = "other",
    title    = "Other cell types",
    # dnT is not plotted (too few cells), but its DE tasks still run.
    columns  = c("MAIT"                = "MAIT",
                 "gdT"                 = "gdT",
                 "NK cells"            = "NK_ILC",
                 "B cells"             = "Bcell",
                 "Myeloid cells"       = "Monocyte",
                 "DCs"                 = "DC",
                 "Hematopoietic cells" = "Hematopoietic"))
)

# Every panel column must be a de_cell_subsets key.
local({
  want <- unlist(lapply(heatmap_figure_panels, function(p) unname(p$columns)))
  miss <- setdiff(want, names(de_cell_subsets))
  if (length(miss))
    stop("heatmap_figure_panels names subset(s) absent from de_cell_subsets: ",
         paste(miss, collapse = ", "))
  dup <- want[duplicated(want)]
  if (length(dup))
    stop("heatmap_figure_panels repeats subset(s) across panels: ",
         paste(unique(dup), collapse = ", "))
})

# --- Azimuth L2 labels and coverage check -------------------------------------
# Every L2 label must be in a figure group, unless listed in
# azimuth_l2_unplotted.
azimuth_l2_labels <- c(
  "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating",
  "Treg",
  "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating",
  "MAIT", "dnT", "gdT",
  "NK", "NK Proliferating", "NK_CD56bright", "ILC",
  "B naive", "B intermediate", "B memory", "Plasmablast",
  "CD14 Mono", "CD16 Mono",
  "cDC1", "cDC2", "pDC", "ASDC",
  "Platelet", "Eryth", "HSPC",
  "Doublet")

# Labels not shown in any figure group.
azimuth_l2_unplotted <- c("Doublet", "dnT")

local({
  covered <- unique(unlist(
    lapply(heatmap_figure_panels,
           function(p) unlist(de_cell_subsets[unname(p$columns)],
                              use.names = FALSE)),
    use.names = FALSE))
  orphan <- setdiff(azimuth_l2_labels, c(covered, azimuth_l2_unplotted))
  if (length(orphan))
    stop("Azimuth L2 label(s) in NO heatmap_figure_panels column - those cells ",
         "would be missing from every heatmap without any error: ",
         paste(orphan, collapse = ", "),
         "\n  Add them to a de_cell_subsets entry that a panel column names, ",
         "or list them in azimuth_l2_unplotted if that is deliberate.")
  # A panel may not name a label that Azimuth never produces.
  unknown <- setdiff(covered, c(azimuth_l2_labels, "CD4 T", "CD8 T"))
  if (length(unknown))
    stop("heatmap_figure_panels column(s) name label(s) that are not Azimuth L2: ",
         paste(unknown, collapse = ", "))
})

# --- Violin cell-type groups (scripts 07 and 07a) -----------------------------
# Taken from heatmap_figure_panels, so violins and heatmaps use the same groups.
# T-cell groups overlap (Total CD4 = Naive CD4 + Activated CD4).
.VIOLIN_SUFFIX <- c(tcell = "Tcells", other = "other")

violin_celltype_groups <- stats::setNames(
  lapply(names(heatmap_figure_panels), function(nm) {
    p <- heatmap_figure_panels[[nm]]
    list(suffix = .VIOLIN_SUFFIX[[nm]],
         label  = p$title,
         #          label -> member L2 labels, in panel order
         groups = lapply(p$columns, function(k) de_cell_subsets[[k]]))
  }),
  names(heatmap_figure_panels))

local({
  bad <- names(which(vapply(
    unlist(lapply(violin_celltype_groups, function(p) p$groups),
           recursive = FALSE),
    function(g) length(g) == 0L, logical(1))))
  if (length(bad))
    stop("violin_celltype_groups: group(s) with no member labels: ",
         paste(bad, collapse = ", "))
})

# --- UMAP T-cell panel labels (umaps_by_cohort.R) -----------------------------
# The T-cell heatmap labels plus MAIT, dnT and gdT. The UMAPs keep these three
# on the T-cell panel, unlike the heatmaps and violins.
umap_tcell_labels <- setdiff(
  unique(c(unlist(violin_celltype_groups$tcell$groups, use.names = FALSE),
           "MAIT", "dnT", "gdT")),
  c("CD4 T", "CD8 T"))

# Script 08's DO GSEA cache, also read by scripts 10 and 11.
dea08_dir       <- file.path(pipeline_root,
                             if (isTRUE(use_cellsweep)) "dea_plots_cellsweep"
                                                       else "dea_plots")
dea08_cache_dir <- file.path(dea08_dir, "cache")

# --- Disease Ontology gene-set size bounds ------------------------------------
# Used by gseDO / enrichDO in script 08. Sets outside these bounds are not tested
# and show as grey tiles in script 11. Size is counted within each subset's
# gene list.
do_gsea_min_size <- 10
do_gsea_max_size <- 2000

# --- Immune pathways for script 06 --------------------------------------------
# Seven Hallmark pathways, named exactly as in the MSigDB .gmt files.
hallmark_immune_pathways <- c(
  "HALLMARK_APOPTOSIS",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

# Axis labels for the pathways above, in the same order.
hallmark_immune_pathway_labels <- c(
  "APOPTOSIS",
  "IL2-STAT5\nSIGNALING",
  "IL6-JAK-STAT3\nSIGNALING",
  "INFLAMMATORY\nRESPONSE",
  "INTERFERON-γ\nRESPONSE",
  "PI3K-AKT-MTOR\nSIGNALING",
  "TNFα SIGNALING\nVIA NFκB"
)

# Dotplot rows, top to bottom, and the subset each row reads from.
dotplot_subset_rows <- data.frame(
  group   = c("CD4+ clusters","CD4+ clusters","CD4+ clusters",
              "CD8+ clusters","CD8+ clusters","CD8+ clusters"),
  cluster = c("Naive","TCM","TEM","Naive","TCM","TEM"),
  subset  = c("CD4_Naive","CD4_TCM","CD4_TEM",
              "CD8_Naive","CD8_TCM","CD8_TEM"),
  stringsAsFactors = FALSE
)

# --- Spot matrices (script 06) ------------------------------------------------
# Cell types shown in the spot matrices.
spot_matrix_cell_types <- c("CD4_CTL", "CD4_Naive", "CD4_TCM", "CD4_TEM",
                            "CD8_Naive", "CD8_TCM", "CD8_TEM")

# Y-axis order. The first level is at the bottom.
spot_matrix_y_order <- c("CD8_TEM", "CD8_TCM", "CD8_Naive", "CD4_CTL",
                         "CD4_TEM", "CD4_TCM", "CD4_Naive")

# Curated DisGeNET terms per contrast. Contrasts not listed use terms found in
# at least 4 cell types.
disgenet_curated_terms <- list(
  HBD_vs_E42K_affected = c("Epstein-Barr Virus Infections", "Juvenile arthritis",
                         "Myasthenia Gravis", "Lymphoma, T-Cell, Cutaneous",
                         "Cytomegalovirus Infections",
                         "Glucocorticoid Receptor Deficiency",
                         "Immune System Diseases", "Nephritis",
                         "Felty Syndrome"),
  HBD_vs_T504S       = c("Juvenile arthritis", "Juvenile rheumatoid arthritis",
                         "Adult Classical Hodgkin Lymphoma", "Nephritis",
                         "Leukemia, T-Cell", "Sarcoidosis"),
  E42K_affected_vs_T504S = c("Enteritis", "Adult Classical Hodgkin Lymphoma",
                         "Adult Hodgkin Lymphoma",
                         "Epstein-Barr Virus Infections", "Nephritis",
                         "Sarcoidosis"),
  GEM108_pre_vs_HBD  = c("Adult Classical Hodgkin Lymphoma",
                         "Pauciarticular juvenile rheumatoid arthritis",
                         "Tumor Immunity", "Felty Syndrome", "Infection"),
  GEM108_post_vs_HBD = c("Adult Classical Hodgkin Lymphoma",
                         "Lymphoid neoplasm", "Myasthenia Gravis",
                         "Juvenile arthritis", "Felty Syndrome"),
  GEM108_pre_vs_post = c("Juvenile arthritis", "Arthritis, Psoriatic",
                         "Glucocorticoid Receptor Deficiency",
                         "Adult Classical Hodgkin Lymphoma",
                         "Autoinflammatory disease",
                         "Epstein-Barr Virus Infections", "Glomerulonephritis")
)

# PNG widths for the DisGeNET plots (height 8). Unlisted contrasts use 15.
spot_matrix_dgn_widths <- list(
  HBD_vs_E42K_affected   = 15,
  HBD_vs_T504S           = 12,
  E42K_affected_vs_T504S = 12,
  GEM108_pre_vs_HBD  = 11,
  GEM108_post_vs_HBD = 11,
  GEM108_pre_vs_post = 13
)

# --- Cohort groups ------------------------------------------------------------
# The sample-sheet cohort column is a partition:
#   HBD              healthy blood donors
#   E42K_affected    PMAI0017, PMAI0018, PMAI0023 and GEM108
#   E42K_unaffected  PMAI0024
#   T504S, D135Y     ITK T504S and D135Y carriers
#   (blank)          PMAI0025, removed in script 03
# E42K_carriers is the union of E42K_affected and E42K_unaffected. In a contrast,
# a member cohort that is also a level of that contrast keeps its own label.
# Wherever E42K_affected or E42K_carriers is used, GEM108 contributes its
# pre-treatment cells only.
cohort_groups <- list(
  E42K_carriers = c("E42K_affected", "E42K_unaffected")
)

# --- Contrasts ----------------------------------------------------------------
# Pseudobulk contrasts between cohorts. Batch is a covariate in the design.
# GEM108 pre-treatment cells count as one E42K_affected donor.
pseudobulk_contrasts <- list(
  HBD_vs_E42K_affected             = c("HBD",  "E42K_affected"),
  HBD_vs_T504S                     = c("HBD",  "T504S"),
  HBD_vs_D135Y                     = c("HBD",  "D135Y"),
  HBD_vs_E42K_unaffected           = c("HBD",  "E42K_unaffected"),
  E42K_affected_vs_T504S           = c("E42K_affected", "T504S"),
  E42K_affected_vs_D135Y           = c("E42K_affected", "D135Y"),
  E42K_affected_vs_E42K_unaffected = c("E42K_affected", "E42K_unaffected"),
  # E42K_carriers contrasts. E42K_carriers_vs_E42K_unaffected gives the same
  # result as E42K_affected_vs_E42K_unaffected.
  HBD_vs_E42K_carriers             = c("HBD", "E42K_carriers"),
  E42K_carriers_vs_T504S           = c("E42K_carriers", "T504S"),
  E42K_carriers_vs_D135Y           = c("E42K_carriers", "D135Y"),
  E42K_carriers_vs_E42K_unaffected = c("E42K_carriers", "E42K_unaffected")
)

# Within-patient: GEM108 pre vs post, single-cell limma blocked on capture.
# Descriptive: GEM108 vs HBD, pseudobulk logFC and GSEA, no per-gene p-values.
within_patient_contrasts <- list(
  GEM108_pre_vs_post = list(
    patient = "GEM108",
    group_col = "treatment",
    levels = c("pre_treatment", "post_treatment")
  )
)

descriptive_contrasts <- list(
  GEM108_vs_HBD = list(
    group_col = "cohort_or_patient",      # GEM108 (pre + post pooled) vs HBD
    levels    = c("GEM108", "HBD")
  ),
  GEM108_pre_vs_HBD = list(
    group_col = "cohort_or_patient_tx",   # GEM108_pre split out
    levels    = c("GEM108_pre", "HBD")
  ),
  GEM108_post_vs_HBD = list(
    group_col = "cohort_or_patient_tx",   # GEM108_post split out
    levels    = c("GEM108_post", "HBD")
  )
)

# MSigDB collections to run for every contrast's ranked list.
gsea_collections <- c("h", "c2", "c5", "c7", "c8")

# Seed for fgsea. run_gsea() uses gsea_seed + sum(utf8ToInt(collection)), so
# results do not depend on which collections are run.
gsea_seed <- 42

# --- Plot palettes ------------------------------------------------------------
# Cell-type colours come from make_celltype_palette() in utils.R.
hto_palette_name     <- "Paired"
capture_palette_name <- "Accent"
