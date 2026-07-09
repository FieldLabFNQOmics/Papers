# Single-cell RNA-seq analysis pipeline

Analysis code for a 10x Genomics single-cell RNA-seq experiment (mouse, Flex /
probe-based chemistry) covering three samples: `SD_BCG`, `SD_PE25SS`, and
`SD_RD1`. The pipeline runs from Cell Ranger output through droplet QC,
annotation, clustering, cell-type identification, differential expression, and
gene set enrichment analysis (GSEA).

## Run order

Scripts are numbered and run in sequence; each consumes the output of the
previous one.

| Step | Script | Purpose |
|------|--------|---------|
| 0 | `cellranger_multi.qsub` + `config_template.csv` | Run `cellranger multi` per sample (one job each) under PBS |
| 1 | `1_create_singleCellExperiment.R` | Build a SingleCellExperiment (SCE) per sample from Cell Ranger output |
| 2 | `2_dropletSelection.R` | Barcode-rank inspection, within-sample doublet annotation, merge samples |
| 3 | `3_preprocessing.R` | Gene annotation (Ensembl + NCBI), per-cell QC, outlier/doublet removal |
| 4 | `4_clustering_and_annotation.R` | Seurat clustering, cell-type annotation, differential expression, GSEA |

## Directory layout

All R scripts resolve paths with the [`here`](https://here.r-lib.org/) package,
so run them from the project root (or open the project's `.Rproj` file). The
scripts read from and write to the following structure:

```
project_root/
├── cellranger/
│   └── <sample>/outs/per_sample_outs/<sample>/count/
│       └── sample_filtered_feature_bc_matrix/     # Cell Ranger input
├── data/
│   ├── <sample>/SCEs/                             # per-sample SCEs (step 1)
│   └── merged/
│       ├── SCEs/                                  # merged + preprocessed SCEs, clustered Seurat object
│       ├── results/                               # QC report, UMAPs, TPM tables
│       ├── plots/<cell_type>/                     # volcano plots (step 4)
│       ├── sanity/<cell_type>/                    # Sanity inputs/outputs (step 4)
│       ├── gsea/<cell_type>/                      # GSEA rank files + reports (step 4)
│       └── visualizations/<cell_type>/            # GSEA dot/network plots (step 4)
└── .cache/R/AnnotationHub/                        # project-local AnnotationHub cache
```

### Key intermediate files

- `data/<sample>/SCEs/<sample>.CellRanger.SCE.rds` — per-sample SCE (step 1)
- `data/merged/SCEs/merged.doublets.rds` — doublet calls (step 2)
- `data/merged/SCEs/merged.emptyDrops.SCE.rds` — merged SCE (step 2)
- `data/merged/SCEs/merged.preprocessed.SCE.rds` — QC'd SCE (step 3)
- `data/merged/SCEs/merged.clustered.SEU.rds` — clustered Seurat object (step 4)

## Configuration to edit before running

**`cellranger_multi.qsub`** — set `#PBS -P PROJECT_CODE`, the `storage` mounts,
and the path to the `cellranger` binary for your environment.

**`config_template.csv`** — copy to `<sample>.csv` per sample and set the
`reference`, `probe-set`, and `fastqs` paths, plus the `fastq_id`.

**`4_clustering_and_annotation.R`** — set the external-tool paths at the top:

```r
SANITY_BIN   <- "/path/to/Sanity"
GSEA_CLI     <- "/path/to/gsea-cli.sh"
GSEA_DB_PATH <- "/path/to/msigdb/mouse"   # directory containing the .gmt files
```

## External tools and reference data

- **Cell Ranger** (10x Genomics) with the mouse mm10-2020-A reference and the
  Chromium Mouse Transcriptome probe set (v1.0.1).
- **Sanity** — expression denoising. https://github.com/jmbreda/Sanity
- **GSEA** (command-line) — https://www.gsea-msigdb.org/gsea/
- **MSigDB mouse collections** as `.gmt` files, placed in `GSEA_DB_PATH`. The
  release used here (edit the filenames in step 4 if you use another):
  `m2.all.v2024.1.Mm.symbols.gmt`, `m3.all.v2024.1.Mm.symbols.gmt`,
  `m5.all.v2024.1.Mm.symbols.gmt`, `mh.all.v2024.1.Mm.symbols.gmt`.
- **AnnotationHub** EnsDb record `AH89211` (mouse, GRCm38/mm10). If unavailable
  in your AnnotationHub snapshot, re-run the query shown in step 3 and update
  the accession.

## R packages

Developed under R (Bioconductor). Core packages, by role:

- **Single-cell infrastructure:** SingleCellExperiment, DropletUtils, scran,
  scater, scuttle, scds, scDblFinder, DropletQC, Seurat, sctransform
- **Annotation:** AnnotationHub, ensembldb, org.Mm.eg.db, Mus.musculus,
  biomaRt, msigdbr
- **Cell-type calling:** SingleR, celldex, scMCA, scRNAseq
- **DE / enrichment:** limma, edgeR, statmod, clusterProfiler, enrichplot, DOSE
- **Clustering / viz:** clustree, igraph, ggraph, patchwork, cowplot, ggrepel,
  UpSetR, RColorBrewer, speckle
- **Utilities:** here, glue, tidyverse (dplyr, stringr, etc.), janitor,
  gridExtra, kableExtra, scales, Matrix, BiocStyle, CellBench

Each script prints `sessionInfo()` at the end; capturing that output (or an
`renv.lock`) alongside a run is recommended for exact version provenance.

## Notes on reproducibility

- The random seed is fixed (`set.seed(42)`) in steps 2 and 4.
- The marker-based gating in step 4 references **hardcoded cluster IDs** (e.g.
  `c(6, 8, 9, 10, 13, 16, 20, 22, 25)`) that correspond to the deposited
  clustered object (`merged.clustered.SEU.rds`). A run re-clustered from scratch
  may assign different cluster numbers; map against the deposited object.
- A DropletQC nuclear-fraction step is intentionally omitted: the Flex chemistry
  is probe-based and yields no intronic reads, so nuclear-fraction scores cannot
  be computed (noted in step 3).
