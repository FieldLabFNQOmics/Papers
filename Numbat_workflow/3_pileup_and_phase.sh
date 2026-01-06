#!/bin/bash

Rscript pileup_and_phase.R \
--label label \
--samples sample \
--bams sample_alignments.bam \
--barcodes sample_filtered_feature_bc_matrix/barcodes.tsv.gz \
--outdir /path/to/outdir \
--gmap /path/to/genetic_map_hg38_withX.txt.gz \
--snpvcf /path/to/genome1k/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /path/to/genome1k/1000G_hg38 \
--ncores 28
