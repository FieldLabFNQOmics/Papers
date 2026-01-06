library(numbat)
library(pagoda2)
library(data.table)
library(dplyr)

# Read in matrix, barcodes and features from cellranger
count_mat = readMM("filtered_feature_bc_matrix/matrix.mtx.gz")
cells = fread('filtered_feature_bc_matrix/barcodes.tsv.gz', header = F)$V1
genes = fread('filtered_feature_bc_matrix/features.tsv.gz', header = F)$V2

colnames(count_mat) = cells
rownames(count_mat) = genes
count_mat = as.matrix(count_mat)
count_mat = rowsum(count_mat, rownames(count_mat))
count_mat = as(count_mat, "dgCMatrix")
count_mat_ref = count_mat

p2 = pagoda2::basicP2proc(count_mat_ref, n.cores = 12)

clusters = p2$clusters$PCA$multilevel

ref_annot = data.frame(
  cell = names(clusters),
  group = unname(clusters)
)

ref = numbat::aggregate_counts(
  count_mat_ref,
  ref_annot %>% group_by(group) %>% filter(n() > 100)
)

saveRDS(ref, file = "ref.rds")