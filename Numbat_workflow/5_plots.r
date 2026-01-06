library(Seurat)
library(here)
library(Azimuth)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(biomaRt)
library(scuttle)
library(glue)
library(scds)
library(scDblFinder)
library(magrittr)
library(UpSetR)
library(dplyr)
library(viridis)

# Read in seurat object from 1_qc.r
load(here("seu.RData"))

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))

ElbowPlot(seu, ndims = 50)

# Run UMAP and clustering
seu <- RunUMAP(seu, dims = 1:25)  
seu <- FindNeighbors(seu, dims = 1:25)
seu <- FindClusters(seu, resolution = 0.5)

# Plot UMAP
p1 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", label = TRUE)
p2 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# Show both plots
p1 + p2

##### Cell type annotation #####
# Convert to seurat v3 object for upload to Azimuth online
options(Seurat.object.assay.version = "v3")

# Create a v3-style Seurat object
counts_matrix <- LayerData(seu, assay = "RNA", layer = "counts")
data_matrix <- LayerData(seu, assay = "RNA", layer = "data")

# Create new v3-style object
seu_v3 <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = seu@meta.data,
  project = "2046t2"
)

# Add normalized data
if (!is.null(data_matrix)) {
  seu_v3[["RNA"]]@data <- data_matrix
}

# Add reductions
if ("umap" %in% Reductions(seu)) {
  seu_v3[["umap"]] <- seu[["umap"]]
}

if ("pca" %in% Reductions(seu)) {
  seu_v3[["pca"]] <- seu[["pca"]]
}

# Save diet as RDS
seu_clean <- DietSeurat(seu, 
                        layers = c("counts", "data"),
                        features = rownames(seu),
                        assays = "RNA",
                        dimreducs = NULL,
                        graphs = NULL)

saveRDS(seu_clean, file = here("seu_diet_v3_compatible.rds"))


# After uploading seu_v3_compatible.rds to azimuth online and then obtaining the tsv and Rds objects
azimuth_umap <- readRDS(here("azimuth_umap.Rds"))
azimuth_predictions <- read.table(here("azimuth_pred.tsv"), 
                                  sep = "\t", header = TRUE, row.names = 1)
# Check what you got
head(azimuth_predictions)
colnames(azimuth_predictions)

# Add the Azimuth UMAP coordinates
seu[["azimuth_umap"]] <- azimuth_umap

# Check it was added
Reductions(seu)

identical(rownames(azimuth_predictions), colnames(seu))

# If they match, add the predictions
seu <- AddMetaData(seu, azimuth_predictions)

# Plot using Azimuth UMAP
p1 <- DimPlot(seu, reduction = "umap",
              group.by = "predicted.celltype.l2", 
              label = TRUE, repel = TRUE,
              cols = viridis(n = length(unique(seu$predicted.celltype.l2)), 
                             option = "H")) +
  ggtitle("Azimuth Cell Type Predictions (Level 2)") +
  NoLegend()

# Plot clones of interest after extracting barcodes from numbat analysis
clone2_barcodes <- readLines("/path/to/numbat/clone_2_barcodes.txt")

seu$clone_2_cells <- ifelse(colnames(seu) %in% clone2_barcodes, 1, 0)
table(seu$clone_2_cells)

p2 <- FeaturePlot(seu, features = "clone_2_cells", 
                  cols = c("lightgray", "red"), 
                  pt.size = 0.5) +
  ggtitle("Clone 2 Cells") +
  NoLegend()

clone3_barcodes <- readLines("/path/to/numbat/clone_3_barcodes.txt")
seu$clone_3_cells <- ifelse(colnames(seu) %in% clone3_barcodes, 1, 0)

table(seu$clone_3_cells)

p3 <- FeaturePlot(seu, features = "clone_3_cells", 
                  cols = c("lightgray", "red"), 
                  pt.size = 0.5) +
  ggtitle("Clone 3 Cells") +
  NoLegend()

# Create a combined clone annotation
seu$clone_status <- case_when(
  seu$clone_2_cells == 1 ~ "Clone 2",
  seu$clone_3_cells == 1 ~ "Clone 3", 
  TRUE ~ "Other"
)

p4 <- DimPlot(seu, group.by = "clone_status", 
              cols = c("Clone 2" = "red", "Clone 3" = "blue", "Other" = "lightgray"),
              pt.size = 0.5) +
  ggtitle("Clones combined")

plot_grid(p1, p2, p3, p4, ncol = 2)

