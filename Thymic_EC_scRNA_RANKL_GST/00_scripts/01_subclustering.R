## Goal 1: Subclustering of the scRNA-seq data to identify thymic EC subsets

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)

# Load the Seurat object containing the EC clusters
EC_data <- LoadH5Seurat(here::here(
  "Thymic_EC_scRNA_RANKL_GST",
  "02_processed_data",
  "GST_RANKL_EC.h5seurat"
))

# Verify
DimPlot(object = EC_data)

# Subcluster
# No need to rerun normalization, start from the find variable feature step

# Find variable features
DefaultAssay(EC_data) <- "RNA"

EC_data <-
  FindVariableFeatures(EC_data,
                       dispersion.function = "vst",
                       nfeatures = 2000, assay = "RNA")

# Scale the data
EC_data <-
  ScaleData(EC_data, features = rownames(EC_data))

# Linear dimensional reduction
EC_data <- RunPCA(EC_data)

# Choose the number of PCs to be used in the subsequent steps
ElbowPlot(EC_data, ndims = 50)

# UMAP embedding
EC_data <-
  RunUMAP(EC_data, reduction = "pca", dims = 1:12)

# Cluster the cells
EC_data <- FindNeighbors(EC_data, dims = 1:12)
EC_data <- FindClusters(EC_data, resolution = 0.5)

# Visualization 
DimPlot(EC_data, label = TRUE)
DimPlot(EC_data, split.by = "sample", label = TRUE)

# Check that all clusters express known EC markers 
FeaturePlot(EC_data, features = c("Pecam1", "Cdh5"))

# Use a DEG analysis to identify cell-types
all_markers <-
  FindAllMarkers(
    object = EC_data,
    assay = "RNA",
    slot = "data",
    verbose = TRUE,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.25,
    test.use = "bimod"
  )

filtered_genes <-
  all_markers[all_markers$p_val_adj < 1e-30 &
                abs(all_markers$avg_log2FC) > log2(2), ]

top_genes <-
  filtered_genes %>% group_by(cluster) %>% top_n(10, avg_log2FC)

heatmap_colors <-
  c(
    "#27408B",
             "#3A5FCD",
             "#3288BD",
             "#66C2A5",
             "#ABDDA4",
             "#E6F598",
             "#FEE08B",
             "#FDAE61",
             "#F46D43",
             "#D53E4F",
             "#8B2323"
  )

avg_exp <-
  AverageExpression(
    EC_data,
    features = top_genes$gene,
    return.seurat = TRUE,
    slot = "data", 
    assays = "RNA"
  )

DoHeatmap(
  avg_exp,
  features = top_genes$gene,
  size = 2,
  label = FALSE,
  draw.lines = FALSE
) +
  ggplot2::scale_fill_gradientn(colors = heatmap_colors, na.value = "white") +
  theme(axis.text.y = element_text(size = 8, face = "bold")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.5, 'cm')
  ) 

# Based on the DEGs and the literature, we can identify the subclusters below

# Cluster 0: Capillary_1
# Cluster 1: Capillary_2
# Cluster 2: Choroid plexus capillary ECs (https://www.cell.com/cell/pdf/S0092-8674(20)30062-3.pdf)
# Cluster 3: Arterial 
# Cluster 4: Venous 2
# Cluster 5 Capillary_2
# Cluster 6: Venous 1 
# Cluster 7 and 8: smooth muscle cell (SMC) - expression of Myl9 and Vim

# Remove contaminants (SMC)
EC_data <-
  subset(EC_data,
         idents = c(7,8),
         invert = TRUE)

DimPlot(EC_data, label = TRUE)

# Rename clusters
EC_data[["cell_types"]] <- NA

subsets <- list(
  "Capillary 1" = 0,
  "Capillary 2" = c(1,5),
  "Choroid plexus like EC" = 2,
  "Arterial" = 3,
  "Venous 2" = 4,
  "Venous 1" = 6
)

for (list_name in names(subsets)) {
  EC_data$cell_types[which(EC_data@active.ident %in% subsets[[list_name]])] <- list_name
}

DimPlot(EC_data, reduction = "umap", label = TRUE, group.by = "cell_types") 

# Export cell-types
Idents(object = EC_data) <- "cell_types"

# Save the Seurat object containing the subclustering information
SaveH5Seurat(
  EC_data,
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "02_processed_data",
    "GST_RANKL_EC_subclusters"
  ),
  overwrite = TRUE
)
