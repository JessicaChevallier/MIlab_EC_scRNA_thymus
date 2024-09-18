## Goal 1: Subclustering of the scRNA-seq data to identify human thymic EC subsets

## Data taken from the following publication:

# 1: Bautista, Jhoanne L., et al. "Single-cell transcriptional profiling 
# of human thymic stroma uncovers novel cellular heterogeneity 
# in the thymic medulla." Nature Communications 12.1 (2021): 1096.

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)

# Load the data
EC_subset <- LoadH5Seurat(here::here(
  "Thymic_EC_scRNA_reanalysis",
  "Bautista_et_al_Nature_Communications_2021",
  "02_processed_data",
  "EC_subset.h5seurat")
)

# No need to rerun normalization, start from the find variable feature step
# Find variable features
EC_subset <-
  FindVariableFeatures(EC_subset,
                       dispersion.function = "vst",
                       nfeatures = 2000, assay = "RNA")

# Run fastMNN
EC_subset <-
  RunFastMNN(
    object.list = SplitObject(EC_subset, split.by = "orig.ident"),
    features = 2000, assay = "RNA")

# UMAP embedding and cluster the cells
EC_subset <-
  RunUMAP(
    EC_subset,
    reduction = "mnn",
    dims = 1:30
  )

EC_subset <-
  FindNeighbors(
    EC_subset,
    reduction = "mnn",
    dims = 1:30
  )

EC_subset <- FindClusters(EC_subset, resolution = 0.2)

# # Check for potential cell-cycle bias
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# EC_subset <-
#   CellCycleScoring(EC_subset,
#                    s.features = s.genes,
#                    g2m.features = g2m.genes,
#                    set.ident = TRUE)
# 
# head(EC_subset[[]])
# 
# EC_subset <- ScaleData(EC_subset, features = rownames(EC_subset))
# 
# EC_subset <- RunPCA(EC_subset, features = c(s.genes, g2m.genes))
# 
# DimPlot(EC_subset) # Cells do not separate by cell-cycle phases.
# 
# FeaturePlot(EC_subset,
#             feature = c("MKI67", "TOP2A", "CDK1", "HMGB2", "PCNA"),
#             label = TRUE) # Cells do not separate by cell-cycle phases.

# Cluster visualization
Idents(object = EC_subset) <- "seurat_clusters"

DimPlot(EC_subset)

DimPlot(
  EC_subset,
  group.by = "seurat_clusters",
  split.by = "orig.ident",
  ncol = 2,
  label = FALSE
)

# Cluster classification using marker genes 
FeaturePlot(
  EC_subset,
  feature = c("PECAM1", "ACKR1", "SELE", "VCAM1", "ICAM1",
              "CXCL12", "BMP4", "CLDN5", "PLVAP", "CA4", "HEY1"),
  label = TRUE
)

# Rename clusters using cell-type names
EC_subset[["cell_types"]] <- NA

subsets <- list(
  "Arterial" = 3,
  "Capillary" = 0,
  "Venous 1" = 2,
  "Venous 2" = 1,
  "Lymphatic" = 4
)

for (list_name in names(subsets)) {
  EC_subset$cell_types[which(EC_subset@active.ident %in% subsets[[list_name]])] <- list_name
}

# Export cell-types
Idents(object = EC_subset) <- "cell_types"

# Verify
DimPlot(EC_subset)

# Reorder clusters
Idents(EC_subset) <- "cell_types"
levels <- c("Arterial",
            "Capillary",
            "Venous 1",
            "Venous 2",
            "Lymphatic")
EC_subset@active.ident <-
  factor(x = EC_subset@active.ident, levels = levels)

DimPlot(EC_subset)

# Save the Seurat object containing the EC subclusters
SaveH5Seurat(
  EC_subset,
  here::here(
    "Thymic_EC_scRNA_reanalysis",
    "Bautista_et_al_Nature_Communications_2021",
    "02_processed_data",
    "Bautista_EC_subclusters"
  ),
  overwrite = TRUE
)