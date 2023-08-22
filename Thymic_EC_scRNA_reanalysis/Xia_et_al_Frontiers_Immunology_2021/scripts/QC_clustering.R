## Goal 1: scRNA-seq quality control 
## Goal 2: Clustering analysis to subset thymic endothelial cell (EC) cluster

## Data taken from the following publication:

# Analysis of published scRNA-seq dataset
# Xia, Huan, et al. "Thymic Egress Is Regulated by T Cell-Derived LTβR Signal 
# and via Distinct Thymic Portal Endothelial Cells." 
# Frontiers in Immunology 12 (2021): 707404.

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(SeuratWrappers)
library(gprofiler2)

# Import datasets
EC <-
  ReadMtx(
    mtx = here::here(
      "scRNA_data_analysis",
      "Xia_Frontiers_Immunology_2021",
      "raw_data",
      "GSM5324978_thymic_portal_endothelial_cells_matrix.mtx.gz"
    ),
    features = here::here(
      "scRNA_data_analysis",
      "Xia_Frontiers_Immunology_2021",
      "raw_data",
      "GSM5324978_thymic_portal_endothelial_cells_genes.tsv.gz"
    ),
    cells = here::here(
      "scRNA_data_analysis",
      "Xia_Frontiers_Immunology_2021",
      "raw_data",
      "GSM5324978_thymic_portal_endothelial_cells_barcodes.tsv.gz"
    )
  )

# Create Seurat object
EC_data <-
  CreateSeuratObject(
    counts = EC,
    project = "Thymic_EC",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )

# Quality control
# Calculate mitochondrial ratio
EC_data$mitoRatio <-
  PercentageFeatureSet(object = EC_data, pattern = "^mt-")

# Filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  EC_data,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

EC_QC <-
  subset(
    x = EC_data,
    subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 &
      nCount_RNA > 0 &
      nCount_RNA < 20000 & mitoRatio < 10
  )

VlnPlot(
  EC_data,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

# -----------------------------------------------------------------------------
# clustering analysis

# Data normalization
EC_QC <-
  NormalizeData(EC_QC,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

# Identify highly variable features
EC_QC <- FindVariableFeatures(EC_QC, nfeatures = 2000)

# Scale data (linear transformation)
EC_QC <-
  ScaleData(EC_QC) # Only perform scaling on previously identified variable features

# # Check for cell-cycle bias
# # Convert cell-cycle human gene list to mouse gene list
# s.genes <- gorth(cc.genes$s.genes,
#                  source_organism = "hsapiens",
#                  target_organism = "mmusculus")$ortholog_name
# 
# g2m.genes <- gorth(cc.genes$g2m.genes,
#                    source_organism = "hsapiens",
#                    target_organism = "mmusculus")$ortholog_name
# 
# # Assign cell-cycle scores
# EC_data <-
#   CellCycleScoring(
#     EC_data,
#     s.features = s.genes,
#     g2m.features = g2m.genes,
#     set.ident = TRUE
#   )
# 
# EC_data <- RunPCA(EC_data)
# 
# # Check for differences due to cell cycle phase
# DimPlot(
#   EC_data,
#   reduction = "pca",
#   group.by = "Phase",
#   split.by = "Phase"
# )
# 
# DimPlot(
#   EC_data,
#   reduction = "pca"
# )
# 
# # Check for differences due to cell cycle phase
# FeaturePlot(
#   EC_data,
#   features = c("S.Score", "G2M.Score"),
#   reduction = "pca",
#   dims = 1:2,
#   split.by = "orig.ident"
# )

# # Scale data (linear transformation)
# EC_QC <-
#   ScaleData(EC_QC,
#             vars.to.regress = c("S.Score", "G2M.Score"),
#             features = rownames(EC_QC))

# Linear dimensional reduction
EC_QC <-
  RunPCA(EC_QC) # Use variable from FindVariableFeatures above

# Determine the dimensionality of the dataset
ElbowPlot(EC_QC, ndims = 50)

# Cluster the cells
EC_QC <- FindNeighbors(EC_QC, dims = 1:20)
EC_QC <- FindClusters(EC_QC, resolution = 0.3)

# Run non-linear dimensional reduction (UMAP)
EC_QC <- RunUMAP(EC_QC, dims = 1:20)

# Cluster visualization
DimPlot(EC_QC, reduction = "umap", label = TRUE) # Plot UMAP

# Remove predicted doublets from the data
# Fetal_19weeks
nExp <- round(ncol(EC_QC) * 0.039)
set.seed(2020)

EC_QC <-
  doubletFinder_v3(
    EC_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

EC_QC_singlet <-
  EC_QC[, EC_QC@meta.data[, "DF.classifications_0.25_0.09_169"] == "Singlet"]

DimPlot(EC_QC_singlet, reduction = "umap", label = TRUE) # Plot UMAP

# Cluster classification using marker genes

# EC markers
FeaturePlot(EC_QC_singlet, features = c("Pecam1", "Ackr1", "Bmp4")) 

# Arterial EC markers
FeaturePlot(EC_QC_singlet, features = c("Hey1", "Cxcl12", "Gja5"))

# Capillary EC markers
FeaturePlot(EC_QC_singlet, features = c("Car4", "Rgcc")) 

# Venous EC markers
FeaturePlot(EC_QC_singlet, features = c("Selp", "Sele")) 

# Tip EC markers (Apln+ ECs)
FeaturePlot(EC_QC_singlet, features = c("Apln", "Pgf")) 

# Smooth Muscle Cell (SMC) markers
FeaturePlot(EC_QC_singlet, features = c("Myl9", "Vim", "Pdgfrb")) 

#  Mesenchymal / Fibroblast markers
FeaturePlot(EC_QC_singlet, features = c("Pdgfra", "Dcn")) 

# B- and T-cell contaminants 
FeaturePlot(EC_QC_singlet, features = c("Cd3g", "Cd74")) 

# Subset to only keep the EC clusters
EC_subset <-
  subset(
    EC_QC_singlet,
    idents = c(0,1,3),
    invert = FALSE
  )

# Verify
DimPlot(EC_subset, reduction = "umap", label = TRUE) 

# Rename clusters using cell-type names
EC_subset[["cell_types"]] <- NA

subsets <- list(
  "Arterial" = 3,
  "Capillary" = 0,
  "Venous" = 1
)

for (list_name in names(subsets)) {
  EC_subset$cell_types[which(EC_subset@active.ident %in% subsets[[list_name]])] <- list_name
}

# Export cell-types
Idents(object = EC_subset) <- "cell_types"

# Verify
DimPlot(EC_subset, reduction = "umap", label = TRUE) 

# Save the Seurat object containing the EC clusters
SaveH5Seurat(
  EC_subset,
  here::here(
    "scRNA_data_analysis",
    "Xia_Frontiers_Immunology_2021",
    "processed_data",
    "Xia_EC_clusters"
  ),
  overwrite = TRUE
)