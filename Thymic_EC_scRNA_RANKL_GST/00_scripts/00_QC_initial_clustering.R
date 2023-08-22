## Goal 1: scRNA-seq quality control 
## Goal 2: Clustering analysis to identify thymic endothelial cell (EC) clusters

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)

# Import datasets
GST_12m <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_RANKL_GST",
      "01_raw_data",
      "GST",
      "matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_RANKL_GST",
      "01_raw_data",
      "GST",
      "features.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_RANKL_GST",
      "01_raw_data",
      "GST",
      "barcodes.tsv.gz"
    )
  )

RANKL_12m <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_RANKL_GST",
      "01_raw_data",
      "RANKL",
      "matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_RANKL_GST",
      "01_raw_data",
      "RANKL",
      "features.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_RANKL_GST",
      "01_raw_data",
      "RANKL",
      "barcodes.tsv.gz"
    )
  )

# Create Seurat object for each dataset and perform quality control
GST_12m_QC <-
  CreateSeuratObject(
    GST_12m,
    min.cells = 3,
    min.genes = 200,
    project = "GST_12m"
  )

RANKL_12m_QC <-
  CreateSeuratObject(
    RANKL_12m,
    min.cells = 3,
    min.genes = 200,
    project = "RANKL_12m"
  )

# Quality control: GST_12m_QC

# Calculate mitochondrial ratio
GST_12m_QC$mitoRatio <-
  PercentageFeatureSet(object = GST_12m_QC, pattern = "^mt-")

# Filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  GST_12m_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), 
  ncol = 3
) 

GST_12m_QC <-
  subset(
    x = GST_12m_QC,
    subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
      nCount_RNA > 0 & nCount_RNA < 15000 & mitoRatio < 10
  )

VlnPlot(
  GST_12m_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
) 

# Quality control: RANKL_12m_QC

# Calculate mitochondrial ratio
RANKL_12m_QC$mitoRatio <-
  PercentageFeatureSet(object = RANKL_12m_QC, pattern = "^mt-")

# Filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  RANKL_12m_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), # mitoRatio = 0 
  ncol = 3
) 

RANKL_12m_QC <-
  subset(
    x = RANKL_12m_QC,
    subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
      nCount_RNA > 0 & nCount_RNA < 15000 & mitoRatio < 10
  )

VlnPlot(
  RANKL_12m_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
) 

# Remove predicted doublets from the data
# GST_12m_QC
nExp <- round(ncol(GST_12m) * 0.04)
set.seed(2020)

GST_12m_QC <-
  doubletFinder_v3(
    GST_12m_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

GST_12m_QC_singlet <-
  GST_12m_QC[, GST_12m_QC@meta.data[, "DF.classifications_0.25_0.09_209"] == "Singlet"]

# Remove predicted doublets from the data
# RANKL_12m_QC
nExp <- round(ncol(RANKL_12m) * 0.04)
set.seed(2020)

RANKL_12m_QC <-
  doubletFinder_v3(
    RANKL_12m_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

RANKL_12m_QC_singlet <-
  RANKL_12m_QC[, RANKL_12m_QC@meta.data[, "DF.classifications_0.25_0.09_204"] == "Singlet"]

# Clustering 

# Set-up a grouping variable
GST_12m_QC_singlet$sample <- "GST"
RANKL_12m_QC_singlet$sample <- "RANKL"

# Merge the datasets
combined_data <-
  merge(
    GST_12m_QC_singlet,
    y = RANKL_12m_QC_singlet,
    add.cell.ids = c("GST_12m", "RANKL_12m"),
    project = "TEC_EC_aging"
  )

# Normalize
combined_data <-
  NormalizeData(combined_data,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

# Find variable features
combined_data <-
  FindVariableFeatures(combined_data,
                       dispersion.function = "vst",
                       nfeatures = 2000)

# Scale the data
combined_data <-
  ScaleData(combined_data, features = rownames(combined_data))

# Linear dimensional reduction
combined_data <- RunPCA(combined_data)

# Choose the number of PCs to be used in the subsequent steps
ElbowPlot(combined_data, ndims = 50)

# UMAP embedding
combined_data <-
  RunUMAP(combined_data, reduction = "pca", dims = 1:15)

# Cluster the cells
combined_data <- FindNeighbors(combined_data, dims = 1:15)
combined_data <- FindClusters(combined_data, resolution = 0.2)

# Cluster visualization
DimPlot(combined_data, label = TRUE)

# Take a look at markers specific to ECs, TECs and splenic B-cells (added as a control)
# Also take a look at T-cell + Fibroblast markers that maybe contamination

# EC markers (Clusters 0,2,3, maybe 9)
FeaturePlot(
  combined_data,
  feature = c("Pecam1", "Ackr1", "Sele", "Bmp4"),
  label = TRUE
)

# TEC markers (Clusters 1,6,7,8, maybe 10)
FeaturePlot(
  combined_data,
  feature = c("Krt8", "Epcam", "Ccl21a", "Aire", "Spink5", "Prss16"),
  label = TRUE
)

# B-cell markers (Cluster 5)
FeaturePlot(combined_data,
            feature = c("Cd79a", "Cd79b"),
            label = TRUE)

# T-cell markers (no T-cell contamination)
FeaturePlot(combined_data,
            feature = c("Dntt", "Cd8b1", "Cd2"),
            label = TRUE)

# Fibroblast markers (Cluster 4: contamination)
FeaturePlot(
  combined_data,
  feature = c("Pdpn", "Cd34", "Apod", "Pdgfra"),
  label = TRUE
)

# Subset to only keep the EC clusters
data_subset <-
  subset(
    combined_data,
    idents = c(0,2,3,9),
    invert = FALSE
  )

# Verify
DimPlot(data_subset, label = TRUE)

# Save the Seurat object containing the EC clusters
SaveH5Seurat(
  data_subset,
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "02_processed_data",
    "GST_RANKL_EC"
  ),
  overwrite = TRUE
)

# # ---------------------------------------------------------------------------
# # Check for potential cell-cycle bias
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes

# Convert cell-cycle human gene list to mouse gene list
# s.genes = gorth(cc.genes$s.genes,
#                 source_organism = "hsapiens",
#                 target_organism = "mmusculus")$ortholog_name
# 
# g2m.genes = gorth(cc.genes$g2m.genes,
#                   source_organism = "hsapiens",
#                   target_organism = "mmusculus")$ortholog_name
# 
# combined_data <-
#   CellCycleScoring(combined_data,
#                    s.features = s.genes,
#                    g2m.features = g2m.genes,
#                    set.ident = TRUE)
# 
# combined_data <- ScaleData(combined_data, features = rownames(EC_subset))
# 
# combined_data <- RunPCA(combined_data, features = c(s.genes, g2m.genes))
# 
# DimPlot(combined_data) # Cells do not separate by cell-cycle phases

# # Check if batch correction is necessary 

# # Normalize
# combined_data <-
#   NormalizeData(combined_data,
#                 normalization.method = "LogNormalize",
#                 scale.factor = 10000)
# 
# 
# # Find variable features
# combined_data <-
#   FindVariableFeatures(combined_data,
#                        dispersion.function = "vst",
#                        nfeatures = 2000)
# 
# # Scale the data
# combined_data <-
#   ScaleData(combined_data, features = rownames(combined_data))
# 
# # Linear dimensional reduction
# combined_data <- RunPCA(combined_data)
# 
# # Choose the number of PCs to be used in the subsequent steps
# ElbowPlot(combined_data, ndims = 50)
# 
# # UMAP embedding
# combined_data <-
#   RunUMAP(combined_data, reduction = "pca", dims = 1:15)
# 
# # Visualization
# umap_total_per_condition <-
#   DimPlot(combined_data, reduction = "umap", group.by = "sample")
# 
# # The data does not split according to the experimental condition (GST versus RANKL),
# # batch correction is not needed