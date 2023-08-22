## Goal 1: scRNA-seq quality control 

## Data taken from the following publication:

# 1. Wells, Kristen L., et al. "Combined transient ablation and 
# single-cell RNA-sequencing reveals the development of 
# medullary thymic epithelial cells." Elife 9 (2020).

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)

# Import datasets
data_wk10 <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "scRNA_data_integration_Wells_and_Michelson",
      "01_raw_data",
      "Wells_et_al",
      "filtered_feature_bc_matrix_controlwk10",
      "matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "scRNA_data_integration_Wells_and_Michelson",
      "01_raw_data",
      "Wells_et_al",
      "filtered_feature_bc_matrix_controlwk10",
      "features.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "scRNA_data_integration_Wells_and_Michelson",
      "01_raw_data",
      "Wells_et_al",
      "filtered_feature_bc_matrix_controlwk10",
      "barcodes.tsv.gz"
    )
  )

data_wk2 <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "scRNA_data_integration_Wells_and_Michelson",
      "01_raw_data",
      "Wells_et_al",
      "filtered_feature_bc_matrix_controlwk2",
      "matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "scRNA_data_integration_Wells_and_Michelson",
      "01_raw_data",
      "Wells_et_al",
      "filtered_feature_bc_matrix_controlwk2",
      "features.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "scRNA_data_integration_Wells_and_Michelson",
      "01_raw_data",
      "Wells_et_al",
      "filtered_feature_bc_matrix_controlwk2",
      "barcodes.tsv.gz"
    )
  )

# Create Seurat object
data_wk2 <-
  CreateSeuratObject(
    counts = data_wk2,
    project = 'Wells',
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )

data_wk10 <-
  CreateSeuratObject(
    counts = data_wk10,
    project = 'Wells',
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )

# Quality control
  # Calculate mitochondrial ratio
data_wk10$mitoRatio <-
  PercentageFeatureSet(object = data_wk10, pattern = "^mt-")
data_wk2$mitoRatio <-
  PercentageFeatureSet(object = data_wk2, pattern = "^mt-")

  # wk10: filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  data_wk10,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

data_wk10_QC <-
  subset(
    x = data_wk10,
    subset = nFeature_RNA > 0 &
      nFeature_RNA < 4000 &
      nCount_RNA > 0 & nCount_RNA < 20000 & mitoRatio < 10
  )

VlnPlot(
  data_wk10_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

  # wk2: filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  data_wk2,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

data_wk2_QC <-
  subset(
    x = data_wk2,
    subset = nFeature_RNA > 0 &
      nFeature_RNA < 5000 &
      nCount_RNA > 0 & nCount_RNA < 30000 & mitoRatio < 10
  )

VlnPlot(
  data_wk2_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

# Remove predicted doublets from the data
  # wk10
nExp <- round(ncol(data_wk10) * 0.016)
set.seed(2020)

data_wk10_QC <-
  doubletFinder_v3(
    data_wk10_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

data_wk10_QC_singlet <-
  data_wk10_QC[, data_wk10_QC@meta.data[, "DF.classifications_0.25_0.09_33"] == "Singlet"]

# wk2
nExp <- round(ncol(data_wk2) * 0.016)
set.seed(2021)

data_wk2_QC <-
  doubletFinder_v3(
    data_wk2_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

data_wk2_QC_singlet <-
  data_wk2_QC[, data_wk2_QC@meta.data[, "DF.classifications_0.25_0.09_29"] == "Singlet"]

# Merge the datasets
combined <-
  merge(
    data_wk10_QC_singlet,
    y = data_wk2_QC_singlet,
    add.cell.ids = c("wk10", "wk2"),
    project = "Wells"
  )

# Quality control on the merged dataset
  # Calculate mitochondrial ratio
combined$mitoRatio <-
  PercentageFeatureSet(object = combined, pattern = "^mt-")

  # Filter merged dataset according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

combined <-
  subset(
    x = combined,
    subset = nFeature_RNA > 0 &
      nFeature_RNA < 4000 & nCount_RNA > 0 & nCount_RNA < 20000
  )

VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

# ** Save the Seurat object containing the merged dataset 
# ** to be used in the script "batch_correction_clustering.R" 

# Save the Seurat object
SaveH5Seurat(
  combined,
  here::here(
    "Thymic_EC_scRNA_reanalysis",
    "scRNA_data_integration_Wells_and_Michelson",
    "02_processed_data",
    "Wells_QC"
  ),
  overwrite = TRUE
)
