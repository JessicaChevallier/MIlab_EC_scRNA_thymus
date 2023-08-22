## Goal 1: scRNA-seq quality control 
## Goal 2: Batch correction using the Fastmnn Seurat wrapper 
## Goal 3: Clustering analysis to subset human thymic endothelial cell (EC) cluster

## Data taken from the following publication:

# 1: Bautista, Jhoanne L., et al. "Single-cell transcriptional profiling 
# of human thymic stroma uncovers novel cellular heterogeneity 
# in the thymic medulla." Nature Communications 12.1 (2021): 1096.

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(SeuratWrappers)

# Import datasets

# NOTE: replicate Postnatal_10days_1 was excluded because few cells were recovered

Fetal_19weeks <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466780_F_19wks_matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466780_F_19wks_genes.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466780_F_19wks_barcodes.tsv.gz"
    )
  )

Postnatal_10m_2 <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466785_P_10m_2_matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466785_P_10m_2_genes.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466785_P_10m_2_barcodes.tsv.gz"
    )
  )

Adult_25years <-
  ReadMtx(
    mtx = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466786_A_25y_matrix.mtx.gz"
    ),
    features = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466786_A_25y_features.tsv.gz"
    ),
    cells = here::here(
      "Thymic_EC_scRNA_reanalysis",
      "Bautista_et_al_Nature_Communications_2021",
      "01_raw_data",
      "GSE147520_RAW",
      "GSM4466786_A_25y_barcodes.tsv.gz"
    )
  )

# Create Seurat object for each dataset and perform quality control
Fetal_19weeks_QC <-
  CreateSeuratObject(
    Fetal_19weeks,
    min.cells = 3,
    min.genes = 200,
    project = "F_19wk"
  )

Postnatal_10m_2_QC <-
  CreateSeuratObject(
    Postnatal_10m_2,
    min.cells = 3,
    min.genes = 200,
    project = "PN_10m_2"
  )

Adult_25years_QC <-
  CreateSeuratObject(
    Adult_25years,
    min.cells = 3,
    min.genes = 200,
    project = "Adult_25y"
  )

# Quality control: Fetal_19weeks_QC

# Calculate mitochondrial ratio
Fetal_19weeks_QC$mitoRatio <-
  PercentageFeatureSet(object = Fetal_19weeks_QC, pattern = "^MT-")

# Filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  Fetal_19weeks_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

Fetal_19weeks_QC <-
  subset(
    x = Fetal_19weeks_QC,
    subset = nFeature_RNA > 500 & nFeature_RNA < 3000 &
      nCount_RNA > 0 & nCount_RNA < 10000 & mitoRatio < 10
  )

VlnPlot(
  Fetal_19weeks_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
) 

# Quality control: Postnatal_10days_2

# Calculate mitochondrial ratio
Postnatal_10m_2_QC$mitoRatio <-
  PercentageFeatureSet(object = Postnatal_10m_2_QC, pattern = "^MT-")

# Filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  Postnatal_10m_2_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

Postnatal_10m_2_QC <-
  subset(
    x = Postnatal_10m_2_QC,
    subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 &
      nCount_RNA > 0 & nCount_RNA < 15000 & mitoRatio < 10
  )

VlnPlot(
  Postnatal_10m_2_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
) 

# Quality control: Adult_25years

# Calculate mitochondrial ratio
Adult_25years_QC$mitoRatio <-
  PercentageFeatureSet(object = Adult_25years_QC, pattern = "^MT-")

# Filter according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  Adult_25years_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

Adult_25years_QC <-
  subset(
    x = Adult_25years_QC,
    subset = nFeature_RNA > 250 & nFeature_RNA < 4000 &
      nCount_RNA > 0 & nCount_RNA < 15000 & mitoRatio < 10
  )

VlnPlot(
  Adult_25years_QC,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
) 

# Remove predicted doublets from the data
# Fetal_19weeks
nExp <- round(ncol(Fetal_19weeks) * 0.08)
set.seed(2020)

Fetal_19weeks_QC <-
  doubletFinder_v3(
    Fetal_19weeks_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

Fetal_19weeks_QC_singlet <-
  Fetal_19weeks_QC[, Fetal_19weeks_QC@meta.data[, "DF.classifications_0.25_0.09_1868"] == "Singlet"]

# Postnatal_10m_2
nExp <- round(ncol(Postnatal_10m_2) * 0.048)
set.seed(2022)

Postnatal_10m_2_QC <-
  doubletFinder_v3(
    Postnatal_10m_2_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

Postnatal_10m_2_QC_singlet <-
  Postnatal_10m_2_QC[, Postnatal_10m_2_QC@meta.data[, "DF.classifications_0.25_0.09_264"] == "Singlet"]

# Adult_25years
nExp <- round(ncol(Adult_25years) * 0.08)
set.seed(2023)

Adult_25years_QC <-
  doubletFinder_v3(
    Adult_25years_QC,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp,
    PCs = 1:35,
    sct = TRUE
  )

Adult_25years_QC_singlet <-
  Adult_25years_QC[, Adult_25years_QC@meta.data[, "DF.classifications_0.25_0.09_1902"] == "Singlet"]

# Merge the datasets
combined_stroma <-
  merge(
    Fetal_19weeks_QC_singlet,
    y = c(Postnatal_10m_2_QC_singlet, Adult_25years_QC_singlet),
    add.cell.ids = c("F_19wk", "PN_10m_2", "Adult_25y"),
    project = "stroma"
  )

# Save the Seurat object
SaveH5Seurat(
  combined_stroma,
  here::here(
    "Thymic_EC_scRNA_reanalysis",
    "Bautista_et_al_Nature_Communications_2021",
    "02_processed_data",
    "Bautista_combined_stroma"
  ),
  overwrite = TRUE
)

# -----------------------------------------------------------------------------
# perform batch correction using the fastMNN Seurat-wrapper

# Load the data
data <- LoadH5Seurat(here::here(
  "Thymic_EC_scRNA_reanalysis",
  "Bautista_et_al_Nature_Communications_2021",
  "02_processed_data",
  "Bautista_combined_stroma.h5seurat")
)

# Normalize
data <-
  NormalizeData(data,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

# Find variable features
data <-
  FindVariableFeatures(data,
                       dispersion.function = "vst",
                       nfeatures = 2000)

# Run fastMNN
data <-
  RunFastMNN(
    object.list = SplitObject(data, split.by = "orig.ident"),
    features = 2000
  )

# -----------------------------------------------------------------------------
# cluster the cells
data <-
  RunUMAP(data, reduction = "mnn", dims = 1:30)

data <-
  FindNeighbors(data, reduction = "mnn", dims = 1:30)

data <- FindClusters(data, resolution = 0.2)

# Cluster visualization

# Check the integration 
DimPlot(
  data,
  group.by = "orig.ident",
  split.by = "orig.ident",
  ncol = 2,
  label = FALSE
)

# Look at the clusters
DimPlot(
  data,
  label = FALSE
)

# Identify cell-types using marker genes from the article

# Markers specific to the mesenchyme and pericytes
FeaturePlot(data,
            feature = c("PDGFRA", "PDGFRB"),
            label = TRUE)

# Markers specific to the immune cluster
FeaturePlot(data,
            feature = c("PTPRC"),
            label = TRUE)

# Markers specific to the mesothelium cluster
FeaturePlot(data,
            feature = c("MSLN"),
            label = TRUE)

# Markers specific to the epithelium clusters
FeaturePlot(data,
            feature = c("KRT8", "EPCAM"),
            label = TRUE)

# Markers specific to the red blood cell cluster
FeaturePlot(data,
            feature = c("GYPA"),
            label = TRUE)

# Markers specific to endo-1,2,3 clusters
FeaturePlot(data,
            feature = c("PECAM1", "ACKR1", "SELE", "VEGFC", "BMP4"),
            label = TRUE)

# Markers specific to endo-4 cluster
FeaturePlot(data,
            feature = c("LYVE1"),
            label = TRUE)

# Endo-1: top 5 markers (logFC) according to the authors
FeaturePlot(data,
            feature = c("ACKR1", "LCN6", "RAMP3", "NOSTRIN", "PLVAP"),
            label = TRUE)

# Endo-2: top 5 markers (logFC) according to the authors
FeaturePlot(data,
            feature = c("SEMA3G", "GJA5", "SRGN", "GJA4", "PLLP"),
            label = TRUE)

# Endo-3: top 5 markers (logFC) according to the authors
FeaturePlot(data,
            feature = c("APLNR", "COL4A1", "FABP4", "COL15A1", "RBP7"),
            label = TRUE)

# Subset to only keep the EC cluster (based on PECAM1 expression)
EC_subset <-
  subset(
    data,
    idents = 0,
    invert = FALSE
  )

# Verify
DimPlot(EC_subset)
DimPlot(EC_subset, split.by = "orig.ident")

# Save the Seurat object containing just the EC cluster
SaveH5Seurat(
  EC_subset,
  here::here(
    "Thymic_EC_scRNA_reanalysis",
    "Bautista_et_al_Nature_Communications_2021",
    "02_processed_data",
    "EC_subset"
  ),
  overwrite = TRUE
)
