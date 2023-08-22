## Goal 1: Demultiplexing using hashtag oligos (HTOs)
## Goal 2: Quality control 

## Data taken from the following publication:

# 1. Michelson, Daniel A., et al. "Thymic epithelial cells co-opt 
# lineage-defining transcription factors to eliminate autoreactive T cells." 
# Cell 185.14 (2022): 2542-2558.

## Tutorial used: https://satijalab.org/seurat/archive/v3.1/hashing_vignette.html

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)

# Import datasets
expression_matrix <-
  ReadMtx(
    mtx = here::here(
      "scRNA_data_analysis",
      "Michelson_Cell_2022",
      "raw_data",
      "GSM5831744_adult_perinate_gex_matrix.mtx.gz"
    ),
    features = here::here(
      "scRNA_data_analysis",
      "Michelson_Cell_2022",
      "raw_data",
      "GSM5831744_adult_perinate_gex_features.tsv.gz"
    ),
    cells = here::here(
      "scRNA_data_analysis",
      "Michelson_Cell_2022",
      "raw_data",
      "GSM5831744_adult_perinate_gex_barcodes.tsv.gz"
    )
  )

htos <-
  ReadMtx(
    mtx = here::here(
      "scRNA_data_analysis",
      "Michelson_Cell_2022",
      "raw_data",
      "GSM5831745_adult_perinate_hash_matrix.mtx.gz"
    ),
    features = here::here(
      "scRNA_data_analysis",
      "Michelson_Cell_2022",
      "raw_data",
      "GSM5831745_adult_perinate_hash_features.tsv.gz"
    ),
    cells = here::here(
      "scRNA_data_analysis",
      "Michelson_Cell_2022",
      "raw_data",
      "GSM5831745_adult_perinate_hash_barcodes.tsv.gz"
    ), 
    feature.column = 1, cell.column = 1
  )

colnames(htos) <- paste0(colnames(htos), "-1")

# Create Seurat object
hashtag <-
  CreateSeuratObject(
    counts = expression_matrix,
    min.cells = 3,
    min.features = 100,
    project = "Mathis"
  )

# Select cell barcodes detected by both RNA and HTOs
joint <- intersect(colnames(hashtag), colnames(htos))
hashtag <- hashtag[, joint]
htos <- as.matrix(htos[, joint])

# Change row names to Hash 1, Hash 2 and so on
rownames(htos) <-
  gsub(pattern = "-.*", replacement = "", rownames(htos))

# Add HTO data as a new assay independent from RNA
hashtag[["HTO"]] <- CreateAssayObject(counts = htos)

# Normalize HTO data using centered log-ratio (CLR) transformation
hashtag <-
  NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")

# Demultiplex cells 
hashtag <- MULTIseqDemux(
  hashtag,
  assay = "HTO",
  autoThresh = TRUE,
  maxiter = 10,
  qrange = seq(from = 0.1, to = 0.9, by = 0.05),
  verbose = TRUE
)

# Visualize demultiplexing results 
table(hashtag$MULTI_ID)

# Create ridgeplot 
dev.off()

png(here::here(
  "scRNA_data_analysis",
  "Michelson_Cell_2022",
  "figures",
  "ridgeplot.png"
), width = 12, height = 9, units = "in", res=300)

RidgePlot(hashtag, features = rownames(hashtag[["HTO"]]), assay = "HTO")

dev.off()

# Group cells based on the max HTO signal
Idents(hashtag) <- "MULTI_ID"

# Remove HTOs with low cell numbers
singlet <-
  subset(
    hashtag,
    idents = c("Doublet", "Negative", "Hash1", "Hash2", "Hash9"),
    invert = T
  ) 

# -----------------------------------------------------------------------------
# Quality control

  # Calculate mitochondrial ratio
singlet$mitoRatio <-
  PercentageFeatureSet(object = singlet,
                       pattern = "^mt-",
                       assay = "RNA")

  # Filter dataset according to nFeature_RNA, nCount_RNA, and mitoRatio
VlnPlot(
  singlet,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

singlet <-
  subset(
    x = singlet,
    subset = nFeature_RNA > 0 &
      nFeature_RNA < 5000 &
      nCount_RNA > 0 &
      nCount_RNA < 20000 &
      mitoRatio < 10 
  )  

VlnPlot(
  singlet,
  features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
  ncol = 3
)

# ** Save the Seurat object containing the merged dataset 
# ** to be used in the script "batch_correction_clustering.R"
# ** Doublets were not removed from the data 

# Save the Seurat object
SaveH5Seurat(
  singlet,
  here::here(
    "scRNA_data_analysis",
    "Mathis_Cell_2022",
    "processed_data",
    "Mathis"
  ),
  overwrite = TRUE
)