## Goal 1: Batch correction using the Fastmnn Seurat wrapper 
## Goal 3: Clustering analysis to identify thymic epithelial cell (TEC) subsets

## Data taken from the following publication:

# 1. Michelson, Daniel A., et al. "Thymic epithelial cells co-opt 
# lineage-defining transcription factors to eliminate autoreactive T cells." 
# Cell 185.14 (2022): 2542-2558.

# 1. Wells, Kristen L., et al. "Combined transient ablation and 
# single-cell RNA-sequencing reveals the development of 
# medullary thymic epithelial cells." Elife 9 (2020).

# Tutorials used:
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2020/scRNAseq/html/featSelec.nb.html
# https://github.com/MarioniLab/FurtherMNN2018
# http://bioconductor.org/packages/devel/bioc/vignettes/batchelor/inst/doc/correction.html
# https://nbisweden.github.io/single-cell_sib_scilifelab/session-batch_correction/batch_correction.html
# https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html
# https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.md

## Author: Jessica Chevallier

# Load packages
library(here)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

# Import datasets
Mathis <- LoadH5Seurat(here::here(
  "scRNA_data_analysis",
  "scRNA-seq_data_integration",
  "raw_data",
  "Mathis.h5seurat"
))

Wells <- LoadH5Seurat(here::here(
  "scRNA_data_analysis",
  "scRNA-seq_data_integration",
  "raw_data",
  "Wells.h5seurat"
))

# Set the seed
set.seed(10403)

# # Save the seed
# oldseed <- .Random.seed
# # Restore the seed
# .Random.seed <- oldseed

# Merge the datasets
TEC_combined <-
  merge(
    Wells,
    y = Mathis,
    add.cell.ids = c("Wells", "Mathis"),
    project = "TEC"
  )

# Normalize
TEC_combined <-
  NormalizeData(TEC_combined,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

# Find variable features
TEC_combined <-
  FindVariableFeatures(TEC_combined,
                       dispersion.function = "vst",
                       nfeatures = 2000)


# Run fastMNN
TEC_combined <-
  RunFastMNN(object.list = SplitObject(TEC_combined, split.by = "orig.ident"),
             features = 2000)

# UMAP embedding and cluster the cells
TEC_combined <-
  RunUMAP(TEC_combined, reduction = "mnn", dims = 1:35)

TEC_combined <-
  FindNeighbors(TEC_combined, reduction = "mnn", dims = 1:35)

TEC_combined <- FindClusters(TEC_combined, resolution = 0.8)

# Check the integration 
DimPlot(
  TEC_combined,
  group.by = c("orig.ident"),
  ncol = 2,
  label = TRUE
)

# Cluster visualization
DimPlot(
  TEC_combined,
  ncol = 2,
  label = TRUE
)

# -----------------------------------------------------------------------------
# Identify cell-types using marker genes of our populations of interest

# Intertypical TECs
FeaturePlot(TEC_combined,
            feature = c("Ccl21a", "Krt5", "Ifitm3"),
            label = TRUE)

# cTECs and Pdpn+ TECs
FeaturePlot(
  TEC_combined,
  feature = c("Prss16", "Psmb11", "Ccl25", "Tbata", "Pdpn"),
  label = TRUE
)

# Aire+ mTECs and Tuft-like TECs
FeaturePlot(
  TEC_combined,
  feature = c("Aire", "Cd52", "Spink5", "Avil", "Trpm5", "Dclk1", "Fezf2"),
  label = TRUE
)

# Proliferating TECs
FeaturePlot(TEC_combined,
            feature = c("Hmgb2", "Hmgn2"),
            label = TRUE)

# Microfold TEC
FeaturePlot(TEC_combined,
            feature = c("Spib", "Ccl9", "Ccl20"),
            label = TRUE) 

# Endothelial TEC
FeaturePlot(
  TEC_combined,
  feature = c("Sod3", "Dpt", "Cd177", "Car8", "Dnajc12", "Ceacam10"),
  label = TRUE
) 

# Keratinocyte 
FeaturePlot(TEC_combined,
            feature = c("Grhl1", "Gata3", "Pou2f3"),
            label = TRUE) 

# Ciliated
FeaturePlot(
  TEC_combined,
  feature = c(
    "Foxj1",
    "Trp73",
    "Dynlrb2",
    "Dnah12",
    "Ccdc153",
    "Ccdc113",
    "Mlf1",
    "Lztfl1"
  ),
  label = TRUE
) 

# Muscle
FeaturePlot(TEC_combined, features = "Myog", label = TRUE)

# T-cell contamination
FeaturePlot(TEC_combined,
            feature = c("Cd3e", "Cd3g", "Cd3d", "Satb1"),
            label = TRUE) 

# -----------------------------------------------------------------------------
# subcluster populations of interest

  # Cluster 10
TEC_combined <-
  FindSubCluster(
    TEC_combined,
    "10",
    subcluster.name = "Ciliated",
    resolution = 1,
    graph.name = "RNA_snn"
  )

Idents(object = TEC_combined) <-
  "Ciliated" # Doing this transfers the subcluster information to the active ident

DimPlot(
  TEC_combined,
  ncol = 2,
  label = TRUE
)

FeaturePlot(
  TEC_combined,
  features = c("Foxj1", "Trp73", "Dynlrb2", "Dnah12"),
  label = TRUE
)

  # Cluster 6
TEC_combined <-
  FindSubCluster(
    TEC_combined,
    "6",
    subcluster.name = "Proliferating_TEC",
    resolution = 0.5,
    graph.name = "RNA_snn"
  )

Idents(object = TEC_combined) <-
  "Proliferating_TEC" # Doing this transfers the subcluster information to the active ident

DimPlot(
  TEC_combined,
  ncol = 2,
  label = TRUE
)

FeaturePlot(TEC_combined,
            feature = c("Hmgb2", "Hmgn2", "Aire"),
            label = TRUE)

  # Cluster 5
TEC_combined <-
  FindSubCluster(
    TEC_combined,
    "5",
    subcluster.name = "Keratinocyte",
    resolution = 0.5,
    graph.name = "RNA_snn"
  )

Idents(object = TEC_combined) <-
  "Keratinocyte" # Doing this transfers the subcluster information to the active ident

DimPlot(
  TEC_combined,
  ncol = 2,
  label = TRUE
)

FeaturePlot(TEC_combined,
            feature = c("Grhl1", "Gata3", "Pou2f3"),
            label = TRUE)

# Rename clusters for downstream analysis
TEC_combined <- RenameIdents(object = TEC_combined, `10_0` = "100")
TEC_combined <- RenameIdents(object = TEC_combined, `10_1` = "101")
TEC_combined <- RenameIdents(object = TEC_combined, `10_2` = "102")
TEC_combined <- RenameIdents(object = TEC_combined, `10_3` = "103")
TEC_combined <- RenameIdents(object = TEC_combined, `6_0` = "60")
TEC_combined <- RenameIdents(object = TEC_combined, `6_1` = "61")
TEC_combined <- RenameIdents(object = TEC_combined, `6_2` = "62")
TEC_combined <- RenameIdents(object = TEC_combined, `6_3` = "63")
TEC_combined <- RenameIdents(object = TEC_combined, `5_0` = "50")
TEC_combined <- RenameIdents(object = TEC_combined, `5_1` = "51")
TEC_combined <- RenameIdents(object = TEC_combined, `5_2` = "52")

# Cluster visualization
DimPlot(
  TEC_combined,
  ncol = 2,
  label = TRUE
)

# Rename clusters using cell-type names
TEC_combined[["cell_types"]] <- NA

subsets <- list(
  "Intertypical_TEC" = c(0, 12),
  "Aire_plus_mTEC" = c(3, 62, 63),
  "Proliferating_TEC" = c(60,61),
  "cTEC" = 8,
  "Pdpn_plus_TEC" = 9,
  "Muscle_mTEC" = 14,
  "Microfold_mTEC" = 11,
  "Ciliated_mTEC" = 102,
  "Contamination_T_cell" = 13,
  "Neuroendocrine_mTEC" = 7,
  "Keratinocyte_mTEC" = c(51,52),
  "Late_Aire_mTEC" = 1,
  "Tuft_like_mTEC" = c(2, 4),
  "Post_Aire_mTEC" = c(50, 100, 101, 103)
)
             
for (list_name in names(subsets)) {
  TEC_combined$cell_types[which(TEC_combined@active.ident %in% subsets[[list_name]])] <-
    list_name
}

# Verify
DimPlot(
  TEC_combined,
  reduction = "umap",
  label = TRUE,
  group.by = "cell_types"
) 

# Export cell-types
Idents(object = TEC_combined) <- "cell_types"

# Remove contaminating T-cells
TEC_combined <-
  subset(
    TEC_combined,
    idents = c("Contamination_T_cell"),
    invert = TRUE
  )

# Verify
DimPlot(TEC_combined)

# Save the Seurat object containing the TEC clusters
SaveH5Seurat(
  TEC_combined,
  here::here(
    "scRNA_data_analysis",
    "scRNA-seq_data_integration",
    "processed_data",
    "Integrated_Wells_Michelson"
  ),
  overwrite = TRUE
)
