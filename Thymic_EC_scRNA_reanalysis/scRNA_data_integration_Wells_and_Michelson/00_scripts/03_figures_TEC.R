## Goal 1: Generate figures

## Author: Jessica Chevallier

# Load packages
library(here)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(RColorBrewer)

# Load the Seurat object containing the TEC clustering information 
# from the Bautista dataset
TEC_data <- LoadH5Seurat(here::here(
  "scRNA_data_analysis",
  "scRNA-seq_data_integration",
  "processed_data",
  "Integrated_Wells_Michelson.h5seurat"
))

# Verify
DimPlot(TEC_data)

# Prepare the data

# Rename cell-types
TEC_data[["cell_types"]] <- NA

subsets <- list(
  "Intertypical TEC" = "Intertypical_TEC",
  "Aire+ mTEC" = "Aire_plus_mTEC",
  "Proliferating TEC" = "Proliferating_TEC",
  "cTEC" = "cTEC",
  "Pdpn+ TEC" = "Pdpn_plus_TEC",
  "Muscle mTEC" = "Muscle_mTEC",
  "Microfold mTEC" = "Microfold_mTEC",
  "Ciliated mTEC" = "Ciliated_mTEC",
  "Neuroendocrine mTEC" = "Neuroendocrine_mTEC",
  "Keratinocyte mTEC" = "Keratinocyte_mTEC",
  "Late-Aire mTEC" = "Late_Aire_mTEC",
  "Tuft-like mTEC" = "Tuft_like_mTEC",
  "Post-Aire mTEC" = "Post_Aire_mTEC"
)

for (list_name in names(subsets)) {
  TEC_data$cell_types[which(TEC_data@active.ident %in% subsets[[list_name]])] <- list_name
}

# Export cell-types
Idents(object = TEC_data) <- "cell_types"

# Reorder the cell-types
cell_levels <-
  c(
    "cTEC",
    "Pdpn+ TEC",
    "Intertypical TEC",
    "Proliferating TEC",
    "Aire+ mTEC",
    "Late-Aire mTEC",
    "Post-Aire mTEC",
    "Keratinocyte mTEC",
    "Ciliated mTEC",
    "Neuroendocrine mTEC",
    "Microfold mTEC",
    "Muscle mTEC",
    "Tuft-like mTEC"
  )

levels(TEC_data) <- cell_levels

cell_colors <-
  c("#8B8682",
             "#4E67C8FF",
             "#F14124FF",
             "#D098EEFF",
             "#5DCEAFFF",
             "#990099FF",
             "#FF8021FF",
             "#579D1CFF",
             "#7E0021FF",
             "#FCC66DFF",
             "#00718BFF",
             "#FC719EFF",
             "#5ECCF3FF"
  )

# UMAP
UMAP_plot <- DimPlot(TEC_data, cols = cell_colors) +
  theme(legend.text = element_text(size=16),
        legend.key.size = unit(0.75, 'cm')) # change legend text font size

ggsave(
  here::here(
    "scRNA_data_analysis",
    "scRNA-seq_data_integration",
    "figures",
    "UMAP.png"
  ),
  UMAP_plot,
  width = 9,
  height = 7,
  units = "in"
)

# Feature plot
genes <- c("Tnfrsf11a")

feature_plot_colors = c(
  "1" = "#CFD8DC",
  "2" = "#3A97FF",
  "3" = "#8816A7",
  "4" = "black"
)

for (gene in genes) {
  plot <- FeaturePlot(TEC_data, features = gene, order = TRUE) +
    ggplot2::theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.key.size = unit(0.75, 'cm')
    ) & 
    scale_colour_gradientn(colours = feature_plot_colors) 
  print(plot)
  ggsave(
    file = here::here(
      "scRNA_data_analysis",
      "scRNA-seq_data_integration",
      "figures",
      paste0(gene, "_feature.png")
    ),
    width = 7,
    height = 6,
    units = "in",
    device = "png",
    dpi = 500
  )
  dev.off()
}