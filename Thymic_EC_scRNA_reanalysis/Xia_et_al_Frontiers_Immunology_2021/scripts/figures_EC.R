## Goal 1: Generate figures

## Author: Jessica Chevallier

# Load packages
library(here)
library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(Nebulosa)

# Load the Seurat object containing the EC clustering information 
# from the Xia dataset
EC_data <- LoadH5Seurat(here::here(
  "scRNA_data_analysis",
  "Xia_Frontiers_Immunology_2021",
  "processed_data",
  "Xia_EC_clusters.h5seurat")
)

# Verify
DimPlot(EC_data)

## EC UMAP figures
# Change the cluster order

my_levels <- c("Venous", "Capillary", "Arterial")
levels(EC_data) <- my_levels

cell_colors <-
  c("#F14124FF",
               "#4E67C8FF",
               "#990099FF")
               
# UMAP 
UMAP_plot <- DimPlot(EC_data, cols = cell_colors) +
  theme(legend.text = element_text(size=16),
        legend.key.size = unit(0.75, 'cm')) # change legend text font size

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Xia_Frontiers_Immunology_2021",
    "figures",
    "UMAP.png"
  ),
  width = 8,
  height = 5,
  units = "in",
  device = "png",
  UMAP_plot,
  dpi = 500
)

# Dotplot
genes <-
  c("Selp",
    "Sele",
    "Vcam1",
    "Icam1",
    "Lrg1",
    "Ackr1",
    "Plaur",
    "Ramp3",
    "Rgcc",
    "Lpl",
    "Fabp5",
    "Car4",
    "Hey1",
    "Gja5",
    "Gja4",
    "Fbln5",
    "Stmn2",
    "Bmx"
  )

# Re order cell-types
Idents(EC_data) <- "cell_types"

levels <-
  c("Arterial",
    "Capillary",
    "Venous")

EC_data@active.ident <-
  factor(x = EC_data@active.ident, levels = levels)

dotplot_cell_types <- DotPlot_scCustom(
  seurat_object = EC_data,
  features = genes,
  x_lab_rotate = TRUE,
  colors_use = cell_colors
) + scale_color_viridis(option = "plasma") + scale_color_viridis(option = "plasma") +
  theme(axis.text.x = element_text(size = 10))

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Xia_Frontiers_Immunology_2021",
    "figures",
    "dotplot.png"
  ),
  width = 8,
  height = 3.5,
  units = "in",
  device = "png",
  dotplot_cell_types,
  dpi = 500
)

# Feature plots
genes <- c(
  "Tnfrsf11a",
  "Selp"
)

feature_plot_colors = c(
  "1" = "#CFD8DC",
  "2" = "#3A97FF",
  "3" = "#8816A7",
  "4" = "black"
)

for (gene in genes) {
  plot <- FeaturePlot(EC_data, features = gene, order = TRUE) +
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
      "Xia_Frontiers_Immunology_2021",
      "figures",
      paste0(gene, "_feature_EC.png")
    ),
    width = 6,
    height = 5,
    units = "in",
    device = "png",
    dpi = 500
  )
  dev.off()
}

# Density plots
for (gene in genes) {
  plot <- plot_density(
    EC_data,
    features = gene,
    joint = FALSE,
    pal = "plasma",
    slot = "data"
  )
  print(plot)
  ggsave(
    file = here::here(
      "scRNA_data_analysis",
      "Xia_Frontiers_Immunology_2021",
      "figures",
      paste0(gene, "_density_EC.png")
    ),
    width = 6,
    height = 5,
    units = "in",
    device = "png",
    dpi = 500
  )
  dev.off()
}
