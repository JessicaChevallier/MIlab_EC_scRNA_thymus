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
library(scuttle)
library(viridis)
library(dittoSeq)
library(Nebulosa)
library(RColorBrewer)

# Load the Seurat object containing the EC subclustering information 
# from the Bautista dataset
EC_data <- LoadH5Seurat(here::here(
  "scRNA_data_analysis",
  "Bautista_Nature_Communications_2021",
  "processed_data",
  "Bautista_EC_subclusters.h5seurat"
))

DimPlot(EC_data)

## EC UMAP figures
cell_colors <-
  c("#990099FF",
               "#5DCEAFFF",
               "#4E67C8FF",
               "#F14124FF",
               "#579D1CFF")

# UMAP total
UMAP_plot <- DimPlot(EC_subset, cols = cell_colors) +
  theme(legend.text = element_text(size=16),
        legend.key.size = unit(0.75, 'cm')) # change legend text font size

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Bautista_Nature_Communications_2021",
    "figures",
    "UMAP_EC_total.png"
  ),
  width = 9,
  height = 7,
  units = "in",
  device = "png",
  UMAP_plot,
  dpi = 500
)

# UMAP split by samples
# Export samples (orig.ident) clusters
Idents(object = EC_subset) <- "orig.ident"

UMAP_cluster_per_sample <- DimPlot_scCustom(
  seurat_object = EC_subset,
  group.by = "cell_types",
  split.by = "orig.ident",
  num_columns = 3,
  repel = TRUE,
  colors_use = cell_colors, 
  label = FALSE
) +
  scale_color_discrete(
    name = "Identity",
    labels = c("Fetal: 19 weeks", "Postnatal: 10 months", "Adult: 25 years"),
  )

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Bautista_Nature_Communications_2021",
    "figures",
    "UMAP_EC_per_sample.png"
  ),
  width = 12,
  height = 4,
  units = "in",
  device = "png",
  UMAP_sample,
  dpi = 500
)

# Dotplot according to the cell-types
# Re order cell-types
Idents(EC_data) <- "cell_types"

levels <-
  c("Lymphatic",
    "Venous 2",
    "Venous 1",
    "Capillary",
    "Arterial")

EC_data@active.ident <-
  factor(x = EC_data@active.ident, levels = levels)

genes <-
  c(
    "PECAM1",
    "CLDN5",
    "CDH5",
    "LTBR",
    "DLL4",
    "KITLG",
    "CXCL12",
    "CCL19",
    "CCL21",
    "ICAM1",
    "SELE",
    "VCAM1",
    "SELP",
    "BST1",
    "IL7",
    "BMP4"
  )

dotplot_cell_types <- DotPlot_scCustom(
  seurat_object = EC_data,
  features = genes,
  x_lab_rotate = TRUE
) + scale_color_viridis(option = "plasma") + scale_color_viridis(option = "plasma") +
  theme(axis.text.x = element_text(size = 8))

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Bautista_Nature_Communications_2021",
    "figures",
    "dotplot_cell_types.png"
  ),
  width = 7,
  height = 3.5,
  units = "in",
  device = "png",
  dotplot_cell_types,
  dpi = 500
)

# Dotplot according to the samples
# Rename samples 
Idents(EC_data) <- "orig.ident"

EC_data <- RenameIdents(EC_data, "F_19wk" = "F19",
                          "PN_10d_2" = "P10",
                          "Adult_25y" = "A25")

EC_data$orig.ident <- Idents(EC_data)

# Export samples
Idents(object = EC_data) <- "orig.ident"

# Re order samples
Idents(EC_data) <- "orig.ident"

levels <-
  c("A25",
    "P10",
    "F19")

EC_data@active.ident <-
  factor(x = EC_data@active.ident, levels = levels)

dotplot_samples <- DotPlot_scCustom(
  seurat_object = EC_data,
  features = genes,
  x_lab_rotate = TRUE
) + scale_color_viridis(option = "plasma") + scale_color_viridis(option = "plasma") +
  theme(axis.text.x = element_text(size = 8))

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Bautista_Nature_Communications_2021",
    "figures",
    "dotplot_samples.png"
  ),
  width = 7,
  height = 3.5,
  units = "in",
  device = "png",
  dotplot_samples,
  dpi = 500
)

# Stacked bar chart displaying cell percentage per sample
cell_colors <-
  c("#579D1CFF",
               "#F14124FF",
               "#4E67C8FF",
               "#5DCEAFFF",
               "#990099FF")
               
stacked_barplot <- dittoBarPlot(
  object = EC_data,
  var = "cell_types",
  group.by = "orig.ident",
  color.panel = cell_colors,
  x.reorder = c(2, 3, 1),
  var.labels.reorder = c(3, 5, 4, 2, 1)
) + scale_y_continuous(
  breaks = seq(0, 1, by = 0.20),
  labels = c("0", "20", "40", "60", "80", "100")
)

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Bautista_Nature_Communications_2021",
    "figures",
    "stacked_barplot.png"
  ),
  width = 3,
  height = 4,
  units = "in",
  device = "png",
  stacked_barplot,
  dpi = 500
)

# Density plots using the Nebulosa package
genes <- c(
  "TNFRSF11A",
  "VCAM1",
  "SELP",
  "ICAM1",
  "SELE",
  "CXCL12",
  "BST1"
)

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
      "Bautista_Nature_Communications_2021",
      "figures",
      paste0(gene, "_density.png")
    ),
    width = 9,
    height = 7,
    units = "in",
    device = "png",
    dpi = 500
  )
  dev.off()
}

# DEG analysis and heatmap visualization 
EC_data <-
  ScaleData(EC_data, features = rownames(EC_data@assays[["RNA"]]))

Idents(EC_data) <- "cell_types"

# Change the cluster order
my_levels <- c("Arterial", "Capillary", "Venous 1", "Venous 2", "Lymphatic")
levels(EC_data) <- my_levels

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
  all_markers[all_markers$p_val_adj < 0.05 &
                abs(all_markers$avg_log2FC) > log2(1),]

top_genes <-
  filtered_genes %>% group_by(cluster) %>% top_n(10, avg_log2FC)

avg_exp <-
  AverageExpression(
    EC_data,
    features = top_genes$gene,
    return.seurat = TRUE,
    slot = "data", 
    assays = "RNA"
  )

# Heatmap
cell_colors <-
  c("#990099FF",
               "#5DCEAFFF",
               "#4E67C8FF",
               "#F14124FF",
               "#579D1CFF")
               
heatmap_colors <- rev(brewer.pal(11, "RdBu"))

heatmap_DEG <-
  DoHeatmap(
    avg_exp,
    features = top_genes$gene,
    size = 2,
    label = FALSE,
    draw.lines = FALSE,
    group.colors = cell_colors
  ) +
  ggplot2::scale_fill_gradientn(colors = heatmap_colors, na.value = "white") +
  theme(axis.text.y = element_text(size = 8, face = "bold")) +
  theme(legend.text = element_text(size=10), legend.title=element_text(size=10),
        legend.key.size = unit(0.5, 'cm')) +
  scale_color_discrete(
    name = "Identity",
    labels = c("Arterial", "Capillary", "Venous 1", "Venous 2", "Lymphatic"),
  ) + 
  scale_color_manual(values = cell_colors)

ggsave(
  here::here(
    "scRNA_data_analysis",
    "Bautista_Nature_Communications_2021",
    "figures",
    "heatmap_DEG.png"
  ),
  width = 4,
  height = 6,
  units = "in",
  device = "png",
  heatmap_DEG,
  dpi = 500
)

# 1 Dotplot for each cluster according to the sample type
genes <-
  c(
    "PECAM1",
    "CLDN5",
    "CDH5",
    "LTBR",
    "DLL4",
    "KITLG",
    "CXCL12",
    "CCL19",
    "CCL21",
    "ICAM1",
    "SELE",
    "VCAM1",
    "SELP",
    "BST1",
    "IL7",
    "BMP4"
  )

Idents(object = EC_data) <- "cell_types"
cell_type_list <- c("Arterial", "Capillary", "Venous 1", "Venous 2", "Lymphatic")

for(cell_type in cell_type_list) {
  subset_obj <-
    subset(
      EC_data,
      idents = cell_type,
      invert = FALSE
    )
  Idents(object = subset_obj) <- "orig.ident"
  plot <- DotPlot_scCustom(
      seurat_object = subset_obj,
      features = genes,
      x_lab_rotate = TRUE
    ) + scale_color_viridis(option = "plasma") +
      theme(axis.text.x = element_text(size = 8))
  print(plot)
  ggsave(
    file = here::here(
      "scRNA_data_analysis",
      "Bautista_Nature_Communications_2021",
      "figures",
      paste0(cell_type, "_dotplot.png")
    ),
    width = 9,
    height = 7,
    units = "in",
    device = "png",
    dpi = 500
  )
  dev.off()
}
