## Goal 1: Generate figures

## Author: Jessica Chevallier

library(here)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(scuttle)
library(dittoSeq)
library(pheatmap)
library(escape)
library(msigdbr)

# Load the Seurat object containing the EC subclustering information
EC_data <- LoadH5Seurat(here::here(
  "Thymic_EC_scRNA_RANKL_GST",
  "02_processed_data",
  "GST_RANKL_EC_subclusters.h5seurat"
))

## EC UMAPs
cell_colors <-
  c("#4E67C8FF",
               "#D098EEFF",
               "#7E0021FF",
               "#5ECCF3FF",
               "#579D1CFF",
               "#F14124FF")

# UMAP total         
UMAP_total <- DimPlot(EC_data, cols = cell_colors)

ggsave(
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "UMAP_EC_total.png"
  ),
  width = 8,
  height = 6,
  units = "in",
  device = "png",
  UMAP_total,
  dpi = 500
)

# Dotplot

# Reorder clusters
Idents(EC_data) <- "cell_types"
levels <-
  c("Capillary 1",
    "Capillary 2",
    "Arterial",
    "Choroid plexus like EC",
    "Venous 1",
    "Venous 2")
EC_data@active.ident <-
  factor(x = EC_data@active.ident, levels = levels)

genes <-
  c("Selp",
    "Sele",
    "Vcam1",
    "Icam1",
    "Lrg1",
    "Ackr1",
    "Egr1",
    "Atf3",
    "Esm1",
    "Cd24a",
    "Plaur",
    "Ramp3",
    "Hey1",
    "Gja5",
    "Gja4",
    "Fbln5",
    "Stmn2",
    "Bmx",
    "Rgcc",
    "Lpl",
    "Fabp5",
    "Car4"
  )

dotplot_cell_types <- DotPlot_scCustom(
  seurat_object = EC_data,
  features = genes,
  x_lab_rotate = TRUE,
  colors_use = cell_colors
) + scale_color_viridis(option = "plasma") + scale_color_viridis(option = "plasma") +
  theme(axis.text.x = element_text(size = 10))

ggsave(
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "dotplot.png"
  ),
  width = 8,
  height = 3.5,
  units = "in",
  device = "png",
  dotplot_cell_types,
  dpi = 500
)

# Heatmap of enrichment pathways

chosen_pathways <-
  c(
    "BIOCARTA_AKT_PATHWAY",
    "BIOCARTA_P38MAPK_PATHWAY",
    "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES",
    "REACTOME_ACTIVATED_TAK1_MEDIATES_P38_MAPK_ACTIVATION",
    "BIOCARTA_NFKB_PATHWAY",
    "REACTOME_TNF_RECEPTOR_SUPERFAMILY_TNFSF_MEMBERS_MEDIATING_NON_CANONICAL_NF_KB_PATHWAY",
    "REACTOME_NF_KB_IS_ACTIVATED_AND_SIGNALS_SURVIVAL",
    "REACTOME_TAK1_ACTIVATES_NFKB_BY_PHOSPHORYLATION_AND_ACTIVATION_OF_IKKS_COMPLEX",
    "REACTOME_TRAF6_MEDIATED_NF_KB_ACTIVATION"
  )

# Retrieve genes for each pathway above
gene_sets <- msigdbr(species = "mouse", category = "C2")
genes_per_pathway <- split(gene_sets$gene_symbol, gene_sets$gs_name)
gene_sets <- gene_sets[gene_sets$gs_name %in% chosen_pathways,]
genes_per_pathway <- split(gene_sets$gene_symbol, gene_sets$gs_name)

# Rename clusters
EC_data[["combined_cell_types"]] <- NA

subsets <- list(
  "Capillary" = c("Capillary 1", "Capillary 2"),
  "Choroid plexus like EC" = "Choroid plexus like EC",
  "Arterial" = "Arterial",
  "Venous" = c("Venous 1", "Venous 2")
)

for (list_name in names(subsets)) {
  EC_data$combined_cell_types[which(EC_data@active.ident %in% subsets[[list_name]])] <-
    list_name
}

# Check the cell types are renamed correctly 
Idents(EC_data) <- "combined_cell_types"
DimPlot(EC_data) 

# Calculate enrichment scores per cell per cluster per condition
ES <- enrichIt(
  obj = EC_data,
  gene.sets = genes_per_pathway,
  method = "ssGSEA",
  groups = 1000,
  cores = 4,
  min.size = 10
)

# Add the info as a metadata slot in the Seurat object
EC_data <- AddMetaData(EC_data, ES)
EC_data@meta.data$active.idents <- EC_data@active.ident

# Prepare the data for the heatmap
EC_data$cluster_cnd <-
  paste0(EC_data$combined_cell_types, " ", EC_data$sample)


Idents(EC_data) <- EC_data$cluster_cnd # Change active ident

meta <- EC_data[["cluster_cnd"]]
meta <- merge(meta, ES, by = "row.names")

heatmap <- meta[, c("cluster_cnd", colnames(ES))]
melted <- reshape2::melt(heatmap, id.vars = c("cluster_cnd"))
medianvalues <- melted %>%
  group_by(cluster_cnd, variable) %>%
  summarise(median(value))

matrix <- reshape2::dcast(medianvalues, cluster_cnd ~ variable)
rownames(matrix) <- matrix[, 1]
matrix <- matrix[, -1]

# Reorder the cell types in the heatmap
col_order <- c(
  "Venous GST",
  "Venous RANKL",
  "Choroid plexus like EC GST",
  "Choroid plexus like EC RANKL",
  "Arterial GST",
  "Arterial RANKL",
  "Capillary GST",
  "Capillary RANKL"
)

# Order rows
matrix_ordered <- matrix %>%
  arrange(factor(rownames(matrix), levels = col_order))

# Order columns
matrix_ordered <- matrix_ordered[, chosen_pathways]

# Build heatmap annotations
x <- c("GST", "RANKL")
samp_annot <- data.frame(sample = rep(x, times = 4))
rownames(samp_annot) <- rownames(matrix_ordered)

rownames(samp_annot) <- rownames(matrix_ordered)
ann_colors = list(sample = c("GST" = "#ADDD8E", "RANKL" = "#006837"))

# Color scale
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

# Save the heatmap 
pdf(
  file = here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "enriched_pathways.pdf"
  ),
  width = 28,
  height = 16
)

pheatmap(
  t(matrix_ordered),
  color = col,
  scale = "row",
  fontsize = 25,
  cluster_rows = F,
  cluster_cols = F, 
  clustering_distance_rows = "none",
  angle_col = 90, 
  annotation_col = samp_annot, 
  annotation_colors = ann_colors,
  border_color = NA 
)

dev.off()

# Heatmap HEV (High endothelial venules): clusters venous 1 and 2 per condition

HEV_gene_list <- c("Glycam1",
                   "Serpina1b",
                   "Enpp2",
                   "Apoe",
                   "C1s1",
                   "Clu",
                   "Ubd",
                   "Lifr",
                   "Il2rg",
                   "Il6st",
                   "Ackr1",
                   "Csf2rb",
                   "Serpinb1a",
                   "Cd34")

# Prepare the data 
EC_data$cluster_cnd <-
  paste0(EC_data$cell_types, " ", EC_data$sample)

Idents(EC_data) <- EC_data$cluster_cnd

# Subset to only keep the venous 1 and venous 2 clusters
Idents(EC_data) <- "cell_types"
venous_subset <- subset(EC_data, idents = c("Venous 1", "Venous 2"))

# Convert to a single-cell experiment
venous_sce <- as.SingleCellExperiment(venous_subset)

# Calculate the mean marker expression per condition for all cells
celltype_mean <- aggregateAcrossCells(venous_sce,
                                      ids = venous_sce$cluster_cnd,
                                      statistics = "mean")

# Save the heatmap 
colors <-
  c("#D098EEFF",
               "#F14124FF",
               "#ADDD8E",
               "#006837")

dev.off()
               
pdf(
  file = here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "heatmap_HEV.pdf"
  ),
  width = 6,
  height = 4
)

dittoHeatmap(
  celltype_mean,
  genes = HEV_gene_list,
  cluster_cols = FALSE, 
  scaled.to.max = TRUE, 
  order.by = "cell_types",
  heatmap.colors.max.scaled = magma(100),
  annot.by = c("cell_types", "sample"),
  annot.colors = colors,
  cluster_rows = FALSE, 
  clustering_distance_rows = "euclidean"
)

dev.off()

# Heatmap functional molecule list 
functional_molecule_list <- c("Selp",
                              "Sele",
                              "Vcam1",
                              "Bst1",
                              "Ltbr",
                              "Dll1",
                              "Kitl",
                              "Cxcl9",
                              "Cldn5")

# Save the heatmap               
pdf(
  file = here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "heatmap_functional_molecule.pdf"
  ),
  width = 6,
  height = 4
)

dittoHeatmap(
  celltype_mean,
  genes = functional_molecule_list,
  cluster_cols = FALSE, 
  scaled.to.max = TRUE, 
  order.by = "cell_types",
  heatmap.colors.max.scaled = magma(100),
  annot.by = c("cell_types", "sample"),
  annot.colors = colors,
  cluster_rows = FALSE, 
  clustering_distance_rows = "euclidean"
)

dev.off()

# Heatmap angiogenesis gene list
angiogenesis_gene_list <- c("Tap1",
                            "Cd74",
                            "Slco2a1",
                            "Lima1",
                            "Cd47",
                            "Alas1",
                            "Bmpr2",
                            "Gbp6")
# Prepare the data
# Convert to a single-cell experiment
EC_sce <- as.SingleCellExperiment(EC_data)

# Calculate the mean marker expression per condition for all cells
celltype_mean <- aggregateAcrossCells(EC_sce,
                                      ids = EC_sce$cluster_cnd,
                                      statistics = "mean")

cell_colors <-
  c("#5ECCF3FF",
               "#4E67C8FF",
               "#7E0021FF",
               "#579D1CFF",
               "#D098EEFF",
               "#F14124FF",
               "#ADDD8E",
               "#006837"
  )

# Save the heatmap               
pdf(
  file = here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "heatmap_angiogenesis.pdf"
  ),
  width = 8,
  height = 4
)

dittoHeatmap(
  celltype_mean,
  genes = angiogenesis_gene_list,
  cluster_cols = FALSE, 
  scaled.to.max = TRUE, 
  order.by = "cell_types",
  heatmap.colors.max.scaled = magma(100),
  annot.by = c("cell_types", "sample"),
  annot.colors = cell_colors,
  cluster_rows = FALSE, 
  clustering_distance_rows = "euclidean"
)

dev.off()

# RANK (Tnfrsf11a) density plot

genes <- c("Tnfrsf11a")

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
      "Thymic_EC_scRNA_RANKL_GST",
      "04_figures",
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
