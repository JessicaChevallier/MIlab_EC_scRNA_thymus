## Goal 1: Compute module scores using published gene sets

## Publications used:

# 1: Glasner, Ariella, et al. "Conserved transcriptional connectivity of 
# regulatory T cells in the tumor microenvironment informs new combination 
# cancer therapy strategies." Nature Immunology (2023): 1-16.

# 2: Camara, Abdouramane, et al. "Lymph node mesenchymal and endothelial 
# stromal cells cooperate via the RANK-RANKL cytokine axis to shape 
# the sinusoidal macrophage niche." Immunity 50.6 (2019): 1467-1481.

## Authors: 
# 1. Jessica Chevallier 
# 2. Lionel Spinelli: research engineer at the CIML 
# (Computational Biology, Biostatistics & Modeling (CB2M) platform)

# Load packages
library(here)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(rstatix)

# Load the data

# Seurat object containing the EC clusters
EC_data <- LoadH5Seurat(here::here(
  "Thymic_EC_scRNA_RANKL_GST",
  "02_processed_data",
  "GST_RANKL_EC_subclusters.h5seurat"
))

# angiogenesis gene set
angiogenesis_gene_set <- read.csv(
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "03_references",
    "gene_factor_Glasner_2023.csv"
  ),
  sep = ';'
)

# Response to RANK (RANK-deficient lymphatic EC) gene set
RANK_response_gene_set <- read.csv(
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "03_references",
    "down_LEC_genes_Camara_2019.csv"
  ),
  sep = ','
)

# Compute the module scores for each gene set

# GENE SET 1 (angiogenesis gene set)

# Select the column gene_endo_factor8 which corresponds to the 
# angiogenesis gene set
gene_endo_factor8 <- as.vector(angiogenesis_gene_set$gene_endo_factor8)

# Keep common genes
total_genes <- list(gene_endo_factor8, rownames(EC_data))
common_genes <- Reduce(f = intersect, x = total_genes)

# Compute the Seurat module score
EC_data <-
    AddModuleScore(EC_data,
                   features = list(common_genes),
                   name = "angiogenesis_score")

# Plot the module scores by cluster and condition
EC_venous <- subset(EC_data, idents = c("Venous 1", "Venous 2"))

feature <- "angiogenesis_score1"
pathway_df = EC_venous@meta.data[, c("cell_types", "sample", feature)]
pathway_df$sample = factor(pathway_df$sample, levels = c("GST", "RANKL"))

sample_color <- c("#ADDD8E", "#006837")


vln_plot_angiogenesis <-  ggplot(pathway_df,
                                 aes_string(
                                   x = "cell_types",
                                   y = feature,
                                   fill = "sample",
                                   color = "sample"
                                 )) +
  geom_violin(position = "dodge") +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "black",
    position = position_dodge(width = 0.9)
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    size = 0.5,
    width = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = sample_color) +
  scale_color_manual(values = sample_color) +
  theme_minimal() + theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Save the violin plot summarizing the angiogenesis score 
# for venous 1 and venous 2 per condition 
ggsave(
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "module_score_vln_plot_angiogenesis.png"
  ),
  width = 6,
  height = 4,
  units = "in",
  device = "png",
  vln_plot_angiogenesis,
  dpi = 500
)

# Compute the module score comparison statistic by cluster, 
# between conditions using a wilcoxon test

pathway_stats_df <- data.frame()
for (clusterid in levels(Idents(EC_venous))) {
  data_cluster = EC_venous@meta.data[which(EC_venous@meta.data$cell_types == clusterid), c("cell_types", "sample", feature)]
  stats_df = wilcox_test(data = data_cluster, formula = as.formula(paste(feature, "~ sample")))
  stats_df$cell_types = clusterid
  pathway_stats_df = rbind(pathway_stats_df,
                           stats_df)
}

# Adjust the p-value for multiple testing 
# using the Benjamini-Hochberg (BH) correction method

pathway_stats_df$FDR = p.adjust(pathway_stats_df$p, method = "BH")
pathway_stats_df$signif = sapply(pathway_stats_df$FDR, function(fdr){if(fdr <= 0.05) "*" else ""})
pathway_stats_df = pathway_stats_df[, c("cell_types",
                                        "group1",
                                        "group2",
                                        "n1",
                                        "n2",
                                        "p",
                                        "FDR",
                                        "signif",
                                        ".y.")]


# GENE SET 2 (RANK_response_gene_set)

# Select the column gene_endo_factor8 which corresponds to the 

# Keep common genes
RANK_response_gene_set$X <- gsub('\\s+', '', RANK_response_gene_set$X) # Remove white space from the gene column
total_genes <- list(RANK_response_gene_set$X, rownames(EC_data))
common_genes <- Reduce(f = intersect, x = total_genes)

# Compute the Seurat module score
EC_data <-
  AddModuleScore(EC_data,
                 features = list(common_genes),
                 name = "RANK_response_score")

# Plot the module scores by cluster and condition
EC_venous <- subset(EC_data, idents = c("Venous 1", "Venous 2"))

feature <- "RANK_response_score1"
pathway_df = EC_venous@meta.data[, c("cell_types", "sample", feature)]
pathway_df$sample = factor(pathway_df$sample, levels = c("GST", "RANKL"))

sample_color <- c("#ADDD8E", "#006837")

vln_plot_RANK_response <-  ggplot(pathway_df,
                                 aes_string(
                                   x = "cell_types",
                                   y = feature,
                                   fill = "sample",
                                   color = "sample"
                                 )) +
  geom_violin(position = "dodge") +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "black",
    position = position_dodge(width = 0.9)
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    color = "black",
    size = 0.5,
    width = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = sample_color) +
  scale_color_manual(values = sample_color) +
  theme_minimal() + theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Save the violin plot summarizing the angiogenesis score 
# for venous 1 and venous 2 per condition 
ggsave(
  here::here(
    "Thymic_EC_scRNA_RANKL_GST",
    "04_figures",
    "module_score_vln_plot_RANK_response.png"
  ),
  width = 6,
  height = 4,
  units = "in",
  device = "png",
  vln_plot_RANK_response,
  dpi = 500
)

# Compute the module score comparison statistic by cluster, 
# between conditions using a wilcoxon test

pathway_stats_df <- data.frame()
for (clusterid in levels(Idents(EC_venous))) {
  data_cluster = EC_venous@meta.data[which(EC_venous@meta.data$cell_types == clusterid), c("cell_types", "sample", feature)]
  stats_df = wilcox_test(data = data_cluster, formula = as.formula(paste(feature, "~ sample")))
  stats_df$cell_types = clusterid
  pathway_stats_df = rbind(pathway_stats_df,
                           stats_df)
}

# Adjust the p-value for multiple testing 
# using the Benjamini-Hochberg (BH) correction method

pathway_stats_df$FDR = p.adjust(pathway_stats_df$p, method = "BH")
pathway_stats_df$signif = sapply(pathway_stats_df$FDR, function(fdr){if(fdr <= 0.05) "*" else ""})
pathway_stats_df = pathway_stats_df[, c("cell_types",
                                        "group1",
                                        "group2",
                                        "n1",
                                        "n2",
                                        "p",
                                        "FDR",
                                        "signif",
                                        ".y.")]
