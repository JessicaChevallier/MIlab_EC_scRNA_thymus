## Overview of each script

### Note: All scripts are ordered numerically and should be run as so.

**00_QC_initial_clustering.R:** Perform quality control, cluster the cells, identify thymic EC and remove contaminants. The Seurat object generated will be save in the folder `Thymic_EC_scRNA_RANKL_GST/02_processed_data` and be used in the `01_subclustering.R` script.

**01_subclustering.R:** Subcluster the data to identify thymic EC subsets. We will generate a new Seurat object containing the EC subclusters that will be save in the folder `Thymic_EC_scRNA_RANKL_GST/02_processed_data` and be used in the remaining scripts.

**02_gene_set_enrichment.R:** After running the first two scripts, use this script to generate **fig. 5E**. The gene sets used to compute the module scores are located in the folder `Thymic_EC_scRNA_RANKL_GST/03_references`. 

**03_figures_EC.R:** Use this script to generate **fig. 5B,C,D,F** and **fig. S6G,H**. 
