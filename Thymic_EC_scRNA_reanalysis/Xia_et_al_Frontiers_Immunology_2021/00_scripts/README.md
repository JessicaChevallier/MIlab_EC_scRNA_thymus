## Overview of each script

### Note: All scripts are ordered numerically and should be run as so.

**00_QC_clustering.R:** Perform quality control, cluster the cells, identify thymic EC and remove contaminants. The Seurat object generated will be save in the folder `Thymic_EC_scRNA_reanalysis/Xia_et_al_Frontiers_Immunology_2021/02_processed_data` and be used in the `01_figures_EC.R` script.

**01_figures_EC.R:** Use this script to generate **fig. 1F**