## Overview of each script

### Note: All scripts are ordered numerically and should be run as so.

**00_QC_batch_correction.R:** Perform quality control, batch correction and cluster the cells to identify human thymic EC. The Seurat object generated will be saved in the folder `Thymic_EC_scRNA_reanalysis/Bautista_et_al_Nature_Communications_2021/02_processed_data` and be used in the `01_subclustering.R` script.

**01_subclustering.R:** Subcluster the data to identify thymic EC subsets. We will generate a new Seurat object containing the EC subclusters that will be saved in the folder mentioned above and be used in the `02_figures_EC.R` script.

**02_figures_EC.R:** Use this script to generate **fig. 8D,E,F,G,H** and **fig. S12,A,B,C**. 

