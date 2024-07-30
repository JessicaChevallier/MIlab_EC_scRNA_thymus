## Overview of each script

### Note: All scripts are ordered numerically and should be run as so.

**00_Wells_QC.R:** Perform quality control on the scRNA-seq data from Wells KL., et al. The Seurat object generated will be save in the folder `Thymic_EC_scRNA_reanalysis/scRNA_data_integration_Wells_and_Michelson/02_processed_data` and be used in the `02_batch_correction_clustering.R` script.

**01_Michelson_QC.R:** Perform quality control on the scRNA-seq data from Michelson DA., et al. The Seurat object generated will be save in the same folder mentioned above and be used in the `02_batch_correction_clustering.R` script.

**02_batch_correction_clustering.R:** Perform batch correction and cluster the cells to identify TEC subsets. The Seurat object containing the TEC subclusters will be save in the same folder mentioned above and be used in the `03_figures_TEC.R` script.

**03_figures_TEC.R:** Use this script to generate **fig. S1D,E**. 
