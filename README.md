# RANKL treatment rejuvenates thymic function and improves T-cell immune responses during aging
## Article information

### Title:
RANKL treatment rejuvenates thymic function and improves T-cell immune responses during aging

### Authors:
Jérémy C Santamaria 1 Jessica Chevallier 1 , Léa Dutour 2,3 , Amandine Picart 4,5 , Renaud
Vincentelli 6 , Emmanuel Clave 2,3 , Arnauld Sergé 7 , Martine Cohen-Solal 4,5 , Antoine Toubert 2,3 ,
Magali Irla 8,*

1 Centre d'Immunologie de Marseille-Luminy, CIML, CNRS, INSERM, Aix-Marseille Université,
Marseille, Turing Centre for Living Systems, Marseille, France  
2 Université de Paris, Institut de Recherche Saint Louis, EMiLy, Inserm U1160, F-75010, Paris, France  
3 Laboratoire d’Immunologie et d’Histocompatibilité, Hôpital Saint‐Louis, AP‐HP, Paris France  
4 Université de Paris, INSERM, UMR-S 1132 BIOSCAR, F-75010 Paris, France  
5 Rheumatology Department, AP-HP, Lariboisière Hospital, F-75010 Paris, France  
6 Architecture et Fonction des Macromolécules Biologiques (AFMB), UMR 7257 CNRS-Aix-Marseille
Université, Marseille, France  
7 Laboratoire Adhesion &amp; Inflammation, LAI, CNRS, INSERM, Aix Marseille Université, Turing Centre
for Living Systems, Marseille, France  
8 Centre d'Immunologie de Marseille-Luminy, CIML, CNRS, INSERM, Aix-Marseille Université,
Marseille, Turing Centre for Living Systems, Marseille, France  

\* For correspondence: Magali.Irla@inserm.fr

### Abstract:
Age-related thymic involution is one of the major causes of immunosenescence,
characterized by reduced production of naive T cells, leading to increased susceptibility to
cancers, infections, autoimmunity and reduced vaccine efficacy. Here, we reveal that the
RANK/RANKL axis is altered in the thymus during aging, which results in a decline in thymic
epithelial cell (TEC) and endothelial cell (EC) function and, consequently, to thymic
involution. Whereas RANKL neutralization in young mice mimics thymic involution, RANKL
treatment in aged individuals restores the thymic architecture, TEC and EC functional
properties, and thymic T-cell production. This cascade of events results in the renewal of
peripheral T cells and effective anti-tumor and anti-viral immune responses. Furthermore, we
provide the proof-of-concept that RANKL stimulates TEC and EC in human thymic organo-
cultures. Collectively, our findings identify RANKL treatment as a potent therapeutic strategy
to rejuvenate thymic function and improve T-cell immune function in the elderly.

### DOI:

### GEO: 
***
## Repository goal 
This github repository contains the scripts and dockerfile necessary to reproduce the analyses / figures described in the article. 

## Description of the repository structure
3 main folders
```
Docker
Thymic_EC_scRNA_RANKL_GST
Thymic_EC_scRNA_reanalysis
```
Docker: contains the dockerfile and instructions on how to build the docker image and run the container.  
Thymic_EC_scRNA_RANKL_GST: code reproduce the analysis / figures using the **RANKL- GST- scRNA-seq dataset**.  
Thymic_EC_scRNA_reanalysis: code to reproduce the analysis / figures using **publicly available thymic EC datasets**. 

The following subfolders are present
```
THYMIC_EC_scRNA_RANKL_GST
  00_scripts
  01_raw_data
  02_processed_data
  03_references
  04_figures

THYMIC_EC_scRNA_reanalysis
  Bautista_et_al_Nature_Communications_2021
   00_scripts
   01_raw_data
   02_processed_data
   03_figures
  Xia_et_al_Frontiers_Immunology_2019
   00_scripts
   01_raw_data
   02_processed_data
   03_figures
  scRNA_data_integration_Wells_and_Michelson
   00_scripts
   01_raw_data
   02_processed_data
   03_figures
```
00_scripts: subfolder containing all scripts.   
01_raw_data: subfolder to place the raw data into.  
02_processed_data: subfolder where the h5Seurat files containing the clustering information will be saved.   
03_references: subfolder containing csv files necessary for the gene set enrichment analysis.  
03_figures / 04_figures: subfolder where the figures will be saved.   

**NOTE:** Each subfolder contains a README.txt describing the subfolder content. The README.txt in all the 01_raw_data subfolders describes which dataset(s) to download from the GEO database to redo the analysis. 
***
## Steps to run the analysis 

**Step 1:**  
Clone the existing github repository to a new Rstudio project in your chosen folder. In doing so, a folder named "MIlab_EC_scRNA_thymus" will be created along with a "MIlab_EC_scRNA_thymus.Rproj" file. It is recommended to create a new Rstudio project because we will be using the ```here``` package to identify the top-level directory using the .Rproj and build paths relative to it throughout the analysis. 
