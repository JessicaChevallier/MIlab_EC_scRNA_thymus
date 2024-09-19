# RANKL treatment rejuvenates thymic function and improves T-cell immune responses during aging
## Article information

### Title:
RANKL treatment rejuvenates thymic function and improves T-cell immune responses during aging

### Authors:
Jérémy C Santamaria $$1$$, Jessica Chevallier $$1$$, Léa Dutour $$2$$, Amandine Picart $${3,4}$$, Camille Kergaravat $$2$$, Agata Cieslak $${5,6}$$, Mourad Amrane $$7$$, Renaud Vincentelli $$8$$, Denis Puthier $$11$$, Emmanuel Clave $$2$$, Arnauld Sergé $$9$$, Martine Cohen-Solal $${3,4}$$, Antoine Toubert $${2,10}$$, Magali Irla $${1,*}$$

$$1$$ Centre d'Immunologie de Marseille-Luminy, CIML, CNRS, INSERM, Aix-Marseille Université, Marseille, Turing Centre for Living Systems, Marseille, France  
$$2$$ Université de Paris Cité, Institut de Recherche Saint Louis, EMiLy, INSERM UMRS 1160, F-75010, Paris, France  
$$3$$ Université de Paris Cité, INSERM, UMR-S 1132 BIOSCAR, F-75010 Paris, France  
$$4$$ Departement de Rhumatologie, Hôpital Lariboisière, AP-HP, F-75010 Paris, France  
$$5$$ Laboratoire d’Onco-Hematologie, Hôpital Necker Enfants Malades, AP-HP, F-75015 Paris, France  
$$6$$ Université Paris Cité, CNRS, INSERM U1151, Institut Necker Enfants Malades (INEM), Paris, France  
$$7$$ Service de Chirurgie Cardiovasculaire, Hôpital Européen Georges Pompidou, AP-HP, F-75015 Paris, France  
$$8$$ Architecture et Fonction des Macromolécules Biologiques (AFMB), UMR 7257 CNRS-Aix-Marseille Université, Marseille, France  
$$9$$ Laboratoire Adhesion & Inflammation, LAI, CNRS, INSERM, Aix Marseille Université, Turing Centre for Living Systems, Marseille, France  
$$10$$ Laboratoire d’Immunologie et d’Histocompatibilité, Hôpital Saint‐Louis, AP‐HP, F-75010 Paris France  
 
\* For correspondence: Magali.Irla@inserm.fr

### Abstract:
<p align="justify"> Age-related thymic involution is one of the major causes of immunosenescence, characterized by a reduced T-cell production, resulting in an increased susceptibility to cancers, infections, autoimmunity and a reduced vaccine efficacy. Here, we identify that the RANK/RANKL axis in the thymus is altered during aging. Using a novel conditional transgenic mouse model, we demonstrate that endothelial cells (EC) depend on RANK for their cellularity and functional maturation. Thus, we uncover that the decreased RANKL
availability during aging results in a decline in cellularity and function of both EC and thymic epithelial cells (TEC), leading to thymic involution. We then show that, whereas RANKL neutralization in young mice mimics thymic involution, RANKL cytokine treatment in aged mice restores thymic architecture, EC and TEC cellularity and functional properties. Consequently, RANKL improves T-cell progenitor homing to the thymus and boosts T-cell production. Importantly, this cascade of events results in peripheral T-cell renewal and effective anti-tumor and vaccine responses. Furthermore, we provide the proof-of-concept that RANKL stimulates EC and TEC in human thymic organo-cultures. Overall, our findings identify this cytokine-based treatment as a potent therapeutic strategy that rejuvenates thymic function and improves T-cell immunity in the elderly. </p>

***
## Repository goal 
This github repository contains the scripts and dockerfile necessary to reproduce the analyses / figures. 

## Description of the repository structure
3 main folders:
```
Docker
Thymic_EC_scRNA_RANKL_GST
Thymic_EC_scRNA_reanalysis
```
```Docker```: contains the dockerfile and instructions on how to build the docker image and run the container.  
```Thymic_EC_scRNA_RANKL_GST```: reproduce the analyses / figures using the **RANKL- GST- scRNA-seq dataset**.  
```Thymic_EC_scRNA_reanalysis```: reproduce the analyses / figures using **publicly available thymic EC datasets**. 

The following subfolders are present:
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
  Xia_et_al_Frontiers_Immunology_2021
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
```00_scripts```: subfolder containing all scripts to reproduce the analyses / figures.   
```01_raw_data```: subfolder to place the raw data into.  
```02_processed_data```: subfolder where the h5Seurat files containing the clustering information will be saved.   
```03_references```: subfolder containing csv files necessary for the gene set enrichment analysis.  
```03_figures``` / ```04_figures```: subfolder where the figures will be saved.   

**NOTE:** Each subfolder contains a README.txt describing the subfolder content. The README.txt in all the ```01_raw_data``` subfolders describes which dataset(s) to download from the [GEO database](https://www.ncbi.nlm.nih.gov/geo/) to redo the analyses. 
***
## Steps to run the analysis 

**Step 1:**    

Clone the github repository into your chosen folder. A folder called "MIlab_EC_scRNA_thymus" will be created.     

**Step 2:**  

Set the variable WORKING_DIR using the path to the "MIlab_EC_scRNA_thymus" folder as your value. 
```
export WORKING_DIR=/home/chevallier/Desktop/projects/MIlab/MIlab_EC_scRNA_thymus
```
**Step 3:**    

Download raw data from the [GEO database](https://www.ncbi.nlm.nih.gov/geo/) and place it in the **```01_raw_data```** subfolders. The **README.txt** files in each **```01_raw_data```** subfolder tells you which raw data needs to be downloaded.  

Raw data generated in this study can be downloaded **here**:  
As a quick summary, we utilized the publicly available datasets below. 

| Author(s) | Year | Dataset title | Datatset URL | Database and Identifier  
| :---: | :---: | :---: | :---: | :---:
| Michelson DA., et al. | 2022 | Thymic epithelial cells co-opt lineage-defining transcription factors to eliminate autoreactive T cells| https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194253| NCBI Gene Expression Omnibus, GSE194253
| Xia, Huan., et al. | 2021 | T cell derived LTR signal regulates thymic egress via distinct thymic portal endothelial cells| https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174732 | NCBI Gene Expression Omnibus, GSE174732
| Bautista JL., et al. | 2021 | Single-cell RNA sequencing of human thymic samples| https://www.ncbi.nlm.nih.gov/geo/query/acc.cgiacc=GSE147520 | NCBI Gene Expression Omnibus, GSE147520 
|Wells KL., et al. | 2020 | Single cell sequencing defines a branched progenitor population of stable medullary thymic epithelial cells|  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137699 | NCBI Gene Expression Omnibus, GSE137699 

**Step 4:**    

Install Docker: https://docs.docker.com/engine/install/  
Build an image from the Dockerfile by running the following in your terminal.  
```
cd $WORKING_DIR/Docker
sudo docker build -t scrna_data_analysis .
```
**Step 5:**    

Run the container using the previously built image.   
Before running the command replace \<PASSWORD\> with a password of your choice that will be used to login to the Rstudio server.  

```
sudo docker run --rm --name cont_scrna_data_analysis -d -p 8888:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) scrna_data_analysis
```
**Step 6:**  

In an Internet browser, use the url : http://127.0.0.1:8888 to connect to the Rstudio server.   
Use the name of the user session your are working with and your chosen password to login. 

**Step 7:**  

Create a new Rstudio project using the path to the "MIlab_EC_scRNA_thymus" folder as your Existing Directory. In doing so, a **"MIlab_EC_scRNA_thymus.Rproj"** file will be created.   

```
In RStudio: File > New Project > Existing Directory > Browse > "MIlab_EC_scRNA_thymus" > Select Folder > Create Project
```

**NOTE:** It's **important** to create a new Rstudio project using the cloned git repository because we will be using the ```here``` package to identify the top-level directory (.Rproj file) and build paths relative to it throughout the analysis. This prevents us from using absolute paths and makes switching from one operating system to another easier.  

You can learn more about the ```here``` package in this post: https://software.cqls.oregonstate.edu/tips/posts/r-tips-here-package/#:~:text=The%20here%20package%20builds%20your,top%2Dlevel%20project%20directory

**Step 8:**  

We're almost ready to run the analyses!  
In the Rstudio session, open the "MIlab_EC_scRNA_thymus.Rproj" file and run the following in the console.  

```
library(here)
```
You should see something similar to this, listing the path to the "MIlab_EC_scRNA_thymus" folder on your computer.  

```
here() starts at /home/chevallier/Desktop/projects/MIlab/MIlab_EC_scRNA_thymus  
```
**Step 9:**  

Everything is now set for you to run the analyses. The container is connected to the WORKING_DIR you defined and you have access to all the files in the Rstudio environment. Go to the folder of your choice and start exploring.  

**NOTE**: All scripts in the ```00_scripts``` subfolders are ordered numerically and should be run as so.    
