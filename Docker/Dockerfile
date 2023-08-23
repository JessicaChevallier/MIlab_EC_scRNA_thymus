# Base image
FROM rocker/tidyverse:4.2.2

MAINTAINER Jessica Chevallier (jesschvllr82@gmail.com)

# Install libraries
RUN apt-get update && apt-get install -y build-essential \
    libhdf5-dev \
    libxt-dev \
    libudunits2-dev \
    libbz2-dev \
    texlive-xetex \
    libgdal-dev \
    apt-utils \
    libglpk40

# Install Various R packages
RUN R -e "install.packages('here', dependencies = T)"
RUN R -e "install.packages('tidyr', dependencies = T)"
RUN R -e "install.packages('Matrix', dependencies = T)"
RUN R -e "install.packages('ggplot2', dependencies = T)"
RUN R -e "install.packages('tibble', dependencies = T)"
RUN R -e "install.packages('devtools', dependencies = T)"
RUN R -e "install.packages('forcats', dependencies = T)"
RUN R -e "install.packages('RColorBrewer', dependencies = T)"
RUN R -e "install.packages('R.utils', dependencies = T)"
RUN R -e "install.packages('gprofiler2', dependencies = T)"
RUN R -e "install.packages('reshape', dependencies = T)"
RUN R -e "install.packages('scCustomize', dependencies = T)"
RUN R -e "install.packages('amap', dependencies = T)"
RUN R -e "install.packages('Hmisc', dependencies = T)"
RUN R -e "install.packages('heatmaply', dependencies = T)"
RUN R -e "install.packages('paletteer', dependencies = T)"
RUN R -e "install.packages('Seurat', version = '4.3.0', dependencies = T)"
RUN R -e "BiocManager::install('dittoSeq', dependencies = T)"
RUN R -e "install.packages('rsvd', dependencies= T)" # Necessary for seurat-wrapper

# Increase timeout for github archive download without "timeout fail"
RUN R -e 'options(download.file.method = "libcurl", url.method = "libcurl", timeout = 10000000000)'

# Install packages using devtools
RUN R -e "devtools::install_github('mojaveazure/seurat-disk')"
RUN R -e "devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')"
RUN R -e "devtools::install_github('powellgenomicslab/Nebulosa')"
RUN R -e "devtools::install_github('jokergoo/ComplexHeatmap')"

# Install 'glmGamPoi' to use advance arguments in the Seurat SCTransform function
RUN R -e "BiocManager::install('glmGamPoi', dependencies = T)"

# Install seurat-wrapper by downloading the package link and installing it manually because "devtools::install_github('satijalab/seurat-wrappers')" does not work
RUN mkdir /tmp/downloads && \ 
	cd /tmp/downloads && \ 
	wget -cO - https://api.github.com/repos/satijalab/seurat-wrappers/tarball/HEAD > satijalab-seurat-wrappers.tar.gz 
	
RUN R -e "install.packages('/tmp/downloads/satijalab-seurat-wrappers.tar.gz', repos = NULL, type = 'source')"

# Install openssh-client to be able to use git in the Rstudio server
RUN apt-get update && apt-get install -y build-essential openssh-client

# Rstudio listens on port 8787 by default
EXPOSE 8787
