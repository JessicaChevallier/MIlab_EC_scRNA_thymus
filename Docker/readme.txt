This image contains:

 - Rstudio server 4.2.0
 - R 4.2.2

# ######################
     BUILD THE IMAGE
# ######################

# cd to the directory containing the docker file

sudo docker build -t scrna_data_analysis .

# ######################
    RUN THE CONTAINER 
# ######################

# 1. Set a working directory variable. 
# Your directory should be where you have chosen to clone the github repository.

export WORKING_DIR=/home/chevallier/Desktop/projects/MIlab

# 2. Start the container from the previously built image.
# Before running the command replace <PASSWORD> with a password of your choice that will be used to login to the Rstudio session.

sudo docker run --rm --name cont_scrna_data_analysis -d -p 8888:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) scrna_data_analysis

# ######################
    CONNECT TO RSTUDIO
# ######################

In an Internet browser, type as a url : http://127.0.0.1:8888 
Use the name of the user session your are working with and your chosen password to login. 

# ######################
    STOP THE CONTAINER
# ######################

sudo docker stop cont_scrna_data_analysis

# ######################
          INFO
# ######################

# Check to see if containers are running
sudo docker ps

# Check to see if images remain
sudo docker images
