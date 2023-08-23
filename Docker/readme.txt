This image contains:

 - Rstudio server 4.2.0
 - R 4.2.2
 
# ############################
    SET A WORKING DIRECTORY
# ############################
 
# Set a working directory variable. 
# Your directory should be where you have chosen to clone the github repository.

export WORKING_DIR=/home/chevallier/Desktop/projects/MIlab

# ############################
        BUILD THE IMAGE
# ############################

# cd to the directory containing the docker file and run the following: 

sudo docker build -t scrna_data_analysis .

# ############################
       RUN THE CONTAINER 
# ############################

# Start the container from the previously built image.
# Before running the command replace <PASSWORD> with a password of your choice that will be used to login to the Rstudio session.

sudo docker run --rm --name cont_scrna_data_analysis -d -p 8888:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) scrna_data_analysis

# ############################
       CONNECT TO RSTUDIO
# ############################

In an Internet browser, type as a url : http://127.0.0.1:8888 
Use the name of the user session your are working with and your chosen password to login. 

# ############################
       STOP THE CONTAINER
# ############################

sudo docker stop cont_scrna_data_analysis

# ############################
            INFO
# ############################

# See what containers are running. 
sudo docker ps

# Get a list of built images. 
sudo docker images
