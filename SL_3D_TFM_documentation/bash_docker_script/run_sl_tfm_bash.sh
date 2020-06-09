#!/bin/bash

#~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
#
#Bash script to initialize the Dockerized FEniCS system developed for this
#package.
#
#This script needs to be set with executable permissions
#(run "chmod +x ./run_sl_tfm.sh" in the shell) and in the same folder as
#the displacment datafile.
#
#This must be run from somewhere in "C:/Users/" in Windows environments
# (a FEniCS on Docker limitation)
#
#--- INPUTS ---
#  $@ : The name of the Python script to run (typically "sl_tfm_call.py")
#
# June, 2019; Alex Landauer
# Franck Lab, Brown Univerisity and University of Wisc - Madison
echo "Running FEniCS"
#create a temporary instant cashe
docker volume create --name instant-cache > /dev/null>&1

#now run the Docker container, using that cache; the -v command shares the
#current directory into the Docker at the shared folder; the -w command sets the
#cd of in the Docker machine; the package is fetched from my repo at quay.io,
#and the parameters from running this shell script (i.e. script to run) are passed
for script_file in ./sl_tfm_call_*.py
do
echo "Running FEniCS for: $script_file"
#docker run --rm -v instant-cache:/home/fenics/.instant -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/alandauer/sltfm_dev "python3 python3 $script_file"
echo "FEniCS complete for: $script_file"
done
echo "FEniCS compete"
