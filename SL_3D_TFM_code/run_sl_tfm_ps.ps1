#~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
#
#PowerShell script to initialize the Dockerized FEniCS system developed for this
#package (Docker Desktop version).
#
#This script needs to be set with executable permissions
#(run "Set-ExecutionPolicy -ExecutionPolicy Unrestricted -Scope process" in the
#shell) and in the same folder as the displacement datafiles.
#
#
#--- INPUTS ---
# *None*
#
# Sept, 2019; Alex Landauer
# Franck Lab, Brown Univerisity and University of Wisc - Madison

#create a temporary instant cashe
(docker volume create --name instant-cache) 2> $null
#now run the Docker container, using that cache; the -v command shares the
#current directory into the Docker at the shared folder; the -w command sets the
#cd of in the Docker machine; the package is fetched from my repo at quay.io,
#and the parameters from running this shell script (i.e. script to run) are passed
$files = Get-ChildItem -Name ./sl_tfm_call_*.py
ForEach($file in $files) {
echo "Running FEniCS for : $file"
docker run --rm -v instant-cache:/home/fenics/.instant -v ${pwd}:/home/fenics/shared:z -w /home/fenics/shared quay.io/alandauer/sltfm_dev "python3 $file"
echo "FEniCS complete for : $file"
}
