#!/bin/bash 
#SBATCH -p batch        # partition name
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=440G      # Memory allocation
#SBATCH -t 70:00:00    # Time limit (hh:mm:ss)
#SBATCH --job-name=We10
#SBATCH -o out-We10.txt
#SBATCH -e error-We10.txt

filename="Droplet_Rough_Surfaces_2D_Axi-We10"
echo "Start $filename"
comsol batch -np 128 -inputfile ${filename}.mph -outputfile ${filename}.mph -batchlog log_${filename}
