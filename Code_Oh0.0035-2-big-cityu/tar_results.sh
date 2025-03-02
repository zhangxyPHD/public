#!/bin/bash 
#SBATCH --partition=cn
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=120:00:00
#SBATCH --account=ehpc-reg-2023r03-178
#SBATCH --qos=ehpc-reg-2023r03-178
#SBATCH --job-name=tar_Results
#SBATCH -o out-tar_Results.txt
#SBATCH -e error-tar_Results.txt

module load gmp/6
module load openmpi/4/gcc/4.1.5
module load gcc/10/10.2.0

# Parameter arrays
tsnap=0.01
Ldomain=8.0
tmax=10.0
Bos=(0.0)
epsilons=(0.001)
MAXlevels=(11)
Wes=("1" "5" "10")
Ohs=("0.0035" "0.01")
Js=(0 0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4)
MAXlevel1=13

for MAXlevel in "${MAXlevels[@]}"; do
  for We in "${Wes[@]}"; do
    for J in "${Js[@]}"; do
      for Oh in "${Ohs[@]}"; do
        for Bo in "${Bos[@]}"; do
          for epsilon in "${epsilons[@]}"; do
            {
              FILENAME="Bo${Bo}-We${We}-J${J}-Oh${Oh}-MAXlevel${MAXlevel}-epsilon${epsilon}"
              # Skip if result already exists
              if [ -e "./Results_Running/video_$FILENAME.mp4" ]; then
                echo "$FILENAME already exists."
                tar -cvf $FILENAME.tar.gz $FILENAME > /dev/null 2>&1
                rm -rf $FILENAME
                continue
              fi
            }&
          done
        done
      done
    done
  done
done
wait