#!/bin/bash 
#SBATCH -J We1-2             # Job name which can be changed
#SBATCH -p batch        # partition name
#SBATCH -o out-We1-2.txt    # Output file name
#SBATCH -e error-We1-2.txt    # Error file name
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=440G      # Memory allocation
#SBATCH -t 70:00:00    # Time limit (hh:mm:ss)


module load gcc/10.2.0
module load openmpi/4.1.5_gcc

# Parameter arrays
tsnap=0.01
Ldomain=6.0
Bos=(0.0)
MAXlevels=(9)
Wes=(1)
# Ohs=(0.001 0.002 0.003 0.004 0.006 0.008 0.01 0.02 0.03 0.04) #10
Ohs=(0.06 0.08 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.5) #10
Js=(0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0)  #12
epsilons=(1e-3)
tmax=20
DT=1e-4
CFL=0.01
###################
# Parallel settings
###################
NPARA=128         # Maximum number of concurrent tasks
set_threads_by_maxlevel() {
  case $1 in
    9) echo 10 ;;  # MAXlevel = 9 -> THREADS = 4
    10) echo 16 ;; # MAXlevel = 10 -> THREADS = 5
    11) echo 16 ;; # MAXlevel = 10 -> THREADS = 5
    12) echo 10 ;; # MAXlevel = 10 -> THREADS = 5
    13) echo 10 ;; # MAXlevel = 10 -> THREADS = 5
    *) echo 1 ;;  # Default
  esac
}

mkdir -p Results_Running

############################################
# Function: Check tasks and remove finished
############################################
# This function iterates over the 'tasks' array,
# removes any PID that has exited, and returns
# how many are still running.
check_and_clean_tasks() {
  local still_running=()
  for pid in "${tasks[@]}"; do
    if kill -0 "$pid" 2>/dev/null; then
      # PID is still alive
      still_running+=( "$pid" )
    fi
  done
  tasks=("${still_running[@]}")
  echo "${#tasks[@]}"  # Return how many remain
}

#################################
# Main loop over all parameters
#################################
declare -a tasks=()  # To store PIDs of launched jobs

for MAXlevel in "${MAXlevels[@]}"; do
  for We in "${Wes[@]}"; do
    for J in "${Js[@]}"; do
      for Oh in "${Ohs[@]}"; do
        for Bo in "${Bos[@]}"; do
          for epsilon in "${epsilons[@]}"; do
            FILENAME="Bo${Bo}-We${We}-J${J}-Oh${Oh}-MAXlevel${MAXlevel}-epsilon${epsilon}"
            THREADS=$(set_threads_by_maxlevel $MAXlevel)
            # Skip if result already exists
            if [ -e "./Results_Running/video_$FILENAME.mp4" ]; then
              echo "$FILENAME already exists."
              continue
            fi

            # (Re)create working directory
            if [ -d "$FILENAME" ]; then
              cp -r basicmodel/* "$FILENAME"
            else
              cp -r basicmodel "$FILENAME"
            fi

            ##############################
            # Wait until concurrency < NPARA
            ##############################
            while true; do
              running_count=$(check_and_clean_tasks)
              if [ "$running_count" -lt "$NPARA" ]; then
                break
              fi
              sleep 1
            done

            ###################################
            # Launch job in background 
            ###################################
            (
              cd "$FILENAME" || exit 1
              echo "[$(date)] Start running $FILENAME"
              # Run the bounce program
              ./bounce "$MAXlevel" "$J" "$We" "$Oh" "$Bo" "$epsilon" "$tmax" "$Ldomain"  "$DT" "$CFL" >log_error 2>&1   # DT=1e-4, CFL=0.1
              cat log_error >> log_run_All

              if grep -q "Success" log_error; then  
                echo "[$(date)] Start postprocess $FILENAME: [$(date)]"
                export OMP_NUM_THREADS=1
                python3 getResults.py --We=$We --Oh=$Oh --J=$J --tMAX=$tmax --tSNAP=$tsnap --CPUs=$THREADS >log_results 2>&1
                python3 getVideo.py --RMAX=3 --ZMAX=6 --tMAX=$tmax --tSNAP=$tsnap --CPUs=$THREADS >log_video 2>&1
                ffmpeg -framerate 30 -pattern_type glob -i 'Video/*.png' -vf scale=850:880 -c:v mpeg4 -r 30 -pix_fmt yuv420p video_$FILENAME.mp4 -y >log_video1 2>&1
                cp log_run_All ../Results_Running/log_run_All_$FILENAME.csv
                cp $FILENAME.csv ../Results_Running/$FILENAME.csv
                cp video_$FILENAME.mp4 ../Results_Running/video_$FILENAME.mp4
                cd ..
                tar -czf $FILENAME.tar.gz $FILENAME > /dev/null 2>&1
                rm -rf $FILENAME
                echo "[$(date)] Finish running $FILENAME"
              fi
            ) &
            # Record the PID of the background job
            tasks+=( "$!" )
          done
        done
      done
    done
  done
done

###############################
# Final wait for all tasks
###############################
# Check if any tasks still running, wait for them
while [ "${#tasks[@]}" -gt 0 ]; do
  # Wait for any job to finish
  wait -n 2>/dev/null || true
  # Clean up finished tasks from array
  check_and_clean_tasks >/dev/null
done

echo "All tasks have completed."