#!/bin/bash 
#SBATCH --partition=cn
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=120:00:00
#SBATCH --account=ehpc-reg-2023r03-178
#SBATCH --qos=ehpc-reg-2023r03-178
#SBATCH --job-name=Oh_OhVALUE_We_WeVALUE_INDEX
#SBATCH -o out-Oh_OhVALUE_We_WeVALUE_INDEX.txt
#SBATCH -e error-Oh_OhVALUE_We_WeVALUE_INDEX.txt

module load gmp/6
module load openmpi/4/gcc/4.1.5
module load gcc/10/10.2.0

# Parameter arrays
tsnap=0.01
Ldomain=8.0
tmax=0.03
Bos=(0.0)
epsilons=(0.001)
MAXlevels=(13)
Wes=(5)
Ohs=(0.0035)
Js=(0.3)
MAXlevel1=11
###################
# Parallel settings
###################
NPARA=1         # Maximum number of concurrent tasks
set_threads_by_maxlevel() {
  case $1 in
    8) echo 16 ;;  # MAXlevel = 9 -> THREADS = 4
    9) echo 16 ;;  # MAXlevel = 9 -> THREADS = 4
    10) echo 16 ;; # MAXlevel = 10 -> THREADS = 5
    11) echo 16 ;; # MAXlevel = 10 -> THREADS = 5
    12) echo 16 ;; # MAXlevel = 10 -> THREADS = 5
    13) echo 16 ;; # MAXlevel = 10 -> THREADS = 5
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
            if [ -e "./Results_Running/log_$FILENAME.csv" ]; then
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
              echo "[$(date)] Start running $FILENAME: [$(date)]"
              # if-else is used for continure.
              while true; do
                # Run the bounce program
                mpirun -n $THREADS ./bounce "$MAXlevel" "$J" "$We" "$Oh" "$Bo" "$epsilon" "$tmax" "$Ldomain" "1e-4" "0.5" "$FILENAME" "$MAXlevel1" >log_error 2>&1
                echo "[$(date)] Log for epsilon=1e-4" | tee -a log_run_All
                cat log_error >> log_run_All

                # Check for errors
                if grep -q "Assertion" log_error || grep -q "Aborted" log_error || grep -q "Ke_Error" log_error; then
                  rm -f "../Results_Running/log_$FILENAME"
                  mpirun -n $THREADS ./bounce "$MAXlevel" "$J" "$We" "$Oh" "$Bo" "$epsilon" "$tmax" "$Ldomain" "1e-5" "0.5" "$FILENAME" "$MAXlevel1" >log_error1 2>&1
                  echo "[$(date)] Log for epsilon=1e-5" | tee -a log_run_All
                  cat log_error1 >> log_run_All
                  if grep -q "Assertion" log_error1 || grep -q "Aborted" log_error1 || grep -q "Ke_Error" log_error1; then
                    rm -f "../Results_Running/log_$FILENAME"
                    mpirun -n $THREADS ./bounce "$MAXlevel" "$J" "$We" "$Oh" "$Bo" "$epsilon" "$tmax" "$Ldomain" "1e-6" "0.5" "$FILENAME" "$MAXlevel1" >log_error2 2>&1
                    echo "[$(date)] Log for epsilon=1e-6" | tee -a log_run_All
                    cat log_error2 >> log_run_All
                  fi              
                else
                  # No errors detected, exit the retry loop
                  break
                fi
              done
              
              if [ -e "../Results_Running/log_$FILENAME.csv" ]; then  
                echo "[$(date)] Start postprocess $FILENAME: [$(date)]"
                export OMP_NUM_THREADS=1
                python3.11 getResults.py --We=$We --Oh=$Oh --J=$J --tMAX=$tmax --tSNAP=$tsnap --CPUs=$THREADS --Level="$MAXlevel" >log_results 2>&1
                python3.11 getVideo.py --RMAX=3 --ZMAX=6 --tMAX=$tmax --tSNAP=$tsnap --CPUs=$THREADS  >log_video 2>&1
                ffmpeg -framerate 30 -pattern_type glob -i 'Video/*.png' -vf scale=850:880 -c:v mpeg4 -r 30 -pix_fmt yuv420p video_$FILENAME.mp4 -y >log_video1 2>&1
                cp log_run_All ../Results_Running/log_run_All_$FILENAME.csv
                cp $FILENAME.csv ../Results_Running/$FILENAME.csv
                cp video_$FILENAME.mp4 ../Results_Running/video_$FILENAME.mp4
                cd ..
                #tar -cvf $FILENAME.tar.gz $FILENAME > /dev/null 2>&1
                #rm -rf $FILENAME
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