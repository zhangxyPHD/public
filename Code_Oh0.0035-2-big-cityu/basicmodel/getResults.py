import os
import subprocess
import multiprocessing as mp
import csv
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
import glob

def compile_executable():
    if os.path.exists("getResults"):
        os.remove("getResults")
    cmd = "qcc -w -Wall -O2 -disable-dimensions getResults.c -o getResults -lm"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError("Compilation of getResults.c failed!")

def run_executable(executable, args):
    args_str = [str(a) for a in args]
    try:
        completed = subprocess.run(
            [executable, *args_str],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Split the stderr output into lines
        lines = completed.stderr.strip().split("\n")

        if lines:
            # Split the first line by whitespace
            parts = lines[0].split()
            # Convert each part to float
            return tuple(float(x) for x in parts)
    except (subprocess.CalledProcessError, ValueError) as e:
        print(f"Error running {executable} with args {args}: {e}")

    # If there's an error or no lines, return an empty tuple
    return ()

def process_Results(ti, tsnap,Oh,We,J,Level):
    t = tsnap * ti
    filepath = f"intermediate/snapshot-{t:.4f}"  
    # print(f"[PID {mp.current_process().pid}] Processing: {filepath}", flush=True)  
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None    
    data = run_executable("./getResults",[filepath,Oh,We,J,Level])
    # If any value in 'data' is None, we consider it a parsing failure.
    if any(v is None for v in data):
        return None 
    print(f"[PID {mp.current_process().pid}] Successfully processed: {filepath}", flush=True)   
    return data

def getResults(tsnap, tmax, results_csv, Oh, We, J,Level, CPUStoUse):
    if os.path.exists(results_csv):
        os.remove(results_csv)
    file_list = sorted(glob.glob("intermediate/snapshot-*"))
    # nsteps = len(file_list)  # This is now the actual number of snapshot files
    nsteps = int(tmax / tsnap)    
    if os.path.exists(results_csv):
        os.remove(results_csv)
    process_func = partial(process_Results, tsnap=tsnap,Oh=Oh,We=We,J=J,Level=Level)
    results = []
    with mp.Pool(processes=CPUStoUse) as pool:
        for res in pool.imap(process_func, range(nsteps + 1)):
            results.append(res)
    # Write results to a CSV file
    header = ["t","n","ke","vcm","vc","vcm1","vc1","Rmax","Zmax","Zmin","vR","uH","pForce","pForce1"]    
    with open(results_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in results:
            if row is not None:
                writer.writerow(row)
    print(f"Parallel processing finished. Results have been written to {results_csv}.")
    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--tMAX', type=float, default=10.0, help='tMAX')
    parser.add_argument('--tSNAP', type=float, default=0.01, help='tSNAP')
    parser.add_argument('--We', type=float, default=10.0, help='We')
    parser.add_argument('--Oh', type=float, default=0.01, help='Oh')
    parser.add_argument('--J', type=float, default=0.0, help='J')
    parser.add_argument('--Level', type=int, default=7, help='Level')
    args = parser.parse_args()
    CPUStoUse = args.CPUs
    tsnap = args.tSNAP
    tmax = args.tMAX
    We = args.We
    Oh = args.Oh
    J = args.J
    Level = args.Level
    #  file dirs
    current_folder = os.path.basename(os.path.dirname(__file__))
    output_results_csv = f"{current_folder}.csv"
    results=getResults(tsnap, tmax, output_results_csv, Oh, We, J, Level, CPUStoUse)
