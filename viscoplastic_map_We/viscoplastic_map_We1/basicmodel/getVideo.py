# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: Oct 27, 2024

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import multiprocessing as mp
from functools import partial
import argparse
import glob
import matplotlib.colors as mcolors

# Remove LaTeX usage:
# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# Optionally, you may still set fonts without LaTeX:
matplotlib.rcParams['font.family'] = 'serif'

# Custom color map (not LaTeX-dependent):
custom_colors = ["white", "#DA8A67", "#A0522D", "#400000"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_hot", custom_colors)

def compile_executable():
    for f in ["getFacet2D", "getData2D-VP"]:
        if os.path.exists(f):
            os.remove(f)
    cmd = (
        "qcc -w -Wall -O2 getFacet2D.c -o getFacet2D -lm && "
        "qcc -w -Wall -O2 getData2D-VP.c -o getData2D-VP -lm"
    )
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError("Compilation of getResults.c failed!")

def gettingFacets(filename, includeCoat='true'):
    exe = ["./getFacet2D", filename, includeCoat]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2, z2)))
                    segs.append(((-r1, z1),(-r2, z2)))
                    skip = True
    return segs

def gettingfield(filename, zmin, zmax, rmax, nr):
    exe = ["./getData2D-VP", filename, str(zmin), str(0), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Rtemp, Ztemp, D2temp, veltemp  = [], [], [], []

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            D2temp.append(float(temp3[2]))
            veltemp.append(float(temp3[3]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    D2 = np.asarray(D2temp)
    vel = np.asarray(veltemp)
    nz = int(len(Z)/nr)

    print("nz is %d" % nz)

    R.resize((nz, nr))
    Z.resize((nz, nr))
    D2.resize((nz, nr))
    vel.resize((nz, nr))

    return R, Z, D2, vel, nz

def process_timestep(ti, tSNAP, folder, GridsPerR, rmin, rmax, zmin, zmax, lw):
    t = ti*tSNAP
    place = f"intermediate/snapshot-{t:.4f}"
    print(f"[PID {mp.current_process().pid}] Processing: {place}", flush=True)  
    name = f"{folder}/{int(t*1000):08d}.png"

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return

    if os.path.exists(name):
        print(f"{name} Image present!")
        return

    segs1 = gettingFacets(place)
    segs2 = gettingFacets(place, 'false')

    if not segs1 and not segs2:
        print(f"Problem in the available file {place}")
        return

    nr = int(GridsPerR * rmax)
    R, Z, taus, vel, nz = gettingfield(place, zmin, zmax, rmax, nr)
    zminp, zmaxp, rminp, rmaxp = Z.min(), Z.max(), R.min(), R.max()

    # Plotting
    AxesLabel, TickLabel = 50, 20
    fig, ax = plt.subplots()
    fig.set_size_inches(19.20, 10.80)

    ax.plot([0, 0], [zmin, zmax], '-.', color='grey', linewidth=lw)
    ax.plot([rmin, rmin], [zmin, zmax], '-', color='black', linewidth=lw)
    ax.plot([rmin, rmax], [zmin, zmin], '-', color='black', linewidth=lw)
    ax.plot([rmin, rmax], [zmax, zmax], '-', color='black', linewidth=lw)
    ax.plot([rmax, rmax], [zmin, zmax], '-', color='black', linewidth=lw)

    line_segments = LineCollection(segs2, linewidths=4, colors='green', linestyle='solid')
    ax.add_collection(line_segments)
    line_segments = LineCollection(segs1, linewidths=4, colors='blue', linestyle='solid')
    ax.add_collection(line_segments)

    cntrl1 = ax.imshow(
        taus, cmap="hot_r", interpolation='none',
        origin='lower', extent=[-rminp, -rmaxp, zminp, zmaxp],
        vmax=2.0, vmin=-2.0
    )
    cntrl2 = ax.imshow(
        vel, interpolation='Bilinear', cmap='Blues',
        origin='lower', extent=[rminp, rmaxp, zminp, zmaxp],
        vmax=2.0, vmin=0.0
    )

    ax.set_aspect('equal')
    ax.set_xlim(rmin, rmax)
    ax.set_ylim(zmin, zmax)
    # Remove ticks and labels
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Simple text instead of LaTeX:
    ax.text(
        0.02, 0.98, f"t = {t:4.2f}",
        transform=ax.transAxes,
        fontsize=60,
        verticalalignment='top',
        horizontalalignment='left',
        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
    )

    plt.savefig(name, bbox_inches="tight")
    plt.close()
    print(f"[PID {mp.current_process().pid}] Successfully processed: {place}", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--tMAX', type=float, default=10.0, help='tMAX')
    parser.add_argument('--tSNAP', type=float, default=10.0, help='tSNAP')
    parser.add_argument('--ZMAX', type=float, default=4.0, help='Maximum Z value')
    parser.add_argument('--RMAX', type=float, default=4.0, help='Maximum R value')
    parser.add_argument('--ZMIN', type=float, default=0.0, help='Minimum Z value')
    args = parser.parse_args()

    CPUStoUse = args.CPUs
    tMAX = args.tMAX
    tSNAP = args.tSNAP
    ZMAX = args.ZMAX
    RMAX = args.RMAX
    ZMIN = args.ZMIN

    file_list = sorted(glob.glob("intermediate/snapshot-*"))
    nsteps = len(file_list)
    num_processes = CPUStoUse
    rmin, rmax, zmin, zmax = [-RMAX, RMAX, ZMIN, ZMAX]
    GridsPerR = 128

    lw = 2
    folder = 'Video'

    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        process_func = partial(
            process_timestep,
            folder=folder, tSNAP=tSNAP,
            GridsPerR=GridsPerR, rmin=rmin, rmax=rmax,
            zmin=zmin, zmax=zmax, lw=lw
        )
        for _ in pool.imap(process_func, range(nsteps)):
            pass

if __name__ == "__main__":
    # If you need to recompile:
    # compile_executable()
    main()
