#!/bin/bash 
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions bounce.c -o bounce -lm
mpirun -n THREADS ./bounce parameters