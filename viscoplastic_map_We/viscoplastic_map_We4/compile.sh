#!/bin/bash

module load gcc/10.2.0
module load openmpi/4.1.5_gcc

cd basicmodel
qcc -Wall -O2 -disable-dimensions bounce.c -o bounce -lm
qcc -Wall -O2 -disable-dimensions getFacet2D.c -o getFacet2D -lm
qcc -Wall -O2 -disable-dimensions getData2D-VP.c -o getData2D-VP -lm
qcc -Wall -O2 -disable-dimensions getResults.c -o getResults -lm
rm -rf .*
cd ..