#!/bin/bash
Wes=("14")  #13
Ohs=("0.0025")
J1s="0 0.05 0.1 0.15 0.2 0.25 0.3 0.35"
J2s="0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1"
tmax=3.5
for We in "${Wes[@]}"; do
    for Oh in "${Ohs[@]}"; do
        cp run_basic.sh run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/WeVALUE/$We/" run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/OhVALUE/$Oh/" run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/JsVALUE/$J1s/" run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/tmaxValue/$tmax/" run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/INDEX/We"$We"-1/" run_Oh"$Oh"_We"$We"-1.sh
        sbatch run_Oh"$Oh"_We"$We"-1.sh
    done
done

for We in "${Wes[@]}"; do
    for Oh in "${Ohs[@]}"; do
        cp run_basic.sh run_Oh"$Oh"_We"$We"-2.sh
        sed -i "s/WeVALUE/$We/" run_Oh"$Oh"_We"$We"-2.sh
        sed -i "s/OhVALUE/$Oh/" run_Oh"$Oh"_We"$We"-2.sh
        sed -i "s/JsVALUE/$J2s/" run_Oh"$Oh"_We"$We"-2.sh
        sed -i "s/tmaxValue/$tmax/" run_Oh"$Oh"_We"$We"-2.sh
        sed -i "s/INDEX/We"$We"-2/" run_Oh"$Oh"_We"$We"-2.sh
        sbatch run_Oh"$Oh"_We"$We"-2.sh
    done
done