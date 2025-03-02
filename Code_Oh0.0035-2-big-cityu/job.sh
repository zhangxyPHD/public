#!/bin/bash
Wes=("1" "5" "10")
Ohs=("0.1")
J1s="0.6 0.8 1.0 1.5 2.0"
J2s="3.0 4.0 6.0 8.0 10.0"

for We in "${Wes[@]}"; do
    for Oh in "${Ohs[@]}"; do
        cp run_basic.sh run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/WeVALUE/$We/" run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/OhVALUE/$Oh/" run_Oh"$Oh"_We"$We"-1.sh
        sed -i "s/JsVALUE/$J1s/" run_Oh"$Oh"_We"$We"-1.sh
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
        sed -i "s/INDEX/We"$We"-2/" run_Oh"$Oh"_We"$We"-2.sh
        sbatch run_Oh"$Oh"_We"$We"-2.sh
    done
done
