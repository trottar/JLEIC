#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-08-10 11:04:28 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

particle="pion"
# particle="kaon"

# eonP="18on275"
eonP="10on135"
# eonP="10on100"
# eonP="5on100"
# eonP="5on41"

xq="x0.001-1.000_q1.0-1000.0"
# xq="x0.100-1.000_q1.0-1000.0"
# xq="x0.001-1.000_q1.0-80.0"
# xq="x0.100-1.000_q1.0-200.0" # here
# xq="x0.010-1.000_q1.0-200.0"
# xq="x0.001-1.000_q1.0-200.0"
# xq="x0.001-0.500_q1.0-200.0" # here
# xq="x0.000-1.000_q1.0-200.0"

if [ $particle = "pion" ]; then
    kinematics="pi_n_${eonP}_${xq}"
elif [ $particle = "kaon" ]; then
    kinematics="k_lambda_$eonP"
else
    echo "Invalid particle $particle"
    exit 2
fi

echo "Running plot_${particle}Structure.py for $kinematics"
python3 plot_${particle}Structure.py "$kinematics"

cd OUTPUTS/

if [[ $particle = "pion" ]]; then
    convert phase_space.png fpivxpi.png fpivxpi_nolog.png fpivt.png fpivt_xbin.png fpivtPrime.png fpivtPrime_xbin.png fpi_$kinematics.pdf
    rm -rf *.png
elif [[ $particle = "kaon" ]]; then
    convert phase_space.png fkvxk.png fkvxk_nolog.png fkvt.png fk_$kinematics.pdf
    rm -rf *.png
else
    echo "Invalid particle $particle"
    exit 2
fi
