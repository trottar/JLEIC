#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2021-01-29 09:45:50 trottar"
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

xq="x0.010-1.000_q1.0-500.0" # table low Q2, high x
# xq="x0.001-0.100_q1.0-100.0" # table low Q2, low x

if [ $particle = "pion" ]; then
    kinematics="pi_n_${eonP}_${xq}"
elif [ $particle = "kaon" ]; then
    kinematics="k_lambda_$eonP"
else
    echo "Invalid particle $particle"
    exit 2
fi

echo "Running theory_table.py for $kinematics"
python3 theory_table.py "$kinematics"
