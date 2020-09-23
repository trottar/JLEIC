#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-09-14 11:08:13 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

# kinematics="pi_n_18on275"
kinematics="pi_n_10on135"
# kinematics="pi_n_10on100"
# kinematics="pi_n_5on100"
# kinematics="pi_n_5on41"
# kinematics="k_lambda_5on100"
# kinematics="k_lambda_18on275"

xq="x0.001-1.000_q1.0-1000.0" # full t

echo "Running plot_MC.py for ${kinematics}_${xq}"
python3 plot_MC.py "${kinematics}_${xq}"

convert hist.png xQ2.png pimomt.png pimomt.png pithetat.png nmomt.png nthetat.png ephaseSpace.png piphaseSpace.png nphaseSpace.png pipolar.png epolar.png npolar.png MC_$kinematics.pdf

rm -rf *.png
