#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-07-17 19:46:35 trottar"
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

echo "Running plot_MC.py for $kinematics"
python3 plot_MC.py $kinematics

convert hist.png xQ2.png pimomt.png pimomt.png pithetat.png nmomt.png nthetat.png ephaseSpace.png piphaseSpace.png nphaseSpace.png pipolar.png epolar.png npolar.png MC_$kinematics.pdf

rm -rf *.png
