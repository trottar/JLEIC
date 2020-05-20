#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-05-18 14:01:31 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#


kinematics="pi_n_18on275"
# kinematics="pi_n_10on135"
# kinematics="pi_n_10on100"
# kinematics="pi_n_5on100"
# kinematics="pi_n_5on41"
# kinematics="k_lambda_5on100"
# kinematics="k_lambda_18on275"


echo "Running plot_mesonStructure.py for $kinematics"
python3 plot_mesonStructure.py $kinematics

cd OUTPUTS/
convert phase_space.png fpivxpi.png fpivxpi_nolog.png fpivt.png fpi_$kinematics.pdf
rm -rf *.png
