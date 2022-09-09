#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2022-05-02 16:36:08 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

# Different versions are required for custom search
python3.8 -m pip install pyqt5==5.15.0

numEvts=$1
IP=$2

python3.8 plot_Fun4all_ecce.py ${numEvts} ${IP}

cd OUTPUTS/DET/
convert *_${IP}.png det_plots_${IP}.pdf

cd ../MOM/
convert *_${IP}.png mom_plots_${IP}.pdf

cd ../KIN/
convert *_${IP}.png kin_plots_${IP}.pdf

cd ../../
evince OUTPUTS/DET/det_plots_${IP}.pdf
evince OUTPUTS/MOM/mom_plots_${IP}.pdf
evince OUTPUTS/KIN/kin_plots_${IP}.pdf
