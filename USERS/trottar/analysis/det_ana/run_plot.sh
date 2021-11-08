#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2021-11-08 09:38:53 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

numEvts=$1
IP=$2

python3.8 plot_Fun4all.py ${numEvts} ${IP}

cd OUTPUTS/DET/

convert *.png det_plots_${IP}.pdf

cd ../MOM/

convert *.png mom_plots_${IP}.pdf
