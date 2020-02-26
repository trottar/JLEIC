#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-02-26 13:20:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

mainDIR="../../eic_SF"

echo "================================================="
echo "Copying to $mainDIR"
echo "        TDISMC_EIC.cpp, TDISMC_EIC.h, run_batch.sh, run_g4e.sh..."
echo "================================================="

eval "cp TDISMC_EIC.cpp \"$mainDIR/TDISMC_EIC.cpp\""

eval "cp TDISMC_EIC.h \"$mainDIR/TDISMC_EIC.h\""

eval "cp run_batch.sh \"$mainDIR/run_batch.sh\""

eval "cp run_g4e.sh \"$mainDIR/run_g4e.sh\""
