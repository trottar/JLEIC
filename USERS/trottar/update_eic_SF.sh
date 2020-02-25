#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-02-25 15:14:30 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

mainDIR="../../eic_SF"

echo "================================================="
echo "Copying to $mainDIR"
echo "        TDISMC_EIC.cpp, TDISMC_EIC.h, batch.cc..."
echo "================================================="

eval "cp TDISMC_EIC.cpp \"$mainDIR/TDISMC_EIC.cpp\""

eval "cp TDISMC_EIC.h \"$mainDIR/TDISMC_EIC.hz\""

eval "cp batch.cc \"$mainDIR/batch.cc\""
