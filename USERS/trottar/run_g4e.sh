#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2020-02-25 16:24:52 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

g4e_FILE="eicMC_g4e.py"

echo "Compiling ejpm enviroment..."
source /home/trottar/.local/share/ejpm/env.csh

cd g4e_files/
python3 $g4e_FILE
