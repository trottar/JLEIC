#! /bin/bash

eval "python plotTDISpionHisto.py"
eval "root -l -q -b run_TDISpion.C"
