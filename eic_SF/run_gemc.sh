#! /bin/bash

input_LUND="TDIS_lund.txt"
output_EVIO="TDIS_lund.evio"
gCARD="~/ResearchNP/gemc/detectors/gcard/det1_hybrid_full.gcard"
#gCARD="~/ResearchNP/gemc/detectors/gcard/det1_dual_full.gcard"

echo "Running input $input_LUND for geometry $gCARD"
eval "gemc -INPUT_GEN_FILE=\"LUND,$input_LUND\" -OUTPUT=\"evio, $output_EVIO\" -USE_GUI=0 $gCARD"

echo "Converting $output_EVIO to root file"
eval "evio2root -INPUTF=$input_EVIO"

