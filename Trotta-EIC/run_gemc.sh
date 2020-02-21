#! /bin/bash

lund_LOC="$HOME/ResearchNP/JLEIC/Trotta-EIC/"

# gCARD_LOC="$HOME/ResearchNP/gemc/detectors_JLab/basic_examples/exampleCentralTOF/"
gCARD_LOC="$HOME/ResearchNP/gemc/detectors_JLab/basic_examples/my_example/"
# gCARD_LOC="/home/trottar/ResearchNP/gemc/detectors_JLab/eic/script/"

# gCARD="ctof.gcard"
gCARD="example.gcard"
# gCARD="detabx1_beamline.gcard"

FOUT="TDIS_lund"
input_LUND="$FOUT.dat"
output_EVIO="$FOUT.ev"

cd "$gCARD_LOC"

pwd

# With geometry input
echo "Running input $input_LUND for geometry $gCARD"
eval "gemc $gCARD -INPUT_GEN_FILE=\"LUND, $lund_LOC$input_LUND\" -OUTPUT=\"evio, $lund_LOC$output_EVIO\" -USE_GUI=1 -N=100"
# eval "gemc $gCARD -USE_GUI=1 -N=100"

cd $lund_LOC

# echo "Converting $output_EVIO to root file"
# eval "evio2root $input_EVIO $FOUT"

