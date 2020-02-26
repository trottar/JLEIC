#! /bin/bash

INPUT=$1

python3 mcInputs.py

inputFile="tmp"

tmp=()
while IFS='' read -r line || [[ -n "$line" ]]; do
    tmp+=("$line")
done < "$inputFile"

XMIN=${tmp[0]}
XMAX=${tmp[1]}
Q2MIN=${tmp[2]}
Q2MAX=${tmp[3]}
RANNUM=$( date '+%H%M%S' )
NEVTS=${tmp[4]}
PBEAM=${tmp[5]}
KBEAM=${tmp[6]}

if [[ $INPUT == "kaon" ]]; then
    echo "Kaon selected"
    SCRIPT="TDISMC_EICk.cpp"
elif [[ $INPUT == "pion" ]]; then
    echo "Pion selected"
    SCRIPT="TDISMC_EIC.cpp"
else
    echo "Invalid input"
    exit 2
fi

root -l<<EOF
.L $SCRIPT+
mainx($XMIN,$XMAX,$Q2MIN,$Q2MAX,$RANNUM,$NEVTS,$PBEAM,$KBEAM)
EOF

rm tmp
