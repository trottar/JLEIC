#! /bin/bash

INPUT=$1

python3 "src/mcInputs.py"

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

cd "src/"

if [[ $INPUT == "kaon" ]]; then
    echo "Kaon selected"
    SCRIPT="TDISMC_EICK.cpp"
elif [[ $INPUT == "pion" ]]; then
    echo "Pion selected"
    SCRIPT="TDISMC_EIC.cpp"
elif [[ $INPUT == "neutron" ]]; then
    echo "Pion with leading neutron selected"
    SCRIPT="TDISMC_EICn.cpp"
else
    echo "Invalid input"
    exit 2
fi

root -l<<EOF
.L $SCRIPT+
mainx($XMIN,$XMAX,$Q2MIN,$Q2MAX,$RANNUM,$NEVTS,$PBEAM,$KBEAM)
EOF

cd "../"

rm tmp
