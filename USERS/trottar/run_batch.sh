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

if [[ $INPUT == "k/lambda" ]]; then
    echo "Kaon with lambda final state selected"
    SCRIPT="TDISMC_EICK.cpp"
elif [[ $INPUT == "pi/p" ]]; then
    echo "Pion with proton final state selected"
    SCRIPT="TDISMC_EIC.cpp"
elif [[ $INPUT == "pi/n" ]]; then
    echo "Pion with neutron final state selected"
    SCRIPT="TDISMC_EICn.cpp"
else
    echo "Please select final states (i.e. pi/p, pi/n, k/lambda)"
    exit 2
fi

root -l<<EOF
.L $SCRIPT+
mainx($XMIN,$XMAX,$Q2MIN,$Q2MAX,$RANNUM,$NEVTS,$PBEAM,$KBEAM)
EOF

cd "../"

rm tmp
