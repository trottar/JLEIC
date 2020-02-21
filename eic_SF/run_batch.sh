#! /bin/bash

INPUT=$1

if [[ $INPUT == "kaon" ]]; then
    echo "Kaon selected"
    SCRIPT="batchk.cc"
elif [[ $INPUT == "pion" ]]; then
    echo "Pion selected"
    SCRIPT="batch.cc"
else
    echo "Invalid input"
    exit 2
fi

eval "root -q -b -l $SCRIPT"

# mv *.png OUTPUT/
