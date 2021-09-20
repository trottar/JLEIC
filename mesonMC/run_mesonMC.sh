#! /bin/bash

while getopts 'hgbp' flag; do
    case "${flag}" in
        h) 
        echo "----------------------------------------"
        echo "./run_mesonMC.sh -{flags} {final states}"
        echo "----------------------------------------"
        echo
        echo "The following flags can be called for the EIC_mesonMC generator..."
        echo "    -h, help"
        echo "    -g, generate events"
        echo "    -b, bin events for analysis"
        echo "    -p, plot binned events with test plot script"
        exit 0
        ;;
        g) g_flag='true' ;;
        b) b_flag='true' ;;
        p) p_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

INPUT=$2

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
RANDNUM=$( date +'%Y%M%S' )
NEVTS=${tmp[4]}
PBEAM=${tmp[5]}
KBEAM=${tmp[6]}
SMEAR=${tmp[7]}

python3 "src/binInputs.py"

inputFile="tmp2"

tmp2=()
while IFS='' read -r line || [[ -n "$line" ]]; do
    tmp2+=("$line")
done < "$inputFile"
xbinwidth=${tmp2[0]}
qbinwidth=${tmp2[1]}
tbinwidth=${tmp2[2]}
xLbinwidth=${tmp2[3]}

rm tmp
rm tmp2

cd "src/"

if [[ $INPUT == "k/lambda" ]]; then
    echo
    echo "Kaon with lambda final state selected"
    echo
    particle="kaon"
    SCRIPT="EIC_mesonMC_Lambda.cpp"
elif [[ $INPUT == "k/sigma" ]]; then
    echo
    echo "Kaon with sigma final state selected"
    echo
    particle="kaon"
    SCRIPT="EIC_mesonMC_Sigma.cpp"
elif [[ $INPUT == "pi/p" ]]; then
    echo
    echo "Pion with proton final state selected"
    echo
    particle="pion"
    SCRIPT="EIC_mesonMC.cpp"
elif [[ $INPUT == "pi/n" ]]; then
    echo
    echo "Pion with neutron final state selected"
    echo
    particle="pion"
    SCRIPT="EIC_mesonMC_n.cpp"
elif [[ $INPUT == "bare" ]]; then
    echo
    echo "Bare run (pion with neutron final state) selected. No physics calculations beyond kinematics."
    echo
    particle="N/A"
    SCRIPT="EIC_mesonMC_bare.cpp"
else
    echo "Please select final states (i.e. pi/p, pi/n, k/lambda)"
    exit 2
fi

LOG=$RANDNUM${SCRIPT:0:-4}.log
rootrun="{
root -l<<EOF
.L $SCRIPT+
mainx($XMIN,$XMAX,$Q2MIN,$Q2MAX,$RANDNUM,$NEVTS,$PBEAM,$KBEAM,$SMEAR)
EOF
} 2> $LOG"

if [[ $g_flag = "true" ]]; then
    eval "$rootrun"
    mv $LOG log/
    mv *.d *.so *.pcm lib/    
fi

cd "../"

eonP="${KBEAM}on${PBEAM}"
xq="x${XMIN}-${XMAX}_q${Q2MIN}-${Q2MAX}"

if [ $particle = "pion" ]; then
    kinematics="pi_n_${eonP}_${xq}"
elif [ $particle = "kaon" ]; then
    kinematics="k_lambda_$eonP"
else
    echo "Invalid particle $particle"
    exit 2
fi

BINNING="./src/process/"
ANALYSIS="./analysis/"
if [[ $b_flag = "true" ]]; then
    echo
    echo
    echo "Binning data for $kinematics"
    python3 ${BINNING}bindata.py "$kinematics" "$xbinwidth" "$qbinwidth" "$tbinwidth" "$xLbinwidth"
fi

if [[ $p_flag = "true" ]]; then
    echo
    echo
    echo "Plotting data for $kinematics"
    python3 ${ANALYSIS}plot_test.py "$kinematics" "$xbinwidth" "$qbinwidth" "$tbinwidth" "$xLbinwidth"
fi