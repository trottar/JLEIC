#! /bin/bash

particle="pion"
eonP="10on135"
xq="x0.001-1.000_q1.0-1000.0"

if [ $particle = "pion" ]; then
    kinematics="pi_n_${eonP}_${xq}"
elif [ $particle = "kaon" ]; then
    kinematics="k_lambda_$eonP"
else
    echo "Invalid particle $particle"
    exit 2
fi

echo "Binning data for $kinematics"
python3 plot_test.py "$kinematics"

