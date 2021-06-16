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

while getopts 'bpv' flag; do
  case "${flag}" in
    b) b_flag='true' ;;
    p) p_flag='true' ;;
    v) verbose='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [[ $b_flag = "true" ]]; then
  echo "Binning data for $kinematics"
  python3 bindata.py "$kinematics"
fi

ANALYSIS="../../analysis/"

if [[ $p_flag = "true" ]]; then
    echo "Plotting data for $kinematics"
    python3 ${ANALYSIS}plot_test.py "$kinematics"
fi
