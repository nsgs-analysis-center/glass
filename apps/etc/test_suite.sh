#!/bin/bash

export GLASS_PREFIX="/Users/tyson/opt/bin"
export GLASS_REPO="/Users/tyson/glass"
export GLASS_TEST="/Users/tyson/glass_test"

mkdir -p ${GLASS_TEST}

bash run_ucb.sh | tee ${GLASS_TEST}/ucb_dft.out 
bash run_ucb_wavelet.sh | tee ${GLASS_TEST}/ucb_dwt.out 
bash run_amcvn.sh | tee ${GLASS_TEST}/amcvn.out 
bash run_noise.sh | tee ${GLASS_TEST}/noise.out
bash run_vgb.sh | tee ${GLASS_TEST}/vgb.out
bash run_mbh.sh | tee ${GLASS_TEST}/mbh.out
bash run_globalfit.sh | tee ${GLASS_TEST}/globalfit.out 
