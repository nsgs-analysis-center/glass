#! /usr/bin/env bash

build/apps/src/noise_wavelet_mcmc \
        --steps 10 \
        --chains 1 \
        --threads 1 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*24)) \
        --CD1L \
        --h5-data ./NOISE_731d_2.5s_L1_source0_0_20251206T220508924302Z.h5
#        --threads 1 \
#        --timeseries \
#        --data ./simulated_data_timeseries.dat
