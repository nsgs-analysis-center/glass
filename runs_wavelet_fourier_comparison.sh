#! /usr/bin/env bash

# let's try an instrument-only injection and recovery, in both WDM and wavelet basis

OUTPUT_DIR=data-inst-wdm

if [[ ! -d "$OUTPUT_DIR" ]]; then
    OMP_NUM_THREADS=1 build/apps/src/noise_wavelet_mcmc --sim-noise --steps 10000 --chains 1 --cheat --fmin 5e-4 --fmax 8e-3 --duration $((60*60*24*30*1)) --threads 1
    mv data data-inst-wdm
    mv chains chains-inst-wdm
fi


OUTPUT_DIR=data-inst-fft

if [[ ! -d "$OUTPUT_DIR" ]]; then
    OMP_NUM_THREADS=1 build/apps/src/noise_mcmc --sim-noise --steps 10000 --chains 1 --cheat --fmin 5e-4 --fmax 8e-3 --duration $((60*60*24*30*1)) --threads 1
    mv data data-inst-fft
    mv chains chains-inst-fft
fi




#OMP_NUM_THREADS=1 build/apps/src/noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --steps 10000 --chains 1 --cheat --fmin 5e-4 --fmax 8e-3 --duration $((60*60*24*30*1)) --threads 1
