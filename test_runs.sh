#! /usr/bin/env bash

rm -rf ./chains
rm -rf ./data
OMP_NUM_THREADS=1 build/apps/src/noise_mcmc --sim-noise --sgwb-template 0 --steps 2000 --chains 1 --cheat --fmin 1e-4 --fmax 8e-3 --duration 1966080
#OMP_NUM_THREADS=6 build/apps/src/noise_mcmc --sim-noise --sgwb-template 0 --steps 2000 --chains 6 --cheat --fmin 1e-4 --fmax 8e-3 --duration 1966080

