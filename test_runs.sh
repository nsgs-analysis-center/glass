#! /usr/bin/env bash

rm -rf ./chains
rm -rf ./data
OMP_NUM_THREADS=24 build/apps/src/noise_mcmc --sim-noise --conf-noise --steps 2000 --chains 24 --cheat --fmin 1e-4 --fmax 1e-2 --duration 1966080

