#! /usr/bin/env bash


OMP_NUM_THREADS=6 build/apps/src/noise_mcmc --sim-noise --conf-noise --steps 1000 --chains 6 --cheat --fmin 1e-4 --fmax 1e-2
