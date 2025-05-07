#! /usr/bin/env bash

OMP_NUM_THREADS=6 build/apps/src/ucb_wavelet_mcmc --inj ucb/etc/sources/precision/PrecisionSource_0.txt --cheat --chains 6 --steps 2000
