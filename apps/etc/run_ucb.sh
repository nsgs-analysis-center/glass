#!/bin/bash

injection_file="${GLASS_REPO}/ucb/etc/sources/precision/PrecisionSource_1.txt"
threads=12

ucb_sampler="${GLASS_PREFIX}/ucb_mcmc"
outdir="${GLASS_TEST}/ucb_dft"
cmd="${ucb_sampler} --rundir ${outdir} --cheat --threads ${threads} --inj ${injection_file}"
echo $cmd
$cmd
