#!/bin/bash

mbh_sampler="${GLASS_PREFIX}/mbh_mcmc"
injection_file="${GLASS_REPO}/apps/etc/mbh_injection.dat"
outdir="${GLASS_TEST}/mbh"
threads=24
Tobs=7864320

cmd="${mbh_sampler} --rundir ${outdir} --cheat --threads ${threads} --inj ${injection_file} --duration ${Tobs}"

echo $cmd

$cmd

