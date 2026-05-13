
# don't exit if a command errors
set +e

NSTEPS=20000
DURATION=$((7680*338*24))
NCH=6
COARSE_Q=13

RUNDIR=wdmruns-coarse2/stat-noise-conf-pl-detection
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc_coarse \
        --steps $NSTEPS \
        --chains $NCH \
        --threads $NCH \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $DURATION \
        --stationary \
        --conf-noise \
        --sgwb-template 0 \
        --coarse-Q $COARSE_Q \
        --rundir $RUNDIR \
        --sim-noise
#        --stationary-conf \

RUNDIR=wdmruns-coarse2/noise-conf-pl-169-detection
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc_coarse \
        --steps $NSTEPS \
        --chains $NCH \
        --threads $NCH \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $DURATION \
        --conf-noise \
        --sgwb-template 0 \
        --coarse-Q 169 \
        --rundir $RUNDIR \
        --sim-noise
