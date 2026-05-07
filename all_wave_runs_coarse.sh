
# don't exit if a command errors
set +e

NSTEPS=20000
DURATION=$((7680*338*12))
NCH=6
COARSE_Q=13


# noise only
#RUNDIR=wdmruns-coarse/stat-noise
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc_coarse \
#        --steps $NSTEPS \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --stationary \
#        --coarse-Q $COARSE_Q \
#        --rundir $RUNDIR \
#        --sim-noise

#rsync -avz --progress wdmruns-coarse desktop:/mnt/storage/research/
#
## noise + conf stationary
#RUNDIR=wdmruns-coarse/stat-noise-conf
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc_coarse \
#        --steps $NSTEPS \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --stationary \
#        --conf-noise \
#        --coarse-Q $COARSE_Q \
#        --rundir $RUNDIR \
#        --sim-noise
#
#rsync -avz --progress wdmruns-coarse desktop:/mnt/storage/research/

## noise + conf non-stationary
#RUNDIR=wdmruns-coarse/noise-conf
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc_coarse \
#        --steps $NSTEPS \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --conf-noise \
#        --coarse-Q $COARSE_Q \
#        --rundir $RUNDIR \
#        --sim-noise
#
#rsync -avz --progress wdmruns-coarse desktop:/mnt/storage/research/

# noise + conf + sgwb stationary
RUNDIR=wdmruns-coarse/stat-noise-conf-pl
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

RUNDIR=wdmruns-coarse/noise-conf-pl-169
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
#        --stationary-conf \
#rsync -avz --progress wdmruns-coarse desktop:/mnt/storage/research/
# noise + conf + sgwb non-stationary
#RUNDIR=wdmruns-coarse/noise-conf-pl-13
#mkdir -p $RUNDIR
#nice -n 15 build/apps/src/noise_wavelet_mcmc_coarse \
#        --steps $NSTEPS \
#        --chains $NCH \
#        --threads $NCH \
#        --cheat \
#        --fmin 5e-4 \
#        --fmax 8e-3 \
#        --duration $DURATION \
#        --conf-noise \
#        --sgwb-template 0 \
#        --coarse-Q $COARSE_Q \
#        --rundir $RUNDIR \
#        --sim-noise
#
#rsync -avz --progress wdmruns-coarse desktop:/mnt/storage/research/
