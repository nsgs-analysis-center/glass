# noise only
RUNDIR=wdmruns/stat-noise
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps 5000 \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --stationary \
        --rundir $RUNDIR \
        --sim-noise

# noise + conf stationary
RUNDIR=wdmruns/stat-noise-conf
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps 5000 \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --stationary \
        --conf-noise \
        --rundir $RUNDIR \
        --sim-noise

# noise + conf non-stationary
RUNDIR=wdmruns/noise-conf
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps 5000 \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --conf-noise \
        --rundir $RUNDIR \
        --sim-noise

# noise + conf + sgwb stationary
RUNDIR=wdmruns/stat-noise-conf-pl
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps 5000 \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --stationary \
        --conf-noise \
        --sgwb-template 0 \
        --rundir $RUNDIR \
        --sim-noise

# noise + conf + sgwb non-stationary
RUNDIR=wdmruns/noise-conf-pl
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps 5000 \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --conf-noise \
        --sgwb-template 0 \
        --rundir $RUNDIR \
        --sim-noise
# missing: upper limit runs!
