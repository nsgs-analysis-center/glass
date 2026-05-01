
# don't exit if a command errors
set +e

NSTEPS=10

# noise only
RUNDIR=wdmruns/stat-noise
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --stationary \
        --rundir $RUNDIR \
        --sim-noise

rsync -avz --progress wdmruns desktop:/mnt/storage/research/

# noise + conf stationary
RUNDIR=wdmruns/stat-noise-conf
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
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

rsync -avz --progress wdmruns desktop:/mnt/storage/research/

# noise + conf non-stationary
RUNDIR=wdmruns/noise-conf
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
        --chains 12 \
        --threads 12 \
        --cheat \
        --fmin 5e-4 \
        --fmax 8e-3 \
        --duration $((7680*338*12)) \
        --conf-noise \
        --rundir $RUNDIR \
        --sim-noise

rsync -avz --progress wdmruns desktop:/mnt/storage/research/
# noise + conf + sgwb stationary
RUNDIR=wdmruns/stat-noise-conf-pl
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
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

rsync -avz --progress wdmruns desktop:/mnt/storage/research/
# noise + conf + sgwb non-stationary
RUNDIR=wdmruns/noise-conf-pl
mkdir -p $RUNDIR
nice -n 15 build/apps/src/noise_wavelet_mcmc \
        --steps $NSTEPS \
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
rsync -avz --progress wdmruns desktop:/mnt/storage/research/
