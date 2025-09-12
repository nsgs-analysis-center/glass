Prototype data analysis software for LISA analysis

# Doxygen code documentation
https://tlittenberg.github.io/glass/html/index.html

# Acknowledgment

```bibtex
@misc{glass,
	author = "Tyson B. Littenberg, Neil J. Cornish",
	title = "GLASS",
	howpublished = "free software (Apache 2.0)"
	year = "2025"
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2026177.svg)](https://doi.org/10.5281/zenodo.2026177)

# C Dependencies 
```bash
openmp
mpi
hdf5
openblas
lapack
```
# Installation

Example build script:
```bash
#!/bin/bash

# set prefix for install directories
export GLASS_PREFIX=<EDIT: path to desired install directories>

# build codes
./install.sh ${GLASS_PREFIX}

# add location of binaries to PATH 
export PATH=${GLASS_PREFIX}/bin:${PATH}
```

Directory tree for the GLASS installation:
```
├── bin
│   ├── gaussian_mixture_model
│   ├── global_fit
│   ├── mbh_mcmc
│   ├── mbh_waveform_benchmark
│   ├── noise_mcmc
│   ├── noise_spline_mcmc
│   ├── ucb_catalog
│   ├── ucb_chirpmass
│   ├── ucb_grid
│   ├── ucb_match
│   ├── ucb_mcmc
│   ├── ucb_waveform_benchmark
│   ├── ucb_wavelet_mcmc
│   └── vgb_mcmc
├── include
│   ├── astrometry
│   │   └── astrometry.h
│   ├── glass
│   │   ├── gitversion.h
│   │   ├── glass_catalog.h
│   │   ├── glass_constants.h
│   │   ├── glass_data.h
│   │   ├── glass_galaxy.h
│   │   ├── glass_gmm.h
│   │   ├── glass_lisa.h
│   │   ├── glass_math.h
│   │   ├── glass_mbh_data.h
│   │   ├── glass_mbh_IMRPhenom.h
│   │   ├── glass_mbh_io.h
│   │   ├── glass_mbh_model.h
│   │   ├── glass_mbh_prior.h
│   │   ├── glass_mbh_proposal.h
│   │   ├── glass_mbh_sampler.h
│   │   ├── glass_mbh_waveform.h
│   │   ├── glass_mbh.h
│   │   ├── glass_model.h
│   │   ├── glass_noise_io.h
│   │   ├── glass_noise_model.h
│   │   ├── glass_noise_sampler.h
│   │   ├── glass_noise.h
│   │   ├── glass_sampler.h
│   │   ├── glass_ucb_data.h
│   │   ├── glass_ucb_fstatistic.h
│   │   ├── glass_ucb_io.h
│   │   ├── glass_ucb_model.h
│   │   ├── glass_ucb_prior.h
│   │   ├── glass_ucb_proposal.h
│   │   ├── glass_ucb_sampler.h
│   │   ├── glass_ucb_waveform.h
│   │   ├── glass_ucb.h
│   │   ├── glass_utils.h
│   │   ├── glass_wavelet.h
│   │   ├── IMRPhenomD_internals.h
│   │   ├── IMRPhenomD.h
│   │   └── IMRPhenomT.h
│   └── kissfft
│       ├── kiss_fft.h
│       ├── kiss_fftr.h
│       └── LICENSES
└── lib
    └── libglass.a

```


# Issue tracker
https://github.com/lisa-analysis-center/glass/issues

# Other resources
 + [Example LISA analysis flow chart](https://www.draw.io/#Hlisa-analysis-center%2Fglass%2Fmaster%2FLISADataFlow.drawio)
 + [Global Fit block diagram](https://app.diagrams.net/#Htlisa-analysis-center%2Fglass%2Fmaster%2FGlobalFit.drawio)
