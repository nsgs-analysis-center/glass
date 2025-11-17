# Noise model refactor
Ideally, the noise model is more modular, represents stochastic contributions on equal footing.
For now, all are probably freq psd outer product with modulation.
breaking changes that would be nice:
  - more stack-allocated stuff
  - covariance allocated in one memory block
  - PSD/Scaleogram struct? Stochastic contribution. members: PSD, modulation, logL, params.
  - model types build off of this so a pointer cast works
  - modular sampler with priors / proposals custom per source but that's it.
