/*
 * Copyright 2025 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef glass_model_h
#define glass_model_h

#include <stdio.h>

/**
\brief Hierarchical structure of signal model
 */
struct Model
{
    ///@name Source parameters
    ///@{
    int NP;     //!<number of parameters in model elements
    int Nmax;   //!<maximum number of signals in model
    int Neff;   //!<effective maximum number of signals during burn in
    int Nlive;  //!<current number of signals in model
    struct Source **source; //!<source structures for each signal in the model
    ///@}
    
    /// Noise parameters
    struct Noise *noise;
    
    /// Calibration parameters
    struct Calibration *calibration;
    
    ///@name TDI structures
    ///@{
    struct TDI *tdi; //!<joint signal model
    struct TDI *residual; //!<joint residual
    ///@}
    
    ///@name Segment start time
    ///@{
    double t0; //!<start time
    ///@}
    
    ///@name Source parameter priors
    ///@{
    double **prior; //!<upper and lower bounds for uniform priors \f$ [\theta_{\rm min},\theta_{\rm max}]\f$
    double *logPriorVolume; //!<prior volume \f$ -\Sum \log(\theta_{\rm max}-\theta_{\rm min})\f$
    ///@}
    
    ///@name Model likelihood
    ///@{
    double logL; //!<unnormalized log likelihood \f$ -(d-h|d-h)/2 \f$
    double logLnorm; //!<normalization of log likelihood \f$ \propto -\log \det C \f$
    ///@}

    ///@name Wavelet bookkeeping
    ///@{
    int *list; //!<list of active wavelet pixels
    int Nlist; //!<number of active wavelet pixels
    ///@}
    
    /// Likelihood function
    double (*log_likelihood)(struct Data *, struct Model *);
    
    /// Signal generator
    void (*generate_signal_model)(struct Orbit *, struct Data *, struct Model *, int);
};

/**
\brief Structure containing parameters and meta data for a single binary source.
*/
struct Source
{
    /// Number of parameters (i.e. size of params array etc.
    int NP;
    
    /// Array containing parameters to be passed to ucb_waveform()
    double *params;

    /// Instrument response to signal with Source::params \f$ h(\vec\theta) \f$
    struct TDI *tdi;
    

    ///@name Phenomenological Parameters
    ///@{
    double f0;    //!< Gravitational wave frequency at start of observations \f$f_0 = 2/P\f$
    double dfdt;  //!< First time derivitive of GW frequemcy \f$ \frac{df}{dt}\f$, see galactic_binary_fdot()
    double d2fdt2;//!< Second time derivitive of GW frequemcy \f$ \frac{d^2f}{dt^2}\f$
    double amp;   //!< Gravitational wave amplitude \f$ \mathcal{A} = 2 \frac{\mathcal{M}^{5/3} (\pi f)^{2/3}}{D_L} \f$, see galactic_binary_Amp()
    ///@}
    
    ///@name Extrinsic Parameters
    ///@{
    double psi;     //!< Polarization angle
    double cosi;    //!< Cosine of orbital inclination angle
    double phiref;  //!< Gravitational wave phase at reference
    double phi;     //!< Ecliptic longitude
    double costheta;//!< Cosine of ecliptic co-latitude
    double dL;      //!< Luminosity Distance \f$D_L\f$, see galactic_binary_dL()
    double tref;    //!< Merger time
    ///@}
    
    ///@name Intrinsic Parameters
    ///See data.c for functions to convert between intrinsic and phenomenological parameters
    ///@{
    double m1;    //!< Primary component mass \f$m_1 \geq m_2\f$
    double m2;    //!< Secondary component mass \f$m_2 \leq m_1\f$
    double Mchirp;//!< Chirp mass \f$\mathcal{M}\f$, see galactic_binary_Mc()
    double Mtotal;//!< Total mass
    double chi1;  //!< Dimensionless spin aligned with for m1
    double chi2;  //!< Dimensionless spin aligned with for m2
    ///@}
    
    ///@name Template waveform alignment
    ///See galactic_binary_alignment()
    ///@{
    int BW;   //!< Signal bandwidth in frequency bins, see galactic_binary_bandwidth()
    int qmin; //!< Minimum frequency bin of signal (relative to 0 Hz)
    int qmax; //!< Maximum frequency bin of signal (relative to 0 Hz)
    int imin; //!< Minimum frequency bin of signal relative to segment start
    int imax; //!< Maximum frequency bin of signal relative to segment start
    ///@}

    ///@name Fisher Information Matrix
    ///See galactic_binary_fisher()
    ///@{
    double **fisher_matrix; //!<Fisher approximation to inverse covariance matrix
    double **fisher_evectr; //!<Eigenvectors of covariance matrix
    double *fisher_evalue;  //!<Eigenvalues of covariance matrix
    int fisher_update_flag; //!<1 if fisher needs update, 0 if not
    ///@}

    ///@name Wavelet bookkeeping
    ///@{
    int *list; //!<list of active wavelet pixels
    int Nlist; //!<number of active wavelet pixels
    ///@}

};

/** @name Allocate memory for structures */
///@{
void alloc_model(struct Data *data, struct Model *model, int NP, int Nmax);
void alloc_source(struct Source *source, int N, int NP, int Nchannel);
///@}

/**
 \brief Shallow copy of Data structure
 */
void copy_data(struct Data *origin, struct Data *copy);

/** @name Deep copy structure contents */
///@{
void copy_source(struct Source *origin, struct Source *copy);
void copy_model(struct Model *origin, struct Model *copy);
void copy_model_lite(struct Model *origin, struct Model *copy);
///@}

/**
 \brief swap model pointers
 */
void swap_model(struct Model **ptr1, struct Model **ptr2);

/** @name Free memory for structures */
///@{
void free_model(struct Model *model);
void free_source(struct Source *source);
///@}

/**
 \brief Deep comparison of Model contents
 
 @return 0 if models are identical
 @return 1 if models are different
 */
int compare_model(struct Model *a, struct Model *b);

/**
 \brief Check for increase in maximum log likelihood
 */
int update_max_log_likelihood(struct Model **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Compute difference in log Likelihood from changing parameters of one source
 
 Updates residuals for model y and computes change in likelihood by summin only over data where waveforms are non-zero.
 @return \f$ -\frac{1}{2}\left[(d-h_{\rm new}|d-h__{\rm new}) - (d-h_{\rm old}|d-h_{\rm old})\right] \f$
 */
double delta_log_likelihood(struct Data *data, struct Model *model, struct Model *trial);

/**
 \brief Compute argument of Gaussian likelihood
 
 Computes residual of data and meta-template from Model and noise weighted inner product.
 @return \f$ -\frac{1}{2}(d-h|d-h) \f$
 */
double gaussian_log_likelihood(struct Data *data, struct Model *model);
double gaussian_log_likelhood_wavelet(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for constant noise level
 
 For noise models that are (approximately) a constant over the band
 of \f$N\f$ Fourier bins, parameterized by a multiplyer \f$\eta\f$
 
 @return \f$ \sum_{\rm TDI} N\log\eta_{\rm TDI} \f$
 */
double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model);

/**
 \brief Compute normalization of Gaussian likelihood for arbitrary noise level
 
 For noise models that are free to vary over the band
 of \f$N\f$ Fourier bins
 
 @return \f$ \sum_{\rm TDI} \sum_f \log S_{n,{\rm TDI}}(f) \f$
 */
double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model);

/**
 \brief Create LISA instrument calibration model
 
 Computes phase and amplitude corrections from Model::calibration parameters.
 */
void generate_calibration_model(struct Data *data, struct Model *model);

/**
\brief Apply amplitude and phase corrections.

Computes new LISA instrument response Model::tdi after applying
 amplitude and phase corrections from calibration parameters.
 */
void apply_calibration_model(struct Data *data, struct Model *model);

/// Print current state of waveform and residuals during run for diagnostics. Disabled when Flags::quiet=`TRUE`.
void print_waveform_draw(struct Data *data, struct Model *model, struct Flags *flags);

/// Print waveform power
void print_waveform(struct Data *data, struct Model *model, FILE *fptr);

#endif /* glass_model_h */
