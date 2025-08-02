/*
 * Copyright 2023 Tyson B. Littenberg & Neil J. Cornish
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


/**
 @file glass_noise_model.h
 \brief Functions and structures defining noise models.

 */
#ifndef noise_model_h
#define noise_model_h

#define MIN_SPLINE_STENCIL 5
#define MIN_SPLINE_SPACING 256.

/*!
 * \brief Data strictire for spline model used in noise model.
 *
 * The SplineModel structure stores metadata and state of spline model,
 * including size of model, current state of spline control points
 * and likelihood, and interpolated model.
 */
struct SplineModel
{
    int Nmin; //!< Minimum number of spline control points
    int Nmax; //!< Maximum number of spline control points
    int Nchannel; //!< Number of TDI channels
    double logL; //!< Log likelhood of spline model
    struct Noise *spline; //!< Spline control points
    struct Noise *psd; //!< Reconstructed noise model
};

struct InstrumentModel
{
    int Nlink;         //!< number of interferometer links (6 for 3-spacecraft constellation)
    double logL;       //!< log Likelihood of model
    double *soms;      //!< optical metrology system noise parameters
    double *sacc;      //!< acceleration noise parameters
    struct Noise *psd; //!< power and cross spectral densities
    
    /** @name Link level noise parameters */
     ///@{
    double sacc12;
    double sacc13;
    double sacc21;
    double sacc23;
    double sacc31;
    double sacc32;
    double soms12;
    double soms13;
    double soms21;
    double soms23;
    double soms31;
    double soms32;
    ///@}
};

struct ForegroundModel
{
    int Nparams;       //!< number of foreground parameters (7)
    double Tobs;       //!< observation time (in seconds)
    double *sgal;      //!< galactic foreground parameters
    double logL;       //!< log Likelihood of model
    struct Noise *psd; //!< power and cross spectral densities

    /** @name galactic foreground parameters from Digman & Cornish 10.3847/1538-4357/ac9139*/
     ///@{
    double Amp;   //!< overall amplitude  (~2.13e-37)
    double f1;
    double alpha; //!< spectral index of low frequency part (~1.58)
    double fk;
    double f2;    //!< some other frequency parameter? (~0.53 mHz)
     ///@}

    struct GalaxyModulation *modulation; //!< time-series of modulation
};

typedef enum {
    SGWB_TEMPLATE_POWERLAW,
    SGWB_TEMPLATE_COUNT // leave this here, counts length of enum
} SGWB_t;
// note that if this array ever becomes extremely large, maybe swap to extern
static const char* SGWB_TEMPLATE_NAMES[] = {"powerlaw"};
static const int SGWB_TEMPLATE_NPARAMS[] = {2};
_Static_assert(sizeof(SGWB_TEMPLATE_NAMES)/sizeof(SGWB_TEMPLATE_NAMES[0]) == SGWB_TEMPLATE_COUNT,
        "Did you add an SGWB template but not its name?");
_Static_assert(sizeof(SGWB_TEMPLATE_NPARAMS)/sizeof(SGWB_TEMPLATE_NPARAMS[0]) == SGWB_TEMPLATE_COUNT,
        "Did you add an SGWB template but not its parameter count?");

struct SGWBModel
{
    int Nparams;       //!< number of SGWB parameters
    double Tobs;       //!< observation time (in seconds)
    double *params;    //!< SGWB parameters
    double logL;       //!< log Likelihood of model
    struct Noise *psd; //!< power and cross spectral densities
    struct SGWBResponse *R; //!< response of the instrument to the isotropic SGWB
    SGWB_t SGWB_type;  //!< which SGWB template is used
};

/**
 \brief Converts physical noise parameters to array expected by InstrumentModel
 */
void map_noise_params_to_array(struct InstrumentModel *model);

/**
 \brief Converts array expected by InstrumentModel to
 physical noise parameters
 */
void map_array_to_noise_params(struct InstrumentModel *model);

/**
 \brief Converts phenomenological foreground parameters to array expected by ForegroundModel
 */
void map_foreground_params_to_array(struct ForegroundModel *model);

/**
 \brief Converts array expected by ForegroundModel to
 phenomenological foreground parameters
 */
void map_array_to_foreground_params(struct ForegroundModel *model);

/**
 \brief Allocates spline model structure and contents.
 */
void alloc_spline_model(struct SplineModel *model, int Ndata, int Nlayer, int Nchannel, int Nspline);

/**
 \brief Allocates instrument model structure and contents.
 */
void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nlayer, int Nchannel);

/**
 \brief Allocates galactic foreground model structure and contents.
 */
void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nlayer, int Nchannel);

/**
 \brief Allocates SGWB model structure and contents.
 */
void alloc_sgwb_model(struct SGWBModel *model, int Ndata, int Nlayer, int Nchannel, SGWB_t SGWB_type);

/**
 \brief Allocates SGWB response.
 */
void alloc_pop_sgwb_response(struct SGWBResponse* sgwbr, char* fname);

/**
 \brief Free allocated spline model.
 */
void free_spline_model(struct SplineModel *model);

/**
 \brief Free allocated spline model.
 */
void free_instrument_model(struct InstrumentModel *model);

/**
 \brief Free allocated foreground model.
 */
void free_foreground_model(struct ForegroundModel *model);

/**
 \brief Free allocated SGWB model.
 */
void free_sgwb_model(struct SGWBModel *model);

/**
 \brief Deep copy of SplineModel structure from `origin` into `copy`
 */
void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy);

/**
 \brief Deep copy of InstrumentModel structure from `origin` into `copy`
 */
void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy);

/**
 \brief Deep copy of ForegroundModel structure from `origin` into `copy`
 */
void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy);

/**
 \brief Deep copy of SGWBModel structure from `origin` into `copy`
 */
void copy_sgwb_model(struct SGWBModel *origin, struct SGWBModel *copy);

/**
 \brief Wrapper to `CubicSplineGSL` functions for generating PSD model based on current state of `model`
 */
void generate_spline_noise_model(struct SplineModel *model);

/**
 \brief Compute instrument model contribution to noise covariance matrix based on current state of `model`
 */
void generate_instrument_noise_model(struct Orbit *orbit, struct InstrumentModel *model);
void generate_instrument_noise_model_wavelet(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model);

/**
 \brief Compute galactic foreground contribution to covariance matrix based on current state of `model`
 */
void generate_galactic_foreground_model(struct ForegroundModel *model);
void generate_galactic_foreground_model_wavelet(struct Wavelets *wdm, struct ForegroundModel *model);

/**
 \brief Compute SGWB contribution to covariance matrix based on current state of `model`
 */
void generate_sgwb_model(struct SGWBModel *model);
void generate_sgwb_model_wavelet(struct Wavelets *wdm, struct SGWBModel *model);

/**
 \brief Add components to `full` noise covariance matrix `C`
 */
void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel);
void generate_full_dynamic_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full);

/**
\brief Compute spline model only where interpolant changes
 
 Interpolates spline points only in the vicinity of `new_knot`
 */
void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot);

/**
 \brief Log likelihood for noise model.
 
 @param data Data structure
 @param model SplineModel structure containing current state
 @return \f$  \ln p({\rm data}|{\rm spline}) \f$
 */
double noise_log_likelihood(struct Data *data, struct Noise *noise);
double noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise);

/**
 \brief Change in log likelihood for noise model.
 
 @param data Data structure
 @param model_x Current SplineModel structure containing current state
 @param model_y Current SplineModel structure containing current state
 @return \f$  \ln p({\rm data}|{\rm spline}_y)-\ln p({\rm data}|{\rm spline}_x) \f$
 */
double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax, int ic);

/**
 \brief Set initial state of spline `model`
 */
void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline);

/**
 \brief Set initial state of instrument noise `model`
 */
void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);
void initialize_instrument_model_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model);
/**
 \brief Set initial state of galactic foreground `model`
 */
void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);
void initialize_foreground_model_wavelet(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model);

/**
 \brief Set initial state of sgwb `model`
 */
void initialize_sgwb_model(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type);
void initialize_sgwb_model_wavelet(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type);


void GetDynamicNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags);
void GetStationaryNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Noise *noise);

/**
 \brief functional form for a powerlaw SGWB
 */
double sgwb_powerlaw(double f, const double* params);

#endif /* noise_model_h */
