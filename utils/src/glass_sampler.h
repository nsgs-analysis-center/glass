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

#ifndef glass_sampler_h
#define glass_sampler_h

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>

/*!
 \brief Prototype structure for prior distributions.
 
 Generic data structure for holding all information needed by prior distributions.
 Structure contains parameters for different supported priors and various book-keeping scalars, vectors, and matrices to
 hold needed metadata.

*/
struct Prior
{
    ///@name Uniform prior
    ///@{
    double **prior; //!<upper and lower bounds for uniform priors \f$ [\theta_{\rm min},\theta_{\rm max}]\f$
    double logPriorVolume; //!<prior volume \f$ -\sum \log(\theta_{\rm max}-\theta_{\rm min})\f$
    ///@}

    ///@name Anisotropic sky prior
    ///@{
    double *skyhist; //!<2D histogram of prior density on sky
    double dcostheta; //!<size of `skyhist` bins in \f$\cos\theta\f$ direction
    double dphi; //!<size of `skyhist` bins in \f$\phi\f$ direction
    double skymaxp; //!<max prior density of `skyhist`
    int ncostheta; //!<number of `skyhist` bins in \f$\cos\theta\f$ direction
    int nphi; //!<number of `skyhist` bins in \f$\phi\f$ direction
    ///@}
    
    ///@name workspace
    ///@{
    double *vector;  //!<utility 1D array for prior metadata
    double **matrix; //!<utility 2D array for prior metadata
    double ***tensor;//!<utility 3D array for prior metadata
    ///@}

    /// Gaussian Mixture Model prior
    struct GMM *gmm;
    
    /**
     \brief Compute prior density given parameters
     @param params parameter vector
     @return logP prior density
     */
    double (*density)(struct Flags*,struct Data*,struct Model*,struct Prior*, double *);
};

/*!
 \brief Prototype structure for proposal distributions.
 
 Generic data structure for holding all information needed by proposal distributions.
 Structure contains function for drawing new parameters, evaluating the proposal density,
 tracking acceptance ratios, and various book-keeping scalars, vectors, and matrices to
 hold needed metadata.

*/
struct Proposal
{
    /**
     \brief Function that generates updated source parameters.
     
     @param  params parameter vector
     @return logQ proposal density
     */
    double (*function)(struct Data*,struct Model*,struct Source*,struct Proposal*,double*,unsigned int*);

    /**
     \brief Compute proposal density given parameters.
     
     @param[in]  params parameter vector
     @param[out] logQ proposal density
     */
    double (*density)(struct Data*, struct Model*, struct Source*,struct Proposal*,double*);
    
    int *trial;      //!<total number of trials for proposal
    int *accept;     //!<total number of accepted trials for proposals*/
    char name[128];  //!<string identifying proposal type
    double norm;     //!<proposal normalization
    double maxp;     //!<max value of proposal density for rejection sampling
    double weight;   //!<proposal weight [0,1] for fixed dimension moves
    double rjweight; //!<proposal weight [0,1] for trans dimensional moves
    int size;        //!<size of proposal arrays
    double *vector;  //!<utility 1D array for proposal metadata
    double **matrix; //!<utility 2D array for proposal metadata
    double ***tensor;//!<utility 3D array for proposal metadata
    
    /** @name Gaussian mixture model
     */
    ///@{
    size_t Ngmm; //!< number of mixture models (1/source)
    struct GMM **gmm; //!<array of individual mixture models
    ///@}
};


/** \brief Parallel tempering exchange
 
 Cycles through all chains and proposes swaps between adjacent pairs.
 */
void ptmcmc(struct Model **model, struct Chain *chain, struct Flags *flags);

/**
 \brief Adaptive temperature spacing
 
 Adjusts temperature spacing between tempered chains according with the goal of having an even acceptance rate between all pairs. The sensitivity of the adjustment asymptotically goes to zero as the sampler approaches the end of the burn-in phase.
 */
void adapt_temperature_ladder(struct Chain *chain, int mcmc);

/**
\brief Fair draw from uniform ranges for each parameter
 
 @param params (updates \f$\vec\theta\f$)
 @return logQ = \f$\ln p(\f$ \c params \f$)\f$

 */
double draw_from_uniform_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed);

/**
 \brief Returns (log) uniform prior density
 
 Returns \f$\sum \log\frac{1}{\Delta V}\f$ for each parameter except the SNR prior for \f$\mathcal{A}\f$.
 */
double uniform_prior_density(struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params);

/**
 \brief Placeholder for symmetric proposal.  Returns 0.0
 
 The Proposal structure requires a proposal density function. This function serves that role for proposals which are symmetric and therefore do not need use any resources computing proposal densities that will just cancel in the Hastings ratio.
 */
double symmetric_density(UNUSED struct Data *data, UNUSED struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, UNUSED double *params);

/**
 \brief Compute and print acceptance ratios for each proposal
 */
void print_acceptance_rates(struct Proposal **proposal, int NProp, int ic, FILE *fptr);

#endif /* glass_sampler_h */
