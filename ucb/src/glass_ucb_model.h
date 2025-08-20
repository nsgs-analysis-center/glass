/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish
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
 @file glass_ucb_model.h
 \brief Functions for defining, manipulating, and evaluating UCB model
 */


#ifndef ucb_model_h
#define ucb_model_h

#include <stdio.h>

#define UCB_MODEL_NP 8 ///< Number of source parameters for UCB model

/**
\brief Create galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and creates meta template of all sources in the model.
 */
void generate_ucb_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);
void generate_ucb_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
\brief Modify galactic binary model waveform

 Computes LISA response to signal with parameters indexed by `source_id`
 in the Model::Source, and updates meta template of all sources in the model.
 */
void update_ucb_model(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);
void update_ucb_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id);

/**
\brief F-statistic maximization of galactic binary parameters

 Wrapper of UCBFstatistic.c functions for maximizing waveform
 over \f$(\mathcal{A},\cos\iota,\psi,\varphi_0)\f$ by filtering on
 original data.
 
 \todo test filtering on residuals
 */
void maximize_ucb_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

/**
 \brief Create LISA instrument noise model
 
 Computes \f$S_n(f)\f$ from Model::noise.
 */
void generate_noise_model(struct Data *data, struct Model *model);
void generate_noise_model_wavelet(struct Data *data, struct Model *model);

/**
 \brief Converts physical UCB parameters to array expected by glass_ucb_waveform.c
 */
void map_ucb_params_to_array(struct Source *source, double *params, double T);

/**
 \brief Converts array expected by glass_mbh_waveform.c to
 physical UCB parameters
 */
void map_array_to_ucb_params(struct Source *source, double *params, double T);

void remove_ucb_model(struct Data *data, struct Model *model, struct Source *source);
void remove_ucb_model_wavelet(struct Data *data, struct Model *model, struct Source *source);

void add_ucb_model(struct Data *data, struct Model *model, struct Source *source);
void add_ucb_model_wavelet(struct Data *data, struct Model *model, struct Source *source);

/**
\brief Adds input `sample` to existing catalog entry and increments counters.
*/
void append_ucb_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel);

/**
 \brief Wrapper for using functions in GMM_with_EM.c to represent posterior samples of `entry` as a Gaussian Mixture Model.
 */
int ucb_gmm_wrapper(double **ranges, struct Flags *flags, struct Entry *entry, char *outdir, size_t NMODE, size_t NTHIN, unsigned int *seed, double *BIC);



#endif /* ucb_model_h */
