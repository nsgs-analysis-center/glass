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

/**
 @file glass_mbh_model.h
 \brief Functions for defining, manipulating, and evaluating MBH model
 */


#ifndef mbh_model_h
#define mbh_model_h

#include <stdio.h>


#define MBH_MODEL_NP 11 ///< Number of source parameters for UCB model

void map_array_to_mbh_params(struct Source *source, double *params);
void map_mbh_params_to_array(struct Source *source, double *params);

void generate_mbh_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id);

#endif /* mbh_model_h */

