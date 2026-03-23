/*
 * Copyright 2025 Tyson B. Littenberg & Neil J. Cornish
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
 @file glass_mbh_prior.h
 \brief Functions supporting MBH prior distributions.
 
 Includes functions for creating and evaluating MBH model priors.
 */


#ifndef mbh_prior_h
#define mbh_prior_h

#define MBH_MIN_MASS_RATIO 0.1
#define MBH_MAX_MASS_RATIO 0.99

int mbh_mass_ratio_check(double m1, double m2);

void set_uniform_mbh_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose);

int check_mbh_prior(double *params, double **uniform_prior);

double evaluate_mbh_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params);

#endif /* mbh_prior_h */

