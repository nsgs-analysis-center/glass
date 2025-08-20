/*
 * Copyright 2025 Tyson B. Littenberg & Neil Cornish
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
@file glass_mbh_proposal.h
\brief Massive Black Hole binary proposal distributions
 */

#ifndef glass_mbh_proposal_h
#define glass_mbh_proposal_h

#define MBH_PROPOSAL_NPROP 3 ///< Number of defined proposal distributions for MBH sampler
#define MBH_DIFFERENTIAL_EVOLUTION_BUFFER 1000 ///< Size of differential buffer for MBH proposal

double draw_from_mbh_prior(struct Data *data, struct Model *model, struct Source * source, struct Proposal *proposal, double *params, unsigned int *seed);

///<
void setup_differential_evolution_proposal(struct Proposal *proposal);

double differential_evolution_jump(struct Data *data, struct Model *model, struct Source *source, struct Proposal *proposal, double *params, unsigned int *seed);

void update_differential_evolution_buffer(struct Proposal **proposal, struct Model *model, unsigned int *seed);


/**
 \brief Sets up different proposals and assigns frequencies with which they are used.
 */
void initialize_mbh_proposal(struct Orbit *orbit, struct Data *data, struct Prior *prior, struct Chain *chain, struct Flags *flags, struct Catalog *catalog, struct Proposal **proposal, int NMAX);

#endif /* glass_mbh_proposal_h */
