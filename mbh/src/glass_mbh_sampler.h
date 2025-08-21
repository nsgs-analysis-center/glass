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

#ifndef glass_mbh_sampler_h
#define glass_mbh_sampler_h

#include <stdio.h>

void mbh_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic);

/**
 \brief Setup UCB sampler
 
 Get all UCB structures into their initial state.
 */
void initialize_mbh_state(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, struct Proposal **proposal, struct Model **model, struct Model **trial, struct Source **inj_vec);

#endif /* glass_mbh_sampler_h */
