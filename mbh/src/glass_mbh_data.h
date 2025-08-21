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
 @file glass_mbh_data.h
 \brief Functions for MBH data handling, including reading external data, simulating internal data, and signal injections.
 
 */

#ifndef glass_mbh_data_h
#define glass_mbh_data_h

#include <stdio.h>

/**
 \brief Injection routine for generic binaries
 
 Unlike verification sources, this code expects the polarization angle \f$\psi\f$ or (currently) initial phase \f$\varphi_0\f$ to be included in the parameter files.
 */
void MBHInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source **inj_vec);


#endif /* glass_mbh_data_h */
