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

#ifndef glass_mbh_io_h
#define glass_mbh_io_h

#include <stdio.h>

/**
 \brief Print command line options for MBH module
 */
void print_mbh_usage();

/**
 \brief Parse part of command line setting MBH module args
*/
void parse_mbh_args(int argc, char **argv, struct Flags *flags);

/**
 \brief Wrapper function that calls all of the chain print functions
 */
void print_mbh_chain_files(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step);

/** @name Massive Black Hole Merger Chain File
 Print/read current state of source model, e.g. to Chain::parameterFile
 */
///@{
void print_mbh_source_params(struct Data *data, struct Source *source, FILE *fptr);
void scan_mbh_source_params(struct Data *data, struct Source *source, FILE *fptr);
///@}


#endif /* glass_mbh_io_h */
