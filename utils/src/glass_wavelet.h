/*
 * Copyright 2024 Neil J. Cornish & Tyson B. Littenberg
 * Copyright 2026 Robert J. Rosati
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

#ifndef glass_wavelet_h
#define glass_wavelet_h

#define WAVELET_FILTER_CONSTANT 6  // filter steepness in frequency
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

struct Wavelets
{
    /** @name define time-frequency-fdot grid */
     ///@{
    int NF; //!<total number of frequency layers
    int NT; //!<total number of time layers

    double cadence; //!<data sample cadence
    double dt;      //!<wavelet pixel duration
    double df;      //!<wavelet pixel bandwidth
    
    int frequency_steps; //!<frequency steps for wavelet filters
    int fdot_steps;      //!<fdot steps
    double d_fdot;       //!<fractional fdot increment
    ///@}

    int oversample;   //!<oversampling factor for basis functions
    
    double Omega;
    double dOmega;
    double domega;
    double inv_root_dOmega;
    double B;
    double A; //!<roll off of Meyer window. Ranges from 0 to dOmega/2
    double BW;
    double deltaf;
    double *fdot;
    double *window;
    double *window_freq_forward;
    double norm;
    
    int N;
    double T;

    int kmin; //!<minimum wavelet index
    int kmax; //!<maximum wavelet index
    
    int imin; //!<minimum time layer
    int imax; //!<maximum time layer
    
    int *n_table; //!< number of terms in the lookup table at each frequency
    double **table; //!< lookup table of wavelet coefficients
};

struct TimeFrequencyTrack
{
    int min_layer;     //!<minimum frequency layer of track
    int max_layer;     //!<maximum frequency layer of track
    int *segment_size; //!<size of time segment in layer
    int *segment_midpt;//!<midpoint pixel of segment
};

struct TimeFrequencyTrack * malloc_time_frequency_track(struct Wavelets *wdm);
void free_time_frequency_track(struct TimeFrequencyTrack *track);

void initialize_wavelet(struct Wavelets *wdm, double T);
void wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k);
void wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k);
void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *reverse_list, int *Nwavelet, int *jmin, int *jmax);

/* // DEPRECATED WAVELET FUNCS, see new ones below
void wavelet_window_frequency(struct Wavelets *wdm, double *window, int Nlayers);
void wavelet_transform(struct Wavelets *wdm, double *data);
void wavelet_transform_inverse_time(struct Wavelets *wdm, double *data);
void wavelet_transform_inverse_fourier(struct Wavelets *wdm, double *data);
void wavelet_transform_by_layers(struct Wavelets *wdm, int jmin, int Nlayers, double *window, double *data);
void wavelet_transform_segment(struct Wavelets *wdm, int N, int layer, double *data);
void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int *list, int *rlist, int Nmax);
*/

// new transform funcs from Robbie (March 2026)
// these are mostly translations / simplifications of WDMWaveletTransforms

// Note that this is called the timefreq transform in other codes:
// FFT first, then use the fast algorithm for FFT->WDM
void wavelet_transform_timefreq(struct Wavelets *wdm, double *timedata);

// NOTE: this implementation matches WDMWaveletTransforms by Matt Digman
// as well as Olitas.jl
// We assume freqdata has ND+2 elements
void wavelet_transform_freq(struct Wavelets *wdm, double *freqdata, double *wdmdata);

// inverse frequency transform, assumes that freqdata has enough room for ND/2+1 complexes
void wavelet_transform_inverse_freq(struct Wavelets *wdm, double *wdmdata, double *freqdata);

// This is a utility function mostly used in the MBHB waveform
// Here we assume that the freqdata only has frequency content in one layer
// N is number of time data points, freqdata is length N/2+1
void wavelet_transform_freq_segment(struct Wavelets *wdm, int N, int layer, double *data);

// Basically the same thing as wavelet_transform_timefreq, but we assume
// that the output fits into the layers with indices jmin...jmin+Nlayers
// used only in the UCB waveform
// expects timedata of length wdm->NT*(Nlayers+1)
void wavelet_transform_timefreq_by_layers(struct Wavelets* wdm, int jmin, int Nlayers, double *window, double* timedata);

// this function puts the Meyer window into phif properly normalized
void build_wdm_filter_freq(double* phif, int Nf, int Nt, double A, bool forward);

#endif /* glass_wavelet_h */
