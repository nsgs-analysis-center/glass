/*
 * Copyright 2024 Neil J. Cornish & Tyson B. Littenberg
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

#include <assert.h>
#include "glass_utils.h"
#include "glass_math.h"

struct TimeFrequencyTrack * malloc_time_frequency_track(struct Wavelets *wdm)
{
    struct TimeFrequencyTrack *track = malloc(sizeof(struct TimeFrequencyTrack));
    track->min_layer = 1;
    track->max_layer = wdm->NF;
    track->segment_size  = int_vector(wdm->NF);
    track->segment_midpt = int_vector(wdm->NF);
    return track;
}

void free_time_frequency_track(struct TimeFrequencyTrack *track)
{
    free_int_vector(track->segment_size);
    free_int_vector(track->segment_midpt);
    free(track);
}

static void setup_wdm_basis(struct Wavelets *wdm, int NF)
{
    wdm->NF = NF;
    wdm->cadence = WAVELET_DURATION/wdm->NF;
    wdm->Omega = M_PI/wdm->cadence;
    wdm->dOmega = wdm->Omega/(double)wdm->NF;
    wdm->inv_root_dOmega = 1.0/sqrt(wdm->dOmega);
    wdm->B = wdm->dOmega;
    wdm->A = (wdm->dOmega-wdm->B)/2.0;
    wdm->BW = (wdm->A+wdm->B)/M_PI;
}

static double phitilde(struct Wavelets *wdm, double omega)
{
    double insDOM = wdm->inv_root_dOmega;
    double A = wdm->A;
    double B = wdm->B;
    
    double x, y, z;
    
    z = 0.0;
    
    if(fabs(omega) >= A && fabs(omega) < A+B)
    {
        x = (fabs(omega)-A)/B;
        y = incomplete_beta_function(WAVELET_FILTER_CONSTANT, WAVELET_FILTER_CONSTANT, x);
        z = insDOM*cos(y*M_PI/2.0);
    }
    
    if(fabs(omega) < A) z = insDOM;
    
    return(z);
    
}

static void wavelet(struct Wavelets *wdm, int m, double *wave)
{
    
    int N = wdm->N;
    double dom = wdm->domega;
    double DOM = wdm->dOmega;
    
    double omega;
    double x, y, z;
    
    double *DE = (double*)malloc(sizeof(double)*(2*N));
    
    // zero and postive frequencies
    for(int i=0; i<=N/2; i++)
    {
        omega = (double)(i)*dom;
        
        y = phitilde(wdm,omega+(double)(m)*DOM);
        z = phitilde(wdm,omega-(double)(m)*DOM);
        
        x = y+z;
        
        REAL(DE,i) = M_SQRT1_2*x;
        IMAG(DE,i) = 0.0;
        
    }
    
    // negative frequencies
    for(int i=1; i< N/2; i++)
    {
        omega = -(double)(i)*dom;
        
        y = phitilde(wdm,omega+(double)(m)*DOM);
        z = phitilde(wdm,omega-(double)(m)*DOM);
        
        x = y+z;
        
        REAL(DE,N-i) = M_SQRT1_2*x;
        IMAG(DE,N-i) = 0.0;
        
    }
    
    glass_inverse_complex_fft(DE,N);

    for(int i=0; i < N/2; i++)
    {
        wave[i] = REAL(DE,N/2+i)/wdm->norm;
        wave[i+N/2] = REAL(DE,i)/wdm->norm;
    }

    free(DE);
}

static void wavelet_window_time(struct Wavelets *wdm)
{
    double *DX = (double*)malloc(sizeof(double)*(2*wdm->N));
    
    //zero frequency
    REAL(DX,0) =  wdm->inv_root_dOmega;
    IMAG(DX,0) =  0.0;
    
    for(int i=1; i<= wdm->N/2; i++)
    {
        int j = wdm->N-i;
        double omega = (double)(i)*wdm->domega;
        
        // postive frequencies
        REAL(DX,i) = phitilde(wdm,omega);
        IMAG(DX,i) =  0.0;
        
        // negative frequencies
        REAL(DX,j) =  phitilde(wdm,-omega);
        IMAG(DX,j) =  0.0;
    }
        
    glass_inverse_complex_fft(DX, wdm->N);

    wdm->window = (double*)malloc(sizeof(double)* (wdm->N));
    for(int i=0; i < wdm->N/2; i++)
    {
        wdm->window[i] = REAL(DX,wdm->N/2+i);
        wdm->window[wdm->N/2+i] = REAL(DX,i);
    }
    
    wdm->norm = sqrt((double)wdm->N * wdm->cadence / wdm->domega);

    free(DX);
}

static void wavelet_lookup_table(struct Wavelets *wdm)
{    
    // it turns out that all the wavelet layers are the same modulo a
    // shift in the reference frequency. Just have to do a single layer
    // we pick one far from the boundaries to avoid edge effects
    
    double *wave = (double*)malloc(sizeof(double)*(wdm->N));
    int ref_layer = wdm->NF/2;
    wavelet(wdm, ref_layer, wave);
    
    // The odd wavelets coefficienst can be obtained from the even.
    // odd cosine = -even sine, odd sine = even cosine
    // each wavelet covers a frequency band of width DW
// execept for the first and last wasvelets
    // there is some overlap. The wavelet pixels are of width
    // DOM/PI, except for the first and last which have width
    // half that
    
    double f0 = ref_layer*wdm->df;
    
    for(int j=0; j<wdm->fdot_steps; j++)  // loop over f-dot slices
    {

        int NT = wdm->n_table[j];
        
        
        for(int n=0; n<NT; n++)  // loop over time slices
        {
            double f = f0 + ((double)(n-NT/2)+0.5)*wdm->deltaf;
            
            double real_coeff = 0.0;
            double imag_coeff = 0.0;
            
            for(int i=0; i<wdm->N; i++)
            {
                double t = ((double)(i-wdm->N/2))*wdm->cadence;
                double phase = PI2*f*t + M_PI*wdm->fdot[j]*t*t;
                real_coeff += wave[i]*cos(phase)*wdm->cadence;
                imag_coeff += wave[i]*sin(phase)*wdm->cadence;
            }
            wdm->table[j][2*n]   = real_coeff;
            wdm->table[j][2*n+1] = imag_coeff;
        }
    }
    
    free(wave);
}
void initialize_wavelet(struct Wavelets *wdm, double T)
{
    fprintf(stdout,"\n========= Initialize Wavelet Basis ==========\n");
    wdm->NT = (int)ceil(T/WAVELET_DURATION);
    wdm->NF = WAVELET_DURATION/LISA_CADENCE;
    wdm->df = WAVELET_BANDWIDTH; // Pixel spacing in f
    wdm->dt = WAVELET_DURATION; // Pixel spacing in t

    //TODO: this needs to be fixed instead of fudged
    //avoid edge effects at start and stop times of segment.
    wdm->imin = WAVELET_EDGE_BUFFER;
    wdm->imax = wdm->NT - WAVELET_EDGE_BUFFER;
    
    setup_wdm_basis(wdm, wdm->NF);

    wdm->frequency_steps = 400;
    wdm->fdot_steps = 50;
    wdm->d_fdot = 0.1;
    wdm->oversample = 8.0;

    wdm->N = wdm->oversample * 2 * wdm->NF;
    wdm->T = wdm->N*wdm->cadence;

    wdm->domega = PI2/wdm->T;
    
    wdm->deltaf = wdm->BW/(double)(wdm->frequency_steps);

    wdm->fdot = malloc(wdm->fdot_steps*sizeof(double));

    /* Only needed when using the lookup table waveform generator */
    wdm->table   = malloc(wdm->fdot_steps*sizeof(double *));
    wdm->n_table = malloc(wdm->fdot_steps*sizeof(int));

    double fdot_step = wdm->df/wdm->T*wdm->d_fdot; // sets the f-dot increment
        
    for(int n=0; n<wdm->fdot_steps; n++)
    {
        wdm->fdot[n] = -fdot_step*wdm->fdot_steps/2 + n*fdot_step;
        size_t N = (int)((wdm->BW+fabs(wdm->fdot[n])*wdm->T)/wdm->deltaf);
        if(N%2 != 0) N++; // makes sure it is an even number
        wdm->n_table[n] = N;
        wdm->table[n] = double_vector(2*N);
    }

    //stores window function and normalization
    wavelet_window_time(wdm);
    wdm->window_freq_forward = double_vector((wdm->NT/2 + 1)); // memory leak
    build_wdm_filter_freq(wdm->window_freq_forward, wdm->NF, wdm->NT, wdm->A, true);

    //stores lookup table of wavelet basis functions
    //wavelet_lookup_table(wdm);

    //set defaults for min and maximum pixels
    wavelet_pixel_to_index(wdm,0,1,&wdm->kmin);         //first pixel of second layer
    wavelet_pixel_to_index(wdm,0,wdm->NF-1,&wdm->kmax); //first pixel of last layer

    fprintf(stdout,"  Number of time pixels:        %i\n", wdm->NT);
    fprintf(stdout,"  Duration of time pixels:      %g [hr]\n", wdm->dt/3600);
    fprintf(stdout,"  Number of frequency layers:   %i\n", wdm->NF);
    fprintf(stdout,"  Bandwidth of frequency layer: %g [uHz]\n", wdm->df*1e6);
    fprintf(stdout,"=============================================\n");
}

void inline wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k)
{
    int NT = wdm->NT;
    
    //which time
    *i = k%NT; 
    
    //which frequency
    *j = (k - (*i))/NT; //scary integer math
}

void inline wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k)
{
    int NT = wdm->NT;
    
    *k = i + j*NT;
}

// FFT first, then use the fast algorithm for FFT->WDM
// to the user, appears to be an in-place operation
void wavelet_transform_timefreq(struct Wavelets *wdm, double *timedata) {
    int ND = wdm->NT*wdm->NF;
    double freqdata[ND + 2]; // extra complex element

    glass_forward_real_fft_outplace(timedata, freqdata, ND);

    // note: below will overwrite timedata
    wavelet_transform_freq(wdm, freqdata, timedata);
}

static inline double my_phitilde(double omega, double DeltaOmega, double A)
{
    double B = DeltaOmega - 2*A;
    double insDOM = 1.0/sqrtf(DeltaOmega);
    
    double x, y, z;
    
    z = 0.0;
    
    if(fabs(omega) >= A && fabs(omega) < A+B)
    {
        x = (fabs(omega)-A)/B;
        y = incomplete_beta_function(WAVELET_FILTER_CONSTANT, WAVELET_FILTER_CONSTANT, x);
        z = insDOM*cos(y*M_PI_2);
    }
    
    if(fabs(omega) < A) z = insDOM;
    
    return(z);
    
}

void build_wdm_filter_freq(double* phif, int Nf, int Nt, double A, bool forward) {
    // this function puts the Meyer window into phif properly normalized
    double DeltaOmega = M_PI / Nf;
    double domega = 2*DeltaOmega / (double)Nt;
    for (int i = 0; i < Nt/2 + 1; i++) {
        phif[i] = my_phitilde(domega*i, DeltaOmega, A);
    }
    // normalize the window
    double norm = 0.0;
    for (int i = 0; i<Nt/2+1; i++) {
        if (i==0)
            norm += phif[i]*phif[i];
        else
            norm += 2*phif[i]*phif[i];
    }
    norm *= domega * M_1_PI;
    norm = sqrtf(norm);
    if (forward)
        norm *= Nf / 2.0;
    for (int i = 0; i<Nt/2+1; i++) {
        phif[i] /= norm;
    }
}
// inverse frequency transform, assumes that freqdata has enough room for ND/2+1 complexes
void wavelet_transform_inverse_freq(struct Wavelets *wdm, double *wdmdata, double *freqdata) {
    int Nf = wdm->NF;
    int Nt = wdm->NT;
    int ND = Nf*Nt;
    memset(freqdata, 0, (ND+2)*sizeof(double));
    
    double prefactors[2*Nt];
    // use filled normalized Meyer window into phif
    double *phif = wdm->window_freq_forward;
    double fft_result[2*Nt];
    for (int m=0; m < Nf+1; m++) {
        int center = m*Nt/2;
        memset(prefactors, 0, 2*Nt*sizeof(double));
        if (m==0) {
            // DC layer, in even rows of column 0
            for (int n=0; n < Nt; n++)
                prefactors[2*n] = wdmdata[(2*n)%Nt] * M_SQRT1_2; // undo normalization
        }
        else if (m==Nf) {
            for (int n=0; n < Nt; n++)
                prefactors[2*n] = wdmdata[(2*n)%Nt + 1] * M_SQRT1_2; // undo normalization
        }
        else {
            for (int n = 0; n < Nt; n++) {
                if ((n+m)%2==0) {
                    prefactors[2*n] = wdmdata[n + Nt*m]; 
                }
                else {
                    prefactors[2*n + 1] = -wdmdata[n + Nt*m]; 
                }
            } //end for n
        }
        glass_forward_complex_fft_outplace(prefactors, fft_result, Nt);

        int imin = center - Nt/2;
        if (imin < 0) imin = 0;
        int imax = center + Nt/2;
        if (imax > ND/2 +1) imax = ND/2+1;
        for (int i = imin; i < imax; i++) {
            int iind = abs(i-center);
            if (iind > Nt/2)
                continue;
            // note that Nf/2 here to adapt window for inverse!
            if (m==0 || m==Nf) {
                freqdata[2*i]     += fft_result[2*((2*i)%Nt)] * phif[iind] * Nf/2;
                freqdata[2*i + 1] += fft_result[2*((2*i)%Nt) + 1] * phif[iind] * Nf/2;
            }
            else {
                freqdata[2*i]     += fft_result[2*(i%Nt)] * phif[iind] * Nf/2;
                freqdata[2*i+ 1]  += fft_result[2*(i%Nt) + 1] * phif[iind] * Nf/2;
            }
        }

    } // end for m
}

// NOTE: this implementation matches WDMWaveletTransforms by Matt Digman
// as well as Olitas.jl
// We assume freqdata has ND+2 elements
void wavelet_transform_freq(struct Wavelets *wdm, double *freqdata, double *wdmdata) {

    int Nf = wdm->NF;
    int Nt = wdm->NT;

    double   DX[2*Nt]; // length NT complex array, we will ifft later
    double DX_t[2*Nt]; // length NT complex array, will be one freq layer
    // fill normalized Meyer window into phif
    double *phif = wdm->window_freq_forward;

    for (int m=0; m <= Nf; m++) {
        int center = m*Nt / 2; // note integer division
        memset(DX, 0, 2*Nt*sizeof(double));
        for (int j=-Nt/2 + 1; j < Nt/2; j++) {
            int idx = center + j;
            // various edge checks:
            if ((m==0) && (idx < 0))
                continue; // m==0  is DC bin , no freqs below 
            if ((m==Nf) && idx > center)
                continue; // m==Nf is Nyquist, no freqs above
            double weight = phif[abs(j)];
            if ( (j==0) && (m==0 || m==Nf) )
                weight /= 2.0; // this accounts for less d.o.f. in boundary layers
            DX[2*(Nt/2 + j)]     = weight * freqdata[2*idx];
            DX[2*(Nt/2 + j) + 1] = weight * freqdata[2*idx+1];
        }
        glass_inverse_complex_fft_outplace(DX, DX_t, Nt);

        // now extract wdm coeffs
        if (m==0) {
            // DC layer, even rows of col 0
            for (int n = 0; n < Nt; n+=2)
                wdmdata[n] = DX_t[2*n] * M_SQRT2;
        }
        else if (m==Nf) {
            // Nyquist layer, odd rows of col 0
            for (int n = 0; n < Nt; n+=2)
                wdmdata[n+1] = DX_t[2*n] * M_SQRT2;
        }
        else {
            // General layer
            int imag_sign = m%2 == 0 ? 1 : -1;
            for (int n = 0; n < Nt; n++) {
                if ((n+m)%2==0)
                    wdmdata[n+Nt*m] = DX_t[2*n];
                else
                    wdmdata[n+Nt*m] = imag_sign * DX_t[2*n+1];
            }
        }
    }
}

// Transform a short segment of FFT data into one WDM layer.
// freqdata is N complex entries (2*N doubles), centered on the layer frequency
// (i.e. the layer's center frequency is at index N/2).
// wdmdata receives N WDM coefficients for this layer.
// Each output pixel has the same ΔT and ΔF as the full transform;
// N > Nt gives more time coverage, N < Nt gives less.
void wavelet_transform_freq_one_layer(struct Wavelets *wdm, int N, double *freqdata, double *wdmdata, int layer) {

    int Nf = wdm->NF;
    int Nt = wdm->NT;
    // TODO: which N are allowed compared to Nt?

    double   DX[2*N]; // length N complex array, we will ifft later
    double DX_t[2*N]; // length N complex array, result of ifft
    memset(DX, 0, 2*N*sizeof(double));
    int m = layer;

    // Evaluate the window at each bin's actual frequency in the short FFT.
    // Bin spacing: domega = 2*DeltaOmega/N, where DeltaOmega = PI/Nf (sample space).
    // The window has support over N/2 bins (goes to zero at bin N/2).
    double DeltaOmega = M_PI / Nf;
    double domega = 2.0 * DeltaOmega / N;

    // Build and normalize the filter on this grid
    double phif[N/2 + 1];
    for (int i = 0; i <= N/2; i++)
        phif[i] = my_phitilde(domega * i, DeltaOmega, wdm->A);
    // normalize: same formula as build_wdm_filter_freq
    double norm = 0.0;
    for (int i = 0; i <= N/2; i++) {
        if (i == 0)
            norm += phif[i] * phif[i];
        else
            norm += 2 * phif[i] * phif[i];
    }
    norm *= domega * M_1_PI;
    norm = sqrtf(norm);
    norm *= Nf / 2.0; // forward transform
    for (int i = 0; i <= N/2; i++)
        phif[i] /= norm;

    for (int j = -N/2 + 1; j < N/2; j++) {
        int local_idx = N/2 + j; // index into the short FFT array
        // edge checks for boundary layers
        if ((m == 0) && (j < 0))
            continue;
        if ((m == Nf) && (j > 0))
            continue;
        if (local_idx < 0 || local_idx >= N)
            continue;
        double weight = phif[abs(j)];
        if ((j == 0) && (m == 0 || m == Nf))
            weight /= 2.0;
        DX[2*(N/2 + j)]     = weight * freqdata[2*local_idx];
        DX[2*(N/2 + j) + 1] = weight * freqdata[2*local_idx + 1];
    }
    glass_inverse_complex_fft_outplace(DX, DX_t, N);

    // extract wdm coefficients
    if (m == 0) {
        for (int n = 0; n < N; n += 2)
            wdmdata[n] = DX_t[2*n] * M_SQRT2;
    }
    else if (m == Nf) {
        for (int n = 0; n < N; n += 2)
            wdmdata[n+1] = DX_t[2*n] * M_SQRT2;
    }
    else {
        int imag_sign = m % 2 == 0 ? 1 : -1;
        for (int n = 0; n < N; n++) {
            // normally this lhs index would be m*Nt+n
            if ((n + m) % 2 == 0)
                wdmdata[n] = DX_t[2*n];
            else
                wdmdata[n] = imag_sign * DX_t[2*n+1];
        }
    }
}

// Basically the same thing as wavelet_transform_timefreq, but we assume
// that the output fits into the layers with indices jmin...jmin+Nlayers
// used only in the UCB waveform
// expects timedata of length wdm->NT*(Nlayers+1)
void wavelet_transform_timefreq_by_layers(struct Wavelets* wdm, int jmin, int Nlayers, double *window, double* timedata) {
    int Nf = wdm->NF;
    int Nt = wdm->NT;
    int N = Nt * (Nlayers+1); // length of input data
    double freqdata[N + 2]; // full-length FFT

    // TODO: tukey window was here, removed for invertibility in tests
    /*
    double alpha = 8.0/Nt;
    tukey(timedata, alpha, N);
    */
    // FFT the full timeseries
    glass_forward_real_fft_outplace(timedata, freqdata, N);

    double   DX[2*Nt];
    double DX_t[2*Nt];
    double *phif = window; // should be length Nt/2+1

    // output is packed as Nlayers * Nt elements
    memset(timedata, 0, Nlayers * Nt * sizeof(double));

    for (int m = jmin; m < jmin + Nlayers; m++) {
        int center = (m - jmin + 1) * Nt / 2;
        memset(DX, 0, 2*Nt*sizeof(double));
        for (int j = -Nt/2 + 1; j < Nt/2; j++) {
            int idx = center + j;
            if (idx < 0 || idx > N/2)
                continue;
            double weight = phif[abs(j)];
            if ((j == 0) && (m == 0 || m == Nf))
                weight /= 2.0;
            DX[2*(Nt/2 + j)]     = weight * freqdata[2*idx];
            DX[2*(Nt/2 + j) + 1] = weight * freqdata[2*idx+1];
        }
        glass_inverse_complex_fft_outplace(DX, DX_t, Nt);

        int layer_offset = (m - jmin) * Nt;
        if (m == 0) {
            for (int n = 0; n < Nt; n += 2)
                timedata[layer_offset + n] = DX_t[2*n] * M_SQRT2;
        }
        else if (m == Nf) {
            for (int n = 0; n < Nt; n += 2)
                timedata[layer_offset + n + 1] = DX_t[2*n] * M_SQRT2;
        }
        else {
            int imag_sign = m % 2 == 0 ? 1 : -1;
            for (int n = 0; n < Nt; n++) {
                if ((n + m) % 2 == 0)
                    timedata[layer_offset + n] = DX_t[2*n];
                else
                    timedata[layer_offset + n] = imag_sign * DX_t[2*n+1];
            }
        }
    }
}

void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int *list, int *rlist, int Nmax)
{
    
    int n, k, jj, kk;
    double dx, dy;
    double f, fdot;
    double fmid,fsam;
    double cos_phase, sin_phase, y, z, yy, zz;

    double df = wdm->deltaf;

    // maximum frequency and frequency derivative
    double f_max    = wdm->df*(wdm->NF-1);
    double fdot_max = wdm->fdot[wdm->fdot_steps-1];
    double fdot_min = wdm->fdot[0];
    double d_fdot   = wdm->fdot[1]-wdm->fdot[0]; // f-dot increment
    
    for(int i=0; i<wdm->NT; i++)
    {
        f     = freq[i];
        fdot  = freqd[i];
        
        //skip this step if f or fdot violate bounds
        if(f>=f_max || fdot>=fdot_max || fdot<=fdot_min)  continue;
        
        cos_phase = amp[i]*cos(phase[i]);
        sin_phase = amp[i]*sin(phase[i]);
        
        n = (int)floor((fdot-fdot_min)/d_fdot);  // lower f-dot layer
        dy = (fdot-fdot_min)/d_fdot - n;         // where in the layer
                                    
        for(int j=jmin[i]; j<=jmax[i]; j++)
        {
        
            // central frequency
            fmid = j*wdm->df;
                
            kk = (int)floor( ( f - (fmid + 0.5*df) )/df );
            fsam = fmid + (kk + 0.5)*df;
            dx = (f - fsam)/df; // used for linear interpolation
                
            // interpolate over frequency
            y = 0.0;
            z = 0.0;
            yy = 0.0;
            zz = 0.0;

            jj = kk + wdm->n_table[n]/2;
            if(jj>=0 && jj< wdm->n_table[n]-1)
            {
                y = (1.0-dx)*wdm->table[n][2*jj]   + dx*wdm->table[n][2*(jj+1)];
                z = (1.0-dx)*wdm->table[n][2*jj+1] + dx*wdm->table[n][2*(jj+1)+1];
            }

            jj = kk + wdm->n_table[n+1]/2;
            if(jj >=0 && jj < wdm->n_table[n+1]-1)
            {
                yy = (1.0-dx)*wdm->table[n+1][2*jj]   + dx*wdm->table[n+1][2*(jj+1)];
                zz = (1.0-dx)*wdm->table[n+1][2*jj+1] + dx*wdm->table[n+1][2*(jj+1)+1];
            }
                
            // interpolate over fdot
            y = (1.0-dy)*y + dy*yy;
            z = (1.0-dy)*z + dy*zz;

            // make sure pixel is in range
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {  
                int n = rlist[k - wdm->kmin];
                if(n<Nmax)
                {
                    if((i+j)%2 == 0) wave[n] =  (cos_phase*y - sin_phase*z);
                    else             wave[n] = -(cos_phase*z + sin_phase*y);
                }
                else
                {
                    //fprintf(stderr,"Warning, wavelet_transform_from_table tried accessing array out of bounds\n");
                    //fflush(stderr);
                }
            }
            
        } //loop over frequency layers
        
    } //loop over time steps
}

void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *reverse_list, int *Nwavelet, int *jmin, int *jmax)
{
    
    int n;
    int k;
    int N;
    double fmx, fdmx, fdmn, dfd, HBW;
    double fmax, fmin;
    double fdotmax, fdotmin;
    int Xflag, Yflag, Zflag;
    
    double df = wdm->deltaf;
    double DF = wdm->df;
    double *fd = wdm->fdot;

    // maximum frequency and frequency derivative
    fmx  = (double)(wdm->NF-1)*DF;
    fdmx = fd[wdm->fdot_steps-1];
    fdmn = fd[0];
    dfd  = fd[1]-fd[0]; // f-dot increment
    
    N = 0;
    for(int i=0; i<wdm->NT; i++)
    {
        // check to see if any of the channels are ok
        Xflag = Yflag = Zflag = 0;
        if(freqX[i] < fmx) Xflag = 1;
        if(freqY[i] < fmx) Yflag = 1;
        if(freqZ[i] < fmx) Zflag = 1;

        // shut off any channel that does not have valid fdots
        if(fdotX[i] < fdmn || fdotX[i] > fdmx) Xflag = 0;
        if(fdotY[i] < fdmn || fdotY[i] > fdmx) Yflag = 0;
        if(fdotZ[i] < fdmn || fdotZ[i] > fdmx) Zflag = 0;

        // skip if no channels have valid values
        if(!Xflag && !Yflag && !Zflag) continue;

        /*  find the largest and smallest frequencies and frequency derivatives
        but only using the valid channels */
        fmin = 1;
        fmax = 0;
        fdotmin =  1;
        fdotmax = -1;

        if(Xflag)
        {
            if(freqX[i]>fmax) fmax=freqX[i];
            if(freqX[i]<fmin) fmin=freqX[i];
            if(fdotX[i]>fdotmax) fdotmax=fdotX[i];
            if(fdotX[i]<fdotmin) fdotmin=fdotX[i];
        }
        
        if(Yflag)
        {
            if(freqY[i]>fmax) fmax=freqY[i];
            if(freqY[i]<fmin) fmin=freqY[i];
            if(fdotY[i]>fdotmax) fdotmax=fdotY[i];
            if(fdotY[i]<fdotmin) fdotmin=fdotY[i];
        }
        
        if(Zflag)
        {
            if(freqZ[i]>fmax) fmax=freqZ[i];
            if(freqZ[i]<fmin) fmin=freqZ[i];
            if(fdotZ[i]>fdotmax) fdotmax=fdotZ[i];
            if(fdotZ[i]<fdotmin) fdotmin=fdotZ[i];
        }
       
        //skip if max/min fdot go out of bounds 
        if(fdotmax >= fdmx || fdotmin <= fdmn) continue;

        // lowest f-dot layer
        n = (int)(floor(fdotmin-fdmn/dfd));
        int NL = wdm->n_table[n];

        // highest f-dot layer
        n = (int)(floor(fdotmax-fdmn/dfd));
        int NH = wdm->n_table[n];

        // find which has the largest number of samples
        if(NL > NH) NH = NL;

        // half bandwidth of layer        
        HBW = 0.5*(NH-1)*df;
        
        // lowest frequency layer
        jmin[i] = (int)ceil((fmin-HBW)/DF);
        
        // highest frequency layer
        jmax[i] = (int)floor((fmax+HBW)/DF);   

        // skip any out-of-bounds layers
        if(jmin[i] < 0) jmin[i] = 0;
        if(jmax[i] > wdm->NF-1) jmax[i] = wdm->NF-1;
        
        for(int j=jmin[i]; j<=jmax[i]; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            
            //check that the pixel is in range
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                wavelet_list[N]=k-wdm->kmin;
                reverse_list[k-wdm->kmin]=N;
                N++;
            }
        }  
    }
    *Nwavelet = N;
}

void wavelet_window_frequency(struct Wavelets *wdm, double *window, int Nlayers)
{
    int i;
    int N;
    double T;
    double domega;
    double omega;
    double norm=0.0;
    
    N = (Nlayers+1);

    //mini wavelet structure for basis covering just N layers
    struct Wavelets *wdm_temp = malloc(sizeof(struct Wavelets));
    setup_wdm_basis(wdm_temp, N);
    
    T = wdm->dt*wdm->NT;
    
    domega = PI2/T;

    //wdm window function
    for(i=0; i<=wdm->NT/2; i++)
    {
        omega = i*domega;
        window[i] = phitilde(wdm_temp,omega);
    }
    
    //normalize
    for(i=-wdm->NT/2; i<= wdm->NT/2; i++) norm += window[abs(i)]*window[abs(i)];
    norm = sqrt(norm/wdm_temp->cadence);
    
    for(i=0; i<=wdm->NT/2; i++) window[i] /= norm;
    
    free(wdm_temp);
    
}

// This is a utility function mostly used in the MBHB waveform
// Here we assume that the freqdata only has frequency content in one layer
// N is number of time data points, freqdata is length N/2+1
void wavelet_transform_freq_segment(struct Wavelets *wdm, int N, int layer, double *freqdata)
{
    double wdmdata[N];
    memset(wdmdata, 0, N*sizeof(double));
    wavelet_transform_freq_one_layer(wdm, N, freqdata, wdmdata, layer);
    for (int i = 0; i < N; i++)
        freqdata[i] = wdmdata[i];
}

