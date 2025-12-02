/*
 * Copyright 2024 Robert J. Rosati
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
 @file noise_likelihood_benchmarks.c
 \brief various benchmarks for noise likelihoods and samplers
 */

#include <glass_utils.h>
#include <glass_noise.h>
#include <time.h>
#include <assert.h>

#define NTOT 10000
#define NBINS 20


void whitening_test(const double* fft_data, const double* psd, int N, double fudge_factor)
{
    // whitening test
    printf("Running whitening test...\n");
    double max_bin = 8.0;
    double min_bin = -max_bin;
    double delta_bin = (max_bin - min_bin)/NBINS;
    double wht;
    int counts[NBINS] = {0};
    double mean = 0.0;
    double var = 0.0;
    for (int i=0; i<N; i++) {
        double asd = fudge_factor*sqrt(psd[i]);
        wht = fft_data[i] / asd;
        //printf("X[%5d] %12lg C[0][0][%5d] %12lg wht: %10lg\n", i, data.tdi->X[i], i, inst.psd->C[0][0][i], wht);
        mean += wht;
        if (wht<min_bin) { counts[0]++; continue; }
        if (wht>max_bin) { counts[NBINS-1]++; continue; }
        int b = (int)((wht - min_bin) / delta_bin);
        assert(b >= 0 && b < NBINS);
        counts[b]++;
    }
    mean /= N;
    for (int i=0; i<N; i++) {
        double asd = fudge_factor*sqrt(psd[i]);
        wht = fft_data[i] / asd;
        var += (wht - mean)*(wht - mean);
    }
    var /= N;

    printf("\n");
    printf("mean: %lg\n", mean);
    printf(" std: %lg\n", sqrt(var));

    // print out histogram
    printf("Histogram:\n");
    for (int i=0; i<NBINS; i++) {
        printf("\t%5.2lf %6d\n", min_bin + delta_bin*i, counts[i]);
    }

}

int main(int argc, char *argv[])
{

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);

    fprintf(stdout, "\n================= FOURIER LIKELIHOOD BENCHMARKS ================\n");

    struct Data data = {0};
    struct Flags flags = {0};
    struct Orbit orbit = {0};
    
    /* Initialize data structures */
    data.fmin = 1e-4;
    data.fmax = 1e-2;
    data.T = 60*60*24*30*1;
    data.sqT = sqrt(data.T);
    data.Nchannel = 3;
    strcpy(data.basis , "fourier");
    strcpy(data.format, "sangria");
    data.NFFT = (int)((data.fmax - data.fmin)*data.T);
    data.N = data.NFFT*2;
    data.Nlayer = 1;
    flags.simNoise = 1;
    alloc_data(&data, &flags);
    
    /* Initialize LISA orbit model */
    initialize_orbit(&data, &orbit, &flags);
    // simulate data
    SimulateData(&data, &orbit, &flags);
    
    struct InstrumentModel inst = {0};
    initialize_instrument_model(&orbit, &data, &inst);

    /* get initial likelihood */
    invert_noise_covariance_matrix(inst.psd);
    inst.logL = my_noise_log_likelihood(&data, inst.psd);

    printf("Running likelihood benchmarks...\n");
    clock_t start = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = my_noise_log_likelihood(&data, inst.psd);
    }
    clock_t stop = clock();
    double logL1 = inst.logL;
    clock_t start2 = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = noise_log_likelihood(&data, inst.psd);
    }
    clock_t stop2 = clock();
    double logL2 = inst.logL;
    
    printf("my_noise_log_likelihood (=%g) :\n\t%d evals in %g seconds\n",logL1,NTOT,(double)(stop-start)/CLOCKS_PER_SEC);
    printf("noise_log_likelihood (=%g) :\n\t%d evals in %g seconds\n",logL2,NTOT,(double)(stop2-start2)/CLOCKS_PER_SEC);
    
    whitening_test(data.tdi->X, inst.psd->C[0][0], data.NFFT, 0.25);

    // repeat for wavelets
    fprintf(stdout, "\n================= WAVELET LIKELIHOOD BENCHMARKS ================\n");

    // memory leak! reset structures
    memset(&data,0,sizeof(data)); 
    memset(&flags,0,sizeof(flags)); 
    memset(&orbit,0,sizeof(orbit)); 
    memset(&inst,0,sizeof(inst)); 
    
    /* Initialize data structures */
    data.fmin = 1e-4;
    data.fmax = 1e-2;
    data.T = 60*60*24*30*1;
    data.t0 = 0.0;
    data.sqT = sqrt(data.T);
    data.Nchannel = 3;
    strcpy(data.basis , "wavelet");
    strcpy(data.format, "sangria");
    data.N = (int) data.T / LISA_CADENCE;
    flags.simNoise = 1;
    flags.confNoise = 0;
    flags.sgwbTemplate = -1;
    alloc_data(&data, &flags);
    
    /* Initialize LISA orbit model */
    initialize_orbit(&data, &orbit, &flags);
    initialize_interpolated_analytic_orbits(&orbit, data.T, data.t0);

    // fix wavelet members of struct
    data.lmin = (int)floor(data.fmin / data.wdm->df);
    data.lmax =  (int)ceil(data.fmax / data.wdm->df);
    data.Nlayer = data.lmax - data.lmin;
    //reset wavelet basis max and min ranges
    wavelet_pixel_to_index(data.wdm,0,data.lmin,&data.wdm->kmin); 
    wavelet_pixel_to_index(data.wdm,0,data.lmax,&data.wdm->kmax); 
    // this points kmin/kmax to first time pixel of first/last freq layer
    
    //recompute fmin and fmax so they align with a bin
    data.fmin = data.lmin*WAVELET_BANDWIDTH;
    data.fmax = data.lmax*WAVELET_BANDWIDTH;


    // simulate data
    struct Noise scaleogram = {0};
    alloc_noise(&scaleogram, data.Nlayer*data.wdm->NT, data.Nlayer, data.Nchannel);
    scaleogram.kmin = data.wdm->kmin;
    int k;
    for (size_t i = 0; i < data.wdm->NT; i++)
        for (size_t j=0; j < data.wdm->NF; j++) {
            wavelet_pixel_to_index(data.wdm,i,j,&k); 
            scaleogram.f[k] = (data.lmin + j)*data.wdm->df;
        }
    initialize_instrument_model_wavelet(&orbit, &data, &inst);
    generate_full_dynamic_covariance_matrix(data.wdm, &inst, NULL, NULL, data.noise);
    alloc_tdi(data.tdi, data.N, 3);
    MyAddNoiseWavelet(&data,data.tdi);

    /* get initial likelihood */
    generate_full_dynamic_covariance_matrix(data.wdm, &inst, NULL, NULL, &scaleogram);
    invert_noise_covariance_matrix(&scaleogram);
    inst.logL = my_noise_log_likelihood_wavelet(&data, &scaleogram);

    printf("Running likelihood benchmarks...\n");
    start = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = my_noise_log_likelihood_wavelet(&data, &scaleogram);
    }
    stop = clock();
    logL1 = inst.logL;
    start2 = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = noise_log_likelihood_wavelet(&data, &scaleogram);
    }
    stop2 = clock();
    logL2 = inst.logL;
    
    printf("my_noise_log_likelihood_wavelet (=%g) :\n\t%d evals in %g seconds\n",logL1,NTOT,(double)(stop-start)/CLOCKS_PER_SEC);
    printf("noise_log_likelihood_wavelet (=%g) :\n\t%d evals in %g seconds\n",logL2,NTOT,(double)(stop2-start2)/CLOCKS_PER_SEC);


    whitening_test(data.tdi->X, data.noise->C[0][0], data.N, 0.5/sqrt(2));
    // 1.0          0.444
    // 1/sqrt(2)    0.628
    // 0.5/sqrt(2)  0.628
    // 0.25         1.778

    return 0;
}

