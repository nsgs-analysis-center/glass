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

#define NTOT 100
#define NBINS 10
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
    SimulateDataWithoutPrinting(&data, &orbit, &flags);
    
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
    clock_t start2 = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = noise_log_likelihood(&data, inst.psd);
    }
    clock_t stop2 = clock();
    
    printf("my_noise_log_likelihood:\n\t%d evals in %g seconds\n",NTOT,(double)(stop-start)/CLOCKS_PER_SEC);
    printf("noise_log_likelihood:\n\t%d evals in %g seconds\n",NTOT,(double)(stop2-start2)/CLOCKS_PER_SEC);
    
    // whitening test
    printf("Running whitening test on X channel...\n");
    double max_bin = 6.0;
    double min_bin = -max_bin;
    double delta_bin = (max_bin - min_bin)/NBINS;
    double wht;
    int counts[NBINS] = {0};
    double mean = 0.0;
    double var = 0.0;
    for (int i=0; i<data.NFFT; i++) {
        wht = data.tdi->X[i] / sqrt(inst.psd->C[0][0][i]);
        //printf("X[%5d] %12lg C[0][0][%5d] %12lg wht: %10lg\n", i, data.tdi->X[i], i, inst.psd->C[0][0][i], wht);
        mean += wht;
        if (wht<min_bin) { counts[0]++; continue; }
        if (wht>max_bin) { counts[NBINS-1]++; continue; }
        int b = (int)((wht - min_bin) / delta_bin);
        assert(b > 0 && b < NBINS);
        counts[b]++;
    }
    mean /= data.NFFT;
    for (int i=0; i<data.NFFT; i++) {
        wht = data.tdi->X[i] / sqrt(inst.psd->C[0][0][i]);
        var += (wht - mean)*(wht-mean);
    }
    var /= data.NFFT;

    printf("\n");
    printf("mean: %lg\n", mean);
    printf(" std: %lg\n", sqrt(var));

    // print out histogram
    printf("Histogram:\n");
    for (int i=0; i<NBINS; i++) {
        printf("\t%2.2lf %5d\n", min_bin + delta_bin*i, counts[i]);
    }

    // repeat for wavelets
    return 0;
}

