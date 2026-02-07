/*
 * Copyright 2025 Robert J. Rosati
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

#define NTOT 10
#define NBINS 20
#define NFFT_TEST (1536*7)

// TODO this wavelet logic is only actually different in the loop. must be a way to combine them?
void whitening_test_wavelet(const struct Data* data)
{
    // whitening test
    double max_bin = 8.0;
    double min_bin = -max_bin;
    double delta_bin = (max_bin - min_bin)/NBINS;
    double wht;
    int counts[NBINS] = {0};
    double mean = 0.0;
    double var = 0.0;
    int k;
    int Nactive = 0;

    const double* wdm_data = data->tdi->X;
    const double* psd = data->noise->C[0][0];
    struct Wavelets* wdm = data->wdm;

    for (int i=0; i<wdm->NT; i++) {
        for (int j=data->lmin; j<data->lmax; j++) {
            wavelet_pixel_to_index(wdm,i,j,&k);
            k -= wdm->kmin;
            double asd = sqrt(psd[k]/2);
            wht = wdm_data[k] / asd;
            //printf("X[%5d] %12lg C[0][0][%5d] %12lg wht: %10lg\n", i, data.tdi->X[i], i, inst.psd->C[0][0][i], wht);
            mean += wht;
            if (wht<min_bin) { counts[0]++; continue; }
            if (wht>max_bin) { counts[NBINS-1]++; continue; }
            int b = (int)((wht - min_bin) / delta_bin);
            assert(b >= 0 && b < NBINS);
            counts[b]++;
            Nactive++;
        }
    }
    mean /= Nactive ;
    for (int i=0; i<wdm->NT; i++) {
        for (int j=data->lmin; j<data->lmax; j++) {
            wavelet_pixel_to_index(wdm,i,j,&k);
            k -= wdm->kmin;
            double asd = sqrt(psd[k]/2);
            wht = wdm_data[k] / asd;
            var += (wht - mean)*(wht - mean);
        }
    }
    var /= Nactive;

    printf("\n");
    printf("mean: %lg\n", mean);
    printf(" std: %lg\n", sqrt(var));

    // print out histogram
    printf("Histogram:\n");
    for (int i=0; i<NBINS; i++) {
        printf("\t%5.2lf %6d\n", min_bin + delta_bin*i, counts[i]);
    }

}

void whitening_test(const double* fft_data, const double* psd, int N)//, double fudge_factor)
{
    // whitening test
    printf("Running whitening test on %d points...\n", N);
    double max_bin = 8.0;
    double min_bin = -max_bin;
    double delta_bin = (max_bin - min_bin)/NBINS;
    double wht;
    int counts[NBINS] = {0};
    double mean = 0.0;
    double var = 0.0;
    for (int i=0; i<N; i++) {
        double asd = sqrt(psd[i]/2);
        wht = fft_data[2*i] / asd; // real part only for now
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
        double asd = sqrt(psd[i]/2);
        wht = fft_data[2*i] / asd;
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

int write_wdm_data(struct Wavelets* wdm, double* data, char* fname) {
    FILE *fptr = NULL;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    int k;
    for (int i=0; i<wdm->NT; i++) {
        double t = i*wdm->dt;
        for (int j=0; j<wdm->NF; j++) {
            double f = j*wdm->df;
            wavelet_pixel_to_index(wdm,i,j,&k);
            k-=wdm->kmin;
            fprintf(fptr,"%lg ",t);
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg",data[k]);
            fprintf(fptr,"\n");
        }
    }
    return 0;
}

int write_fft_data(double df, int NFFT, double* data, char* fname) {
    FILE *fptr = NULL;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    for (int i=0; i<NFFT; i++) {
            double f = i*df;
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg ",data[2*i]);
            fprintf(fptr,"%lg" ,data[2*i+1]);
            fprintf(fptr,"\n");
    }
    return 0;
}

int main(int argc, char *argv[])
{

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    fprintf(stdout, "\n================= CONVENTION CHECKS ================\n");
    // test what FFT coeffs are
    double test_data[NFFT_TEST] = {0};
    double test_data_copy[NFFT_TEST] = {0};
    test_data[0] = 1.0;
    glass_forward_real_fft(&test_data, NFFT_TEST);
    printf("FFT of impulse:\n");
    double real_val = test_data[0];
    double imag_val = test_data[1];
    for (int i=0; i<NFFT_TEST/2; i++) {
        test_data_copy[i] = test_data[i];
        if ((fabs(test_data[2*i] - real_val) > 1e-10) || (fabs(test_data[2*i+1] - imag_val) > 1e-10)) {
            printf("\tNot all FFT coeffs equal!");
            break;
        }
    }
    printf("\t%g + j*%g\n", real_val, imag_val);
    printf("numpy has:\n");
    printf("\t%g + j*%g\n", 1.0, 0.0);

    printf("IFFT of FFT(impulse):\n");
    bool perfect=true;
    double atol=1e-10;
    glass_inverse_real_fft(&test_data, NFFT_TEST);
    if (fabs(test_data[0] - 1.0) > atol) {
        printf("\tMismatch, impulse should be 1.0, got %lf\n", test_data[0]);
        perfect=false;
    }
    for (int i=1; i<NFFT_TEST/2; i++) {
        if ((abs(test_data[i]) > atol)) {
            printf("\tFirst non-impulse result at index %d : %lf\n", i, test_data[i]);
            perfect=false;
            break;
        }
    }
    if (perfect)
        printf("\tMatches within atol=%g\n", atol);

    // timeseries
    memset(test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    struct Wavelets wdm = {0};
    // NFFT_TEST == NF*NT
    // T == NFFT_TEST*LISA_CADENCE
    double T = NFFT_TEST*LISA_CADENCE;
    initialize_wavelet(&wdm, T);
    wavelet_transform(&wdm, &test_data);
    write_wdm_data(&wdm, &test_data, "./dbg_wdm_impulse.dat");
    wavelet_transform_inverse_fourier(&wdm, &test_data);
    write_fft_data(1./T, NFFT_TEST/2, &test_data, "./dbg_wdmfft_impulse.dat");
    write_fft_data(1./T, NFFT_TEST/2, &test_data_copy, "./dbg_fft_impulse.dat");

    perfect = true;
    printf("WDM(impulse)->FFT vs FFT(impulse):\n");
    for (int i=0; i<NFFT_TEST; i++) {
        if (fabs(test_data[i] - test_data_copy[i]) > atol) {
            printf("\tFirst mismatch at index %d\n", i);
            printf("\t\tExpected: %g\n",test_data_copy[i]);
            printf("\t\tGot:      %g\n",test_data[i]);
            perfect = false;
            break;
        }
    }
    if (perfect)
        printf("\tMatches within atol=%g\n", atol);
    
    


    /*
    memset(&test_data,0,sizeof(test_data));
    test_data[0] = 1.0;
    test_data[1] = -1.0;
    glass_forward_real_fft(&test_data, NFFT_TEST);
    FILE* fptr = fopen("./fft_test.dat","w");
    for (int i=0; i<NFFT_TEST/2; i++) {
        fprintf(fptr, "(%g+%gj)\n", test_data[2*i], test_data[2*i+1]);
    }
    fclose(fptr);
    */
    // okay, so, fft conventions with KissFFT seem to match numpy!



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
    // Write data
    FILE* fptr = NULL;
    fptr = fopen("./simulated_data_fft.dat", "w");
    if (fptr) {
        for (int i=0; i<data.NFFT; i++) {
            double f = (double)(i+data.qmin)/data.T;
            //             f   rex imx rey imy rez imz
            fprintf(fptr, "%lg %lg %lg %lg %lg %lg %lg\n", f, data.tdi->X[2*i], data.tdi->X[2*i+1], data.tdi->Y[2*i], data.tdi->Y[2*i+1], data.tdi->Z[2*i], data.tdi->Z[2*i+1]);
        }
        fclose(fptr);
    } else {
        printf("Couldn't open data output file for writing!\n");
    }
    // Write data
    // need to zero pad! we only have a band-limited part FFT, not full freq content
    int N = (int)floor(data.T / LISA_CADENCE);
    double df = 1./data.T;
    double *Xtime = malloc(sizeof(double)*N);
    double *Ytime = malloc(sizeof(double)*N);
    double *Ztime = malloc(sizeof(double)*N);
    int emin = 2 * (int)floor(data.fmin / df);
    int emax = 2 * (int)floor(data.fmax / df);
    printf("Did we get the indices right?\n");
    printf("\temax - emin: %d\n", emax-emin);
    printf("\tdata.N: %d\n", data.N);
    for (int i=0; i<N; i++) {
        double f = (i/2)*df;
        if ((data.fmin < f) && (f < data.fmax)) {
            Xtime[i] = data.tdi->X[i-emin];
            Ytime[i] = data.tdi->Y[i-emin];
            Ztime[i] = data.tdi->Z[i-emin];
        } else {
            Xtime[i] = 0;
            Ytime[i] = 0;
            Ztime[i] = 0;
        }
    }
    glass_inverse_real_fft(Xtime, N);
    glass_inverse_real_fft(Ytime, N);
    glass_inverse_real_fft(Ztime, N);
    fptr = fopen("./simulated_data_timeseries.dat", "w");
    if (fptr) {
        for (int i=0; i<N; i++) {
            double t = i*LISA_CADENCE;
            //             t   x y z
            fprintf(fptr, "%lg %lg %lg %lg\n", t, Xtime[i], Ytime[i], Ztime[i]);
        }
        fclose(fptr);
    } else {
        printf("Couldn't open data output file for writing!\n");
    }
    free(Xtime);
    free(Ytime);
    free(Ztime);
     
    // TODO: GetNoiseModel has C[i][j][k] /= 4. Currently skipping it __entirely__
    // AddNoise has (correct) variance PSD/2
    whitening_test(data.tdi->X, data.noise->C[0][0], data.NFFT);
    
    struct InstrumentModel inst = {0};
    initialize_instrument_model(&orbit, &data, &inst);
    // TODO: generate_instrument_noise_model has C[i][j][k] /= 2.
    // whether correct or not will depend on data validation

    /* get initial likelihood */
    invert_noise_covariance_matrix(inst.psd);
    inst.logL = my_noise_log_likelihood(&data, inst.psd);
    // TODO: my_noise_log_likelihood just uses invC directly. no 2?

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
    /*
    struct Noise data.noisescaleogram = {0};
    alloc_noise(&scaleogram, data.Nlayer*data.wdm->NT, data.Nlayer, data.Nchannel);
    scaleogram.kmin = data.wdm->kmin;
    int k;
    for (size_t i = 0; i < data.wdm->NT; i++)
        for (size_t j=0; j < data.wdm->NF; j++) {
            wavelet_pixel_to_index(data.wdm,i,j,&k); 
            scaleogram.f[k] = (data.lmin + j)*data.wdm->df;
        }
    */
    // TODO: note normalization of cov is /= 8. here
    initialize_instrument_model_wavelet(&orbit, &data, &inst);
    generate_full_dynamic_covariance_matrix(data.wdm, &inst, NULL, NULL, data.noise);
    alloc_tdi(data.tdi, data.N, 3);
    // TODO: are we handling real/imag parts of wavelets correctly?
    // see https://arxiv.org/pdf/2511.10632
    MyAddNoiseWavelet(&data,data.tdi);
    whitening_test_wavelet(&data);

    /* get initial likelihood */
    invert_noise_covariance_matrix(data.noise);
    inst.logL = my_noise_log_likelihood_wavelet(&data, data.noise); 

    printf("Running likelihood benchmarks...\n");
    start = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = my_noise_log_likelihood_wavelet(&data, data.noise);
    }
    stop = clock();
    logL1 = inst.logL;
    start2 = clock();
    for (int i=0; i<NTOT; i++) {
        inst.logL = noise_log_likelihood_wavelet(&data, data.noise);
    }
    stop2 = clock();
    logL2 = inst.logL;

    print_noise_model_dynamic(&data, data.noise, "./debug-data/dbg_wavelet_input_scaleogram.dat");
    strcpy(data.dataDir, "./debug-data/");
    data.qmin = data.lmin;
    if(!flags.strainData) wavelet_layer_to_fourier_transform(&data);
    print_data(&data, &flags);

    printf("my_noise_log_likelihood_wavelet (=%g) :\n\t%d evals in %g seconds\n",logL1,NTOT,(double)(stop-start)/CLOCKS_PER_SEC);
    printf("noise_log_likelihood_wavelet (=%g) :\n\t%d evals in %g seconds\n",logL2,NTOT,(double)(stop2-start2)/CLOCKS_PER_SEC);



    return 0;
}

