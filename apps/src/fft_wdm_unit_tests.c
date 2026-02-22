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
 @file fft_wdm_unit_tests.c
 \brief functions to check the FFT and WDM transform logic and ensure they pass some basic tests
 */

#include <glass_utils.h>
#include <glass_noise.h>
#include <time.h>
#include <assert.h>

#define NFFT_TEST (1536*10)
#define CREATE_DEBUG_FILES (false)
// Note that for us, NF is essentially harcoded by the choice of WAVELET_DURATION
// these tests were written for an NF=1536

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

int write_time_data(double dt, int N, double* data, char* fname) {
    FILE *fptr = NULL;
    fptr = fopen(fname,"w");
    if (!fptr) {
        fprintf(stderr, "Couldn't open %s for writing!\n", fname);
        return 1;
    }
    for (int i=0; i<N; i++) {
            double f = i*dt;
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg" ,data[i]);
            fprintf(fptr,"\n");
    }
    return 0;
}


bool test_array_equality(double* ta, double* tb, int N, double atol, double rtol, int* fail_counter, int* test_counter, char* test_name) {
    printf("Test %d: %s\n", ++(*test_counter), test_name);
    for (int i=0; i<N; i++) {
        if (fabs(ta[i] - tb[i]) > atol + rtol*fabs(tb[i])) {
            printf("\tarrays unequal within tolerance at index %d, got %lg != %lg\n", i, ta[i], tb[i]);
            *fail_counter++;
            return false;
        }
    }
    printf("\tpass\n");
    return true;
}

/* unit test conventions!
 * I'll try to keep output minimal. Let's print name of test and either pass or fail with error
 */

int main(int argc, char *argv[])
{

    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    // start test block
    int test_counter_block = 0;
    int test_counter_total = 0;
    int num_failed_block = 0;
    int num_failed_total= 0;
    double atol = 1e-10;
    double rtol = 0.0;


    fprintf(stdout, "\n================= FFT CONVENTION CHECKS ================\n");

    // test what FFT coeffs are
    double test_data[NFFT_TEST] = {0};
    double ref_data[NFFT_TEST] = {0};
    double ref_fft_data[NFFT_TEST] = {0};
    double test_data_copy[NFFT_TEST] = {0};
    test_data[0] = 1.0;
    // these copies are just for reference, won't be touched
    ref_data[0]  = 1.0;
    for (int i=0; i<NFFT_TEST; i++) ref_fft_data[i] = 1.0;

    bool ok = false;
    glass_forward_real_fft(test_data, NFFT_TEST);
    ok = test_array_equality(test_data,
            ref_fft_data,
            NFFT_TEST,
            atol,
            rtol,
            &num_failed_block,
            &test_counter_block,
            "Real FFT of impulse matches analytic");

    glass_inverse_real_fft(test_data, NFFT_TEST);
    ok = test_array_equality(test_data,
            ref_data,
            NFFT_TEST,
            atol,
            rtol,
            &num_failed_block,
            &test_counter_block,
            "IFFTR of FFTR(impulse)");

    // reset timeseries
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    printf("Test %d: complex IFFT of complex FFT(impulse)\n", ++test_counter_block);
    failed = false;
    glass_forward_complex_fft(test_data, NFFT_TEST/2-1);
    // TODO: check cx fft coefficients match reference here
    glass_inverse_complex_fft(test_data, NFFT_TEST/2-1);
    ok = test_array_equality(test_data,
            ref_data,
            NFFT_TEST,
            atol,
            rtol,
            &num_failed_block,
            &test_counter_block,
            "complex IFFT of complex FFT(impulse)");

    printf("\n\nTest block finished. Failed %d/%d\n\n", num_failed_block, test_counter_block);
    test_counter_total += test_counter_block;
    num_failed_total   += num_failed_block;

    fprintf(stdout, "\n================= WDM CONVENTION CHECKS ================\n");
    // reset test counters for new block
    num_failed_block = 0;
    test_counter_block = 0;

    // timeseries, impulse again
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    struct Wavelets wdm = {0};
    // NFFT_TEST == NF*NT
    // T == NFFT_TEST*LISA_CADENCE
    double T = NFFT_TEST*LISA_CADENCE;
    initialize_wavelet(&wdm, T);
    wavelet_transform(&wdm, test_data);
    // TODO: test window definitions
    // Limits on WD windows are only from Necula et. al eq (7)
    // TODO: test WDM coeffs comparison to reference
    // TODO: test WDM inversion to be correct

    if (CREATE_DEBUG_FILES)
        write_wdm_data(&wdm, test_data, "./dbg_wdm_impulse.dat");

    wavelet_transform_inverse_fourier(&wdm, test_data);
    if (CREATE_DEBUG_FILES) {
        write_fft_data(1./T, NFFT_TEST/2, test_data, "./dbg_wdmfft_impulse.dat");
        write_fft_data(1./T, NFFT_TEST/2, test_data_copy, "./dbg_fft_impulse.dat");
    }
    glass_inverse_real_fft(test_data, NFFT_TEST);
    if (CREATE_DEBUG_FILES)
        write_time_data(T/NFFT_TEST, NFFT_TEST, test_data, "./dbg_invwdm_impulse.dat");
    // try to fourier directly
    memset(test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    // could window here, but let's not for reversibility
    /*
    glass_forward_real_fft(test_data, NFFT_TEST);
    wavelet_transform_fourier(&wdm, test_data);
    write_wdm_data(&wdm, test_data, "./dbg_fftwdm_impulse.dat");
    wavelet_transform_inverse_fourier(&wdm, test_data);
    write_fft_data(1./T, NFFT_TEST/2, test_data, "./dbg_invfftwdm_impulse.dat");
    */
    /*
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
    */
    printf("\n\nTest block finished. Failed %d/%d\n\n", num_failed_block, test_counter_block);
    test_counter_total += test_counter_block;
    num_failed_total   += num_failed_block;

    printf("\n\nAll tests finished! Failed %d/%d\n\n", num_failed_total, test_counter_total);
    return 0;
}

