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

#define N_FREQ_LAYERS ((int)WAVELET_DURATION / LISA_CADENCE)
#define NFFT_TEST (N_FREQ_LAYERS*10)
#define CREATE_DEBUG_FILES (true)
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
            k = j*wdm->NT + i;
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

struct UnitTestBlockInfo {
    int test_counter;
    int fail_counter;
    double atol;
    double rtol;
    char* block_name;
};

void print_test_block_stats(struct UnitTestBlockInfo* testinfo) {
    printf("\n\n%s finished. Failed %d/%d\n\n", testinfo->block_name, testinfo->fail_counter, testinfo->test_counter);
}

void print_test_blocks_summary_stats(struct UnitTestBlockInfo** testinfos, int N) {
    int total_tests = 0;
    int total_fails = 0;
    for (int i=0; i<N; i++) {
        total_tests += testinfos[i]->test_counter;
        total_fails += testinfos[i]->fail_counter;
    }
    printf("\n\nAll tests finished. Failed %d/%d\n\n", total_fails, total_tests);
}
bool test_array_equality_from_index(double* ta, double* tb, int N, int istart, struct UnitTestBlockInfo* testinfo, char* test_name) {
    printf("Test %d: %s\n", ++(testinfo->test_counter), test_name);
    for (int i=istart; i<N; i++) {
        if (fabs(ta[i] - tb[i]) > testinfo->atol + testinfo->rtol*fabs(tb[i])) {
            printf("\tarrays unequal within tolerance at index %d, got %lg != %lg\n", i, ta[i], tb[i]);
            (testinfo->fail_counter)++;
            return false;
        }
    }
    printf("\tpass\n");
    return true;
}

bool test_array_equality(double* ta, double* tb, int N, struct UnitTestBlockInfo* testinfo, char* test_name) {
    return test_array_equality_from_index(ta, tb, N, 0, testinfo, test_name);
}

/* unit test conventions!
 * I'll try to keep output minimal. Let's print name of test and either pass or fail with error
 */

void fill_with_impulse(double* data, int N) {
    memset(data,0,sizeof(double)*N);
    data[0] = 1.0;
}
void fill_with_sines(double* data, int N, double* freqs, int Nfreqs) {
    memset(data,0,sizeof(double)*N);
    for (int i = 0; i < N; i++) {
        double t = LISA_CADENCE*i;
        for (int j = 0; j<Nfreqs; j++) {
            data[i] += sin(2*M_PI*freqs[j]*t);
        }
    }
}

int main(int argc, char *argv[])
{
    _Static_assert(N_FREQ_LAYERS % 2 == 0, "Number of frequency layers must be even for our WDM implementation!");
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    // start test block
    struct UnitTestBlockInfo total_tests = {0};
    struct UnitTestBlockInfo fft_tests = {0};
    struct UnitTestBlockInfo fft_outplace_tests = {0};
    struct UnitTestBlockInfo wdm_tests = {0};
    struct UnitTestBlockInfo wdm_sine_tests = {0};
    total_tests.block_name = "All tests";
    fft_tests.block_name = "FFT inplace tests";
    fft_outplace_tests.block_name = "FFT outplace tests";
    wdm_tests.block_name = "WDM tests";
    wdm_sine_tests.block_name = "WDM sine tests";
    fft_tests.atol = 1e-10;
    fft_tests.rtol = 1e-5;
    fft_outplace_tests.atol = 1e-10;
    fft_outplace_tests.rtol = 1e-5;
    wdm_tests.atol = 1e-10;
    wdm_tests.rtol = 1e-5;
    wdm_sine_tests.atol = 1e-8;
    wdm_sine_tests.rtol = 1e-5;


    // test what FFT coeffs are
    double test_data[NFFT_TEST] = {0};
    double test_data_cx[2*NFFT_TEST] = {0};
    test_data[0]    = 1.0;
    test_data_cx[0] = 1.0;
    double test_fft_data[NFFT_TEST+2] = {0};
    double test_cx_fft_data[2*NFFT_TEST] = {0};
    double ref_data[NFFT_TEST] = {0};
    double ref_data_cx[2*NFFT_TEST] = {0};
    double ref_fft_data[NFFT_TEST+2] = {0};
    double ref_cx_fft_data[2*NFFT_TEST] = {0};
    double test_data_copy[NFFT_TEST] = {0};
    double ref_wdm[NFFT_TEST];
    // these copies are just for reference, won't be touched
    ref_data[0]  = 1.0;
    ref_data_cx[0]  = 1.0;
    for (int i=0; i<NFFT_TEST/2+1; i++) ref_fft_data[2*i] = 1.0; // real part is 1, imag is 0
    for (int i=0; i<NFFT_TEST; i++) ref_cx_fft_data[2*i] = 1.0; // real part is 1, imag is 0

    fprintf(stdout, "\n================= FFT CONVENTION CHECKS - INPLACE ================\n");
    bool ok = false;
    glass_forward_real_fft(test_data, NFFT_TEST);
    ok = test_array_equality(test_data,
            ref_fft_data,
            NFFT_TEST,
            &fft_tests,
            "Real FFT of impulse matches analytic");

    glass_inverse_real_fft(test_data, NFFT_TEST);
    ok = test_array_equality_from_index(test_data,
            ref_data,
            NFFT_TEST,
            1,
            &fft_tests,
            "IFFTR of FFTR(impulse)");

    // reset timeseries
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;

    glass_forward_complex_fft(test_data_cx, NFFT_TEST);
    // check cx fft coefficients match reference here
    ok = test_array_equality(test_data_cx,
            ref_cx_fft_data,
            2*NFFT_TEST,
            &fft_tests,
            "complex FFT(impulse) matches analytic");
    glass_inverse_complex_fft(test_data_cx, NFFT_TEST);
    ok = test_array_equality(test_data_cx,
            ref_data_cx,
            2*NFFT_TEST,
            &fft_tests,
            "complex IFFT of complex FFT(impulse)");

    print_test_block_stats(&fft_tests);

    fprintf(stdout, "\n================= FFT CONVENTION CHECKS - OUTPLACE ================\n");
    // reset timeseries
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    glass_forward_real_fft_outplace(test_data, test_fft_data, NFFT_TEST);
    ok = test_array_equality(test_fft_data,
            ref_fft_data,
            NFFT_TEST+2,
            &fft_outplace_tests,
            "Outplace real FFT of impulse matches analytic");

    glass_inverse_real_fft_outplace(test_fft_data, test_data, NFFT_TEST);
    ok = test_array_equality(test_data,
            ref_data,
            NFFT_TEST,
            &fft_outplace_tests,
            "Outplace IFFTR of FFTR(impulse)");

    // reset timeseries
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;

    glass_forward_complex_fft_outplace(test_data_cx, test_cx_fft_data, NFFT_TEST);
    // check cx fft coefficients match reference here
    ok = test_array_equality(test_cx_fft_data,
            ref_cx_fft_data,
            2*NFFT_TEST,
            &fft_outplace_tests,
            "Outplace complex FFT(impulse) matches analytic");
    glass_inverse_complex_fft_outplace(test_cx_fft_data, test_data_cx, NFFT_TEST);
    ok = test_array_equality(test_data_cx,
            ref_data_cx,
            2*NFFT_TEST,
            &fft_outplace_tests,
            "Outplace complex IFFT of complex FFT(impulse)");

    print_test_block_stats(&fft_outplace_tests);
    fprintf(stdout, "\n================= WDM CONVENTION CHECKS ================\n");

    // timeseries, impulse again
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    struct Wavelets wdm = {0};
    // NFFT_TEST == NF*NT
    // T == NFFT_TEST*LISA_CADENCE
    double T = NFFT_TEST*LISA_CADENCE;
    initialize_wavelet(&wdm, T);
    wavelet_transform_timefreq(&wdm, test_data);

    // Could test window definitions
    // Limits on WD windows are only from Necula et. al eq (7)
    // skip for now, seems hard

    // WDM comparison with olitas
    double olitas_wdm[NFFT_TEST] = {0};
    FILE* f = fopen("./olitas_wdm_impulse.dat","r");
    if (!f) {
        printf("No Olitas.jl output to compare against, skipping WDM code comparison...");
    } else {
        for (int j=0; j < wdm.NT; j++)
            for (int i=0; i< wdm.NF; i++)
                fscanf(f, "%lf ", &olitas_wdm[i*wdm.NT+j]); // note transpose during read
        ok = test_array_equality(test_data,
                olitas_wdm,
                NFFT_TEST,
                &wdm_tests,
                "Olitas.jl WDM(impulse) vs our WDM(impulse)");
    }
    fclose(f);

    if (CREATE_DEBUG_FILES)
        write_wdm_data(&wdm, test_data, "./dbg_wdm_impulse.dat");

    wavelet_transform_inverse_freq(&wdm, test_data, test_fft_data);
    
    ok = test_array_equality(test_fft_data,
            ref_fft_data,
            NFFT_TEST+2,
            &wdm_tests,
            "WDM(impulse) to FFT vs reference FFT");
    if (CREATE_DEBUG_FILES) {
        write_fft_data(1./T, NFFT_TEST/2+1, test_fft_data, "./dbg_wdmfft_impulse.dat");
    }
    glass_inverse_real_fft_outplace(test_fft_data, test_data, NFFT_TEST);
    if (CREATE_DEBUG_FILES)
        write_time_data(T/NFFT_TEST, NFFT_TEST, test_data, "./dbg_invwdm_impulse.dat");
    ok = test_array_equality(test_data,
            ref_data,
            NFFT_TEST,
            &wdm_tests,
            "IWDM of WDM(impulse)");


    // Need to test
    // wavelet_transform_from_table
    // active_wavelet_list
    // wavelet_window_frequency -- probably remove

    // wavelet_transform_segment
        // this one takes a short freq segment and transforms to wdm
        // assumes it contains only one freq layer
        // is length N. what is N? WDM transform wants length Nt...
        // is only used in MBH waveform
        // for now force length Nt
    // take one freq layer of olitas wdm, turn it into FFT
    double test_wdm_full[NFFT_TEST] = {0};
    double test_wdm_layer[wdm.NT];
    int test_layer = 5;
    for (int i=0; i<wdm.NT; i++) {
        int k = test_layer*wdm.NT + i;
        test_wdm_full[k] = olitas_wdm[k];
        test_wdm_layer[i] = olitas_wdm[k];
    }
    wavelet_transform_inverse_freq(&wdm, test_wdm_full, test_fft_data);
    if (CREATE_DEBUG_FILES) {
        write_fft_data(1./T, NFFT_TEST/2+1, test_fft_data, "./dbg_wdmfft_onelayer.dat");
    }
    // extract a short FFT segment centered on the layer
    int center = test_layer * wdm.NT / 2;
    int Nseg = wdm.NT; // segment length in complex bins
    double short_fft[2*Nseg];
    for (int i = 0; i < Nseg; i++) {
        int full_idx = center - Nseg/2 + i;
        if (full_idx >= 0 && full_idx <= NFFT_TEST/2) {
            short_fft[2*i]     = test_fft_data[2*full_idx];
            short_fft[2*i + 1] = test_fft_data[2*full_idx + 1];
        } else {
            short_fft[2*i]     = 0.0;
            short_fft[2*i + 1] = 0.0;
        }
    }
    // now forward it with wavelet_transform_segment
    wavelet_transform_freq_segment(&wdm, Nseg, test_layer, short_fft);

    ok = test_array_equality(short_fft,
            test_wdm_layer,
            wdm.NT,
            &wdm_tests,
            "wavelet_transform_segment (layer 5) vs WDM(impulse) (layer 5)");

    // Build full scalogram by calling my_wavelet_transform_segment on every layer
    // and compare against the full wavelet_transform_freq result
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    // get reference full WDM
    memcpy(ref_wdm, test_data, NFFT_TEST*sizeof(double));
    wavelet_transform_timefreq(&wdm, ref_wdm);
    // get full FFT for extracting segments
    double full_fft[NFFT_TEST+2];
    glass_forward_real_fft_outplace(test_data, full_fft, NFFT_TEST);
    // assemble scalogram layer by layer
    double assembled_wdm[NFFT_TEST];
    memset(assembled_wdm, 0, NFFT_TEST*sizeof(double));
    for (int m = 0; m <= wdm.NF; m++) {
        int seg_center = m * wdm.NT / 2;
        int Nseg = wdm.NT;
        double seg_fft[2*Nseg];
        for (int i = 0; i < Nseg; i++) {
            int full_idx = seg_center - Nseg/2 + i;
            if (full_idx >= 0 && full_idx <= NFFT_TEST/2) {
                seg_fft[2*i]     = full_fft[2*full_idx];
                seg_fft[2*i + 1] = full_fft[2*full_idx + 1];
            } else {
                seg_fft[2*i]     = 0.0;
                seg_fft[2*i + 1] = 0.0;
            }
        }
        wavelet_transform_freq_segment(&wdm, Nseg, m, seg_fft);
        // DC and Nyquist layers both map into column 0
        if (m == 0) {
            for (int n = 0; n < wdm.NT; n += 2)
                assembled_wdm[n] = seg_fft[n];
        } else if (m == wdm.NF) {
            for (int n = 0; n < wdm.NT; n += 2)
                assembled_wdm[n + 1] = seg_fft[n + 1];
        } else {
            for (int n = 0; n < wdm.NT; n++)
                assembled_wdm[m*wdm.NT + n] = seg_fft[n];
        }
    }
    ok = test_array_equality(assembled_wdm,
            ref_wdm,
            NFFT_TEST,
            &wdm_tests,
            "layer-by-layer assembly vs full wavelet_transform_freq");


    print_test_block_stats(&wdm_tests);

    fprintf(stdout, "\n================= WDM CONVENTION CHECKS -- SINE WAVES ================\n");
    double DeltaF = wdm.df;
    double sine_freqs[] = {
        10      *DeltaF,
        20.5    *DeltaF,
        1000.0  *DeltaF,
        1000.5  *DeltaF,
    };
    int Nfreqs = sizeof(sine_freqs) / sizeof(sine_freqs[0]);
    // timeseries, sine waves now
    fill_with_sines(test_data, NFFT_TEST, sine_freqs, Nfreqs);
    fill_with_sines(ref_wdm,   NFFT_TEST, sine_freqs, Nfreqs);
    fill_with_sines(ref_data,  NFFT_TEST, sine_freqs, Nfreqs);

    wavelet_transform_timefreq(&wdm, test_data);
    wavelet_transform_timefreq(&wdm, ref_wdm);
    // TODO: compare to reference?
    
    glass_forward_real_fft_outplace(ref_data, ref_fft_data, NFFT_TEST);
    wavelet_transform_inverse_freq(&wdm, test_data, test_fft_data);
    
    ok = test_array_equality(test_fft_data,
            ref_fft_data,
            NFFT_TEST+2,
            &wdm_sine_tests,
            "WDM(sines) to FFT vs reference FFT");
    // wavelet_transform_by_layers
    // note that this is only used in the UCB wavefrom, which
    // gets heterodyned and downsampled before calling this function
    //
    // We do an analagous thing here to match! 

    test_layer = 999;
    int test_Nlayer = 5;

    double by_layer_sine_freqs[] = {
        1001.0 * DeltaF, 
       // 1001.5 * DeltaF, 
       // note only integer layers work at the moment since we're not heterodyning the test correctly
    };
    /*
     *
     *
      for (int i = 0; i < N_ds; i++) {
          double t = dt_ds * i;
          double re = 0, im = 0;
          for (int j = 0; j < by_layer_Nfreqs; j++) {
              // sin(2π f t) · e^{-2π i f₀ t} kept as complex
              double ang_sig = 2*M_PI*by_layer_sine_freqs[j]*t;
              double ang_het = 2*M_PI*f0*t;
              // sin(a)·e^{-ib} = (cos(b)sin(a)) + i(-sin(b)sin(a))
              re += sin(ang_sig)*cos(ang_het);
              im += -sin(ang_sig)*sin(ang_het);
          }
          het_data_cx[2*i]   = re;
          het_data_cx[2*i+1] = im;
      }
    */
    int by_layer_Nfreqs = sizeof(by_layer_sine_freqs) / sizeof(by_layer_sine_freqs[0]);

    double ref_wdm_sines[NFFT_TEST];
    fill_with_sines(ref_wdm_sines, NFFT_TEST, by_layer_sine_freqs, by_layer_Nfreqs);
    wavelet_transform_timefreq(&wdm, ref_wdm_sines);

    // Heterodyned + downsampled input like UCB: carrier f0 shifts layer
    // (test_layer-1) to DC, new cadence dt_ds = wdm.dt/(Nlayers+1)
    int N_ds = wdm.NT * (test_Nlayer + 1);
    int decimation = wdm.NF / (test_Nlayer + 1);
    if (wdm.NF % (test_Nlayer + 1)) {
        fprintf(stderr, "Couldn't evenly divide the number of frequencies by Nlayers+1 ! Can't decimate test signal cleanly");
        exit(-1);
    }
    double dt_ds = decimation * LISA_CADENCE;
    double f0 = (test_layer - 1) * wdm.df;
    double het_data[N_ds];
    memset(het_data, 0, N_ds*sizeof(double));
    for (int i = 0; i < N_ds; i++) {
        double t = dt_ds * i;
        for (int j = 0; j < by_layer_Nfreqs; j++)
            het_data[i] += sin(2*M_PI*(by_layer_sine_freqs[j] - f0) * t);
    }

    double window[wdm.NT/2 + 1];
    build_wdm_filter_freq(window, wdm.NF, wdm.NT, wdm.A, true);
    wavelet_transform_timefreq_by_layers(&wdm, test_layer, test_Nlayer, window, het_data);

    double wdm_crop_sines[test_Nlayer*wdm.NT];
    for (int n = 0; n < wdm.NT; n++)
        for (int l = 0; l < test_Nlayer; l++)
            wdm_crop_sines[n*test_Nlayer + l] = ref_wdm_sines[(test_layer + l)*wdm.NT + n];
    ok = test_array_equality(het_data,
            wdm_crop_sines,
            test_Nlayer*wdm.NT,
            &wdm_sine_tests,
            "wavelet_transform_timefreq_by_layers(heterodyned, downsampled sine) vs WDM(sine)");

    print_test_block_stats(&wdm_sine_tests);

    struct UnitTestBlockInfo* all_test_blocks[] = {
        &fft_tests,
        &fft_outplace_tests,
        &wdm_tests,
        &wdm_sine_tests,
    };
    int num_blocks = sizeof(all_test_blocks)/sizeof(all_test_blocks[0]);
    print_test_blocks_summary_stats(all_test_blocks, num_blocks);

    return 0;
}

