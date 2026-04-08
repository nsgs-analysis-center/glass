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
            wavelet_pixel_to_index(wdm, i, j, &k);
            fprintf(fptr,"%lg ",t);
            fprintf(fptr,"%lg ",f);
            fprintf(fptr,"%lg",data[k]);
            fprintf(fptr,"\n");
        }
    }
    fclose(fptr);
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
    fclose(fptr);
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
    fclose(fptr);
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

bool test_array_equality(double* ta, double* tb, int N, struct UnitTestBlockInfo* testinfo, char* test_name) {
    printf("Test %d: %s\n", ++(testinfo->test_counter), test_name);
    for (int i=0; i<N; i++) {
        if (fabs(ta[i] - tb[i]) > testinfo->atol + testinfo->rtol*fabs(tb[i])) {
            printf("\tarrays unequal within tolerance at index %d, got %lg != %lg\n", i, ta[i], tb[i]);
            (testinfo->fail_counter)++;
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
    struct UnitTestBlockInfo total_tests = {0};
    struct UnitTestBlockInfo fft_tests = {0};
    struct UnitTestBlockInfo fft_outplace_tests = {0};
    struct UnitTestBlockInfo wdm_tests = {0};
    total_tests.block_name = "All tests";
    fft_tests.block_name = "FFT inplace tests";
    fft_outplace_tests.block_name = "FFT outplace tests";
    wdm_tests.block_name = "WDM tests";
    fft_tests.atol = 1e-10;
    fft_tests.rtol = 1e-5;
    fft_outplace_tests.atol = 1e-10;
    fft_outplace_tests.rtol = 1e-5;
    wdm_tests.atol = 1e-10;
    wdm_tests.rtol = 1e-5;


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
    ok = test_array_equality(test_data,
            ref_data,
            NFFT_TEST,
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
    // wavelet_transform_by_layers
        // this one appears to take time data, tukey window
        // perform wdm essentially assuming content only goes from j=jmin to Nlayers+jmin
        // used in UCB waveform
    test_layer = 5;
    int test_Nlayer = 3;
    // timeseries, impulse again
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    double window[wdm.NT/2 + 1];
    build_wdm_filter_freq(window, wdm.NF, wdm.NT, wdm.A, true);
    wavelet_transform_timefreq_by_layers(&wdm, test_layer, test_Nlayer, window, test_data);
    // output is now the first test_Nlayer*NT elements of test_data
    double wdm_crop[test_Nlayer*wdm.NT];
    for (int i=0; i < test_Nlayer*wdm.NT; i++)
        wdm_crop[i] = olitas_wdm[test_layer*wdm.NT + i];
    ok = test_array_equality(test_data,
            wdm_crop,
            test_Nlayer*wdm.NT,
            &wdm_tests,
            "wavelet_transform_timefreq_by_layers (5-7) vs WDM(impulse) (5-7)");


    // Build full scalogram by calling my_wavelet_transform_segment on every layer
    // and compare against the full wavelet_transform_freq result
    memset(&test_data, 0, NFFT_TEST*sizeof(double));
    test_data[0] = 1.0;
    // get reference full WDM
    double ref_wdm[NFFT_TEST];
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

    fprintf(stdout, "\n================= WDM VARIANCE / PSD CHECKS ================\n");
    struct UnitTestBlockInfo var_tests = {0};
    var_tests.block_name = "WDM variance tests";
    var_tests.atol = 0.0;
    var_tests.rtol = 0.15; // 15% tolerance for statistical test with Nreal realizations

    /*
     * Empirical WDM variance test:
     * Generate white noise in FFT domain (unit variance per complex bin),
     * transform to WDM, accumulate variance per layer.
     * This measures the actual transfer function from FFT PSD to WDM variance.
     */
    {
        int Nreal = 200; // number of realizations for variance estimation
        unsigned int seed = 12345;
        int Nt = wdm.NT;
        int Nf = wdm.NF;
        int ND = Nf * Nt;
        int NFFT = ND / 2 + 1;

        // accumulate variance per WDM layer (NF+1 layers: 0..NF)
        double *layer_var = calloc(Nf + 1, sizeof(double));
        int    *layer_count = calloc(Nf + 1, sizeof(int));
        double *fft_data = calloc(2 * NFFT, sizeof(double));
        double *wdm_data = calloc(ND, sizeof(double));

        for (int r = 0; r < Nreal; r++) {
            // Generate unit-variance white noise in FFT domain
            // Each complex bin: re ~ N(0, 1/sqrt(2)), im ~ N(0, 1/sqrt(2))
            // so E[|X|^2] = 1
            memset(fft_data, 0, 2 * NFFT * sizeof(double));
            for (int k = 1; k < NFFT - 1; k++) {
                fft_data[2*k]     = rand_r_N_0_1(&seed) / M_SQRT2;
                fft_data[2*k + 1] = rand_r_N_0_1(&seed) / M_SQRT2;
            }
            // DC and Nyquist are real
            fft_data[0] = rand_r_N_0_1(&seed);
            fft_data[2*(NFFT-1)] = rand_r_N_0_1(&seed);

            // Forward WDM transform
            memset(wdm_data, 0, ND * sizeof(double));
            wavelet_transform_freq(&wdm, fft_data, wdm_data);

            // Accumulate per-layer variance
            // DC layer (m=0): even rows of column 0
            for (int n = 0; n < Nt; n += 2) {
                layer_var[0] += wdm_data[n] * wdm_data[n];
                layer_count[0]++;
            }
            // Nyquist layer (m=NF): odd rows of column 0
            for (int n = 0; n < Nt; n += 2) {
                layer_var[Nf] += wdm_data[n + 1] * wdm_data[n + 1];
                layer_count[Nf]++;
            }
            // General layers (m=1..NF-1)
            for (int m = 1; m < Nf; m++) {
                for (int n = 0; n < Nt; n++) {
                    int k = n + m * Nt;
                    layer_var[m] += wdm_data[k] * wdm_data[k];
                    layer_count[m]++;
                }
            }
        }

        // Normalize to get empirical variance
        for (int m = 0; m <= Nf; m++) {
            if (layer_count[m] > 0)
                layer_var[m] /= layer_count[m];
        }

        // Compute analytical prediction: Σ phif_fwd[|j|]^2 / (2*Nt^2) for general layers
        // For DC/Nyquist: factor of sqrt(2) in transform, and only Nt/2 coefficients
        double *phif_fwd = wdm.window_freq_forward;
        double window_sum = 0.0;
        for (int j = 0; j <= Nt/2; j++) {
            double w = (j == 0) ? 1.0 : 2.0;
            window_sum += w * phif_fwd[j] * phif_fwd[j];
        }
        double predicted_var_general = window_sum / (2.0 * (double)Nt * (double)Nt);
        // DC/Nyquist: sqrt(2) factor in transform, so variance is 2x, but only re part
        // and only Nt/2 coefficients (even rows). The factor works out differently.
        double predicted_var_dc = window_sum / ((double)Nt * (double)Nt);

        if (CREATE_DEBUG_FILES) {
            FILE *fvar = fopen("./dbg_wdm_layer_variance.dat", "w");
            if (fvar) {
                fprintf(fvar, "# layer  empirical_var  predicted_var_general  predicted_var_dc  Nsamples\n");
                for (int m = 0; m <= Nf; m++) {
                    double predicted = (m == 0 || m == Nf) ? predicted_var_dc : predicted_var_general;
                    fprintf(fvar, "%d  %lg  %lg  %lg  %d\n", m, layer_var[m], predicted, predicted_var_dc, layer_count[m]);
                }
                fclose(fvar);
            }
        }

        // Check general layers (skip a few near edges)
        printf("Test %d: WDM empirical variance vs analytical (general layers)\n", ++(var_tests.test_counter));
        printf("\tpredicted_var_general = %lg (window_sum=%lg, Nt=%d)\n", predicted_var_general, window_sum, Nt);
        printf("\tpredicted_var_dc      = %lg\n", predicted_var_dc);
        int var_fail = 0;
        for (int m = 2; m < Nf - 1; m++) {
            double rel_err = fabs(layer_var[m] - predicted_var_general) / predicted_var_general;
            if (rel_err > var_tests.rtol) {
                if (var_fail < 5) // limit output
                    printf("\tlayer %d: empirical=%lg predicted=%lg rel_err=%lg\n",
                           m, layer_var[m], predicted_var_general, rel_err);
                var_fail++;
            }
        }
        if (var_fail > 0) {
            printf("\t%d/%d general layers failed (>%.0f%% relative error)\n", var_fail, Nf-3, var_tests.rtol*100);
            var_tests.fail_counter++;
        } else {
            printf("\tpass (all general layers within %.0f%%)\n", var_tests.rtol*100);
        }

        // Print the ratio C_fft/C_wdm for reference (what the normalization factor should be)
        printf("\n\tReference: for unit FFT PSD (C_fft=1), WDM variance = %lg\n", predicted_var_general);
        printf("\tThis means C_wdm = C_fft * %lg\n", predicted_var_general);
        printf("\tOr equivalently, C_wdm = C_fft / %lg\n", 1.0/predicted_var_general);
        printf("\tCompare: code currently uses C_wdm = C_fft / 8\n");
        printf("\tNf*Nt = %d, 1/(Nf*Nt) = %lg\n", Nf*Nt, 1.0/(Nf*Nt));

        free(layer_var);
        free(layer_count);
        free(fft_data);
        free(wdm_data);
    }

    print_test_block_stats(&var_tests);

    struct UnitTestBlockInfo* all_test_blocks[] = {&fft_tests, &fft_outplace_tests, &wdm_tests, &var_tests};
    int num_blocks = sizeof(all_test_blocks)/sizeof(all_test_blocks[0]);
    print_test_blocks_summary_stats(all_test_blocks, num_blocks);

    return 0;
}

