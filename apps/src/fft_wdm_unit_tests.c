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

#define NFFT_TEST (1536*30)
#define NBINS 100
#define CREATE_DEBUG_FILES (true)
// Note that for us, NF is essentially harcoded by the choice of WAVELET_DURATION
// these tests were written for an NF=1536

struct UnitTestBlockInfo {
    int test_counter;
    int fail_counter;
    double atol;
    double rtol;
    char* block_name;
};

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

/* Wilson-Hilferty normal approximation: works for any dof >= 1 */
double chi2_sf(double x, int dof) {
    if (x <= 0.0) return 1.0;
    double nu = (double)dof;
    double z = (pow(x / nu, 1.0/3.0) - (1.0 - 2.0/(9.0*nu))) / sqrt(2.0/(9.0*nu));
    return 0.5 * erfc(z / M_SQRT2);
}
bool test_whitening_wdm_3ch(const struct TDI* tdi, const struct Noise* noise, int N, int Nskip, struct UnitTestBlockInfo* testinfo, char* testname)
{
    // whitening test
    double max_bin = 8.0;
    double min_bin = -max_bin;
    double delta_bin = (max_bin - min_bin)/NBINS;
    int k;

    double* psds[] = {
        noise->C[0][0],
        noise->C[1][1],
        noise->C[2][2],
    };
    double* datas[] = {
        tdi->X,
        tdi->Y,
        tdi->Z,
    };
    char* ch_names[] = {
        "X",
        "Y",
        "Z",
    };
    printf("Test %d: %s\n", ++(testinfo->test_counter), testname);
    bool all_ch_good = true;
    for (int ich=0; ich<3; ich++) {
        double* psd = psds[ich];
        double* data = datas[ich];
        double mean = 0.0;
        double var = 0.0;
        double wht;
        mean = 0.0;
        var = 0.0;
        int Nactive = 0;
        int counts[NBINS] = {0};

        // note: we skip DC bins to avoid division by 0
        for (int i=Nskip; i<N; i++) {
            double asd = sqrt(psd[i]);
            wht = data[i] / asd;
            if (isnan(wht)) {
                fprintf(stderr, "Encountered NaN at i=%d, stopping\n", i);
                break;
            }
            mean += wht;
            if (wht<min_bin) { counts[0]++; continue; }
            if (wht>max_bin) { counts[NBINS-1]++; continue; }
            int b = (int)((wht - min_bin) / delta_bin);
            assert(b >= 0 && b < NBINS);
            counts[b]++;
            Nactive++;
        }
        mean /= Nactive ;
        for (int i=Nskip; i<N; i++) {
            double asd = sqrt(psd[i]);
            wht = data[i] / asd;
            if (isnan(wht)) {
                fprintf(stderr, "Encountered NaN at i=%d, stopping\n", i);
                break;
            }
            var += (wht - mean)*(wht - mean);
        }
        var /= Nactive;

        printf("\n");
        // strictly speaking , let's just test 3sigma CI
        double sample_mean_sigma = sqrtf(1.0/Nactive);
        if (fabs(mean) < 3*sample_mean_sigma ) {
            printf("\t%s (mean) pass: |%lg| < %lg\n", ch_names[ich], mean, 3*sample_mean_sigma);
        } else {
            printf("\t%s (mean) fail, inconsistent within 3sigma. Tested |%lf| < %lf\n", ch_names[ich], mean, 3*sample_mean_sigma);
            all_ch_good = false;
        }

        double corrected_var = (Nactive) * var; // sample variance
        int dof = Nactive - 1;
        double chi2_upper = chi2_sf(corrected_var, dof);
        double chi2_lower = 1.0 - chi2_upper;
        double p = 2.0 * fmin(chi2_upper, chi2_lower);
        if (p >  0.0027) { // equivalent of 3sigma, basically 99.7%
            printf("\t%s (var) pass: var=%lg p=%lg\n", ch_names[ich], var, p);
        } else {
            printf("\t%s (var) fail, inconsistent with chi2 distribution at 99.7%% CI. Got var %lf p=%lg\n", ch_names[ich], var, p);
            all_ch_good = false;
        }

        /*
        // print out histogram
        printf("\t\tHistogram:\n");
        for (int i=0; i<NBINS; i++) {
            printf("\t\t\t%5.2lf %6d\n", min_bin + delta_bin*i, counts[i]);
        }
        */
    }
    if (!all_ch_good) {
        testinfo->fail_counter++;
        return false;
    } else {
        return true;
    }
}

bool test_whitening_fft_3ch(const struct TDI* tdi, const struct Noise* noise, int N, struct UnitTestBlockInfo* testinfo, char* testname)
{
    // whitening test
    double max_bin = 8.0;
    double min_bin = -max_bin;
    double delta_bin = (max_bin - min_bin)/NBINS;
    int k;

    double* psds[] = {
        noise->C[0][0],
        noise->C[1][1],
        noise->C[2][2],
    };
    double* datas[] = {
        tdi->X,
        tdi->Y,
        tdi->Z,
    };
    char* ch_names[] = {
        "X",
        "Y",
        "Z",
    };
    printf("Test %d: %s\n", ++(testinfo->test_counter), testname);
    bool all_ch_good = true;
    for (int ich=0; ich<3; ich++) {
        double* psd = psds[ich];
        double* data = datas[ich];
        double mean = 0.0;
        double var = 0.0;
        double wht;
        mean = 0.0;
        var = 0.0;
        int Nactive = 0;
        int counts[NBINS] = {0};

        // note: we skip DC bin to avoid division by 0
        for (int i=1; i<N; i++) {
            double asd = sqrt(psd[i]/2);
            wht = data[2*i] / asd;
            mean += wht;
            if (wht<min_bin) { counts[0]++; continue; }
            if (wht>max_bin) { counts[NBINS-1]++; continue; }
            int b = (int)((wht - min_bin) / delta_bin);
            assert(b >= 0 && b < NBINS);
            counts[b]++;
            Nactive++;
            wht = data[2*i+1] / asd;
            mean += wht;
            if (wht<min_bin) { counts[0]++; continue; }
            if (wht>max_bin) { counts[NBINS-1]++; continue; }
            b = (int)((wht - min_bin) / delta_bin);
            assert(b >= 0 && b < NBINS);
            counts[b]++;
            Nactive++;
        }
        mean /= Nactive ;
        for (int i=1; i<N; i++) {
            double asd = sqrt(psd[i]/2);
            wht = data[2*i] / asd;
            var += (wht - mean)*(wht - mean);
            wht = data[2*i+1] / asd;
            var += (wht - mean)*(wht - mean);
        }
        var /= Nactive;

        printf("\n");
        // strictly speaking , let's just test 3sigma CI
        double sample_mean_sigma = sqrtf(1.0/Nactive);
        if (fabs(mean) < 3*sample_mean_sigma ) {
            printf("\t%s (mean) pass: |%lg| < %lg\n", ch_names[ich], mean, 3*sample_mean_sigma);
        } else {
            printf("\t%s (mean) fail, inconsistent within 3sigma. Tested |%lf| < %lf\n", ch_names[ich], mean, 3*sample_mean_sigma);
            all_ch_good = false;
        }

        double corrected_var = (Nactive) * var; // sample variance
        int dof = Nactive - 1;
        double chi2_upper = chi2_sf(corrected_var, dof);
        double chi2_lower = 1.0 - chi2_upper;
        double p = 2.0 * fmin(chi2_upper, chi2_lower);
        if (p >  0.0027) { // equivalent of 3sigma, basically 99.7%
            printf("\t%s (var) pass: var=%lg p=%lg\n", ch_names[ich], var, p);
        } else {
            printf("\t%s (var) fail, inconsistent with chi2 distribution at 99.7%% CI. Got var %lf p=%lg\n", ch_names[ich], var, p);
            all_ch_good = false;
        }

        /*
        // print out histogram
        printf("\t\tHistogram:\n");
        for (int i=0; i<NBINS; i++) {
            printf("\t\t\t%5.2lf %6d\n", min_bin + delta_bin*i, counts[i]);
        }
        */
    }
    if (!all_ch_good) {
        testinfo->fail_counter++;
        return false;
    } else {
        return true;
    }
}

int main(int argc, char *argv[])
{
    _Static_assert(WAVELET_DURATION == 7680, "currently we have hardcoded the WAVELET_DURATION. Either fix it or go back");

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
    // These arrays scale with NFFT_TEST. Stash them in BSS so the test can
    // grow NT without blowing the stack.
    static double test_data[NFFT_TEST];
    static double test_data_cx[2*NFFT_TEST];
    test_data[0]    = 1.0;
    test_data_cx[0] = 1.0;
    static double test_fft_data[NFFT_TEST+2];
    static double test_cx_fft_data[2*NFFT_TEST];
    static double ref_data[NFFT_TEST];
    static double ref_data_cx[2*NFFT_TEST];
    static double ref_fft_data[NFFT_TEST+2];
    static double ref_cx_fft_data[2*NFFT_TEST];
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

    // save WDM result before inverse chain overwrites test_data
    static double our_wdm_impulse[NFFT_TEST];
    memcpy(our_wdm_impulse, test_data, NFFT_TEST*sizeof(double));

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
    static double olitas_wdm[NFFT_TEST];
    FILE* f = fopen("./olitas_wdm_impulse.dat","r");
    if (!f) {
        printf("No Olitas.jl output to compare against, skipping WDM code comparison...");
    } else {
        for (int j=0; j < wdm.NT; j++)
            for (int i=0; i< wdm.NF; i++)
                fscanf(f, "%lf ", &olitas_wdm[i*wdm.NT+j]); // note transpose during read
        ok = test_array_equality(our_wdm_impulse,
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
    static double test_wdm_full[NFFT_TEST];
    memset(test_wdm_full, 0, NFFT_TEST*sizeof(double));
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
    static double ref_wdm[NFFT_TEST];
    memcpy(ref_wdm, test_data, NFFT_TEST*sizeof(double));
    wavelet_transform_timefreq(&wdm, ref_wdm);
    // get full FFT for extracting segments
    static double full_fft[NFFT_TEST+2];
    glass_forward_real_fft_outplace(test_data, full_fft, NFFT_TEST);
    // assemble scalogram layer by layer
    static double assembled_wdm[NFFT_TEST];
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

    fprintf(stdout, "\n================= FFT/WDM PSD CHECKS ================\n");
    struct UnitTestBlockInfo var_tests = {0};
    var_tests.block_name = "WDM variance tests";
    var_tests.atol = 0.0;
    var_tests.rtol = 1e-2;


    // test FFT PSD definitions
    struct Data data2 = {0};
    data2.NFFT = NFFT_TEST/2 + 1;
    data2.Nchannel = 3;
    data2.Nlayer = 0;
    data2.T = wdm.NT*wdm.dt;
    data2.fmin=0;
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    initialize_interpolated_analytic_orbits(orbit, data2.T, data2.t0);
    struct InstrumentModel *inst_fft = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model(orbit, &data2, inst_fft);

    data2.noise = inst_fft->psd;

    struct TDI* testtdi_fft = malloc(sizeof(struct TDI));
    alloc_tdi(testtdi_fft, 2*data2.NFFT, data2.Nchannel);
    MyAddNoise(&data2, testtdi_fft);

    ok = test_whitening_fft_3ch(testtdi_fft,
            data2.noise,
            data2.NFFT,
            &var_tests,
            "FFT PSD test: InstrumentModel PSD whitens FFT-generated noise");

    // test WDM PSD definitions
    struct Data data3 = {0};
    data3.NFFT = NFFT_TEST/2 + 1;
    data3.Nchannel = 3;
    data3.Nlayer = wdm.NF;
    data3.T = wdm.NT*wdm.dt;
    data3.fmin=0;
    data3.wdm = &wdm;
    data3.lmax = data3.Nlayer;
    data3.lmin = 0;
    data3.noise = malloc(sizeof(struct Noise));
    alloc_noise(data3.noise, data3.Nlayer*wdm.NT, data3.Nlayer, data3.Nchannel);

    // Build the WDM PSD that whitens wavelet_transform_freq output by passing
    // each FFT-domain channel PSD through the new convolution helper.
    stationary_dft_psd_to_wdm_psd(&wdm, inst_fft->psd->C[0][0], data3.noise->C[0][0]);
    stationary_dft_psd_to_wdm_psd(&wdm, inst_fft->psd->C[1][1], data3.noise->C[1][1]);
    stationary_dft_psd_to_wdm_psd(&wdm, inst_fft->psd->C[2][2], data3.noise->C[2][2]);

    struct TDI* testtdi = malloc(sizeof(struct TDI));
    alloc_tdi(testtdi, data3.Nlayer*wdm.NT, data3.Nchannel);
    MyAddNoiseWavelet(&data3, testtdi);

    ok = test_whitening_wdm_3ch(testtdi,
            data3.noise,
            data3.Nlayer*wdm.NT,
            wdm.NT,
            &var_tests,
            "WDM PSD test: InstrumentModel PSD whitens WDM-generated noise");


    // now try to whiten FFT data with WDM PSD
    wavelet_transform_freq(&wdm, testtdi_fft->X, testtdi->X);
    wavelet_transform_freq(&wdm, testtdi_fft->Y, testtdi->Y);
    wavelet_transform_freq(&wdm, testtdi_fft->Z, testtdi->Z);


    ok = test_whitening_wdm_3ch(testtdi,
            data3.noise,
            data3.Nlayer*wdm.NT,
            wdm.NT,
            &var_tests,
            "WDM PSD test: InstrumentModel PSD whitens FFT-generated noise");

    // Same FFT->WDM whitening test, but build data3.noise via the production
    // initialize_instrument_model_wavelet + generate_full_dynamic_covariance_matrix
    // path. This exercises the retooled generate_instrument_noise_model_wavelet
    // end to end.
    {
        struct InstrumentModel *inst_w = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model_wavelet(orbit, &data3, inst_w);
        generate_full_dynamic_covariance_matrix(data3.wdm, inst_w, NULL, NULL, data3.noise);

        ok = test_whitening_wdm_3ch(testtdi,
                data3.noise,
                data3.Nlayer*wdm.NT,
                wdm.NT,
                &var_tests,
                "WDM PSD test (model path): InstrumentModel PSD whitens FFT-generated noise");

        free_instrument_model(inst_w);
    }

    // Repeat the FFT->WDM whitening test with the approximate WDM PSD
    // (single-bin lookup at layer center). Looser tolerance because the
    // approximation drops the phif convolution.
    stationary_dft_psd_to_wdm_psd_approx(&wdm, inst_fft->psd->C[0][0], data3.noise->C[0][0]);
    stationary_dft_psd_to_wdm_psd_approx(&wdm, inst_fft->psd->C[1][1], data3.noise->C[1][1]);
    stationary_dft_psd_to_wdm_psd_approx(&wdm, inst_fft->psd->C[2][2], data3.noise->C[2][2]);

    ok = test_whitening_wdm_3ch(testtdi,
            data3.noise,
            data3.Nlayer*wdm.NT,
            wdm.NT,
            &var_tests,
            "WDM PSD test (approx): InstrumentModel PSD whitens FFT-generated noise");

    free_tdi(testtdi_fft);
    free_tdi(testtdi);
    free_instrument_model(inst_fft);
    print_test_block_stats(&var_tests);

    struct UnitTestBlockInfo* all_test_blocks[] = {&fft_tests, &fft_outplace_tests, &wdm_tests, &var_tests};
    int num_blocks = sizeof(all_test_blocks)/sizeof(all_test_blocks[0]);
    print_test_blocks_summary_stats(all_test_blocks, num_blocks);

    return 0;
}

