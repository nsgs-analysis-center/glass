#include <stdio.h>
#include <glass_utils.h>

#define TOBS 400*WAVELET_DURATION
#define N (int)((TOBS)/(LISA_CADENCE))

int main(int argc, char** argv) {
    // make fake data

    double data[N];
    for (int i = 0; i < N; i++) {
        data[i] = sin(i*2*M_PI*0.1);
    }

    struct Wavelets wdm = {0};
    initialize_wavelet(&wdm, TOBS);

    wavelet_transform_timefreq(&wdm, data);
    
    for (int i=0; i<N; i++) {
        printf("%lg\n", data[i]);
    }

    return 0;
}
