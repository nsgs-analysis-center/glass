#include <glass_utils.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define N 70000

int main(int argc, char** argv) {
    double x[2*N] = {0};
    double y[2*N] = {0};
    double yp[N] = {0};
    double yp2[N] = {0};

    for (int i=0; i<2*N; i++) {
        x[i] = i;
        y[i] = sin(2*M_PI*0.1*i);
    }

    struct CubicSpline *spline = alloc_cubic_spline(2*N);

    initialize_cubic_spline(spline, x, y, SPLINE_EVEN_SAMPLED);

    clock_t start = clock();
    for (int i=0; i<N; i++) {
        yp[i] = spline_interpolation(spline, i+0.1);
    }
    clock_t end = clock();
    spline->interp_type = SPLINE_BINARY_SEARCH;
    spline->index_lookup = &binary_search;
    clock_t start2 = clock();
    for (int i=0; i<N; i++) {
        yp2[i] = spline_interpolation(spline, i+0.1);
    }
    clock_t end2 = clock();
    double errorCS = 0.0;
    double errorCSES = 0.0;
    for (int i=0; i<N; i++) {
        double ytrue = sin(2*M_PI*0.1*(i+0.1));
        errorCS += pow(ytrue - yp2[i], 2); 
        errorCSES += pow(ytrue - yp[i], 2); 
    }
    errorCS = sqrtf(errorCS);
    errorCSES = sqrtf(errorCSES);
    printf("CubicSpline with (%s) interpolation:\n\ttime to interpolate %d times (in input data of size %d): %lf\n", SPLINE_INTERPOLATION_NAMES[SPLINE_EVEN_SAMPLED], N, 2*N, (double)(end-start)/CLOCKS_PER_SEC);
    printf("\terror: %lg (%lg per point)\n", errorCSES, errorCSES/N);
    printf("CubicSpline with %s interpolation:\n\ttime to interpolate %d times (in input data of size %d): %lf\n",SPLINE_INTERPOLATION_NAMES[SPLINE_BINARY_SEARCH], N, 2*N, (double)(end2-start2)/CLOCKS_PER_SEC);
    printf("\terror: %lg (%lg per point)\n", errorCS, errorCS/N);


    free_cubic_spline(spline);
    return 0;
};
