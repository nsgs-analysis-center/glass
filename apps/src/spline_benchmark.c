#include <glass_math.h>
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

    struct CubicSplineEvenSampling *splineES = alloc_cubic_spline_even_sampling(2*N);

    initialize_cubic_spline_even_sampling(splineES, x, y, 1.0);

    clock_t start = clock();
    for (int i=0; i<N; i++) {
        yp[i] = spline_interpolation_even_sampling(splineES, i+0.1);
    }
    clock_t end = clock();
    clock_t start2 = clock();
    for (int i=0; i<N; i++) {
        yp2[i] = spline_interpolation(splineES->cspline, i+0.1);
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
    printf("CubicSplineEvenSampling\n\ttime to interpolate %d times (in input data of size %d): %lf\n", N, 2*N, (double)(end-start)/CLOCKS_PER_SEC);
    printf("\terror: %lg (%lg per point)\n", errorCSES, errorCSES/N);
    printf("CubicSpline\n\ttime to interpolate %d times (in input data of size %d): %lf\n", N, 2*N, (double)(end2-start2)/CLOCKS_PER_SEC);
    printf("\terror: %lg (%lg per point)\n", errorCS, errorCS/N);


    free_cubic_spline_even_sampling(splineES);
    return 0;
};
