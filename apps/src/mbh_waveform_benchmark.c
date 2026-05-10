
#include <time.h>

#include <glass_utils.h>
#include <glass_mbh.h>
#include "IMRPhenomD.h"
#include "IMRPhenomD_internals.h"

#define NTEST 10000

int main(int argc, char* argv[])
{
    if(argc!=1)
    {
        printf("Usage: ./mbh_waveform_benchmark\n");
        return 0;
    }
    
    clock_t start, end;
    double cpu_time_used;
    
    double t0   = 0;
    double Tobs = 31457280;
    int N = (int)(Tobs/LISA_CADENCE);

    /*
     Wavelet basis
     */
    struct Wavelets *wdm = malloc(sizeof(struct Wavelets));
    initialize_wavelet(wdm, Tobs);
    
    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    initialize_analytic_orbit(orbit);

    // get spacecraft ephemerides on the spline time grid
    initialize_interpolated_analytic_orbits(orbit, Tobs, t0);

    /*
     Set injection parameters for test signal
    */
    double *params = double_vector(11);

    double m1   = 2.0e6;
    double m2   = 1.0e6;
    double chi1 = 0.85;
    double chi2 = 0.42;
    double dL   = 1.0; // Gpc
    double tc   = 3.0e7; //seconds
    double phi0 = 0.0;
    double Mtot = m1 + m2;        // total mass in solar masses
    double Mc   = chirpmass(m1,m2);  // chirp mass in solar masses


    params[0]  = log(Mc);   // chirp mass
    params[1]  = log(Mtot); // total mass
    params[2]  = chi1;      // spin of BH1
    params[3]  = chi2;      // spin of BH2
    params[4]  = phi0;      // reference phase
    params[5]  = tc;        // coalesence time
    params[6]  = log(dL);   // luminosity distance
    params[7]  = 2.31;      // ecliptic latitude (theta)
    params[8]  = 0.57;      // ecliptic longitude (phi)
    params[9]  = 0.4;       // polarization
    params[10] = 0.3;  // cos inclination

    /*
     TDI structure to hold simulated waveform(s)
     */
    struct TDI *tdi = malloc(sizeof(struct TDI));
    alloc_tdi(tdi,N,3);
    
    /*
     Keep track of which wavelets are used for template
     */
    int *wavelet_list = int_vector(N);
    int Nwavelet;
    
    /*
     Test frequency domain waveform
     */
    start = clock();
    for(int i=0; i<NTEST; i++) mbh_fd_waveform(orbit, wdm, Tobs, t0, params, wavelet_list, &Nwavelet, tdi->X, tdi->Y, tdi->Z);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC / NTEST;
    printf("PhenomD calculation took %f seconds\n", cpu_time_used);

    /*
     Test time domain waveform
     
    start = clock();
    for(int i=0; i<NTEST; i++) mbh_td_waveform(orbit, wdm, Tobs, t0, params, wavelet_list, &Nwavelet, tdi->X, tdi->Y, tdi->Z);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC / NTEST;
    printf("PhenomT calculation took %f seconds\n", cpu_time_used);*/

    return 0;
}
