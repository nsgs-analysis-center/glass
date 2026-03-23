/*
 * Copyright 2023 Tyson B. Littenberg, Neil J. Cornish, Robert Rosati
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


#include <glass_utils.h>
#include <glass_math.h>

#include "glass_noise.h"


void map_array_to_noise_params(struct InstrumentModel *model)
{
    model->sacc12 = model->sacc[0];
    model->sacc21 = model->sacc[1];

    model->sacc23 = model->sacc[2];
    model->sacc32 = model->sacc[3];

    model->sacc13 = model->sacc[4];
    model->sacc31 = model->sacc[5];

    model->soms12 = model->soms[0];
    model->soms21 = model->soms[1];
        
    model->soms23 = model->soms[2];
    model->soms32 = model->soms[3];
    
    model->soms13 = model->soms[4];
    model->soms31 = model->soms[5];

}

void map_noise_params_to_array(struct InstrumentModel *model)
{
    model->sacc[0] = model->sacc12;
    model->sacc[1] = model->sacc21;
        
    model->sacc[2] = model->sacc23;
    model->sacc[3] = model->sacc32;
    
    model->sacc[4] = model->sacc13;
    model->sacc[5] = model->sacc31;

    model->soms[0] = model->soms12;
    model->soms[1] = model->soms21;
    
    model->soms[2] = model->soms23;
    model->soms[3] = model->soms32;

    model->soms[4] = model->soms13;
    model->soms[5] = model->soms31;
}

void map_array_to_foreground_params(struct ForegroundModel *model)
{
    model->Amp   = exp(model->sgal[0]);
    model->f1    = exp(model->sgal[1]);
    model->alpha = model->sgal[2];
    model->fk    = exp(model->sgal[3]);
    model->f2    = exp(model->sgal[4]);
}

void map_foreground_params_to_array(struct ForegroundModel *model)
{
    model->sgal[0] = log(model->Amp);
    model->sgal[1] = log(model->f1);
    model->sgal[2] = model->alpha;
    model->sgal[3] = log(model->fk);
    model->sgal[4] = log(model->f2);
}

void alloc_spline_model(struct SplineModel *model, int Ndata, int Nlayer, int Nchannel, int Nspline)
{
    model->Nchannel = Nchannel;
    
    model->psd = malloc(sizeof(struct Noise));
    model->spline=malloc(sizeof(struct Noise));
    
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    alloc_noise(model->spline, Nspline, Nlayer, Nchannel);
}

void alloc_instrument_model(struct InstrumentModel *model, int Ndata, int Nlayer, int Nchannel)
{
    model->Nlink=6;
    model->soms = malloc(model->Nlink*sizeof(double));
    model->sacc = malloc(model->Nlink*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
}

void alloc_foreground_model(struct ForegroundModel *model, int Ndata, int Nlayer, int Nchannel)
{
    model->Nparams=5;
    model->sgal = malloc(model->Nparams*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
}

struct SGWBResponse* global_SGWBResponse = NULL;

void alloc_sgwb_model(struct SGWBModel *model, int Ndata, int Nlayer, int Nchannel, SGWB_t SGWB_type)
{
    model->SGWB_type = SGWB_type;
    model->Nparams = SGWB_TEMPLATE_NPARAMS[SGWB_type];
    model->params = malloc(model->Nparams*sizeof(double));
    model->psd = malloc(sizeof(struct Noise));
    alloc_noise(model->psd, Ndata, Nlayer, Nchannel);
    if (global_SGWBResponse == NULL) {
        model->R = malloc(sizeof(struct SGWBResponse)); // memory leak
        alloc_pop_sgwb_response(model->R,"./sgwb_response_xyz1.dat");
        global_SGWBResponse = model->R;
    } else {
        model->R = global_SGWBResponse;
    }
}

void alloc_pop_sgwb_response(struct SGWBResponse* sgwbr, char* fname) {
    FILE* ff = fopen(fname,"r");
    if (ff == NULL) {
        printf("Could not open file %s ! SGWBResponse generation failed.",fname);
        exit(-2);
    }
    int *N = &sgwbr->N;
    if (fscanf(ff,"%d\n",N) != 1) {
        fclose(ff);
        printf("Failed to read first line of %s ! SGWBResponse generation failed.",fname);
        exit(-2);
    }
    if (*N > 3000 || *N < 10) {
        fclose(ff);
        printf("That seems wrong? I read N=%d . SGWBResponse generation failed.",*N);
        exit(-2);
    }

    // memory leak
    sgwbr->logf = malloc(sizeof(double) * *N); 
    sgwbr->logXX = malloc(sizeof(double) * *N);
    sgwbr->logYY = malloc(sizeof(double) * *N);
    sgwbr->logZZ = malloc(sizeof(double) * *N);
    sgwbr->asinhXY = malloc(sizeof(double) * *N);
    sgwbr->asinhXZ = malloc(sizeof(double) * *N);
    sgwbr->asinhYZ = malloc(sizeof(double) * *N);
    sgwbr->asinh_scale = 1e-7;
    for (int i=0;i<*N;i++) {
        double f,XX,YY,ZZ,XY,XZ,YZ;
        fscanf(ff,"%lf %lf %lf %lf %lf %lf %lf\n",&f,&XX,&YY,&ZZ,&XY,&XZ,&YZ);
        sgwbr->logf[i] = log(f);
        sgwbr->logXX[i] = log(XX);
        sgwbr->logYY[i] = log(YY);
        sgwbr->logZZ[i] = log(ZZ);
        sgwbr->asinhXY[i] = asinh(XY / sgwbr->asinh_scale);
        sgwbr->asinhXZ[i] = asinh(XZ / sgwbr->asinh_scale);
        sgwbr->asinhYZ[i] = asinh(YZ / sgwbr->asinh_scale);
    }
    fclose(ff);


    // memory leak
    sgwbr->spline_logRXX   = alloc_cubic_spline_even_sampling(*N);
    sgwbr->spline_logRYY   = alloc_cubic_spline_even_sampling(*N);
    sgwbr->spline_logRZZ   = alloc_cubic_spline_even_sampling(*N);
    sgwbr->spline_asinhRXY = alloc_cubic_spline_even_sampling(*N);
    sgwbr->spline_asinhRXZ = alloc_cubic_spline_even_sampling(*N);
    sgwbr->spline_asinhRYZ = alloc_cubic_spline_even_sampling(*N);
    double dlogf = sgwbr->logf[1] - sgwbr->logf[0];
    initialize_cubic_spline_even_sampling(sgwbr->spline_logRXX, sgwbr->logf, sgwbr->logXX, dlogf);
    initialize_cubic_spline_even_sampling(sgwbr->spline_logRYY, sgwbr->logf, sgwbr->logYY, dlogf);
    initialize_cubic_spline_even_sampling(sgwbr->spline_logRZZ, sgwbr->logf, sgwbr->logZZ, dlogf);
    initialize_cubic_spline_even_sampling(sgwbr->spline_asinhRXY, sgwbr->logf, sgwbr->asinhXY, dlogf);
    initialize_cubic_spline_even_sampling(sgwbr->spline_asinhRXZ, sgwbr->logf, sgwbr->asinhXZ, dlogf);
    initialize_cubic_spline_even_sampling(sgwbr->spline_asinhRYZ, sgwbr->logf, sgwbr->asinhYZ, dlogf);

}

void free_sgwb_response(struct SGWBResponse *sgwbr) {
    free(sgwbr->logf);
    free(sgwbr->logXX);
    free(sgwbr->logYY);
    free(sgwbr->logZZ);
    free(sgwbr->asinhXY);
    free(sgwbr->asinhXZ);
    free(sgwbr->asinhYZ);
    free_cubic_spline_even_sampling(sgwbr->spline_logRXX);
    free_cubic_spline_even_sampling(sgwbr->spline_logRYY);
    free_cubic_spline_even_sampling(sgwbr->spline_logRZZ);
    free_cubic_spline_even_sampling(sgwbr->spline_asinhRXY);
    free_cubic_spline_even_sampling(sgwbr->spline_asinhRXZ);
    free_cubic_spline_even_sampling(sgwbr->spline_asinhRYZ);
    free(sgwbr);
}

void free_spline_model(struct SplineModel *model)
{
    free_noise(model->psd);
    free_noise(model->spline);
    free(model);
}

void free_instrument_model(struct InstrumentModel *model)
{
    free(model->soms);
    free(model->sacc);
    free_noise(model->psd);
    free(model);
}

void free_foreground_model(struct ForegroundModel *model)
{
    free_noise(model->psd);
    free(model->sgal);
    free(model);
}

void free_sgwb_model(struct SGWBModel *model)
{
    free_noise(model->psd);
    // TODO should actually only do this if we're the last SGWBModel instance. 
    // for now, leak this memory
    // free_sgwb_response(model->R);
    free(model->params);
    free(model);
}

void copy_spline_model(struct SplineModel *origin, struct SplineModel *copy)
{
    copy->Nchannel = origin->Nchannel;

    //Spline parameters
    copy_noise(origin->spline,copy->spline);
    
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;
    
    //Priors
    copy->Nmin = origin->Nmin;
    copy->Nmax = origin->Nmax;
}

void copy_instrument_model(struct InstrumentModel *origin, struct InstrumentModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);

    //Model likelihood
    copy->logL = origin->logL;

    //Instrument model Parameters
    copy->Nlink = origin->Nlink;
    memcpy(copy->soms, origin->soms, origin->Nlink*sizeof(double));
    memcpy(copy->sacc, origin->sacc, origin->Nlink*sizeof(double));
}

void copy_foreground_model(struct ForegroundModel *origin, struct ForegroundModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;

    //Foreground model Parameters
    copy->Tobs = origin->Tobs;
    copy->Nparams = origin->Nparams;
    memcpy(copy->sgal, origin->sgal, origin->Nparams*sizeof(double));
    
}
void copy_sgwb_model(struct SGWBModel *origin, struct SGWBModel *copy)
{
    //Noise model parameters
    copy_noise(origin->psd,copy->psd);
    
    //Model likelihood
    copy->logL = origin->logL;

    //SGWB model Parameters
    copy->Tobs = origin->Tobs;
    copy->Nparams = origin->Nparams;
    copy->SGWB_type = origin->SGWB_type;
    memcpy(copy->params, origin->params, origin->Nparams*sizeof(double));
    
}

void update_spline_noise_model(struct SplineModel *model, int new_knot, int min_knot, int max_knot)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    /* find location in data vector of knot */
    double T = 1./(psd->f[1] - psd->f[0]);

    struct CubicSpline **cspline = malloc(model->Nchannel*sizeof(struct CubicSpline *));
    
    /* have to recompute the spline everywhere (derivatives on boundary) */
    for(int n=0; n<model->Nchannel; n++)
    {
        cspline[n] = alloc_cubic_spline(spline->N);
        initialize_cubic_spline(cspline[n], spline->f, spline->C[n][n], SPLINE_BINARY_SEARCH);
    }
    
    int imin = (int)((spline->f[min_knot]-psd->f[0])*T);
    int imax = (int)((spline->f[max_knot]-psd->f[0])*T);
    
    
    for(int i=imin; i<imax; i++)
    {
        for(int n=0; n<model->Nchannel; n++)
        {
            psd->C[n][n][i]=spline_interpolation(cspline[n],psd->f[i]);
            
            /*
             apply transfer function
             -this catches the sharp features in the spectrum from f/fstar
             -without needing to interpolate
             */
            psd->C[n][n][i]+=psd->transfer[i];
        }
    }
    invert_noise_covariance_matrix(psd);

    for(int n=0; n<model->Nchannel; n++)
    {
        free_cubic_spline(cspline[n]);
    }
    free(cspline);
}


void generate_spline_noise_model(struct SplineModel *model)
{
    struct Noise *psd = model->psd;
    struct Noise *spline = model->spline;
    
    for(int i=0; i<model->Nchannel; i++)
        CubicSplineGLASS(spline->N, spline->f, spline->C[i][i], psd->N, psd->f, psd->C[i][i]);

    //set a floor on Sn so likelihood doesn't go crazy
    for(int n=0; n<psd->N; n++)
    {
        /*
         apply transfer function
         -this catches the sharp features in the spectrum from f/fstar
         -without needing to interpolate
         */
        for(int i=0; i<model->Nchannel; i++)
            psd->C[i][i][n]+=psd->transfer[n];
        
    }
    invert_noise_covariance_matrix(psd);
}

void generate_instrument_noise_model(struct Orbit *orbit, struct InstrumentModel *model)
{
    double f,f2;
    double x;
    double cosx;
    double oms_transfer_function;
    double acc_transfer_function;
    double tdi_transfer_function;
    double Sacc;
    double Soms;
    double Sacc12,Sacc21,Sacc13,Sacc31,Sacc23,Sacc32;
    double Soms12,Soms21,Soms13,Soms31,Soms23,Soms32;

    map_array_to_noise_params(model);

    double acc_units = 2.81837551648e-19;//1./(PI2*CLIGHT)/(PI2*CLIGHT)
    double oms_units = 4.39256635604e-16;//PI2*PI2/CLIGHT/CLIGHT
    double acc_fmax  = 0.01; //Hz (maximum frequency for computing acc contribution
    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        x = f/orbit->fstar;
        f2= f*f;
        
        //at low frequency use linear approximation for trig functions
        if(x<0.1)
        {
            cosx = 1.0;
            tdi_transfer_function = x*x;
        }
        else
        {
            cosx = cos(x);
            tdi_transfer_function = noise_transfer_function(x);
        }

        cosx = cos(x);
        tdi_transfer_function = noise_transfer_function(x);

        if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
            acc_transfer_function = acc_units / f2 * (1.0 + pow(0.4e-3/f,2)) * (1.0 + pow(f/8.0e-3,4));
        else
            acc_transfer_function = 0.0;

        oms_transfer_function = oms_units * f2 * (1.0 + pow(2.0e-3/f,4));

        //absorb factor of -8 into tdi transfer function to cut down on multiplications
        tdi_transfer_function *= -8.0;
        
        switch(model->psd->Nchannel)
        {
            case 1:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                break;
            case 2:
                Sacc = model->sacc12 * acc_transfer_function;
                Soms = model->soms12 * oms_transfer_function;
                model->psd->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Sacc, Soms);
                model->psd->C[0][1][n] = 0.0;
                break;
            case 3:
                
                if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
                {
                    Sacc12 = model->sacc12 * acc_transfer_function;
                    Sacc21 = model->sacc21 * acc_transfer_function;
                    Sacc23 = model->sacc23 * acc_transfer_function;
                    Sacc32 = model->sacc32 * acc_transfer_function;
                    Sacc13 = model->sacc13 * acc_transfer_function;
                    Sacc31 = model->sacc31 * acc_transfer_function;
                }
                
                Soms12 = model->soms12 * oms_transfer_function;
                Soms21 = model->soms21 * oms_transfer_function;
                Soms23 = model->soms23 * oms_transfer_function;
                Soms32 = model->soms32 * oms_transfer_function;
                Soms13 = model->soms13 * oms_transfer_function;
                Soms31 = model->soms31 * oms_transfer_function;
            
                //Start from scratch
                model->psd->C[0][0][n] = 0.;
                model->psd->C[1][1][n] = 0.;
                model->psd->C[2][2][n] = 0.;

                model->psd->C[0][1][n] = 0.;
                model->psd->C[0][2][n] = 0.;
                model->psd->C[1][2][n] = 0.;

                //OMS noise
                model->psd->C[0][0][n] += (Soms12+Soms21+Soms13+Soms31)*.25;
                model->psd->C[1][1][n] += (Soms23+Soms32+Soms21+Soms12)*.25;
                model->psd->C[2][2][n] += (Soms31+Soms13+Soms32+Soms23)*.25;

                model->psd->C[0][1][n] += (Soms12 + Soms21)*.5;
                model->psd->C[0][2][n] += (Soms13 + Soms31)*.5;
                model->psd->C[1][2][n] += (Soms32 + Soms23)*.5;
                
                //Acceleration noise
                if(f<acc_fmax) //at f<~.1mHz acceleration noise is < 1% of the total noise budget
                {
                    //distribute through factors of 2
                    /*
                    model->psd->C[0][0][n] += 2. * ( (Sacc12 + Sacc13)*.5 + ((Sacc21 + Sacc31)*.5)*cosx*cosx );
                    model->psd->C[1][1][n] += 2. * ( (Sacc23 + Sacc21)*.5 + ((Sacc32 + Sacc12)*.5)*cosx*cosx );
                    model->psd->C[2][2][n] += 2. * ( (Sacc31 + Sacc32)*.5 + ((Sacc13 + Sacc23)*.5)*cosx*cosx );

                    model->psd->C[0][1][n] += 4. * ( (Sacc12 + Sacc21)*.5 );
                    model->psd->C[0][2][n] += 4. * ( (Sacc13 + Sacc31)*.5 );
                    model->psd->C[1][2][n] += 4. * ( (Sacc32 + Sacc23)*.5 );
                    */
                    model->psd->C[0][0][n] += ( (Sacc12 + Sacc13) + ((Sacc21 + Sacc31))*cosx*cosx );
                    model->psd->C[1][1][n] += ( (Sacc23 + Sacc21) + ((Sacc32 + Sacc12))*cosx*cosx );
                    model->psd->C[2][2][n] += ( (Sacc31 + Sacc32) + ((Sacc13 + Sacc23))*cosx*cosx );

                    model->psd->C[0][1][n] += 2.*(Sacc12 + Sacc21);
                    model->psd->C[0][2][n] += 2.*(Sacc13 + Sacc31);
                    model->psd->C[1][2][n] += 2.*(Sacc32 + Sacc23);
                }
                

                //TDI transfer functions (note tdi_transfer_function has abosrbed a factor of -8)
                model->psd->C[0][0][n] *= -2. * tdi_transfer_function;
                model->psd->C[1][1][n] *= -2. * tdi_transfer_function;
                model->psd->C[2][2][n] *= -2. * tdi_transfer_function;

                model->psd->C[0][1][n] *= tdi_transfer_function * cosx ;
                model->psd->C[0][2][n] *= tdi_transfer_function * cosx ;
                model->psd->C[1][2][n] *= tdi_transfer_function * cosx ;
                
                //Symmetry
                model->psd->C[1][0][n] = model->psd->C[0][1][n];
                model->psd->C[2][0][n] = model->psd->C[0][2][n];
                model->psd->C[2][1][n] = model->psd->C[1][2][n];
                
                break;
        }
        
        //Normalization
        for(int i=0; i<model->psd->Nchannel; i++)
            for(int j=0; j<model->psd->Nchannel; j++)
                model->psd->C[i][j][n] /= 4.0;
    }
}
//TODO: is it worth oversampling these very smooth functional forms...?
/*
void generate_instrument_noise_model_wavelet_coarse(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model) {
    // active layers
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;
    generate_instrument_noise_model(orbit,model);
    //NOTE: normalization fudge factor
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                C[n][m][i]/=8.;
}
*/

void generate_instrument_noise_model_wavelet(struct Wavelets *wdm, struct Orbit *orbit, struct InstrumentModel *model)
{

    /* 
    oversampled frequency grid 
    */
    struct InstrumentModel *grid = malloc(sizeof(struct InstrumentModel));

    // active layers
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;

    // initialize data models
    // NOTE: originally the grid went over every frequency layer, even ones we don't use
    alloc_instrument_model(grid, 2*(imax-imin) + 2, imax-imin, 3);
    /* Basic idea here is to resample at df/2, with one extra sample on each end
     *    model:         0 ------- 1 ------- 2 ------- 3 ------- ... --------- N
     *    grid:   -0.5 - 0 - 0.5 - 1 - 1.5 - 2 - 2.5 - 3 - 3.5 - ... - N-0.5 - N - N+0.5 
     */

    // set up psd frequency grid
    double grid_fmin = imin*wdm->df - wdm->df/2.0;
    for(size_t n=0; n<grid->psd->N; n++)
        grid->psd->f[n] = grid_fmin + wdm->df/2.0*n;

    // initialize noise levels
    for(int i=0; i<grid->Nlink; i++)
    {
        grid->soms[i] = model->soms[i];
        grid->sacc[i] = model->sacc[i];
    }

    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model(orbit,grid);
    
    /*
    integrate instrument noise over each frequency layer
    */
    double ***C     = model->psd->C;
    double ***Cgrid = grid->psd->C;
    for(int i=imin; i<imax; i++)
    {
        int j = 2*(i-imin) + 1;
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++) {
                // TODO: fix formula? spacing is df!
                C[n][m][i-imin] = simpson_integration_3(Cgrid[n][m][j-1],Cgrid[n][m][j],Cgrid[n][m][j+1],1.0);
            }
    }

    //NOTE: normalization fudge factor (and symmetrization)
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++) {
                // TODO: normalization check here
                C[n][m][i]/=8.;
                C[m][n][i] = C[n][m][i];
            }

    free_instrument_model(grid);
}

void generate_galactic_foreground_model(struct ForegroundModel *model)
{
    double f;
    double Sgal;
    
    map_array_to_foreground_params(model);
    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        double fstar=CLIGHT/(PI2*LARM);
        double x = f/fstar;
        double t = 4.*x*x*sin(x)*sin(x); //transfer function

        //skip confusion noise calculation at high frequencies where contribution is negligible
        int high_f_check = 0;
        if(f>model->f1*10) high_f_check = 1;
        
        if(high_f_check) Sgal = 0.0;
        else Sgal = galaxy_foreground(f,model->Amp, model->f1, model->alpha, model->fk, model->f2);
        
        switch(model->psd->Nchannel)
        {
            case 1:
                model->psd->C[0][0][n] = Sgal;
                break;
            case 2:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = 1.5*t*Sgal;
                model->psd->C[0][1][n] = model->psd->C[1][0][n] = 0;
                break;
            case 3:
                model->psd->C[0][0][n] = model->psd->C[1][1][n] = model->psd->C[2][2][n] = t*Sgal;
                model->psd->C[0][1][n] = model->psd->C[0][2][n] = model->psd->C[1][2][n] = -0.5*t*Sgal;
                model->psd->C[1][0][n] = model->psd->C[2][0][n] = model->psd->C[2][1][n] = -0.5*t*Sgal;
                break;
        }
        
        //Normalization
        for(int i=0; i<model->psd->Nchannel; i++)
            for(int j=0; j<model->psd->Nchannel; j++)
                model->psd->C[i][j][n] /= 4.0;
    }
}

void generate_galactic_foreground_model_wavelet(struct Wavelets *wdm, struct ForegroundModel *model)
{
    /* 
    oversampled frequency grid 
    */
    struct ForegroundModel *grid = malloc(sizeof(struct ForegroundModel));

    // active layers
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;

    // initialize data models
    alloc_foreground_model(grid, 2*(imax-imin) + 2, imax-imin, 3);
    /* Basic idea here is to resample at df/2, with one extra sample on each end
     *    model:         0 ------- 1 ------- 2 ------- 3 ------- ... --------- N
     *    grid:   -0.5 - 0 - 0.5 - 1 - 1.5 - 2 - 2.5 - 3 - 3.5 - ... - N-0.5 - N - N+0.5 
     */

    // set up psd frequency grid
    double grid_fmin = imin*wdm->df - wdm->df/2.0;
    for(size_t n=0; n<grid->psd->N; n++)
        grid->psd->f[n] = grid_fmin + wdm->df/2.0*n;

    // initialize foreground parameters levels
    map_array_to_foreground_params(model);
    grid->Tobs  = model->Tobs;
    grid->Amp   = model->Amp;
    grid->f1    = model->f1;
    grid->alpha = model->alpha;
    grid->fk    = model->fk;
    grid->f2    = model->f2;
    map_foreground_params_to_array(grid);

    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(grid);

    /*
    integrate foreground over each frequency layer
    */
    double ***C     = model->psd->C;
    double ***Cgrid = grid->psd->C;

    for(int i=imin; i<imax; i++)
    {
        int j = 2*(i-imin) + 1;
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                C[n][m][i-imin] = simpson_integration_3(Cgrid[n][m][j-1],Cgrid[n][m][j],Cgrid[n][m][j+1],1.0);
    }

    //NOTE: undo isotropc -1/2 on covariance hardcoded in generate_galactic_foreground_model()
    // modulation(t) will be applied to PSDs later
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++)
                if(n!=m) C[n][m][i]*=-2.;


    //NOTE: normalization fudge factor
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++) {
                C[n][m][i]/=8.;
                C[m][n][i] = C[n][m][i];
            }

    free_foreground_model(grid);
}


_Static_assert(SGWB_TEMPLATE_COUNT == 1, "Did you add an SGWB template? Implement its form in frequency in a new function here");
inline double sgwb_powerlaw(double f, const double* params) {
    static double fref = 25.0;
    double A     = pow(10.0,params[0]);
    double alpha = params[1];
    double prefactor = Hscale / (f*f*f);
    return prefactor*A*pow(f/fref,alpha);
}

void generate_sgwb_model(struct SGWBModel *model)
{
    double f;
    double Sgw;

    
    for(int n=0; n<model->psd->N; n++)
    {
        f = model->psd->f[n];
        
        _Static_assert(SGWB_TEMPLATE_COUNT == 1, "Did you add an SGWB template? Edit this switch case, it needs to be exhaustive.");
        switch(model->SGWB_type) {
            case SGWB_TEMPLATE_POWERLAW:
                Sgw = sgwb_powerlaw(f,model->params);
                break;
            default:
                fprintf(stderr,"Unimplemented SGWB type?\n\tTemplate index: %d\n\tTemplate name: %s\n",model->SGWB_type,SGWB_TEMPLATE_NAMES[model->SGWB_type]);
                exit(1);
                break;
        }
        
        switch(model->psd->Nchannel)
        {
            double Rxx, Ryy, Rzz, Rxy, Rxz, Ryz;
            case 1:
            case 2:
                fprintf(stderr, "Unimplemented error! SGWB covariance generation with Nchannels = %d\n",model->psd->Nchannel);
                exit(-3);
                break;
            case 3:
                double logf = log(f);
                Rxx =  exp(spline_interpolation_even_sampling(model->R->spline_logRXX, logf));
                Ryy =  exp(spline_interpolation_even_sampling(model->R->spline_logRYY, logf));
                Rzz =  exp(spline_interpolation_even_sampling(model->R->spline_logRZZ, logf));
                Rxy = model->R->asinh_scale * sinh(spline_interpolation_even_sampling(model->R->spline_asinhRXY, logf));
                Rxz = model->R->asinh_scale * sinh(spline_interpolation_even_sampling(model->R->spline_asinhRXZ, logf));
                Ryz = model->R->asinh_scale * sinh(spline_interpolation_even_sampling(model->R->spline_asinhRYZ, logf));
                // note that, in equal-arm LISA, XYZ: Rxx = Ryy = Rzz, and Rxy = Rxz = Ryz = -0.5*Rxx
                model->psd->C[0][0][n] = Rxx*Sgw;
                model->psd->C[1][1][n] = Ryy*Sgw;
                model->psd->C[2][2][n] = Rzz*Sgw;
                model->psd->C[0][1][n] = model->psd->C[1][0][n] = Rxy*Sgw;
                model->psd->C[0][2][n] = model->psd->C[2][0][n] = Rxz*Sgw;
                model->psd->C[1][2][n] = model->psd->C[2][1][n] = Ryz*Sgw;
                break;
        }
    }
}
void generate_sgwb_model_wavelet(struct Wavelets* wdm, struct SGWBModel *model)
{
    /* 
    oversampled frequency grid 
    */
    struct SGWBModel *grid = malloc(sizeof(struct SGWBModel));
    // active layers
    int imin = (int)round(model->psd->f[0]/wdm->df);
    int imax = (int)round(model->psd->f[model->psd->N-1]/wdm->df)+1;

    // initialize data models
    alloc_sgwb_model(grid, 2*(imax-imin) + 2, imax-imin, 3, model->SGWB_type);
    /* Basic idea here is to resample at df/2, with one extra sample on each end
     *    model:         0 ------- 1 ------- 2 ------- 3 ------- ... --------- N
     *    grid:   -0.5 - 0 - 0.5 - 1 - 1.5 - 2 - 2.5 - 3 - 3.5 - ... - N-0.5 - N - N+0.5 
     */

    // set up psd frequency grid
    double grid_fmin = imin*wdm->df - wdm->df/2.0;
    for(size_t n=0; n<grid->psd->N; n++)
        grid->psd->f[n] = grid_fmin + wdm->df/2.0*n;

    // intialize grid parameter levels
    grid->Nparams = model->Nparams;
    grid->Tobs = model->Tobs;
    for (int i=0; i<model->Nparams; i++)
        grid->params[i] = model->params[i]; 

    // get noise covariance matrix for initial parameters
    generate_sgwb_model(grid);

    /*
    integrate SGWB over each frequency layer
    */
    double ***C     = model->psd->C;
    double ***Cgrid = grid->psd->C;
    for(int i=imin; i<imax; i++)
    {
        int j = 2*(i-imin) + 1;
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++) {
                C[n][m][i-imin] = simpson_integration_3(Cgrid[n][m][j-1],Cgrid[n][m][j],Cgrid[n][m][j+1],1.0);
            }
    }

    // TODO: don't think we need this???
    //NOTE: normalization fudge factor
    for(int i=0; i<model->psd->N; i++)
        for(int n=0; n<3; n++)
            for(int m=n; m<3; m++) {
                C[n][m][i]/=8.;
                C[m][n][i] = C[n][m][i];
            }
    free_sgwb_model(grid);
}

void generate_full_covariance_matrix(struct Noise *full, struct Noise *component, int Nchannel)
{
    for(int n=0; n<full->N; n++)
    {
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=0; j<Nchannel; j++)
            {
                full->C[i][j][n] += component->C[i][j][n];
            }
        }
    }
}


void generate_full_dynamic_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full)
{
    int k;
    int jmin=(int)round(inst->psd->f[0]/wdm->df);
    int jmax=(int)round(inst->psd->f[inst->psd->N-1]/wdm->df)+1; 

    for(int i=0; i<wdm->NT; i++)
    {
        double t = i*wdm->dt;
        for(int j=jmin; j<jmax; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);

            k-=wdm->kmin;

            //stationary instrument noise 
            full->C[0][0][k] = inst->psd->C[0][0][j-jmin];
            full->C[1][1][k] = inst->psd->C[1][1][j-jmin];
            full->C[2][2][k] = inst->psd->C[2][2][j-jmin];
            full->C[0][1][k] = inst->psd->C[0][1][j-jmin];
            full->C[0][2][k] = inst->psd->C[0][2][j-jmin];
            full->C[1][2][k] = inst->psd->C[1][2][j-jmin];

            //modulated galactic foreground
            if (conf) { // this is NULL when confusion is disabled
                full->C[0][0][k] += conf->psd->C[0][0][j-jmin]*spline_interpolation_even_sampling(conf->modulation->XX_spline, t);
                full->C[1][1][k] += conf->psd->C[1][1][j-jmin]*spline_interpolation_even_sampling(conf->modulation->YY_spline, t);
                full->C[2][2][k] += conf->psd->C[2][2][j-jmin]*spline_interpolation_even_sampling(conf->modulation->ZZ_spline, t);
                full->C[0][1][k] += conf->psd->C[0][1][j-jmin]*spline_interpolation_even_sampling(conf->modulation->XY_spline, t);
                full->C[0][2][k] += conf->psd->C[0][2][j-jmin]*spline_interpolation_even_sampling(conf->modulation->XZ_spline, t);
                full->C[1][2][k] += conf->psd->C[1][2][j-jmin]*spline_interpolation_even_sampling(conf->modulation->YZ_spline, t);
            }


            //stationary stochastic background
            if (sgwb) { // this is NULL when sgwb is disabled
                full->C[0][0][k] += sgwb->psd->C[0][0][j-jmin];
                full->C[1][1][k] += sgwb->psd->C[1][1][j-jmin];
                full->C[2][2][k] += sgwb->psd->C[2][2][j-jmin];
                full->C[0][1][k] += sgwb->psd->C[0][1][j-jmin];
                full->C[0][2][k] += sgwb->psd->C[0][2][j-jmin];
                full->C[1][2][k] += sgwb->psd->C[1][2][j-jmin];
            }

            //noise covariance matrix is symmetric
            full->C[1][0][k] = full->C[0][1][k]; 
            full->C[2][0][k] = full->C[0][2][k]; 
            full->C[2][1][k] = full->C[1][2][k]; 
        } //loop over frequency layers
    } //loop over time slices
}

static void generate_full_stationary_covariance_matrix(struct Wavelets *wdm, struct InstrumentModel *inst, struct ForegroundModel *conf, struct SGWBModel *sgwb, struct Noise *full)
{
    int k;
    int jmin=(int)round(inst->psd->f[0]/wdm->df);
    int jmax=(int)round(inst->psd->f[inst->psd->N-1]/wdm->df)+1;

    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=jmin; j<jmax; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);

            k-=wdm->kmin;

            //stationary instrument noise
            full->C[0][0][k] = inst->psd->C[0][0][j-jmin];
            full->C[1][1][k] = inst->psd->C[1][1][j-jmin];
            full->C[2][2][k] = inst->psd->C[2][2][j-jmin];
            full->C[0][1][k] = inst->psd->C[0][1][j-jmin];
            full->C[0][2][k] = inst->psd->C[0][2][j-jmin];
            full->C[1][2][k] = inst->psd->C[1][2][j-jmin];

            //modulated galactic foreground
            full->C[0][0][k] += conf->psd->C[0][0][j-jmin];
            full->C[1][1][k] += conf->psd->C[1][1][j-jmin];
            full->C[2][2][k] += conf->psd->C[2][2][j-jmin];
            full->C[0][1][k] -= conf->psd->C[0][1][j-jmin]/2.;
            full->C[0][2][k] -= conf->psd->C[0][2][j-jmin]/2.;
            full->C[1][2][k] -= conf->psd->C[1][2][j-jmin]/2.; // 2 here is baked in response

            //stationary stochastic background
            full->C[0][0][k] += sgwb->psd->C[0][0][j-jmin];
            full->C[1][1][k] += sgwb->psd->C[1][1][j-jmin];
            full->C[2][2][k] += sgwb->psd->C[2][2][j-jmin];
            full->C[0][1][k] += sgwb->psd->C[0][1][j-jmin];
            full->C[0][2][k] += sgwb->psd->C[0][2][j-jmin];
            full->C[1][2][k] += sgwb->psd->C[1][2][j-jmin];

            //noise covariance matrix is symmetric
            full->C[1][0][k] = full->C[0][1][k];
            full->C[2][0][k] = full->C[0][2][k];
            full->C[2][1][k] = full->C[1][2][k];
        }// loop over frequency layers
    }// loop over time slices
}

double noise_log_likelihood(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    
    struct TDI *tdi = data->tdi;
    
    int N = data->NFFT;
    
    switch(data->Nchannel)
    {
        case 1:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            break;
        case 2:
            logL += -0.5*fourier_nwip(tdi->A, tdi->A, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->E, tdi->E, noise->invC[1][1], N);
            break;
        case 3:
            logL += -0.5*fourier_nwip(tdi->X, tdi->X, noise->invC[0][0], N);
            logL += -0.5*fourier_nwip(tdi->Y, tdi->Y, noise->invC[1][1], N);
            logL += -0.5*fourier_nwip(tdi->Z, tdi->Z, noise->invC[2][2], N);
            logL += -fourier_nwip(tdi->X, tdi->Y, noise->invC[0][1], N);
            logL += -fourier_nwip(tdi->X, tdi->Z, noise->invC[0][2], N);
            logL += -fourier_nwip(tdi->Y, tdi->Z, noise->invC[1][2], N);
            break;
    }
    for(int n=0; n<N; n++)
        logL -= noise->logdetC[n];
    
    return logL;
}

static inline double wavelet_nwip_linear(const double* __restrict a, const double* __restrict b, const double* __restrict invC, size_t N)
{
    // TODO: use lapack/blas here? this is just a.b.invC
    double arg = 0.0;
    #pragma omp simd reduction(+:arg)
    for(size_t k=0; k<N; k++)
    {
        arg += a[k]*b[k]*invC[k];
    }
    return arg;
}

double my_noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise) {
    struct TDI *tdi = data->tdi;
    const double log2pi = log(2*M_PI);
    double logL = 0.0;
    #pragma omp simd reduction(+:logL)
    for (int i=0; i<data->wdm->NT; i++)
        for (int j=data->lmin; j<data->lmax; j++) {
            int k;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k -= data->wdm->kmin;
            logL -= 0.5*tdi->X[k]*tdi->X[k]*noise->invC[0][0][k];
            logL -= 0.5*tdi->Y[k]*tdi->Y[k]*noise->invC[1][1][k];
            logL -= 0.5*tdi->Z[k]*tdi->Z[k]*noise->invC[2][2][k];
            logL -=     tdi->X[k]*tdi->Y[k]*noise->invC[0][1][k];
            logL -=     tdi->X[k]*tdi->Z[k]*noise->invC[0][2][k];
            logL -=     tdi->Y[k]*tdi->Z[k]*noise->invC[1][2][k];
            logL -= noise->logdetC[k];
    }
    logL -= 3 * data->N * log2pi;

    return logL;
}

double my_noise_log_likelihood(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    const double log2pi = log(2*M_PI);
    
    struct TDI *tdi = data->tdi;
    
    int N = data->NFFT;
    #pragma omp simd reduction(+:logL)
    for (int n=0; n<N; n++) {
        double rex = tdi->X[2*n], imx = tdi->X[2*n+1];
        double rey = tdi->Y[2*n], imy = tdi->Y[2*n+1];
        double rez = tdi->Z[2*n], imz = tdi->Z[2*n+1];
        double xx = rex*rex + imx*imx;
        double yy = rey*rey + imy*imy;
        double zz = rez*rez + imz*imz;
        double xy = rex*rey + imx*imy; // note this is really (x*conj(y) + y*conj(x))/2
        double xz = rex*rez + imx*imz;
        double yz = rey*rez + imy*imz;
        logL -= 0.5*xx*noise->invC[0][0][n]*2;
        logL -= 0.5*yy*noise->invC[1][1][n]*2;
        logL -= 0.5*zz*noise->invC[2][2][n]*2;
        logL -=     xy*noise->invC[0][1][n]*2;
        logL -=     xz*noise->invC[0][2][n]*2;
        logL -=     yz*noise->invC[1][2][n]*2;
        logL -= noise->logdetC[n];
    }
    logL -= 3 * N * log2pi;
    
    return logL;
}

double noise_log_likelihood_wavelet(struct Data *data, struct Noise *noise)
{
    double logL = 0.0;
    
    struct TDI *tdi = data->tdi;
    
    int N = data->N;
    
    switch(data->Nchannel)
    {
        case 1:
            logL += -0.5*wavelet_nwip_linear(tdi->X, tdi->X, noise->invC[0][0], N);
            break;
        case 2:
            logL += -0.5*wavelet_nwip_linear(tdi->A, tdi->A, noise->invC[0][0], N);
            logL += -0.5*wavelet_nwip_linear(tdi->E, tdi->E, noise->invC[1][1], N);
            break;
        case 3:
            logL += -0.5*wavelet_nwip_linear(tdi->X, tdi->X, noise->invC[0][0], N);
            logL += -0.5*wavelet_nwip_linear(tdi->Y, tdi->Y, noise->invC[1][1], N);
            logL += -0.5*wavelet_nwip_linear(tdi->Z, tdi->Z, noise->invC[2][2], N);
            logL += -wavelet_nwip_linear(tdi->X, tdi->Y, noise->invC[0][1], N);
            logL += -wavelet_nwip_linear(tdi->X, tdi->Z, noise->invC[0][2], N);
            logL += -wavelet_nwip_linear(tdi->Y, tdi->Z, noise->invC[1][2], N);
            break;
    }
    for(int n=0; n<N; n++)
        logL -= noise->logdetC[n];
    
    return logL;
    // this comment constitutes an offering to the diety responsible for the correct factors of 2
}

double noise_delta_log_likelihood(struct Data *data, struct SplineModel *model_x, struct SplineModel *model_y, double fmin, double fmax,int ic)
{
    double dlogL = 0.0;
    
    struct TDI *tdi = data->tdi;
    struct Noise *psd_x = model_x->psd;
    struct Noise *psd_y = model_y->psd;

    int N = (int)floor((fmax-fmin)*data->T);
    int imin = (int)floor((fmin-psd_x->f[0])*data->T);
    if(imin<0)imin=0;
    
    /* remove contribution for current state x */
    dlogL -= -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_x->invC[0][0]+imin, N);
    dlogL -= -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_x->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL += psd_x->logdetC[n];

    /* add contribution for proposed state y */
    dlogL += -0.5*fourier_nwip(tdi->A+2*imin, tdi->A+2*imin, psd_y->invC[0][0]+imin, N);
    dlogL += -0.5*fourier_nwip(tdi->E+2*imin, tdi->E+2*imin, psd_y->invC[1][1]+imin, N);
    for(int n=imin; n<imin+N; n++)
        dlogL -= psd_y->logdetC[n];
    
    return dlogL;
}

void initialize_spline_model(struct Orbit *orbit, struct Data *data, struct SplineModel *model, int Nspline)
{
    
    // Initialize data models
    alloc_spline_model(model, data->NFFT, data->Nlayer, data->Nchannel, Nspline);
    
    //set max and min spline points
    model->Nmin = MIN_SPLINE_STENCIL;
    model->Nmax = Nspline;
    model->Nchannel = data->Nchannel;
    
    //set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
    {
        double f = data->fmin + (double)n/data->T;
        double Spm, Sop;
        get_noise_levels("sangria", f, &Spm, &Sop);
        model->psd->f[n] = f;
        if(model->Nchannel==2) model->psd->transfer[n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);//noise_transfer_function(f/orbit->fstar);
        else model->psd->transfer[n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
    }
    
    //divide into Nspline control points
    double logdf = (log(data->fmax) - log(data->fmin))/(Nspline-1);
    for(int i=0; i<Nspline; i++)
    {
        double f = exp(log(data->fmin) + (double)i*logdf);
        model->spline->f[i] = f;
        
        for(int n=0; n<model->Nchannel; n++)
        {
            /* initialize model to theoretical level without transfer function applied */
            for(int m=0; m<model->Nchannel; m++)
                model->spline->C[n][m][i] = 0.0;
        }
    }
    //shift first spline control point by half a bin to avoid rounding problems
    model->spline->f[0] -= 0.5/data->T;
    
    generate_spline_noise_model(model);
    model->logL = noise_log_likelihood(data, model->psd);
}

void initialize_instrument_model(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // initialize data models
    alloc_instrument_model(model, data->NFFT, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize noise levels
    /* NOTE: changed to match wavelet */
    /*
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    */
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model(orbit,model);
}

void initialize_instrument_model_wavelet(struct Orbit *orbit, struct Data *data, struct InstrumentModel *model)
{
    // wavelet basis
    struct Wavelets *wdm = data->wdm;

    // initialize data models
    alloc_instrument_model(model, data->lmax - data->lmin, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->lmin+n)*wdm->df;

    // initialize noise levels
    for(int i=0; i<model->Nlink; i++)
    {
        model->soms[i] = 1.28e-22;
        model->sacc[i] = 5.76e-30;
    }
    
    // get noise covariance matrix for initial parameters
    generate_instrument_noise_model_wavelet(wdm, orbit, model);
}

void initialize_foreground_model(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    // initialize data models
    alloc_foreground_model(model, data->NFFT, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    // initialize constant foreground parameters levels 
    model->Tobs  =  data->T;
    model->Amp   =  1.2826e-44;
    model->alpha =  1.629667;
    model->f2    =  4.810781e-4;

    // initialize time dependent foreground parameter elves
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;

    model->f1  =  pow(10., af1*log10(model->Tobs/YEAR) + bf1);
    model->fk  =  pow(10., afk*log10(model->Tobs/YEAR) + bfk);

    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model(model);

}

void default_sgwb_injection(double* params, const SGWB_t SGWB_type) {
    _Static_assert(SGWB_TEMPLATE_COUNT == 1, "Did you add an SGWB template? Edit this switch case, it needs to be exhaustive.");
    // set default values
    switch (SGWB_type) {
        case SGWB_TEMPLATE_POWERLAW:
            //params[0] = -21.0;
            //params[0] = -8.45;
            params[0] = -16.0;
            params[1] = 2./3.;
            //params[1] = 0.0;
            break;
        default:
            fprintf(stderr,"need default values for SGWB type: %s", SGWB_TEMPLATE_NAMES[SGWB_type]);
            exit(1);
            break;
    }
}


void initialize_sgwb_model(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type)
{
    // initialize data models
    alloc_sgwb_model(model, data->NFFT, data->Nlayer, data->Nchannel, SGWB_type);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = data->fmin + (double)n/data->T;

    default_sgwb_injection(model->params, SGWB_type);
    
    model->Tobs  =  data->T;
    // get covariance matrix for initial parameters
    generate_sgwb_model(model);

}

void initialize_sgwb_model_wavelet(struct Orbit *orbit, struct Data *data, struct SGWBModel *model, SGWB_t SGWB_type)
{
    struct Wavelets* wdm = data->wdm;

    // initialize data models
    alloc_sgwb_model(model, data->lmax - data->lmin, data->Nlayer, data->Nchannel, SGWB_type);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->lmin+n)*wdm->df;

    default_sgwb_injection(model->params, SGWB_type);

    
    model->Tobs  =  data->T;
    // get covariance matrix for initial parameters
    generate_sgwb_model_wavelet(wdm, model);

}

void initialize_foreground_model_wavelet(struct Orbit *orbit, struct Data *data, struct ForegroundModel *model)
{
    //wavelet basis
    struct Wavelets *wdm = data->wdm;

    // initialize data models
    alloc_foreground_model(model, data->lmax - data->lmin, data->Nlayer, data->Nchannel);
    
    // set up psd frequency grid
    for(int n=0; n<model->psd->N; n++)
        model->psd->f[n] = (data->lmin+n)*wdm->df;

    // initialize constant foreground parameters levels 
    model->Tobs  =  data->T;
    model->Amp   =  1.2826e-44;
    //model->Amp   =  1.2826e-53;
    model->alpha =  1.629667;
    model->f2    =  4.810781e-4;

    // initialize time dependent foreground parameter elves
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;

    model->f1  =  pow(10., af1*log10(model->Tobs/YEAR) + bf1);
    model->fk  =  pow(10., afk*log10(model->Tobs/YEAR) + bfk);
    map_foreground_params_to_array(model);
    
    // get noise covariance matrix for initial parameters
    generate_galactic_foreground_model_wavelet(wdm,model);

    // get galaxy modulation
    model->modulation = malloc(sizeof(struct GalaxyModulation));
    // TODO: discretize and save this in Galaxy class so we don't have to recalcuate it all the time
    initialize_galaxy_modulation(model->modulation, data->wdm, orbit, data->T, data->t0);
    
    /**************************************************
     * Compute galaxy modulation
    **************************************************/
    // TODO: use Pozzoli et al's analytic response for modulation
    
    double *galaxy_params = double_vector(6); // defines galaxy shape
    galaxy_params[0] = 0.25; // A 0.25    bulge fraction
    galaxy_params[1] = 0.8;  // Rb 0.8    bulge radius (kpc)
    galaxy_params[2] = 2.5;  // Rd 2.5    disk radius (kpc)
    galaxy_params[3] = 0.4;  // Zd 0.4    disk height (kpc)
    galaxy_params[4] = 7.2;  // RGC 7.2   distance from solar BC to GC (kpc)
    galaxy_params[5] = 3.5;  // Rcut 3.5  radius out to which all sources are found (kpc)

    //computes the modulation of the confusion noise
    galaxy_modulation(model->modulation, galaxy_params);
    free_double_vector(galaxy_params);
}

void GetDynamicNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    /**************************************************
     * Compute instrument noise levels
    **************************************************/
    struct InstrumentModel *inst_noise = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model_wavelet(orbit, data, inst_noise);

    /**************************************************
     * Compute galactic foreground noise levels
    **************************************************/
    struct ForegroundModel *conf_noise = malloc(sizeof(struct ForegroundModel));
    initialize_foreground_model_wavelet(orbit, data, conf_noise);

    /**************************************************
     * Compute SGWB noise levels
    **************************************************/
    struct SGWBModel *sgwb = malloc(sizeof(struct SGWBModel));
    initialize_sgwb_model_wavelet(orbit, data, sgwb, flags->sgwbTemplate);

    /**************************************************
     * Combine noise components
    **************************************************/
    generate_full_dynamic_covariance_matrix(data->wdm, inst_noise, conf_noise, sgwb, data->noise);
    invert_noise_covariance_matrix(data->noise);

    char filename[128];
    sprintf(filename,"%s/power_noise.dat",data->dataDir);
    FILE *fptr=fopen(filename,"w");
    int k;
    for(int j=data->lmin; j<data->lmax; j++)
    {
        double f = j*data->wdm->df;
        for(int i=0; i<data->wdm->NT; i++)
        {
            double t = i*data->wdm->dt;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k-=data->wdm->kmin;
            fprintf(fptr,"%lg %lg %.14e\n", t, f, data->noise->C[0][0][k]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);

    free_instrument_model(inst_noise);
    free_foreground_model(conf_noise);
    free_sgwb_model(sgwb);
}

void GetStationaryNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Noise *noise)
{
    /**************************************************
     * Compute instrument noise levels
    **************************************************/
    struct InstrumentModel *inst_noise = malloc(sizeof(struct InstrumentModel));
    initialize_instrument_model_wavelet(orbit, data, inst_noise);

    /**************************************************
     * Compute galactic foreground noise levels
    **************************************************/
    struct ForegroundModel *conf_noise = malloc(sizeof(struct ForegroundModel));
    initialize_foreground_model_wavelet(orbit, data, conf_noise);

    /**************************************************
     * Compute SGWB noise levels
    **************************************************/
    struct SGWBModel *sgwb = malloc(sizeof(struct SGWBModel));
    initialize_sgwb_model_wavelet(orbit, data, sgwb, flags->sgwbTemplate);


    /**************************************************
     * Combine noise components
    **************************************************/
    generate_full_stationary_covariance_matrix(data->wdm, inst_noise, conf_noise, sgwb, noise);
    invert_noise_covariance_matrix(noise);

    char filename[128];
    sprintf(filename,"%s/power_stationary_noise.dat",data->dataDir);
    FILE *fptr=fopen(filename,"w");
    int k;
    // TODO: check bounds here
    for(int j=data->lmin; j<data->lmax; j++)
    {
        double f = j*data->wdm->df;
        for(int i=0; i<data->wdm->NT; i++)
        {
            double t = i*data->wdm->dt;
            wavelet_pixel_to_index(data->wdm,i,j,&k);
            k-=data->wdm->kmin;
            fprintf(fptr,"%lg %lg %.14e\n", t, f,noise->C[0][0][k]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);

    free_instrument_model(inst_noise);
    free_foreground_model(conf_noise);
}
 
