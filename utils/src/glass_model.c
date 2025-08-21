/*
 * Copyright 2025 Tyson B. Littenberg 
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

#include "glass_utils.h"


void alloc_model(struct Data *data, struct Model *model, int NP, int Nmax)
{
    int n;
    model->NP = NP;
    model->Nlive = 1;
    model->Nmax  = Nmax;
    model->Neff  = 2;
    model->t0 = 0.0;

    //Wavelet bookkeeping
    model->Nlist = 0;
    model->list = calloc(data->N,sizeof(int));

    model->source = malloc(model->Nmax*sizeof(struct Source *));
    
    model->calibration = malloc(sizeof(struct Calibration) );
    model->noise       = malloc(sizeof(struct Noise)       );
    model->tdi         = malloc(sizeof(struct TDI)         );
    model->residual    = malloc(sizeof(struct TDI)         );
    
    if(!strcmp(data->basis,"fourier")) alloc_noise(model->noise,data->NFFT, data->Nlayer, data->Nchannel);
    if(!strcmp(data->basis,"wavelet")) alloc_noise(model->noise,data->N, data->Nlayer, data->Nchannel);

    alloc_tdi(model->tdi, data->N, data->Nchannel);
    alloc_tdi(model->residual, data->N, data->Nchannel);
    alloc_calibration(model->calibration);
    
    for(n=0; n<model->Nmax; n++)
    {
        model->source[n] = malloc(sizeof(struct Source));
        alloc_source(model->source[n],data->N,NP,data->Nchannel);
    }
    
    //TODO: move this out of Model structure into Prior structure, create a alloc_prior() and initialize_prior() function
    model->logPriorVolume = calloc(NP,sizeof(double));
    model->prior = malloc(NP*sizeof(double *));
    for(n=0; n<NP; n++) model->prior[n] = calloc(2,sizeof(double));

}

void copy_model(struct Model *origin, struct Model *copy)
{
    //Source parameters
    copy->NP    = origin->NP;
    copy->Nmax  = origin->Nmax;
    copy->Neff  = origin->Neff;
    copy->Nlive = origin->Nlive;
    for(int n=0; n<origin->Nlive; n++)
        copy_source(origin->source[n],copy->source[n]);
    
    //Noise parameters
    copy_noise(origin->noise,copy->noise);
    
    //TDI
    copy_tdi(origin->tdi,copy->tdi);
    
    //Calibration parameters
    copy_calibration(origin->calibration,copy->calibration);
    
    //Residual
    copy_tdi(origin->residual,copy->residual);
    
    //Start time for segment for model
    copy->t0 = origin->t0;
    
    
    //Source parameter priors
    for(int n=0; n<copy->NP; n++)
    {
        for(int j=0; j<2; j++) copy->prior[n][j] = origin->prior[n][j];
        copy->logPriorVolume[n] = origin->logPriorVolume[n];
    }
    //Model likelihood
    copy->logL           = origin->logL;
    copy->logLnorm       = origin->logLnorm;

    //Wavelet bookkeeping
    copy->Nlist = origin->Nlist;
    memcpy(copy->list,origin->list,origin->Nlist*sizeof(int));

}

void copy_model_lite(struct Model *origin, struct Model *copy)
{
    copy->Nlive = origin->Nlive;

    //Source parameters
    for(int n=0; n<origin->Nlive; n++)
        copy_source(origin->source[n],copy->source[n]);
    
    //Source waveforms
    copy_tdi(origin->tdi,copy->tdi);
    copy_tdi(origin->residual,copy->residual);
    
    //Model likelihood
    copy->logL = origin->logL;

    //Wavelet bookkeeping
    copy->Nlist = origin->Nlist;
    memcpy(copy->list,origin->list,origin->Nlist*sizeof(int));
}

int compare_model(struct Model *a, struct Model *b)
{
    
    //Source parameters
    if(a->NP    != b->NP)    return 1;  //number of parameters in singal model
    if(a->Nmax  != b->Nmax)  return 1;  //maximum number of signals in model
    if(a->Nlive != b->Nlive) return 1;  //current number of signals in model
    
    //struct Source **source;
    for(int i=0; i<a->Nlive; i++)
    {
        struct Source *sa = a->source[i];
        struct Source *sb = b->source[i];
        
        //Intrinsic
        if(sa->m1 != sb->m1) return 1;
        if(sa->m2 != sb->m2) return 1;
        if(sa->f0 != sb->f0) return 1;
        
        //Extrinisic
        if(sa->psi    != sb->psi)  return 1;
        if(sa->cosi   != sb->cosi) return 1;
        if(sa->phiref != sb->phiref) return 1;
        
        if(sa->dL       != sb->dL)       return 1;
        if(sa->phi      != sb->phi)      return 1;
        if(sa->costheta != sb->costheta) return 1;
        
        //Derived
        if(sa->amp    != sb->amp)    return 1;
        if(sa->dfdt   != sb->dfdt)   return 1;
        if(sa->d2fdt2 != sb->d2fdt2) return 1;
        if(sa->Mchirp != sb->Mchirp) return 1;
        
        //Book-keeping
        if(sa->BW   != sb->BW)   return 1;
        if(sa->qmin != sb->qmin) return 1;
        if(sa->qmax != sb->qmax) return 1;
        if(sa->imin != sb->imin) return 1;
        if(sa->imax != sb->imax) return 1;
        
        //Response
        //TDI
        struct TDI *tsa = sa->tdi;
        struct TDI *tsb = sb->tdi;
        
        if(tsa->N != tsb->N) return 1;
        if(tsa->Nchannel != tsb->Nchannel) return 1;
        
        for(int i=0; i<2*tsa->N; i++)
        {
            
            //Michelson
            if(tsa->X[i] != tsb->X[i]) return 1;
            if(tsa->Y[i] != tsb->Y[i]) return 1;
            if(tsa->Z[i] != tsb->Z[i]) return 1;
            
            //Noise-orthogonal
            if(tsa->A[i] != tsb->A[i]) return 1;
            if(tsa->E[i] != tsb->E[i]) return 1;
            if(tsa->T[i] != tsb->T[i]) return 1;
        }
        
        //Package parameters for waveform generator
        for(int j=0; j<a->NP; j++) if(sa->params[j] != sb->params[j]) return 1;
        
    }
    
    //Noise parameters
    struct Noise *na = a->noise;
    struct Noise *nb = b->noise;
    
    if(na->N != nb->N) return 1;
    if(na->Nchannel != nb->Nchannel) return 1;

    for(int n=0; n<na->Nchannel; n++)
    {
        if(na->eta[n] != nb->eta[n]) return 1;
        for(int i=0; i<na->N; i++) if(na->detC[i] != nb->detC[i]) return 1;
        for(int m=n; m<na->Nchannel; m++)
        {
            for(int i=0; i<na->N; i++)
            {
                if(na->C[n][m][i] != nb->C[n][m][i]) return 1;
                if(na->invC[n][m][i] != nb->invC[n][m][i]) return 1;
            }
        }
    }
    
    //TDI
    struct TDI *ta = a->tdi;
    struct TDI *tb = b->tdi;
    
    if(ta->N != tb->N) return 1;
    if(ta->Nchannel != tb->Nchannel) return 1;
    
    for(int i=0; i<2*ta->N; i++)
    {
        
        //Michelson
        if(ta->X[i] != tb->X[i]) return 1;
        if(ta->Y[i] != tb->Y[i]) return 1;
        if(ta->Z[i] != tb->Z[i]) return 1;
        
        //Noise-orthogonal
        if(ta->A[i] != tb->A[i]) return 1;
        if(ta->E[i] != tb->E[i]) return 1;
        if(ta->T[i] != tb->T[i]) return 1;
    }
    
    //Start time for segment for model
    if(a->t0 != b->t0)     return 1;
    
    //Source parameter priors
    //double **prior;
    if(a->logPriorVolume != b->logPriorVolume) return 1;
    
    //Model likelihood
    if(a->logL     != b->logL)     return 1;
    if(a->logLnorm != b->logLnorm) return 1;
    
    return 0;
}



void free_model(struct Model *model)
{
    int n;
    for(n=0; n<model->Nmax; n++)
    {
        free_source(model->source[n]);
    }
    free(model->source);

    for(n=0; n<model->NP; n++) free(model->prior[n]);
    free(model->prior);
    free(model->logPriorVolume);

    free_tdi(model->tdi);
    free_tdi(model->residual);
    free_noise(model->noise);
    free_calibration(model->calibration);

    free(model->list);
    
    free(model);
}

void alloc_source(struct Source *source, int N, int NP, int Nchannel)
{
    source->NP = NP;
    
    //Intrinsic
    source->m1=1.;
    source->m2=1.;
    source->f0=0.;
    
    //Extrinisic
    source->psi=0.0;
    source->cosi=0.0;
    source->phiref=0.0;
    
    source->dL=1.0;
    source->phi=0.0;
    source->costheta=0.0;
    
    //Derived
    source->amp=1.;
    source->Mchirp=1.;
    source->dfdt=0.;
    source->d2fdt2=0.;
    
    //Book-keeping
    source->BW   = N/2;
    source->qmin = 0;
    source->qmax = N/2;
    source->imin = 0;
    source->imax = N/2;
    
    
    //Package parameters for waveform generator
    source->params=calloc(NP,sizeof(double));
    
    //Response
    source->tdi = malloc(sizeof(struct TDI));
    alloc_tdi(source->tdi,N, Nchannel);
    
    //Fisher
    source->fisher_matrix = malloc(NP*sizeof(double *));
    source->fisher_evectr = malloc(NP*sizeof(double *));
    source->fisher_evalue = calloc(NP,sizeof(double));
    for(int i=0; i<NP; i++)
    {
        source->fisher_matrix[i] = calloc(NP,sizeof(double));
        source->fisher_evectr[i] = calloc(NP,sizeof(double));
    }

    //Wavelet bookkeeping
    source->Nlist = 0;
    source->list = calloc(N,sizeof(int));
};

void copy_source(struct Source *origin, struct Source *copy)
{
    copy->NP = origin->NP;
    
    //Intrinsic
    copy->m1     = origin->m1;
    copy->m2     = origin->m2;
    copy->Mchirp = origin->Mchirp;
    copy->Mtotal = origin->Mtotal;
    copy->chi1   = origin->chi1;
    copy->chi2   = origin->chi2;

    //Phenomenological
    copy->f0 = origin->f0;
    
    //Extrinisic
    copy->psi      = origin->psi;
    copy->cosi     = origin->cosi;
    copy->phiref   = origin->phiref;
    copy->phi      = origin->phi;
    copy->costheta = origin->costheta;
    copy->dL       = origin->dL;
    copy->tref     = origin->tref;

    //Derived
    copy->f0     = origin->f0;
    copy->dfdt   = origin->dfdt;
    copy->d2fdt2 = origin->d2fdt2;
    copy->amp    = origin->amp;

    //Book-keeping
    copy->BW   = origin->BW;
    copy->qmin = origin->qmin;
    copy->qmax = origin->qmax;
    copy->imin = origin->imin;
    copy->imax = origin->imax;
    
    //Response
    copy_tdi(origin->tdi,copy->tdi);

    
    //Fisher
    memcpy(copy->fisher_evalue, origin->fisher_evalue, copy->NP*sizeof(double));
    memcpy(copy->params, origin->params, copy->NP*sizeof(double));
    copy->fisher_update_flag = origin->fisher_update_flag;
    
    for(int i=0; i<copy->NP; i++)
    {
        memcpy(copy->fisher_matrix[i], origin->fisher_matrix[i], copy->NP*sizeof(double));
        memcpy(copy->fisher_evectr[i], origin->fisher_evectr[i], copy->NP*sizeof(double));
    }

    copy->Nlist = origin->Nlist;
    memcpy(copy->list, origin->list, origin->Nlist*sizeof(int));
    
}

void free_source(struct Source *source)
{
    for(int i=0; i<source->NP; i++)
    {
        free(source->fisher_matrix[i]);
        free(source->fisher_evectr[i]);
    }
    free(source->fisher_matrix);
    free(source->fisher_evectr);
    free(source->fisher_evalue);
    free(source->params);
    
    free_tdi(source->tdi);

    free(source->list);
    
    free(source);
}

double delta_log_likelihood(struct Data *data, struct Model *model)
{
    /*
    *
    * Update residual and only sum over affected bins
    *
    */
      
    double dlogL=0.0;
    int N = data->N;
    struct TDI *d = data->tdi;
    struct TDI *h = model->tdi;
    struct Noise *n = model->noise;
    //int *list = model->list;
    int *list = int_vector(data->N);
    for(int n=0; n<data->N; n++) list[n]=n;

    // add the -2(d|h) contribution to the likelihood over the pixel list
    dlogL += -2.0*wavelet_nwip(d->X, h->X, n->invC[0][0], list, N);
    dlogL += -2.0*wavelet_nwip(d->Y, h->Y, n->invC[1][1], list, N);
    dlogL += -2.0*wavelet_nwip(d->Z, h->Z, n->invC[2][2], list, N);
    dlogL += -2.0*wavelet_nwip(d->X, h->Y, n->invC[0][1], list, N)*2;
    dlogL += -2.0*wavelet_nwip(d->X, h->Z, n->invC[0][2], list, N)*2;
    dlogL += -2.0*wavelet_nwip(d->Y, h->Z, n->invC[1][2], list, N)*2;

    // add the (h|h) conribution to the likelihoood over the pixel list
    dlogL += wavelet_nwip(h->X, h->X, n->invC[0][0], list, N);
    dlogL += wavelet_nwip(h->Y, h->Y, n->invC[1][1], list, N);
    dlogL += wavelet_nwip(h->Z, h->Z, n->invC[2][2], list, N);
    dlogL += wavelet_nwip(h->X, h->Y, n->invC[0][1], list, N)*2;
    dlogL += wavelet_nwip(h->X, h->Z, n->invC[0][2], list, N)*2;
    dlogL += wavelet_nwip(h->Y, h->Z, n->invC[1][2], list, N)*2;
    
    // the (d|d) terms cancel
    // add the (h|h) conribution to the likelihoood over the pixel list
    dlogL += wavelet_nwip(d->X, d->X, n->invC[0][0], list, N);
    dlogL += wavelet_nwip(d->Y, d->Y, n->invC[1][1], list, N);
    dlogL += wavelet_nwip(d->Z, d->Z, n->invC[2][2], list, N);
    dlogL += wavelet_nwip(d->X, d->Y, n->invC[0][1], list, N)*2;
    dlogL += wavelet_nwip(d->X, d->Z, n->invC[0][2], list, N)*2;
    dlogL += wavelet_nwip(d->Y, d->Z, n->invC[1][2], list, N)*2;

    dlogL *= -0.5;
    free_int_vector(list);

    return dlogL;
}

int update_max_log_likelihood(struct Model **model, struct Chain *chain, struct Flags *flags)
{
    int n = chain->index[0];
    int N = model[n]->Nlive;
    
    double logL = 0.0;
    double dlogL= 0.0;
    
    // get full likelihood
    logL = model[n]->logL + model[n]->logLnorm;
    
    // update max
    if(logL > chain->logLmax)
    {
        dlogL = logL - chain->logLmax;
        
        //clone chains if new mode is found (dlogL > D/2)
        if( dlogL > (double)(8*N/2) )
        {
            chain->logLmax = logL;
            
            for(int ic=1; ic<chain->NC; ic++)
            {
                int m = chain->index[ic];
                copy_model(model[n],model[m]);
            }
            if(flags->burnin)return 1;
        }
    }
    
    return 0;
}

double gaussian_log_likelihood(struct Data *data, struct Model *model)
{
    
    /*
    *
    * Form residual and sum
    *
    */
    
    double chi2 = 0.0; //chi^2, where logL = -chi^2 / 2
    
    //loop over time segments
    struct TDI *residual = model->residual;
    
    for(int i=0; i<data->N; i++)
    {
        residual->X[i] = data->tdi->X[i] - model->tdi->X[i];
        residual->Y[i] = data->tdi->Y[i] - model->tdi->Y[i];
        residual->Z[i] = data->tdi->Z[i] - model->tdi->Z[i];
        residual->A[i] = data->tdi->A[i] - model->tdi->A[i];
        residual->E[i] = data->tdi->E[i] - model->tdi->E[i];
    }
            
    switch(data->Nchannel)
    {
        case 1:
            chi2 += fourier_nwip(residual->X, residual->X, model->noise->invC[0][0], data->NFFT);
            break;
        case 2:
            chi2 += fourier_nwip(residual->A, residual->A, model->noise->invC[0][0], data->NFFT);
            chi2 += fourier_nwip(residual->E, residual->E, model->noise->invC[1][1], data->NFFT);
            break;
        case 3:
            chi2 += fourier_nwip(residual->X, residual->X, model->noise->invC[0][0], data->NFFT);
            chi2 += fourier_nwip(residual->Y, residual->Y, model->noise->invC[1][1], data->NFFT);
            chi2 += fourier_nwip(residual->Z, residual->Z, model->noise->invC[2][2], data->NFFT);

            chi2 += 2.0*fourier_nwip(residual->X, residual->Y, model->noise->invC[0][1], data->NFFT);
            chi2 += 2.0*fourier_nwip(residual->X, residual->Z, model->noise->invC[0][2], data->NFFT);
            chi2 += 2.0*fourier_nwip(residual->Y, residual->Z, model->noise->invC[1][2], data->NFFT);

            break;
        default:
            fprintf(stderr,"Unsupported number of channels in gaussian_log_likelihood()\n");
            exit(1);
    }

    return -0.5*chi2;
}

double gaussian_log_likelihood_constant_norm(struct Data *data, struct Model *model)
{
    
    double logLnorm = 0.0;
    
    //loop over time segments
    if(!strcmp(data->basis,"fourier")) logLnorm -= (double)data->NFFT*log(model->noise->detC[0]);
    if(!strcmp(data->basis,"wavelet")) logLnorm -= 0.5*(double)data->N*log(model->noise->detC[0]);
    
    return logLnorm;
}

double gaussian_log_likelihood_model_norm(struct Data *data, struct Model *model)
{
    
    double logLnorm = 0.0;
    int N;
    if(!strcmp(data->basis,"fourier")) N=data->NFFT;
    if(!strcmp(data->basis,"wavelet")) N=data->N;

    //loop over time segments
    for(int n=0; n<N; n++)
        logLnorm -= log(model->noise->detC[n]);

    if(!strcmp(data->basis,"wavelet")) logLnorm *= 0.5;  //normalization of 1/2 for wavelet domain (sum over N, not N/2)
    
    return logLnorm;
}

double gaussian_log_likelhood_wavelet(struct Data *data, struct Model *model)
{
    /*
    Form residual and sum
    */

    double chi2 = 0.0;

    struct TDI *residual = model->residual;

    
    for(int n=0; n<data->N; n++)
    {
        residual->X[n] = data->tdi->X[n];
        residual->Y[n] = data->tdi->Y[n];
        residual->Z[n] = data->tdi->Z[n];
    }
    
    for(int n=0; n<model->Nlist; n++)
    {
        int k = model->list[n];
        if(k>=0 && k<data->N)
        {
            residual->X[k] -= model->tdi->X[k];
            residual->Y[k] -= model->tdi->Y[k];
            residual->Z[k] -= model->tdi->Z[k];
        }
    }
    
    int *list = int_vector(data->N);
    for(int n=0; n<data->N; n++) list[n]=n;

    chi2 += wavelet_nwip(residual->X, residual->X, model->noise->invC[0][0], list, data->N);
    chi2 += wavelet_nwip(residual->Y, residual->Y, model->noise->invC[1][1], list, data->N);
    chi2 += wavelet_nwip(residual->Z, residual->Z, model->noise->invC[2][2], list, data->N);
    chi2 += wavelet_nwip(residual->X, residual->Y, model->noise->invC[0][1], list, data->N)*2;
    chi2 += wavelet_nwip(residual->X, residual->Z, model->noise->invC[0][2], list, data->N)*2;
    chi2 += wavelet_nwip(residual->Y, residual->Z, model->noise->invC[1][2], list, data->N)*2;

    free_int_vector(list);
    
    return -0.5*chi2;

}

void generate_calibration_model(struct Data *data, struct Model *model)
{
    struct Calibration *calibration = model->calibration;
    switch(data->Nchannel)
    {
        case 1:
            calibration->real_dphiX = cos(calibration->dphiX);
            calibration->imag_dphiX = sin(calibration->dphiX);
            break;
        case 2:
            calibration->real_dphiA = cos(calibration->dphiA);
            calibration->imag_dphiA = sin(calibration->dphiA);
            
            calibration->real_dphiE = cos(calibration->dphiE);
            calibration->imag_dphiE = sin(calibration->dphiE);
            break;
        case 3:
            calibration->real_dphiX = cos(calibration->dphiX);
            calibration->imag_dphiX = sin(calibration->dphiX);
            
            calibration->real_dphiY = cos(calibration->dphiY);
            calibration->imag_dphiY = sin(calibration->dphiY);

            calibration->real_dphiZ = cos(calibration->dphiZ);
            calibration->imag_dphiZ = sin(calibration->dphiZ);
            break;
        default:
            break;
    }
}

void apply_calibration_model(struct Data *data, struct Model *model)
{
    double dA;
    double cal_re;
    double cal_im;
    double h_re;
    double h_im;
    int i_re;
    int i_im;
    int i;
    
//apply calibration error to full signal model
    for(i=0; i<data->N; i++)
    {
        i_re = 2*i;
        i_im = i_re+1;
        
        switch(data->Nchannel)
        {
            case 1:
                h_re = model->tdi->X[i_re];
                h_im = model->tdi->X[i_im];
                
                dA     = (1.0 + model->calibration->dampX);
                cal_re = model->calibration->real_dphiX;
                cal_im = model->calibration->imag_dphiX;
                
                model->tdi->X[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->X[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                break;
            case 2:
                h_re = model->tdi->A[i_re];
                h_im = model->tdi->A[i_im];
                
                dA     = (1.0 + model->calibration->dampA);
                cal_re = model->calibration->real_dphiA;
                cal_im = model->calibration->imag_dphiA;
                
                model->tdi->A[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->A[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                
                
                h_re = model->tdi->E[i_re];
                h_im = model->tdi->E[i_im];
                
                dA     = (1.0 + model->calibration->dampE);
                cal_re = model->calibration->real_dphiE;
                cal_im = model->calibration->imag_dphiE;
                
                model->tdi->E[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->E[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                break;
            case 3:
                h_re = model->tdi->X[i_re];
                h_im = model->tdi->X[i_im];
                
                dA     = (1.0 + model->calibration->dampX);
                cal_re = model->calibration->real_dphiX;
                cal_im = model->calibration->imag_dphiX;
                
                model->tdi->X[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->X[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                
                
                h_re = model->tdi->Y[i_re];
                h_im = model->tdi->Y[i_im];
                
                dA     = (1.0 + model->calibration->dampY);
                cal_re = model->calibration->real_dphiY;
                cal_im = model->calibration->imag_dphiY;
                
                model->tdi->Y[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->Y[i_im] = dA*(h_re*cal_im + h_im*cal_re);

                h_re = model->tdi->Z[i_re];
                h_im = model->tdi->Z[i_im];
                
                dA     = (1.0 + model->calibration->dampZ);
                cal_re = model->calibration->real_dphiZ;
                cal_im = model->calibration->imag_dphiZ;
                
                model->tdi->Z[i_re] = dA*(h_re*cal_re - h_im*cal_im);
                model->tdi->Z[i_im] = dA*(h_re*cal_im + h_im*cal_re);
                break;
            default:
                break;
        }//end switch
    }//end loop over data
}
