/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish 
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

#include "glass_ucb.h"

#define FIXME 0
#define find_max(x,y) (((x) >= (y)) ? (x) : (y))
#define find_min(x,y) (((x) <= (y)) ? (x) : (y))

void map_array_to_ucb_params(struct Source *source, double *params, double T)
{
    source->f0       = params[0]/T;
    source->costheta = params[1];
    source->phi      = params[2];
    source->amp      = exp(params[3]);
    source->cosi     = params[4];
    source->psi      = params[5];
    source->phiref   = params[6];
    if(UCB_MODEL_NP>7)
        source->dfdt   = params[7]/(T*T);
    if(UCB_MODEL_NP>8)
        source->d2fdt2 = params[8]/(T*T*T);
}

void map_ucb_params_to_array(struct Source *source, double *params, double T)
{
    params[0] = source->f0*T;
    params[1] = source->costheta;
    params[2] = source->phi;
    params[3] = log(source->amp);
    params[4] = source->cosi;
    params[5] = source->psi;
    params[6] = source->phiref;
    if(UCB_MODEL_NP>7)
        params[7] = source->dfdt*T*T;
    if(UCB_MODEL_NP>8)
        params[8] = source->d2fdt2*T*T*T;
}

void generate_ucb_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    int i,n;

    for(n=0; n<data->N; n++)
    {
        model->tdi->X[n]=0.0;
        model->tdi->Y[n]=0.0;
        model->tdi->Z[n]=0.0;
        model->tdi->A[n]=0.0;
        model->tdi->E[n]=0.0;
    }
    
    //Loop over signals in model
    for(n=0; n<model->Nlive; n++)
    {
        struct Source *source = model->source[n];
        
        if(source_id==-1 || source_id==n)
        {

            for(i=0; i<data->N; i++)
            {
                source->tdi->X[i]=0.0;
                source->tdi->Y[i]=0.0;
                source->tdi->Z[i]=0.0;
                source->tdi->A[i]=0.0;
                source->tdi->E[i]=0.0;
            }

            map_array_to_ucb_params(source, source->params, data->T);
            

            //Book-keeping of injection time-frequency volume
            ucb_alignment(orbit, data, source);

        }
        
        //Simulate gravitational wave signal
        /* the source_id = -1 condition is redundent if the model->tdi structure is up to date...*/

        if(source_id==-1 || source_id==n) ucb_waveform(orbit, data->format, data->T, model->t0, source->params, UCB_MODEL_NP, source->tdi->X, source->tdi->Y, source->tdi->Z, source->tdi->A, source->tdi->E, source->BW, source->tdi->Nchannel);
        
        //Add waveform to model TDI channels
        add_ucb_model(data,model,source);

        
    }//loop over sources
}

void generate_ucb_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    int i,n;
    
    for(n=0; n<data->N; n++)
    {
        model->tdi->X[n]=0.0;
        model->tdi->Y[n]=0.0;
        model->tdi->Z[n]=0.0;
        model->tdi->A[n]=0.0;
        model->tdi->E[n]=0.0;
        model->list[n]=0;
    }
    
    //Loop over signals in model
    for(n=0; n<model->Nlive; n++)
    {
        struct Source *source = model->source[n];
        
        if(source_id==-1 || source_id==n)
        {
            for(i=0; i<data->N; i++)
            {
                source->tdi->X[i]=0.0;
                source->tdi->Y[i]=0.0;
                source->tdi->Z[i]=0.0;
                source->tdi->A[i]=0.0;
                source->tdi->E[i]=0.0;
                source->list[i]=0;
            }
            
            map_array_to_ucb_params(source, source->params, data->T);
        }
        
        //Simulate gravitational wave signal
        /* the source_id = -1 condition is redundent if the model->tdi structure is up to date...*/
        if(source_id==-1 || source_id==n) ucb_waveform_wavelet(orbit, data->wdm, data->T, model->t0, source->params, source->list, &source->Nlist, source->tdi->X, source->tdi->Y, source->tdi->Z);

        //Add waveform to model TDI channels
        add_ucb_model_wavelet(data,model,source);

    }//loop over sources
}

void remove_ucb_model(struct Data *data, struct Model *model, struct Source *source)
{
    //subtract current nth source from  model
    for(int i=0; i<source->BW; i++)
    {
        int j = i+source->imin;
        
        if(j>-1 && j<data->NFFT)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            switch(data->Nchannel)
            {
                case 1:
                    model->tdi->X[j_re] -= source->tdi->X[i_re];
                    model->tdi->X[j_im] -= source->tdi->X[i_im];
                    break;
                case 2:
                    model->tdi->A[j_re] -= source->tdi->A[i_re];
                    model->tdi->A[j_im] -= source->tdi->A[i_im];
                    
                    model->tdi->E[j_re] -= source->tdi->E[i_re];
                    model->tdi->E[j_im] -= source->tdi->E[i_im];
                    break;
                case 3:
                    model->tdi->X[j_re] -= source->tdi->X[i_re];
                    model->tdi->X[j_im] -= source->tdi->X[i_im];
                    
                    model->tdi->Y[j_re] -= source->tdi->Y[i_re];
                    model->tdi->Y[j_im] -= source->tdi->Y[i_im];

                    model->tdi->Z[j_re] -= source->tdi->Z[i_re];
                    model->tdi->Z[j_im] -= source->tdi->Z[i_im];
                    break;
            }
        }//check that source_id is in range
    }//loop over waveform bins
}


void add_ucb_model(struct Data *data, struct Model *model, struct Source *source)
{
    //add current nth source to  model
    for(int i=0; i<source->BW; i++)
    {
        int j = i+source->imin;
        
        if(j>-1 && j<data->NFFT)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            switch(data->Nchannel)
            {
                case 1:
                    model->tdi->X[j_re] += source->tdi->X[i_re];
                    model->tdi->X[j_im] += source->tdi->X[i_im];
                    break;
                case 2:
                    model->tdi->A[j_re] += source->tdi->A[i_re];
                    model->tdi->A[j_im] += source->tdi->A[i_im];
                    
                    model->tdi->E[j_re] += source->tdi->E[i_re];
                    model->tdi->E[j_im] += source->tdi->E[i_im];
                    break;
                case 3:
                    model->tdi->X[j_re] += source->tdi->X[i_re];
                    model->tdi->X[j_im] += source->tdi->X[i_im];
                    
                    model->tdi->Y[j_re] += source->tdi->Y[i_re];
                    model->tdi->Y[j_im] += source->tdi->Y[i_im];

                    model->tdi->Z[j_re] += source->tdi->Z[i_re];
                    model->tdi->Z[j_im] += source->tdi->Z[i_im];
                    break;
            }
        }//check that source_id is in range
    }//loop over waveform bins
}
void add_ucb_model_wavelet(struct Data *data, struct Model *model, struct Source *source)
{
    //insert source into model 
    for(int n=0; n<source->Nlist; n++)
    {
        int k = source->list[n];
        model->tdi->X[k] += source->tdi->X[k];
        model->tdi->Y[k] += source->tdi->Y[k];
        model->tdi->Z[k] += source->tdi->Z[k];
    }

    //get union of list
    list_union(model->list, source->list, model->Nlist, source->Nlist, model->list, &model->Nlist); 
}

void remove_ucb_model_wavelet(struct Data *data, struct Model *model, struct Source *source)
{
    //insert source into model 
    for(int n=0; n<source->Nlist; n++)
    {
        int k=source->list[n];
        if(k>=0 && k<data->N)
        {
            model->tdi->X[k] -= source->tdi->X[k];
            model->tdi->Y[k] -= source->tdi->Y[k];
            model->tdi->Z[k] -= source->tdi->Z[k];
        }
    }

    //get union of list
    if(model->Nlive == 0) model->Nlist = 0;
    else
    {
        for(int n=0; n<model->Nlive; n++)
            list_union(model->list, model->source[n]->list, model->Nlist, model->source[n]->Nlist, model->list, &model->Nlist); 
    }
}

void update_ucb_model(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id)
{
    int i;
    struct Source *source_x = model_x->source[source_id];
    struct Source *source_y = model_y->source[source_id];
    
    //subtract current nth source from  model
    remove_ucb_model(data,model_y,source_x);

    //generate proposed signal model
    for(i=0; i<data->N; i++)
    {
        switch(data->Nchannel)
        {
            case 1:
                source_y->tdi->X[i]=0.0;
                break;
            case 2:
                source_y->tdi->A[i]=0.0;
                source_y->tdi->E[i]=0.0;
                break;
            case 3:
                source_y->tdi->X[i]=0.0;
                source_y->tdi->Y[i]=0.0;
                source_y->tdi->Z[i]=0.0;
                break;
        }
    }

    map_array_to_ucb_params(source_y, source_y->params, data->T);
    ucb_alignment(orbit, data, source_y);

    ucb_waveform(orbit, data->format, data->T, model_y->t0, source_y->params, UCB_MODEL_NP, source_y->tdi->X,source_y->tdi->Y,source_y->tdi->Z, source_y->tdi->A, source_y->tdi->E, source_y->BW, source_y->tdi->Nchannel);

    //add proposed nth source to model
    add_ucb_model(data,model_y,source_y);
    
}

void update_ucb_model_wavelet(struct Orbit *orbit, struct Data *data, struct Model *model_x, struct Model *model_y, int source_id)
{
    int i;
    int N=data->N;
    struct Source *source_x = model_x->source[source_id];
    struct Source *source_y = model_y->source[source_id];
    
    //subtract current nth source from  model
    remove_ucb_model_wavelet(data,model_y,source_x);

    //generate proposed signal model
    for(i=0; i<N; i++)
    {
        source_y->tdi->X[i]=0.0;
        source_y->tdi->Y[i]=0.0;
        source_y->tdi->Z[i]=0.0;
    }

    map_array_to_ucb_params(source_y, source_y->params, data->T);
    ucb_waveform_wavelet(orbit, data->wdm, data->T, model_y->t0, source_y->params, source_y->list, &source_y->Nlist, source_y->tdi->X, source_y->tdi->Y, source_y->tdi->Z);

    //add proposed nth source to model
    add_ucb_model_wavelet(data,model_y,source_y);
    
}

void generate_noise_model(struct Data *data, struct Model *model)
{
    for(int n=0; n<data->NFFT; n++)
    {
        for(int i=0; i<data->Nchannel; i++)
            for(int j=i; j<data->Nchannel; j++)
                model->noise->C[i][j][n] = model->noise->C[j][i][n] = data->noise->C[i][j][n]*sqrt(model->noise->eta[i]*model->noise->eta[j]);

    }
    invert_noise_covariance_matrix(model->noise);
}

void generate_noise_model_wavelet(struct Data *data, struct Model *model)
{
    int Nlayers = data->Nlayer;         //number of frequency layers
    int Nslices = data->N/data->Nlayer; //number of time slices
    for(int n=0; n<Nlayers; n++)
    {
        for(int m=0; m<Nslices; m++)
        {
            int k = n*Nslices+m;
            for(int i=0; i<data->Nchannel; i++)
            {
                for(int j=i; j<data->Nchannel; j++)
                {
                    model->noise->C[i][j][k] = model->noise->C[j][i][k] = data->noise->C[i][j][k]*sqrt(model->noise->eta[i*Nlayers+n]*model->noise->eta[j*Nlayers+n]);
                }
            }
        }

    }
    invert_noise_covariance_matrix(model->noise);
}

void maximize_ucb_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    if(source_id < model->Nlive)
    {
        double *Fparams = calloc(UCB_MODEL_NP,sizeof(double));
        
        struct Source *source = model->source[source_id];
        
        /* save original data */
        //    struct TDI *data_save = malloc(sizeof(struct TDI));
        //    alloc_tdi(data_save, data->N, data->Nchannel);
        //    copy_tdi(data->tdi[FIXME],data_save);
        //
        //    /* save original noise */
        //    struct Noise *noise_save = malloc(sizeof(struct Noise));
        //    alloc_noise(noise_save, data->N);
        //    copy_noise(data->noise[FIXME], noise_save);
        //
        //    /* put current noise model in data structure for get_Fstat_logL() */
        //    copy_noise(model->noise[FIXME],data->noise[FIXME]);
        
        
        /* create residual of all sources but n for F-statistic */
        //    for(int m=0; m<model->Nlive; m++)
        //    {
        //        if(m!=source_id)
        //        {
        //            for(int i=0; i<data->N*2; i++)
        //            {
        //                data->tdi[FIXME]->A[i] -= model->source[m]->tdi->A[i];
        //                data->tdi[FIXME]->E[i] -= model->source[m]->tdi->E[i];
        //                data->tdi[FIXME]->X[i] -= model->source[m]->tdi->X[i];
        //            }
        //        }
        //    }
        
        
        if(!check_range(source->params, model->prior))
        {
            /* maximize parameters w/ F-statistic */
            get_Fstat_xmax(orbit, data, source->params, Fparams);
            
            /* unpack maximized parameters */
            source->amp    = exp(Fparams[3]);
            source->cosi   = Fparams[4];
            source->psi    = Fparams[5];
            source->phiref = Fparams[6];
            map_ucb_params_to_array(source, source->params, data->T);
        }
        
        /* restore original data */
        //    copy_tdi(data_save,data->tdi[FIXME]);
        //    free_tdi(data_save);
        //
        //    /* restore original noise */
        //    copy_noise(noise_save,data->noise[FIXME]);
        //    free_noise(noise_save);
        free(Fparams);
    }
}

void append_ucb_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel)
{
    //malloc source structure for this sample's entry
    entry->source[entry->Nchain] = malloc(sizeof(struct Source));
    
    //only copy source parameters
    entry->source[entry->Nchain]->params=calloc(UCB_MODEL_NP,sizeof(double));
    memcpy(entry->source[entry->Nchain]->params, sample->params, UCB_MODEL_NP*sizeof(double));
    
    //need Tobs which isn't stored in entries
    double T = sample->params[0]/sample->f0;
    
    //get physical parameters for waveform calculations later
    map_array_to_ucb_params(entry->source[entry->Nchain], entry->source[entry->Nchain]->params, T);

    //increment number of stored samples for this entry
    entry->Nchain++;
}

int ucb_gmm_wrapper(double **ranges, struct Flags *flags, struct Entry *entry, char *outdir, size_t NMODE, size_t NTHIN, unsigned int *seed, double *BIC)
{
    if(flags->verbose)fprintf(stdout,"Event %s, NMODE=%i\n",entry->name,(int)NMODE);
    
    // number of samples
    size_t NMCMC = entry->Nchain;
    
    // number of EM iterations
    size_t NSTEP = 100;
    
    // thin chain
    NMCMC /= NTHIN;
    
    struct Sample **samples = malloc(NMCMC*sizeof(struct Sample*));
    for(size_t n=0; n<NMCMC; n++)
    {
        samples[n] = malloc(sizeof(struct Sample));
        samples[n]->x = double_vector(UCB_MODEL_NP);
        samples[n]->p = double_vector(NMODE);
        samples[n]->w = double_vector(NMODE);
    }
    
    // covariance matrices for different modes
    struct MVG **modes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t n=0; n<NMODE; n++)
    {
        modes[n] = malloc(sizeof(struct MVG));
        alloc_MVG(modes[n],UCB_MODEL_NP);
    }
    
    // Logistic mapping of samples onto R
    double y;
    double pmin,pmax;
    double *y_vec = double_vector(NMCMC);
    
    /* parse chain file */
    double **params = malloc(UCB_MODEL_NP*sizeof(double *));
    for(size_t n=0; n<UCB_MODEL_NP; n++) params[n] = double_vector(NMCMC);
    double value[UCB_MODEL_NP];
    for(size_t i=0; i<NMCMC; i++)
    {
        value[0] = entry->source[i*NTHIN]->f0;
        value[1] = entry->source[i*NTHIN]->costheta;
        value[2] = entry->source[i*NTHIN]->phi;
        value[3] = log(entry->source[i*NTHIN]->amp);
        value[4] = entry->source[i*NTHIN]->cosi;
        value[5] = entry->source[i*NTHIN]->psi;
        value[6] = entry->source[i*NTHIN]->phiref;
        if(UCB_MODEL_NP>7)
            value[7] = entry->source[i*NTHIN]->dfdt;
        if(UCB_MODEL_NP>8)
            value[8] = entry->source[i*NTHIN]->d2fdt2;
        
        for(size_t n=0; n<UCB_MODEL_NP; n++)
        {
            params[n][i] = value[n];
        }
    }
    
    /* Use priors to set min and max of each parameter*/
    for(size_t n=0; n<UCB_MODEL_NP; n++)
    {
        // copy max and min into each MVG structure
        for(size_t k=0; k<NMODE; k++)
        {
            modes[k]->minmax[n][0] = ranges[n][0];
            modes[k]->minmax[n][1] = ranges[n][1];
        }
    }
    
    
    /* map params to R with logit function */
    for(size_t n=0; n<UCB_MODEL_NP; n++)
    {
        pmin = modes[0]->minmax[n][0];;
        pmax = modes[0]->minmax[n][1];
        logit_mapping(params[n], y_vec, pmin, pmax, NMCMC);
        
        for(size_t i=0; i<NMCMC; i++)
        {
            y = y_vec[i];
            samples[i]->x[n] = y;
        }
    }
    
    /* The main Gaussian Mixture Model with Expectation Maximization function */
    double logL;
    if(GMM_with_EM(modes,samples,UCB_MODEL_NP,NMODE,NMCMC,NSTEP,seed,&logL,BIC)) return 1;
    
    /* Write GMM results to binary for pick up by other processes */
    char filename[BUFFER_SIZE];
    sprintf(filename,"%s/%s_gmm.bin",outdir,entry->name);
    FILE *fptr = fopen(filename,"wb");
    fwrite(&NMODE, sizeof NMODE, 1, fptr);
    for(size_t n=0; n<NMODE; n++) write_MVG(modes[n],fptr);
    fclose(fptr);
    
    /* print 1D PDFs and 2D contours of GMM model */
    if(flags->verbose) print_model(modes, samples, UCB_MODEL_NP, NMODE, NMCMC, logL, *BIC, NMODE);
    
    /* clean up */
    for(size_t n=0; n<NMCMC; n++)
    {
        free_double_vector(samples[n]->x);
        free_double_vector(samples[n]->p);
        free_double_vector(samples[n]->w);
        free(samples[n]);
    }
    free(samples);
    
    // covariance matrices for different modes
    for(size_t n=0; n<NMODE; n++) free_MVG(modes[n]);
    free(modes);
    
    return 0;
}
